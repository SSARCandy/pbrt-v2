
#include "stdafx.h"
#include "vdb.h"
#include "cameras/realistic.h"
#include "shapes/sphere.h"
#include <fstream>

using namespace std;

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
				 float hither, float yon, 
				 float sopen, float sclose, 
				 float filmdistance, float aperture_diameter, string specfile, 
				 float filmdiag, Film *f)
	: Camera(cam2world, sopen, sclose, f) // pbrt-v2 doesnot specify hither and yon
{
	// YOUR CODE HERE -- build and store datastructures representing the given lens
	// and film placement.
	Vector v;
	//SnellsLaw(Vector(-3, -4, 0), Vector(0, 1, 0), 0.5, 1.5, &v);
	char line[100];
	fstream fin;
	fin.open(specfile, ios::in);
	while (fin.getline(line, sizeof(line), '\n')) {
		if (line[0] == '#') continue;

		Lens l;
		sscanf(line, "%f %f %f %f", &(l.radius), &(l.axpos), &(l.n), &(l.aperture));
		lens.push_back(l);
	}

	assert(lens.size()); // no lens been defined !

	this->yon = yon;
	this->hither = hither;
	this->film = f;
	this->film_diagonal = filmdiag;
	this->film_distance = filmdistance;
	this->film_position = lens.back().axpos - filmdistance;

	this->raster_diagonal = sqrt(f->xResolution*f->xResolution +
								 f->yResolution*f->yResolution);

	this->Raster2Camera = Translate(Vector(f->xResolution / 2, f->yResolution / 2, 0))
		*Scale(film_diagonal / raster_diagonal, film_diagonal / raster_diagonal, 1)
		*Scale(-1, -1, 1);

	this->Camera2Raster = Inverse(Raster2Camera);
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
	// Sample point on lens, scale to aperture radius
	
	float lensU, lensV;
	ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
	lensU *= lens.back().aperture*0.5f;
	lensV *= lens.back().aperture*0.5f;

	// Generate raster and camera samples
	Point Praster(sample.imageX, sample.imageY, 0);
	Point Pcamera;
	Raster2Camera(Praster, &Pcamera);

	Vector rayDirection = Point(lensU, lensV, lens.back().axpos) - Pcamera;
	Ray r(Pcamera, rayDirection, INFINITY, 0);

	for (int i = lens.size() - 1; i >= 0; i--) {
		// Check intersection with this lens
		Point pHit;
		Vector normal;
		if (!LensIntersect(lens[i], r, &pHit, &normal)) {
			return 0;
		}

		// Use SnellsLaw to get new Ray direcion
		Vector rayDirection2;
		if (!SnellsLaw(rayDirection, normal, lens[i].n, (i == 0 ? 1 : lens[i - 1].n), &rayDirection2)) {
			return 0;
		}

		r = Ray(pHit, Normalize(rayDirection2), INFINITY, 0);
	}

	CameraToWorld(r, ray);
	ray->d = Normalize(ray->d);

	// Return the weight of generated Ray
	// E = A*cos^4(theta)/Z^2
	float A = M_PI* pow(lens.back().aperture / 2, 2);
	float Z = film_distance;
	float costheta = Dot(r.d, Vector(0, 0, 1));

	float E = A*pow(costheta, 4) / (Z*Z);

    return E;
}

bool RealisticCamera::LensIntersect(const Lens l, const Ray & r, Point * pHit, Vector * normal) const {
	if (l.radius != 0) {
		float th = 0;
		Vector w2oV(0.f, 0.f, l.radius - l.axpos);
		Transform o2w = Translate(-1 * w2oV);
		Transform w2o = Translate(1 * w2oV);

		Sphere sphere(&o2w, &w2o, false, abs(l.radius), -1 * l.radius, l. radius, 360);

		float rayEpsilon;
		DifferentialGeometry dg;
		if (!sphere.Intersect(r, &th, &rayEpsilon, &dg)) return false;
		if (th > r.maxt || th < r.mint) return false;

		*pHit = r(th);

	} else {
		// using absolute value because ray can come from either side
		float scale = fabs((l.axpos - r.o.z) / r.d.z);
		*pHit = r.o + scale*r.d;
	}

	// check aperture
	if ((pHit->x * pHit->x + pHit->y * pHit->y) >= (l.aperture*l.aperture) / 4) {
		return false;
	}

	// returned normal always points from the sphere surface outward
	*normal = Normalize(*pHit - Point(0, 0, l.axpos - l.radius));

	return true;
}

// https://www.wikiwand.com/en/Snell's_law#/Vector_form
bool RealisticCamera::SnellsLaw(const Vector vin, const Vector normal, const float N1, const float N2, Vector * vout) const {

	Vector l = Normalize(vin);
	Vector n = Normalize(normal);
	// r = n1/n2
	// c = dot(-normal, vin)
	// Vrefract = r x l + (r x c - (1 - mu^2 x (1 - c^2))^0.5)n
	float r = N1 / N2;
	float c = Dot(-n, l);
	int sign = c < 0 ? -1 : 1; // cos(theta) must > 0

	// total reflaction
	float tmp = 1 - r*r*(1 - c*c);
	if (tmp < 0) {
		return false;
	}

	Vector refract = r*l + (r*c - sign*sqrt(tmp))*n;

	*vout = refract;

	//Vector t = Cross(n, l);
	//t = Cross(t, n);
	//vdb_line(-t.x, -t.y, -t.z, t.x, t.y, t.z);
	//vdb_color(0, 1, 0);
	//vdb_line(l.x, l.y, l.z, -l.x, -l.y, -l.z);
	//vdb_color(1, 0, 0);
	//vdb_line(-n.x, -n.y, -n.z, n.x, n.y, n.z);
	//vdb_color(0, 0, 1);
	//vdb_line(refract.x, refract.y, refract.z, -refract.x, -refract.y, -refract.z);

	return true;
}

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	// Extract common camera parameters from \use{ParamSet}
	float hither = params.FindOneFloat("hither", -1);
	float yon = params.FindOneFloat("yon", -1);
	float shutteropen = params.FindOneFloat("shutteropen", -1);
	float shutterclose = params.FindOneFloat("shutterclose", -1);

	// Realistic camera-specific parameters
	string specfile = params.FindOneString("specfile", "");
	float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
 	float fstop = params.FindOneFloat("aperture_diameter", 1.0);	
	float filmdiag = params.FindOneFloat("filmdiag", 35.0);

	Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
		shutterclose != -1 && filmdistance!= -1);
	if (specfile == "") {
	    Severe( "No lens spec file supplied!\n" );
	}
	return new RealisticCamera(cam2world, hither, yon,
				   shutteropen, shutterclose, filmdistance, fstop, 
				   specfile, filmdiag, film);
}

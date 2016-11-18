
#include "stdafx.h"
#include "vdb.h"
#include "cameras/realistic.h"
#include "shapes/sphere.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include <fstream>
#include <string>

using namespace std;

//int vdb_c = 0;
//bool draw = false;

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
				 float hither, float yon, 
				 float sopen, float sclose, 
				 float filmdistance, float aperture_diameter, string specfile, 
				 float filmdiag, Film *f)
	: Camera(cam2world, sopen, sclose, f) // pbrt-v2 doesnot specify hither and yon
{
	string line;
	ifstream fin(specfile.c_str());
	while (getline(fin, line)) {
		if (line[0] == '#') continue;

		Lens l;
		sscanf(line.c_str(), "%f %f %f %f", &(l.radius), &(l.axpos), &(l.n), &(l.aperture));

		if (l.radius == 0) {
			l.n = 1;
			l.aperture = aperture_diameter;
		}

		if (lens.size()) {
			l.abs_axpos =  + lens.back().abs_axpos-lens.back().axpos;
		} else {
			l.abs_axpos = 0;
		}

		lens.push_back(l);
	}

	fin.close();

	assert(lens.size()); // no lens been defined !

	this->yon = yon;
	this->hither = hither;
	this->film = f;
	this->film_diagonal = filmdiag;
	this->film_distance = filmdistance;
	this->film_position = lens.back().abs_axpos - filmdistance;

	this->raster_diagonal = sqrt(f->xResolution*f->xResolution +
								 f->yResolution*f->yResolution);

	this->Raster2Camera = Translate(Vector((f->xResolution *0.5)*film_diagonal / raster_diagonal, -f->yResolution *0.5*film_diagonal / raster_diagonal, 0))
		*Scale(film_diagonal / raster_diagonal, film_diagonal / raster_diagonal, 1)
		*Scale(-1.f, 1.f, 1.f);
		
	this->Camera2Raster = Inverse(Raster2Camera);

	// vdb draw film
	//vdb_point(0,0,0);

	//vdb_color(1.f, 1.f, 1.f);
	//vdb_line(0, 15, film_position, 0, -15, film_position);
	//vdb_line(15, 0, film_position, -15, 0, film_position);
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
	// Sample point on lens, scale to aperture radius
	float lensU, lensV;
	ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
	lensU *= lens.back().aperture*0.5f;
	lensV *= lens.back().aperture*0.5f;

	// Generate raster and camera samples
	Point Praster(sample.imageX, sample.imageY, film_position);
	Point Pcamera;
	Raster2Camera(Praster, &Pcamera);

	Vector rayDirection = Normalize(Point(lensU, lensV, lens.back().abs_axpos) - Pcamera);
	ray->d = rayDirection;
	ray->o = Pcamera;
	ray->mint = 0;
	ray->maxt = INFINITY;
	//ray->time = 0;

	//draw = !((vdb_c++) % 10000);

	for (int i = lens.size() - 1; i >= 0; i--) {
		// Check intersection with this lens
		Point pHit;
		Vector normal;
		if (!LensIntersect(lens[i], *ray, &pHit, &normal)) return 0;

		// Keep normal pointed to the film side
		if (lens[i].radius*ray->d.z > 0) {
			normal = -normal;
		}		

		// Use SnellsLaw to get new Ray direcion
		Vector rayDirection2;
		float n1 = lens[i].n;
		float n2 = i == 0 ? 1 : lens[i - 1].n;
		if (!SnellsLaw(ray->d, &rayDirection2, normal, n1, n2)) return 0;

		//vdb_color(1, 1, 0);
		//vdb_point(pHit.x, pHit.y, pHit.z);

		ray->o = pHit;
		ray->d = Normalize(rayDirection2);
	}

	CameraToWorld(*ray, ray);

	// Return the weight of generated Ray
	// E = A*cos^4(theta)/Z^2
	float r = lens.back().aperture*0.5;
	float A = M_PI* r*r;
	float Z = film_distance;
	float costheta = Dot(rayDirection, Vector(0, 0, 1));
	float E = A*pow(costheta, 4) / (Z*Z);

    return E;
}

bool RealisticCamera::LensIntersect(const Lens l, const Ray & r, Point * pHit, Vector * normal) const {
	if (l.radius == 0) {
		// Use absolute value because ray can come from either side
		float scale = fabs((l.abs_axpos - r.o.z) / r.d.z);
		*pHit = r.o + scale*r.d;
	} else {
		float thit = 0;
		Transform w2o = Translate(Vector(0.f, 0.f, l.radius - l.abs_axpos));

		if (!SphereIntersect(r, &thit, &w2o, l.radius)) return false;
		if (thit > r.maxt || thit < r.mint) return false;

		*pHit = r(thit);
	}

	// Check aperture x^2 + y^2 < (aperture/2)^2
	if ((pHit->x * pHit->x + pHit->y * pHit->y) >= (l.aperture*l.aperture) / 4) {
		return false;
	}

	// Returned normal always points from the sphere surface outward
	*normal = Normalize(*pHit - Point(0, 0, l.abs_axpos - l.radius));

	//if (draw) {
	//	vdb_color(1, 0, 0);
	//	vdb_line(r.o.x, r.o.y, r.o.z, pHit->x, pHit->y, pHit->z);
	//}

	return true;
}

// https://www.wikiwand.com/en/Snell's_law#/Vector_form
bool RealisticCamera::SnellsLaw(const Vector l, Vector * refract, const Vector n, const float N1, const float N2) const {
	// Refract = r x l + (r x c - sqrt(1 - r^2 x (1 - c^2)))n
	float r = N1 / N2;
	float c = Dot(-n, l);
	int sign = c < 0 ? -1 : 1; // cos(theta) must > 0
	float tmp = 1 - r*r*(1 - c*c);

	// Total reflection
	if (tmp < 0) {
		return false;
	}

	*refract = r*l + (r*c - sign*sqrt(tmp))*n;

	return true;
}

bool RealisticCamera::SphereIntersect(const Ray &r, float *tHit, Transform *w2o, float radius) const {
	float phi;
	Point phit;
	// Transform _Ray_ to object space
	Ray ray;
	(*w2o)(r, &ray);

	// Compute quadratic sphere coefficients
	float A = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z;
	float B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y + ray.d.z*ray.o.z);
	float C = ray.o.x*ray.o.x + ray.o.y*ray.o.y + ray.o.z*ray.o.z - radius*radius;

	// Solve quadratic equation for _t_ values
	float t0, t1;
	if (!Quadratic(A, B, C, &t0, &t1))
		return false;

	// Compute intersection distance along ray
	if (t0 > ray.maxt || t1 < ray.mint)
		return false;
	float thit = t0;
	if (t0 < ray.mint) {
		thit = t1;
		if (thit > ray.maxt) return false;
	}

	*tHit = thit;

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

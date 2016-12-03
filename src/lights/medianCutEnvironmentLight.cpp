// lights/medianCutEnvironmentLight.cpp*
#include "stdafx.h"
#include "lights/medianCutEnvironmentLight.h"
#include "sh.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"

// MedianCutEnvironmentLight Utility Classes
struct MedianCutEnvironmentLightCube {
    // MedianCutEnvironmentLightCube Public Methods
	MedianCutEnvironmentLightCube(const MedianCutEnvironmentLight *l, const Scene *s,
                     float t, bool cv, float pe)
        : light(l), scene(s), time(t), pEpsilon(pe), computeVis(cv) { }
    Spectrum operator()(int, int, const Point &p, const Vector &w) {
        Ray ray(p, w, pEpsilon, INFINITY, time);
        if (!computeVis || !scene->IntersectP(ray))
            return light->Le(RayDifferential(ray));
        return 0.f;
    }
    const MedianCutEnvironmentLight *light;
    const Scene *scene;
    float time, pEpsilon;
    bool computeVis;
};



// MedianCutEnvironmentLight Method Definitions
MedianCutEnvironmentLight::~MedianCutEnvironmentLight() {
    delete distribution;
    delete radianceMap;
}


MedianCutEnvironmentLight::MedianCutEnvironmentLight(const Transform &light2world,
        const Spectrum &L, int ns, const string &texmap)
    : Light(light2world, ns) {
	int width = 0, height = 0;
	RGBSpectrum *texels = NULL;
	// Read texel data from _texmap_ into _texels_
	if (texmap != "") {
		texels = ReadImage(texmap, &width, &height);
		if (texels)
			for (int i = 0; i < width * height; ++i)
				texels[i] *= L.ToRGBSpectrum();
	}
	if (!texels) {
		width = height = 1;
		texels = new RGBSpectrum[1];
		texels[0] = L.ToRGBSpectrum();
	}
	radianceMap = new MIPMap<RGBSpectrum>(width, height, texels);
	// Initialize sampling PDFs for infinite area light


	// Compute scalar-valued image _img_ from environment map
	float solidAngle = ((2.f * M_PI) / (width - 1)) * ((M_PI) / (height - 1));
	float filter = 1.f / max(width, height);
	float *img = new float[width*height];
	for (int v = 0; v < height; ++v) {
		float vp = (float)v / (float)height;
		float sinTheta = sinf(M_PI * float(v + .5f) / float(height));
		for (int u = 0; u < width; ++u) {
			float up = (float)u / (float)width;
			img[u + v*width] = radianceMap->Lookup(up, vp, filter).y();
			img[u + v*width] *= sinTheta;

			// Scale the light intensity accroding to the areas
			texels[u + v*width] = texels[u + v*width] * solidAngle * sinTheta;
		}
	}

	// Compute SAT
	ConstructSAT(width, height, img, solidAngle);


	// Subdivided into 64 equal energy regions
	vector<Area> Areas;
	Areas.push_back(Area(0, 0, width - 1, height - 1));
	MedianCut(Areas, width, height, 64);


	// Put virtual light into each Area's centroid
	this->pdf = 1.f / Areas.size();
	for (Area a: Areas) {
		RGBSpectrum spectrum = RGBSpectrum(0.f);
		float cv = 0.f;
		float cu = 0.f;
		float sumf = 0;

#define IDX(u, v) ((u)+(v)*width)
		for (int v = a.miny; v <= a.maxy; v++) {
			for (int u = a.minx; u <= a.maxx; u++) {
				spectrum += texels[IDX(u, v)];
				float f = img[IDX(u, v)];
				cv += v * f;
				cu += u * f;
				sumf += f;
			}
		}
#undef IDX

		this->lights.push_back(virtualLight(cu / sumf / height, cv / sumf / width, spectrum));
	}

	// Compute sampling distributions for rows and columns of image
	distribution = new Distribution2D(img, width, height);
	delete[] texels;
	delete[] img;
}


Spectrum MedianCutEnvironmentLight::Power(const Scene *scene) const {
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    return M_PI * worldRadius * worldRadius *
        Spectrum(radianceMap->Lookup(.5f, .5f, .5f), SPECTRUM_ILLUMINANT);
}


Spectrum MedianCutEnvironmentLight::Le(const RayDifferential &r) const {
    Vector wh = Normalize(WorldToLight(r.d));
    float s = SphericalPhi(wh) * INV_TWOPI;
    float t = SphericalTheta(wh) * INV_PI;
    return Spectrum(radianceMap->Lookup(s, t), SPECTRUM_ILLUMINANT);
}


void MedianCutEnvironmentLight::SHProject(const Point &p, float pEpsilon,
        int lmax, const Scene *scene, bool computeLightVis,
        float time, RNG &rng, Spectrum *coeffs) const {
    // Project _InfiniteAreaLight_ to SH using Monte Carlo if visibility needed
    if (computeLightVis) {
        Light::SHProject(p, pEpsilon, lmax, scene, computeLightVis,
                         time, rng, coeffs);
        return;
    }
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = 0.f;
    int ntheta = radianceMap->Height(), nphi = radianceMap->Width();
    if (min(ntheta, nphi) > 50) {
        // Project _InfiniteAreaLight_ to SH from lat-long representation

        // Precompute $\theta$ and $\phi$ values for lat-long map projection
        float *buf = new float[2*ntheta + 2*nphi];
        float *bufp = buf;
        float *sintheta = bufp;  bufp += ntheta;
        float *costheta = bufp;  bufp += ntheta;
        float *sinphi = bufp;    bufp += nphi;
        float *cosphi = bufp;
        for (int theta = 0; theta < ntheta; ++theta) {
            sintheta[theta] = sinf((theta + .5f)/ntheta * M_PI);
            costheta[theta] = cosf((theta + .5f)/ntheta * M_PI);
        }
        for (int phi = 0; phi < nphi; ++phi) {
            sinphi[phi] = sinf((phi + .5f)/nphi * 2.f * M_PI);
            cosphi[phi] = cosf((phi + .5f)/nphi * 2.f * M_PI);
        }
        float *Ylm = ALLOCA(float, SHTerms(lmax));
        for (int theta = 0; theta < ntheta; ++theta) {
            for (int phi = 0; phi < nphi; ++phi) {
                // Add _InfiniteAreaLight_ texel's contribution to SH coefficients
                Vector w = Vector(sintheta[theta] * cosphi[phi],
                                  sintheta[theta] * sinphi[phi],
                                  costheta[theta]);
                w = Normalize(LightToWorld(w));
                Spectrum Le = Spectrum(radianceMap->Texel(0, phi, theta),
                                       SPECTRUM_ILLUMINANT);
                SHEvaluate(w, lmax, Ylm);
                for (int i = 0; i < SHTerms(lmax); ++i)
                    coeffs[i] += Le * Ylm[i] * sintheta[theta] *
                        (M_PI / ntheta) * (2.f * M_PI / nphi);
            }
        }

        // Free memory used for lat-long theta and phi values
        delete[] buf;
    }
    else {
        // Project _InfiniteAreaLight_ to SH from cube map sampling
        SHProjectCube(MedianCutEnvironmentLightCube(this, scene, time, computeLightVis,
                                       pEpsilon),
                      p, 200, lmax, coeffs);
    }
}

void MedianCutEnvironmentLight::ConstructSAT(const int w, const int h, const float* img, const float solidAngle) {
	float *sat = new float[w*h];

	for (int u = 0; u < w; u++) {
		for (int v = 0; v < h; v++) {
			float currentPixel = img[u + v*w] * solidAngle;
			if (!u && !v) {
				sat[0] = currentPixel;
				continue;
			}
			if (u && !v) {
				sat[u] = sat[u - 1] + currentPixel;
				continue;
			}
			if (!u && v) {
				sat[v*w] = sat[(v - 1)*w] + currentPixel;
				continue;
			}

			sat[u + v*w] =
				sat[(u - 1) + v*w] +
				sat[u + (v - 1)*w] -
				sat[(u - 1) + (v - 1)*w] +
				currentPixel;
		}
	}

	this->SAT = sat;
}

void MedianCutEnvironmentLight::MedianCut(vector<Area>& areas, const int w, const int h, const int partitions) {
	for (int it = 0; it < 6 && areas.size() < partitions; it++) {
		vector<Area> newArea;
		for (Area a : areas) {
			float halfE = a.getEnergy(SAT, w, h) * 0.5f;

			int ret = FindMedianCut(a, halfE, a.getLongestAxis(), w, h);
			if (a.getLongestAxis() == 0) {
				newArea.push_back(Area(a.minx, a.miny, ret, a.maxy));
				if (ret + 1 <= a.maxx) {
					newArea.push_back(Area(ret + 1, a.miny, a.maxx, a.maxy));
				}
			} else {
				newArea.push_back(Area(a.minx, a.miny, a.maxx, ret));
				if (ret + 1 <= a.maxy) {
					newArea.push_back(Area(a.minx, ret + 1, a.maxx, a.maxy));
				}
			}
		}
		areas = newArea;
	}
}

int MedianCutEnvironmentLight::FindMedianCut(const Area a, const float halfEnergy, const int axis, const int w, const int h) {
	int l = axis ? a.miny : a.minx;
	int r = axis ? a.maxy : a.maxx;
	int middle;
	int ret = l;

	while (l <= r) {
		middle = (l + r) / 2;
		float energy = axis
			? Area(a.minx, a.miny, a.maxx, middle).getEnergy(SAT, w, h) 
			: Area(a.minx, a.miny, middle, a.maxy).getEnergy(SAT, w, h);

		if (energy < halfEnergy) {
			l = middle + 1;
			ret = middle;
		} else {
			r = middle - 1;
		}
	}

	return ret;
}

MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
        const ParamSet &paramSet) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    string texmap = paramSet.FindOneFilename("mapname", "");
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
    return new MedianCutEnvironmentLight(light2world, L * sc, nSamples, texmap);
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
	PBRT_INFINITE_LIGHT_STARTED_SAMPLE();

	// random select virtual light
	virtualLight l = lights[Floor2Int(ls.uComponent*lights.size())];

	// Convert infinite light sample point to direction
	float theta = l.v * M_PI, phi = l.u * 2.f * M_PI;
	float costheta = cosf(theta), sintheta = sinf(theta);
	float sinphi = sinf(phi), cosphi = cosf(phi);
	*wi = LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
		costheta));

	// Compute PDF for sampled infinite light direction
	*pdf = this->pdf;
	if (sintheta == 0.f) *pdf = 0.f;

	// Return radiance value for infinite light direction
	visibility->SetRay(p, pEpsilon, *wi, time);
	Spectrum Ls = Spectrum(l.spectrum, SPECTRUM_ILLUMINANT);
	PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
	return Ls;
}


float MedianCutEnvironmentLight::Pdf(const Point &, const Vector &w) const {
    PBRT_INFINITE_LIGHT_STARTED_PDF();
    Vector wi = WorldToLight(w);
    float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
    float sintheta = sinf(theta);
    if (sintheta == 0.f) return 0.f;
    float p = distribution->Pdf(phi * INV_TWOPI, theta * INV_PI) /
           (2.f * M_PI * M_PI * sintheta);
    PBRT_INFINITE_LIGHT_FINISHED_PDF();
    return p;
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, float time,
        Ray *ray, Normal *Ns, float *pdf) const {

	// Shouldn't go here
	assert(false);

	return Spectrum(0);
}

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_LIGHTS_MEDIANCUTENV_H
#define PBRT_LIGHTS_MEDIANCUTENV_H

// lights/medianCutEnvironmentLight.h*
#include "pbrt.h"
#include "light.h"
#include "texture.h"
#include "shape.h"
#include "scene.h"
#include "mipmap.h"

struct virtualLight {
	float u;
	float v;
	Spectrum spectrum;
	virtualLight() {};
	virtualLight(float u, float v, Spectrum s) {
		this->u = u;
		this->v = v;
		this->spectrum = s;
	}
};

struct Area {
	int minx;
	int miny;
	int maxx;
	int maxy;

	Area() {};
	Area(int x0, int y0, int x1, int y1) {
		minx = x0;
		maxx = x1;
		miny = y0;
		maxy = y1;
	}

	float getEnergy(float SAT[], int width, int height) {
#define IDX(u, v) ((u)+(v)*width)
		float e = SAT[IDX(maxx, maxy)];
		if (minx > 0) e -= SAT[IDX(minx - 1, maxy)];
		if (miny > 0) e -= SAT[IDX(maxx, miny - 1)];
		if (miny > 0 && minx > 0) e += SAT[IDX(minx - 1, miny - 1)];
#undef IDX
		return e;
	}

	int getLongestAxis() {
		return maxx - minx > maxy - miny ? 0 : 1;
	}
};

// MedianCutEnvironmentLight Declarations
class MedianCutEnvironmentLight : public Light {
public:
    // MedianCutEnvironmentLight Public Methods
	MedianCutEnvironmentLight(const Transform &light2world, const Spectrum &power, int ns,
        const string &texmap);
    ~MedianCutEnvironmentLight();
    Spectrum Power(const Scene *) const;
    bool IsDeltaLight() const { return true; }
    Spectrum Le(const RayDifferential &r) const;
    Spectrum Sample_L(const Point &p, float pEpsilon, const LightSample &ls,
        float time, Vector *wi, float *pdf, VisibilityTester *visibility) const;
    Spectrum Sample_L(const Scene *scene, const LightSample &ls, float u1, float u2,
        float time, Ray *ray, Normal *Ns, float *pdf) const;
    float Pdf(const Point &, const Vector &) const;
    void SHProject(const Point &p, float pEpsilon, int lmax, const Scene *scene,
        bool computeLightVis, float time, RNG &rng, Spectrum *coeffs) const;
private:
	void ConstructSAT(const int w, const int h, const float* img, const float solidAngle);
	void MedianCut(vector<Area>& areas, const int w, const int h, const int partitions);
	int FindMedianCut(const Area a, const float halfEnergy, const int axis, const int w, const int h);

    // MedianCutEnvironmentLight Private Data
    MIPMap<RGBSpectrum> *radianceMap;
    Distribution2D *distribution;
	vector<virtualLight> lights;
	float pdf;
	float* SAT;
};


MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
        const ParamSet &paramSet);

#endif // PBRT_LIGHTS_MEDIANCUTENV_H

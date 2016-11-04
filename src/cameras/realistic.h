
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "camera.h"
#include "paramset.h"
#include "film.h"

struct Lens {
	float radius;
	float thickness;
	float axpos;
	float n;
	float aperture;
	//bool Intersect(const Ray &r, Point *pHit, Vector *normal) const;
};

// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
	// RealisticCamera Public Methods
	RealisticCamera(const AnimatedTransform &cam2world,
						float hither, float yon, float sopen,
						float sclose, float filmdistance, float aperture_diameter, string specfile,
						float filmdiag, Film *film);
	float GenerateRay(const CameraSample &sample, Ray *) const;
  
private:
	// RealisticCamera Private Data
	vector<Lens> lens;
	float raster_diagonal;
	float film_diagonal;
	float film_distance;
	float film_position;
	float hither;
	float yon;

	Transform Raster2Camera, Camera2Raster;
};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);


#endif	// PBRT_CAMERAS_REALISTIC_H
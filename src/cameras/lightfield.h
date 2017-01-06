#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_LIGHTFIELD_H
#define PBRT_CAMERAS_LIGHTFIELD_H

#include "pbrt.h"
#include "camera.h"
#include "film.h"

struct SurfaceSpec {
    float radius;    // Radius of lens sphere
    float thickness; // Thickness of lens
    float n;         // Index of refraction corresponding the current lens element.
    float aperture;  // Diameter of the lens element for the current interface definition.
    float z_pos;     // The z position of the image-side of the element's center
};

struct LensletStruct {
    int lensesX;        // The number of lenses in x direction.
    int lensesY;        // The number of lenses in y direction.
    float lensletWidth; // in mm
    float focalLength;  // in mm
    float fNumber;      // Should be focalLength/lensletWidth
    
    float filmdiag;
    int xRes;
    int yRes;
    float lensletDistance;
    float apertureDiameter;
};


class LFCamera : public Camera {
public:
	LFCamera(const AnimatedTransform &cam2world,
		float hither, float yon, float sopen,
		float sclose, float lensletDistance, float aperture_diameter,
		const string &specfile,
		const string &lensletSpec,
		float filmdiag,
		Film *film);
	~LFCamera();
	float GenerateRay(const CameraSample &sample, Ray *) const;
	bool SurfaceIntercept(SurfaceSpec ss, Ray r, DifferentialGeometry &dg) const;
	LensletStruct lensletArray;
	float stored_back_lens_a;

private:
	float ShutterOpen;
	float ShutterClose;
	Film * film;

	vector<SurfaceSpec> lensSpec;
	float lensletdistance;
	float aperture_diameter;
	float filmdiag;
	float filmy;
	float filmx;
	int failures;
	int successes;
};

LFCamera *CreateLFCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);

#endif

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_LIGHTFIELD_H
#define PBRT_CAMERAS_LIGHTFIELD_H

#include "pbrt.h"
#include "camera.h"
#include "film.h"

struct SurfaceSpec {
    float radius;
    float thickness;
    float n;
    float aperture;
    float z_pos; // The z position of the image-side of the element's center
};

struct LensletStruct {
    int lensesX; //The number of lenses in x direction.  e.g. 296
    int lensesY;
    float lensletWidth; //in mm
    float focalLength; // in mm
    float fNumber;  // Should simply be focalLength/lensletWidth
    
    //extra stuff
    float filmdiag;
    int xRes;
    int yRes;
    float lensletDistance;
    float apertureDiameter;
    
};

// Example representation of an autofocus zone.
class AfZone {
	public:
	  // from the data file
	  float left, right;
	  float top, bottom;
	  int xres,yres;
};
class LFCamera : public Camera {
public:
   LFCamera(const AnimatedTransform &cam2world,
      float hither, float yon, float sopen,
      float sclose, float lensletDistance, float aperture_diameter,
      const string &specfile,
      const string &lensletSpec,
	  const string &autofocusfile,
      float filmdiag,
	  Film *film);
   ~LFCamera();
   float GenerateRay(const CameraSample &sample, Ray *) const;
    bool surfaceIntercept(SurfaceSpec ss, Ray r, DifferentialGeometry &dg) const;
   void  AutoFocus(Renderer * renderer, const Scene * scene, Sample * origSample);
   void  ParseAfZones(const string& filename);
    float TestAutoFocusDepth(Renderer * renderer, const Scene * scene, Sample * origSample, char* output);
    LensletStruct lensletArray;
    float stored_back_lens_a;
    float calculatePrincipalPlane();

private:
   bool  autofocus;
   vector<AfZone> afZones;
   float ShutterOpen;
   float ShutterClose;
   Film * film;
    
    vector<SurfaceSpec> lensSpec;
    float mylensletdistance;
    float myaperture_diameter;
    float myfilmdiag;
    float myfilmy;
    float myfilmx;
    int failures;
    int successes;
    
    
};

LFCamera *CreateLFCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);

#endif

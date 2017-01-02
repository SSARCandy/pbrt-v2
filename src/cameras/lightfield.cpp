// cameras/lightfield.cpp*
// Extension of RealisticCamera from Project 3
#include "stdafx.h"
#include "cameras/lightfield.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include "filters/box.h"
#include "film/image.h"
#include "shapes/sphere.h"
#include "samplers/stratified.h"
#include "intersection.h"
#include "renderer.h"

#include <stdio.h>
#include <stdlib.h>
#include <tiffio.h>
#include <assert.h>
#include <ImfInputFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <half.h>
#include <algorithm>



#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>

using namespace std;
using namespace Imf;
using namespace Imath;



LFCamera *CreateLFCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	   // Extract common camera parameters from \use{ParamSet}
	   float hither = params.FindOneFloat("hither", -1);
	   float yon = params.FindOneFloat("yon", -1);
	   float shutteropen = params.FindOneFloat("shutteropen", -1);
	   float shutterclose = params.FindOneFloat("shutterclose", -1);

	   // LF camera-specific parameters
	   string specfile = params.FindOneString("specfile", "");
       string lensletSpec = params.FindOneString("lensletSpec", "");
	   float lensletDistance = params.FindOneFloat("lensletDistance", 70.0); // about 70 mm default to film
	   float fstop = params.FindOneFloat("aperture_diameter", 1.0);
	   float filmdiag = params.FindOneFloat("filmdiag", 35.0);
	   string autofocusfile = params.FindOneString("af_zones", "");
	   assert(hither != -1 && yon != -1 && shutteropen != -1 &&
	      shutterclose != -1 && lensletDistance!= -1);
	   if (specfile == "") {
	       Severe( "No lens spec file supplied!\n" );
	   }
    
	   return new LFCamera(cam2world, hither, yon,
	      shutteropen, shutterclose, lensletDistance, fstop,
	      specfile, lensletSpec, autofocusfile, filmdiag, film);
}

static void printLensletArray(LensletStruct lensletArray) {
    fprintf(stderr,"Lenslet Array Structure:\n");
    fprintf(stderr,"\tlensesX: %d\n",lensletArray.lensesX);
    fprintf(stderr,"\tlensesY: %d\n",lensletArray.lensesY);
    fprintf(stderr,"\tlensletWidth: %f\n",lensletArray.lensletWidth);
    fprintf(stderr,"\tfocalLength: %f\n",lensletArray.focalLength);
    fprintf(stderr,"\tfNumber: %f\n",lensletArray.fNumber);
    fprintf(stderr,"\tfilmdiag: %f\n",lensletArray.filmdiag);
    fprintf(stderr,"\txRes: %d\n",lensletArray.xRes);
    fprintf(stderr,"\tyRes: %d\n",lensletArray.yRes);
    fprintf(stderr,"\tlensletDistance: %f\n",lensletArray.lensletDistance);
    fprintf(stderr,"\tapertureDiameter: %f\n",lensletArray.apertureDiameter);
    fprintf(stderr,"\n");
}

static void printLensSpec(vector<SurfaceSpec> lensSpec) {
    for (unsigned int i = 0; i < lensSpec.size(); i++) {
        SurfaceSpec ss = lensSpec[i];
        fprintf(stderr,
                "Surface %d - radius:%f, thickness:%f,n:%f,aperture:%f,z_pos:%f\n",
                i+1,ss.radius,ss.thickness,ss.n,ss.aperture,ss.z_pos);
    }
    fprintf(stderr,"\n");
}

#define RADIUS_INDEX 0
#define THICKNESS_INDEX 1
#define N_INDEX 2
#define APERTURE_INDEX 3
#define Z_INDEX 4

// Tests for intercept with the given surface.  Updates a DifferentialGeometry
// structure if there is an intercept.
bool LFCamera::surfaceIntercept(SurfaceSpec ss, Ray r, DifferentialGeometry &dg) const {
    float z_pos = ss.z_pos;
    float signed_radius = ss.radius;
    if (signed_radius != 0) {
        // is a curved surface
        Transform o2w = Translate(Vector(0.f,0.f,z_pos - signed_radius));
        Transform w2o = Translate(Vector(0.f,0.f,signed_radius - z_pos));
        float radius = fabs(signed_radius);
        Sphere sphere(&o2w,&w2o,false,
                               radius,-radius,radius,360.0f);
        float tHit;
        float rayEpsilon;
        bool intersect = sphere.Intersect(r,&tHit,&rayEpsilon,&dg);
        if (intersect) {
            // Now need to make sure it is within aperture
            Point p = dg.p;
            float aperture = ss.aperture;
            float neededAperture = sqrt(p.x*p.x + p.y*p.y);
            if (neededAperture <= aperture/2.0) {
                // This is a valid surface intercept
                return true;
            } else {
                // Not only is this intercept not valid, but I don't think
                // there should be any valid intercepts.
                //fprintf(stderr,"Needed aperture was %f, but was %f\n",neededAperture,aperture/2.0);
                return false;
            }
        }
        return false;
    } else {
        // planar
        // CITATION: May have looked up implementation for this?
        float t = float(z_pos - r.o.z) / (r.d.z);
        if (t >= 0) {
            // possible intersection, make sure is within aperture
            Point p = r.o + t * (r.d);
            if (sqrt(p.x * p.x + p.y * p.y) <= myaperture_diameter/2.0) {
                //valid intercept
                dg.p = p;
                Normal normal(0.0f,0.0f,-1.0f);
                dg.nn = normal;
                return true;
            } else {
                //not valid
                return false;
            }
        }
        return false;
    }
}

// Goes through the specfile and builds up the lensSpec
static vector<SurfaceSpec> buildLensSpec(const string &specfile) {
    vector<SurfaceSpec> lensSpec;
    const char *specfilename = specfile.c_str();
    ifstream myspecfile (specfilename);
    string line;
    float z_pos = 0.0;
    if (myspecfile.is_open()) {
        while (myspecfile.good()) {
            getline (myspecfile,line);
            cout << line << endl;
            // Order of each line is:
            // Radius / Thickness / N / Aperture
            char *linechars = new char [line.length()+1];
            strcpy(linechars,line.c_str());
            char *hash_char = strchr(linechars,'#');
            int index = 0;
            if (hash_char == NULL && strcmp(linechars,"") != 0) {
                SurfaceSpec ss;
                char *token = strtok(linechars,"\t");
                while (token != NULL) {
                    float value = atof(token);
                    if (index == RADIUS_INDEX) {
                        ss.radius = value;
                    } else if (index == THICKNESS_INDEX) {
                        ss.thickness = value;
                        ss.z_pos = z_pos;
                        z_pos -= value;
                    } else if (index == N_INDEX) {
                        ss.n = value;
                        if (value == 0.0) {
                            ss.n = 1.0;
                        }
                    } else if (index == APERTURE_INDEX) {
                        ss.aperture = value;
                    }
                    index++;
                    token = strtok(NULL,"\t");
                }
                lensSpec.push_back(ss);
            }
        }
        myspecfile.close();
    }
    
    return lensSpec;
}

LensletStruct buildLensletArray(const string &lensletSpec) {
	fprintf(stderr, "Building lensletArray with spec %s\n", lensletSpec.c_str());
	LensletStruct lensletArray;


	string line;
	ifstream fin(lensletSpec.c_str());
	if (fin.fail()) perror("open failed ");

	while (getline(fin, line)) {
		if (line[0] == '#') continue;

		char key[30] = { '\0' };
		int value;
		sscanf(line.c_str(), "%s %d \n", key, &value);

		if (strcmp(key, "lensesX") == 0) {
			lensletArray.lensesX = int(value);
		} else if (strcmp(key, "lensesY") == 0) {
			lensletArray.lensesY = int(value);
		} else if (strcmp(key, "fNumber") == 0) {
			lensletArray.fNumber = value;
		} else if (strcmp(key, "filmdiag") == 0) {
			lensletArray.filmdiag = value;
		} else if (strcmp(key, "xRes") == 0) {
			lensletArray.xRes = int(value);
		} else if (strcmp(key, "yRes") == 0) {
			lensletArray.yRes = int(value);
		} else if (strcmp(key, "lensletDistance") == 0) {
			lensletArray.lensletDistance = value;
		} else if (strcmp(key, "apertureDiameter") == 0) {
			lensletArray.apertureDiameter = value;
		}
	}

	fin.close();


	assert(lensletArray.lensesX == lensletArray.lensesY);


	float filmx = sqrt(pow(lensletArray.filmdiag, 2) / (1 + float(lensletArray.yRes) / lensletArray.xRes));
	lensletArray.lensletWidth = filmx / lensletArray.lensesX;
	lensletArray.focalLength = lensletArray.lensletWidth * lensletArray.fNumber;

	fprintf(stderr, "Calculated focalLength as %f\n", lensletArray.focalLength);

	return lensletArray;
}


LFCamera::LFCamera(const AnimatedTransform &cam2world,
                                 float hither, float yon,
                                 float sopen, float sclose,
                                 float lensletDistance, float aperture_diameter_,
                                 const string &specfile,
                                 const string &lensletSpec,
								 const string &autofocusfile,
                                 float filmdiag,
								 Film *f)
                                 : Camera(cam2world, sopen, sclose, f),
								   ShutterOpen(sopen),
								   ShutterClose(sclose),
								   film(f)
{

	// YOUR CODE HERE -- build and store datastructures representing the given lens
	// and film placement.
    lensSpec = buildLensSpec(specfile);
    printLensSpec(lensSpec);
    
    
    
    lensletArray = buildLensletArray(lensletSpec);
    printLensletArray(lensletArray);
    
    myfilmx = sqrt(pow(filmdiag,2)/(1+float(film->yResolution)/film->xResolution));
    myfilmy = myfilmx * float(film->yResolution)/film->xResolution;
    
    
    
    float back_lens_a = lensSpec[lensSpec.size()-1].aperture;
    stored_back_lens_a = back_lens_a;
    
    mylensletdistance = lensletDistance;
    myaperture_diameter = aperture_diameter_;
    myfilmdiag = filmdiag;
    // calculate myfilmx myfilmy
    // myfilmx^2 + (yres/xres)*myfilmx^2 = myfilmdiag^2
    // (1 + yres/xres)*myfilmx^2 = myfilmdiag^2
    // myfilmx^2 = (myfilmdiag^2)/(1+yres/xres)
    
    float z = calculatePrincipalPlane();
    fprintf(stderr,"PRINCIPAL PLANE AT Z=%f\n",z);
    
    successes = 0;
    failures = 0;
    
    
    
	// If 'autofocusfile' is the empty string, then you should do
	// nothing in any subsequent call to AutoFocus()
	autofocus = false;

	if (autofocusfile.compare("") != 0)  {
		ParseAfZones(autofocusfile);
		autofocus = true;
	}
}


// parses the AF zone file
void LFCamera::ParseAfZones(const string& filename)
{
  ifstream specfile(filename.c_str());
   if (!specfile) {
      fprintf(stderr, "Cannot open file %s\n", filename.c_str());
      exit (-1);
   }

   char line[512];

   while (!specfile.eof()) {
      specfile.getline(line, 512);
      if (line[0] != '\0' && line[0] != '#' &&
         line[0] != ' ' && line[0] != '\t' && line[0] != '\n')
      {
		afZones.resize(afZones.size()+1);
		AfZone& zone = afZones[afZones.size()-1];
		sscanf(line, "%f %f %f %f\n", &zone.left, &zone.right, &zone.top, &zone.bottom);
      }
   }

	printf("Read in %zu AF zones from %s\n", afZones.size(), filename.c_str());
}

LFCamera::~LFCamera()
{

}

// Given a ray, the DifferentialGeometry of the surface intercept, and
// the two indices of refraction, returns the resulting ray.
// CITATION: This implementation comes directly from Section 5 of
// http://www-physics.ucsd.edu/~tmurphy/astr597/exercises/raytrace-3d.pdf
static Ray getNewRay(Ray r, DifferentialGeometry &dg, float n1, float n2) {
    Vector n(dg.nn);
    Vector curDir(r.d);
    if (Dot(n,curDir) < 0) {
        Vector newNormal(-1 * dg.nn);
        n = newNormal;
    }
    Vector ki = Normalize(curDir);
    Vector lhat = Normalize(n - Dot(ki,n)*ki);
    float first = Dot(n,ki);
    float second = pow(first,2);
    float third = 1 - second;
    //http://stackoverflow.com/questions/4453372/sqrt1-0-pow1-0-2-returns-nan
    float sintheta = sqrt(third < 0 ? 0 : third);
    
	float sintheta2 = min(max((n1 * sintheta) / n2, -1.0f), 1.0f);
    
    float theta = asin(sintheta);
    float theta2 = asin(sintheta2);
    float deltatheta = theta2 - theta;
    
    float a = (n1*n2) / fabs(n1*n2);
    Vector b = ki * cos(deltatheta);
    float c = (Dot(ki,n)/fabs(Dot(ki,n)));
    Vector d = lhat * sin(deltatheta);
    
    Vector k2 = a*(b - c * d);
    Vector k2n = Normalize(k2);
    
    Ray newRay(dg.p,k2n,0.0);
    return newRay;
}

// Figures out the starting ray and weight, and then takes this ray through
// the lens system.
float LFCamera::GenerateRay(const CameraSample &sample, Ray *ray) const
{
    // YOUR CODE HERE -- make that ray!

    // use sample->imageX and sample->imageY to get raster-space coordinates
    // of the sample point on the film.
    // use sample->lensU and sample->lensV to get a sample position on the lens
    float dx;
    float dy;
    //ConcentricSampleDisk(sample.lensU,sample.lensV,&dx,&dy);
    dx = 2 * sample.lensU - 1;
    dy = 2 * sample.lensV - 1;
    
    float back_lens_z = lensSpec[lensSpec.size()-1].z_pos;
    float lenslet_z = back_lens_z - mylensletdistance; //middle of lenslet
    float film_z = lenslet_z - lensletArray.focalLength;
    
    //calculate start position for sample
    float film_xfrac = float(sample.imageX) / film->xResolution;
    float film_yfrac = float(sample.imageY) / film->yResolution;
    float startx = -(-myfilmx/2.0 + film_xfrac*myfilmx);
    float starty = (-myfilmy/2.0 + film_yfrac*myfilmy);
    fflush(stderr);
    Point start(startx,starty,film_z);
    fflush(stderr);
    
    //calculate which lenslet this comes from
    int lensletX = lensletArray.lensesX - floor(film_xfrac*lensletArray.lensesX) - 1;
    int lensletY = floor(film_yfrac*lensletArray.lensesY);
    
    //calculate where it's headed on the lenslet
    float back_lens_a = lensletArray.lensletWidth;
    
    //TODO figure out if these should be backward or not
    float lensx = dx * back_lens_a/2.0;
    float lensy = dy * back_lens_a/2.0;
    //fprintf(stderr,"dx:%f,dy:%f\n",dx,dy);
    //fprintf(stderr,"lensesX:%d,lensesY:%d,lensletWidth:%f\n",lensletArray.lensesX,
        //    lensletArray.lensesY,lensletArray.lensletWidth);
    
    
    Point lensPoint(-(lensletArray.lensletWidth*lensletArray.lensesX/2.0)+lensletArray.lensletWidth/2.0+lensx+lensletX*lensletArray.lensletWidth,
                    -(lensletArray.lensletWidth*lensletArray.lensesY/2.0)+lensletArray.lensletWidth/2.0+lensy+lensletY*lensletArray.lensletWidth,lenslet_z);
    fflush(stderr);
    fflush(stderr);
    Vector direc = Normalize(lensPoint-start);
    Ray r(start,direc,0.0f);
    fflush(stderr);
    
    Vector normalized = Normalize(Vector(r.d));
    float costhetahere = Dot(normalized,Vector(0.f,0.f,1.0f));
    
    
    
    float weight = pow(costhetahere,4.0) / float(pow(lensletArray.focalLength,2)) * ((M_PI) * pow(back_lens_a/2.0,2));
    
    Point lensletCenter(-(lensletArray.lensletWidth*lensletArray.lensesX/2.0)+lensletArray.lensletWidth/2.0+lensletX*lensletArray.lensletWidth,
                        -(lensletArray.lensletWidth*lensletArray.lensesY/2.0)+lensletArray.lensletWidth/2.0+lensletY*lensletArray.lensletWidth,lenslet_z);
    Ray straightThrough(start,Normalize(lensletCenter-start),0.0f);
    //trace straightThrough
    float beyondZ = lenslet_z + (lenslet_z - start.z);
    float t = float(beyondZ - straightThrough.o.z) / (straightThrough.d.z);
    Point straightPoint;
    if (t >= 0) {
        // possible intersection, make sure is within aperture
        straightPoint = straightThrough.o + t * (straightThrough.d);
    }
    Vector newDir = Normalize(straightPoint-lensPoint);
    Ray newR(lensPoint,newDir,0.0f);
        
    r = newR;
    
    
    
    
    
    for (int s = lensSpec.size()-1; s >= 0; s--) {
        DifferentialGeometry dg;
        bool validIntercept = surfaceIntercept(lensSpec[s],r,dg);
        if (validIntercept) {
            //fprintf(stderr,"At least one!\n");
            fflush(stderr);
            float n1 = lensSpec[s].n;
            float n2;
            if (s > 0) {
                n2 = lensSpec[s-1].n;
            } else {
                n2 = 1.0;
            }
            r = getNewRay(r,dg,n1,n2);
        } else {
            *ray = r;
            //failures++;
            //exit(-1);
            //fprintf(stderr,"Failure\n");
            return 0.f;
        }
        
    }
    
    CameraToWorld(r,ray);
    ray->d = Normalize(ray->d);
    // GenerateRay() should return the weight of the generated ray
    fflush(stderr);
    
    return weight;
}

int getImageIndex(int x, int y) {
    return 3*(x + y*x);
}

float getIntensity(float* rgb,int start_index) {
    float r = rgb[start_index];
    float g = rgb[start_index+1];
    float b = rgb[start_index+2];
    return r;
}

// Uses a SML-like algorithm for focus calculation
float LFCamera::TestAutoFocusDepth(Renderer * renderer, const Scene * scene, Sample * origSample, char* output) {
    float bestFij = 0.0;
    for (size_t i=0; i<afZones.size(); i++) {
		AfZone & zone = afZones[i];
        
		RNG rng;
		MemoryArena arena;
		Filter * filter = new BoxFilter(.5f,.5f);
		const float crop[] = {zone.left,zone.right,zone.top,zone.bottom};
		ImageFilm sensor(film->xResolution, film->yResolution, filter, crop,output,false);
		int xstart,xend,ystart,yend;
		sensor.GetSampleExtent(&xstart,&xend,&ystart,&yend);
        
		StratifiedSampler sampler(xstart, xend, ystart, yend,
		                          16, 16, true, ShutterOpen, ShutterClose);
        
		// Allocate space for samples and intersections
		int maxSamples = sampler.MaximumSampleCount();
		Sample *samples = origSample->Duplicate(maxSamples);
		RayDifferential *rays = new RayDifferential[maxSamples];
		Spectrum *Ls = new Spectrum[maxSamples];
		Spectrum *Ts = new Spectrum[maxSamples];
		Intersection *isects = new Intersection[maxSamples];
        
		// Get samples from _Sampler_ and update image
		int sampleCount;
		while ((sampleCount = sampler.GetMoreSamples(samples, rng)) > 0) {
			// Generate camera rays and compute radiance along rays
			for (int i = 0; i < sampleCount; ++i) {
				// Find camera ray for _sample[i]_
                
				float rayWeight = this->GenerateRayDifferential(samples[i], &rays[i]);
				rays[i].ScaleDifferentials(1.f / sqrtf(sampler.samplesPerPixel));
                
                
				// Evaluate radiance along camera ray
                
				if (rayWeight > 0.f)
					Ls[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
													 arena, &isects[i], &Ts[i]);
				else {
					Ls[i] = 0.f;
					Ts[i] = 1.f;
				}
                
				// Issue warning if unexpected radiance value returned
				if (Ls[i].HasNaNs()) {
					Error("Not-a-number radiance value returned "
						  "for image sample.  Setting to black.");
					Ls[i] = Spectrum(0.f);
				}
				else if (Ls[i].y() < -1e-5) {
					Error("Negative luminance value, %f, returned"
						  "for image sample.  Setting to black.", Ls[i].y());
					Ls[i] = Spectrum(0.f);
				}
				else if (isinf(Ls[i].y())) {
					Error("Infinite luminance value returned"
						  "for image sample.  Setting to black.");
					Ls[i] = Spectrum(0.f);
				}
                
			}
            
			// Report sample results to _Sampler_, add contributions to image
			if (sampler.ReportResults(samples, rays, Ls, isects, sampleCount))
			{
				for (int i = 0; i < sampleCount; ++i)
				{
                    
					sensor.AddSample(samples[i], Ls[i]);
                    
				}
			}
            
			// Free _MemoryArena_ memory from computing image sample values
			arena.FreeAll();
		}
        
		float * rgb;
		int width;
		int height;
		sensor.WriteRGB(&rgb,&width,&height,1.f);
		// YOUR CODE HERE! The rbg contents of the image for this zone
		// are now stored in the array 'rgb'.  You can now do whatever
		// processing you wish
        fprintf(stderr,"AF Zone width:%d, height:%d\n",width,height);
        float Fij = 0.0;
        int step = 1;
        for (int x = 1; x < width-1; x++) {
            for (int y = 1; y < height-1; y++) {
                //first absolute value
                float f_one = 2 * getIntensity(rgb,getImageIndex(x,y));
                float f_two = getIntensity(rgb,getImageIndex(x-step,y));
                float f_three = getIntensity(rgb,getImageIndex(x+step,y));
                //second absolute value
                float s_one = 2 * getIntensity(rgb,getImageIndex(x,y));
                float s_two = getIntensity(rgb,getImageIndex(x,y-step));
                float s_three = getIntensity(rgb,getImageIndex(x,y+step));
                float MLxy = fabs(f_one-f_two-f_three)+fabs(s_one-s_two-s_three);
                if (MLxy > 0.0) {
                    Fij += MLxy;
                }
            }
        }
        
		//you own rgb  now so make sure to delete it:
		delete [] rgb;
		//if you want to see the output rendered from your sensor, uncomment this line (it will write a file called foo.exr)
		//sensor.WriteImage(1.f);
        
        
		delete[] samples;
		delete[] rays;
		delete[] Ls;
		delete[] Ts;
		delete[] isects;
        
        if (Fij > bestFij) {
            bestFij = Fij;
        }
    }
    return bestFij;
}






float LFCamera::calculatePrincipalPlane() {
    Ray orig(Point(-10,5,1),Vector(0,0,-1),0.0);
    Ray r(Point(-10,5,1),Vector(0,0,-1),0.0);
    fprintf(stderr,"My aperture diameter: %f\n",myaperture_diameter);
    // Figure out where what ray comes out of back of lens
    for (int s = 0; s < lensSpec.size(); s++) {
        DifferentialGeometry dg;
        bool validIntercept = surfaceIntercept(lensSpec[s],r,dg);
        if (validIntercept) {
            //fprintf(stderr,"At least one!\n");
            fflush(stderr);
            float n2 = lensSpec[s].n;
            float n1;
            if (s > 0) {
                n1 = lensSpec[s-1].n;
            } else {
                n1 = 1.0;
            }
            r = getNewRay(r,dg,n1,n2);
            //fprintf(stderr,"Did one surface\n");
        } else {
            fprintf(stderr,"Didn't work\n");
            exit(-1);
        }
        
    }
    fprintf(stderr,"New ray o(%f,%f,%f), d(%f,%f,%f)\n",r.o.x,r.o.y,r.o.z,r.d.x,r.d.y,r.d.z);
    //Figure out where new ray crosses x=0,y=0
    float t = (0 - r.o.x) / r.d.x;
    fprintf(stderr,"First t is %f\n",t);
    float z = r.o.z + t * r.d.z;
    fprintf(stderr,"\n\n\n\nCalculated z to be at %f\n",z);
    float y = r.o.y + t * r.d.y;
    fprintf(stderr,"Y should also be 0.  Is: %f\n",y);
    
    //Now figure out where this ray intersects backward with other ray
    // should intersect plane y=4
    // was number here
    float t2 = (5 - r.o.y) / r.d.y;
    fprintf(stderr,"T2 is %f\n",t2);
    z = r.o.z + t2 * r.d.z;
    float x = r.o.x + t2 * r.d.x;
    fprintf(stderr,"x should also be at 4, is: %f\n",x);
    
    return z;
    
    
}








// Sets film distance and tests for focus.  First traverses by 1mm intervals,
// then by 0.1mm intervals around the best scoring distance.
void  LFCamera::AutoFocus(Renderer * renderer, const Scene * scene, Sample * origSample) {
	// YOUR CODE HERE:
	// The current code shows how to create a new Sampler, and Film cropped to the size of the auto focus zone.
	// It then renders the film, producing rgb values.  You need to:
	//
	// 1. Modify this code so that it can adjust film plane of the camera
	// 2. Use the results of raytracing to evaluate whether the image is in focus
	// 3. Search over the space of film planes to find the best-focused plane.

	if(!autofocus) {
        fprintf(stderr,"NOT DOING AUTOFOCUS\n");
		return;
    }
    
    fprintf(stderr,"DOING AUTOFOCUS\n");
    float nearest = mylensletdistance/2.0;
    float farthest = nearest+60.0;
    
    float bestFij = 0.0;
    float bestDist;
    
    float step1 = 1.0;
    for (float cur_dist = nearest; cur_dist < farthest; cur_dist += step1) {
        fprintf(stderr,"mylensletdistance -> nearest:%f, farthest:%f, current:%f\n",
                nearest,farthest,cur_dist);
        mylensletdistance = cur_dist;
        char* output = new char[100];
        sprintf(output,"littleout%f.exr",cur_dist);

        float Fij = TestAutoFocusDepth(renderer,scene,origSample,output);
        Fij *= (pow(mylensletdistance,2.0)+pow(stored_back_lens_a/2.0,2.0)) / (pow(stored_back_lens_a/2.0,2.0)*M_PI);
        fprintf(stderr,"Fij = %f\n", Fij);
        if (Fij > bestFij) {
            bestFij = Fij;
            bestDist = cur_dist;
        }
    }
    
    float step2 = 0.1;
    for (float cur_dist = bestDist-step1+step2; cur_dist < bestDist+step1;
         cur_dist += step2) {
        fprintf(stderr,"mylensletdistance -> nearest:%f, farthest:%f, current:%f\n",
                nearest,farthest,cur_dist);
        mylensletdistance = cur_dist;
        char* output = new char[100];
        sprintf(output,"littleout%f.exr",cur_dist);

        float Fij = TestAutoFocusDepth(renderer,scene,origSample,output);
        Fij *= (pow(mylensletdistance,2.0)+pow(stored_back_lens_a/2.0,2.0)) / (pow(stored_back_lens_a/2.0,2.0)*M_PI);
        fprintf(stderr,"Fij = %f\n", Fij);
        if (Fij > bestFij) {
            bestFij = Fij;
            bestDist = cur_dist;
        }
    }
    
    mylensletdistance = bestDist;
    fprintf(stderr,"mylensletdistance set to %f\n",bestDist);
		
}

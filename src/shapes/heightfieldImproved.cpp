
// shapes/heightfieldImproved.cpp*
#include "stdafx.h"
#include "shapes/heightfieldImproved.h"
#include "shapes/trianglemesh.h"
#include "paramset.h"

// heightfieldImproved Method Definitions
HeightfieldImproved::HeightfieldImproved(const Transform *o2w, const Transform *w2o,
        bool ro, int x, int y, const float *zs)
    : Shape(o2w, w2o, ro) {
    nx = x;
    ny = y;
    z = new float[nx*ny];
    memcpy(z, zs, nx*ny*sizeof(float));
}


HeightfieldImproved::~HeightfieldImproved() {
    delete[] z;
}


BBox HeightfieldImproved::ObjectBound() const {
    float minz = z[0], maxz = z[0];
    for (int i = 1; i < nx*ny; ++i) {
        if (z[i] < minz) minz = z[i];
        if (z[i] > maxz) maxz = z[i];
    }
    return BBox(Point(0,0,minz), Point(1,1,maxz));
}


bool HeightfieldImproved::CanIntersect() const {
    return true;
}

HeightfieldImproved *CreateHeightfieldImprovedShape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params) {
    int nu = params.FindOneInt("nu", -1);
    int nv = params.FindOneInt("nv", -1);
    int nitems;
    const float *Pz = params.FindFloat("Pz", &nitems);
    Assert(nitems == nu*nv);
    Assert(nu != -1 && nv != -1 && Pz != NULL);
    return new HeightfieldImproved(o2w, w2o, reverseOrientation, nu, nv, Pz);
}

bool HeightfieldImproved::Intersect(const Ray &ray, float *tHit, float *rayEpsilon, DifferentialGeometry *dg) const {
	return false;
}

bool HeightfieldImproved::IntersectP(const Ray &ray) const {
	return false;
}
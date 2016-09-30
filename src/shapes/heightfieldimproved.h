
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_SHAPES_HEIGHTFIELDIMPROVED_H
#define PBRT_SHAPES_HEIGHTFIELDIMPROVED_H

// shapes/heightfieldImproved.h*
#include "shape.h"

// heightfieldImproved Declarations
class HeightfieldImproved : public Shape {
public:
    // heightfieldImproved Public Methods
    HeightfieldImproved(const Transform *o2, const Transform *w2o, bool ro, int nu, int nv, const float *zs);
    ~HeightfieldImproved();
    bool CanIntersect() const;
    void Refine(vector<Reference<Shape> > &refined) const;
	bool Intersect(const Ray & ray, float * tHit, float * rayEpsilon, DifferentialGeometry * dg) const;
	bool IntersectP(const Ray & ray) const;
    BBox ObjectBound() const;
private:
    // heightfieldImproved Private Data
    float *z;
    int nx, ny;
};


HeightfieldImproved *CreateHeightfieldImprovedShape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params);

#endif // PBRT_SHAPES_HEIGHTFIELDIMPROVED_H

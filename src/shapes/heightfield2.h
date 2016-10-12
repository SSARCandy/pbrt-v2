
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_SHAPES_HEIGHTFIELD2_H
#define PBRT_SHAPES_HEIGHTFIELD2_H

// shapes/heightfield2.h*
#include "shape.h"

// heightfield2 Declarations
class Heightfield2 : public Shape {
public:
    // heightfield2 Public Methods
    Heightfield2(const Transform *o2, const Transform *w2o, bool ro, int nu, int nv, const float *zs);
    ~Heightfield2();
    bool CanIntersect() const;
	bool Intersect(const Ray &r, float *tHit, float *rayEpsilon, DifferentialGeometry *dg) const;
	bool IntersectP(const Ray &r) const;
    BBox ObjectBound() const;
	void GetShadingGeometry(const Transform &obj2world, const DifferentialGeometry &dg, DifferentialGeometry *dgShading) const;

private:
	bool VoxelIntersector(const Ray &r, int x, int y, Intersection *in, float *tHit, BBox &bbox) const;
	void ComputeVertexNormal();
	inline bool OutOfBoundary(int x, int y) {
		return x < 0 || y < 0 || x >= nx || y >= ny;
	}
	inline int posToVoxel(const Point &P, int axis) const {
		if (axis == 2) return 0;
		int v = Float2Int(P[axis] * nVoxels[axis]);
		return Clamp(v, 0, nVoxels[axis]);
	}
	inline float voxelToPos(int p, int axis) const {
		return p / (float)(nVoxels[axis]);
	}
	inline float getZ(int x, int y) const {
		return z[y*nx + x];
	}

    // heightfield2 Private Data
    float *z;
    int nx, ny;
	int nVoxels[3]; // [nx-1, ny-1, 1]
	Normal *vertexNormals;
	float minz, maxz; // save it for BBox
};


Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params);

#endif // PBRT_SHAPES_HEIGHTFIELD2_H

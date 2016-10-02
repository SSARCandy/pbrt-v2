
// shapes/heightfieldImproved.cpp*

#include "stdafx.h"
#include "shapes/heightfieldImproved.h"
#include "shapes/trianglemesh.h"
#include "paramset.h"
//#include "core/primitive.h"

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
	return BBox(Point(0, 0, minz), Point(1, 1, maxz));
}

bool HeightfieldImproved::VoxelIntersector(const Ray &ray, int i, int j, Intersection *in) const {
	Point TL(voxelToPos(i, 0),     voxelToPos(j, 1),     getZ(i, j));
	Point TR(voxelToPos(i + 1, 0), voxelToPos(j, 1),     getZ(i + 1, j));
	Point BR(voxelToPos(i + 1, 0), voxelToPos(j + 1, 1), getZ(i + 1, j + 1));
	Point BL(voxelToPos(i, 0),     voxelToPos(j + 1, 1), getZ(i, j + 1));

	return false;
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

bool HeightfieldImproved::Intersect(const Ray &r, float *tHit, float *rayEpsilon, DifferentialGeometry *dg) const {

	return false;
	// Transform _Ray_ to object space
	Ray ray;
	(*WorldToObject)(r, &ray);

	// Check if ray intersect with BBox, and get intersection between ray & bbox
	float rayT;
	BBox bbox = ObjectBound();
	if (bbox.Inside(ray(ray.mint))) {
		rayT = ray.mint;
	} else if (!bbox.IntersectP(ray, &rayT)) {
		return false;
	}

	// move t (unit) from ray.o
	Point gridIntersect = ray(rayT);


	// Set up 3D DDA for ray
	float NextCrossingT[3], DeltaT[3];
	int Step[3], Out[3], Pos[3];
	for (int axis = 0; axis < 3; ++axis) {
		// Compute current voxel for axis
		Pos[axis] = posToVoxel(gridIntersect, axis);
		if (ray.d[axis] >= 0) {
			// Handle ray with positive direction for voxel stepping
			NextCrossingT[axis] = rayT + (voxelToPos(Pos[axis] + 1, axis) - gridIntersect[axis]) / ray.d[axis];
			DeltaT[axis] = 1.0 / ray.d[axis];
			Step[axis] = 1;
			Out[axis] = nVoxels[axis];
		} else {
			// Handle ray with negative direction for voxel stepping
			NextCrossingT[axis] = rayT + (voxelToPos(Pos[axis], axis) - gridIntersect[axis]) / ray.d[axis];
			DeltaT[axis] = -1.0 / ray.d[axis];
			Step[axis] = -1;
			Out[axis] = -1;
		}
	}

	// Walk ray through voxel grid
	bool hitSomething = false;
	Intersection intersection;
	for (;;) {
		int i = Pos[0], j = Pos[1];
		hitSomething = VoxelIntersector(ray, i, j, &intersection);

		// Advance to next voxel

		// Find _stepAxis_ for stepping to next voxel
		int bits = ((NextCrossingT[0] < NextCrossingT[1]) << 2) +
			       ((NextCrossingT[0] < NextCrossingT[2]) << 1) +
			       ((NextCrossingT[1] < NextCrossingT[2]));
		const int cmpToAxis[8] = { 2, 1, 2, 1, 2, 2, 0, 0 };
		int stepAxis = cmpToAxis[bits];
		if (ray.maxt < NextCrossingT[stepAxis])
			break;
		Pos[stepAxis] += Step[stepAxis];
		if (Pos[stepAxis] == Out[stepAxis])
			break;
		NextCrossingT[stepAxis] += DeltaT[stepAxis];
	}

	if (!hitSomething) return false;

	// fill-in the intersection properities


	return false;
}

bool HeightfieldImproved::IntersectP(const Ray &r) const {
	return false;
}
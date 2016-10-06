// shapes/heightfieldImproved.cpp*

#include "stdafx.h"
#include "vdb.h"
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

	nVoxels[0] = nx - 1;
	nVoxels[1] = ny - 1;
	nVoxels[2] = 1;

	ComputeVertexNormal();
	DebugerDrawHeightField(true);

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

bool HeightfieldImproved::VoxelIntersector(const Ray &r, int x, int y, Intersection *in, float *tHit) const {
	Point TL(voxelToPos(x, 0),     voxelToPos(y, 1),     getZ(x, y));
	Point TR(voxelToPos(x + 1, 0), voxelToPos(y, 1),     getZ(x + 1, y));
	Point BR(voxelToPos(x + 1, 0), voxelToPos(y + 1, 1), getZ(x + 1, y + 1));
	Point BL(voxelToPos(x, 0),     voxelToPos(y + 1, 1), getZ(x, y + 1));

	// ray should intersect with voxel
	BBox bbox(Union(BBox(TL, TR), BBox(BR, BL)));
	if (!bbox.IntersectP((*WorldToObject)(r))) return false;

	int vptr[6] = { 0,1,2,0,2,3 };
	Point pts[4] = { TL,TR,BR,BL };
	float uvs[8] = { TL.x, TL.y, TR.x, TR.y, BR.x, BR.y, BL.x, BL.y };

	DebugerDrawTriangle(TL, TR, BR, true);
	DebugerDrawTriangle(TL, BR, BL, true);

#define INDEX(x,y) ((x)+(y)*nx)
	Normal normals[4] = {
		vertexNormals[INDEX(x,y)],
		vertexNormals[INDEX(x+1,y)],
		vertexNormals[INDEX(x+1,y+1)],
		vertexNormals[INDEX(x,y+1)]
	};
#undef INDEX

	// build-up two triangles
	TriangleMesh *triMesh =
		new TriangleMesh(ObjectToWorld, WorldToObject, ReverseOrientation, 2, 4, vptr, pts, normals, NULL, uvs, NULL);
	Triangle *tri1 = new Triangle(ObjectToWorld, WorldToObject, ReverseOrientation, triMesh, 0);
	Triangle *tri2 = new Triangle(ObjectToWorld, WorldToObject, ReverseOrientation, triMesh, 1);

	// test intersection with each trangles, and get the intersection info.
	Intersection i1, i2;
	float tHit1, tHit2;
	bool haveIntersect1 = tri1->Intersect(r, &tHit1, &(i1.rayEpsilon), &(i1.dg));
	bool haveIntersect2 = tri2->Intersect(r, &tHit2, &(i2.rayEpsilon), &(i2.dg));

	// no intersection with both triangles
	if (!haveIntersect1 && !haveIntersect2) return false;

	// both triangles intersect with ray, return the nearest(smaller tHit)
	if (haveIntersect1 && haveIntersect2) {
		if (tHit1 < tHit2) {
			*in = i1;
			*tHit = tHit1;
		} else {
			*in = i2;
			*tHit = tHit2;
		}
		return true;
	}

	// only one of triangles have intersect with ray
	if (haveIntersect1) {
		*in = i1;
		*tHit = tHit1;
	} else {
		*in = i2;
		*tHit = tHit2;
	}

	return true;
}

void HeightfieldImproved::ComputeVertexNormal() {
	vertexNormals = new Normal[ny*nx];

	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x++) {
			// each vertex have 6 neighbors, {T, TL, L, B, BR, R}
			Point T, TL, L, B, BR, R;
			Point M(voxelToPos(x, 0), voxelToPos(y, 1), getZ(x, y)); // Middle point

			// if the point doesn't exist, then equal to M
			T  = OutOfBoundary(x,     y - 1) ? M : Point(voxelToPos(x,     0), voxelToPos(y - 1, 1), getZ(x,     y - 1));
			TL = OutOfBoundary(x - 1, y - 1) ? M : Point(voxelToPos(x - 1, 0), voxelToPos(y - 1, 1), getZ(x - 1, y - 1));
			L  = OutOfBoundary(x - 1,     y) ? M : Point(voxelToPos(x - 1, 0), voxelToPos(y,     1), getZ(x - 1,     y));
			B  = OutOfBoundary(x,     y + 1) ? M : Point(voxelToPos(x,     0), voxelToPos(y + 1, 1), getZ(x,     y + 1));
			BR = OutOfBoundary(x + 1, y + 1) ? M : Point(voxelToPos(x + 1, 0), voxelToPos(y + 1, 1), getZ(x + 1, y + 1));
			R  = OutOfBoundary(x + 1,     y) ? M : Point(voxelToPos(x + 1, 0), voxelToPos(y,     1), getZ(x + 1,     y));

#define X(p1,p2) Cross((p1-M),(p2-M))
			// the true normal = average of all 6 normals
			Vector aggerateNormal = X(TL, T) + X(T, L) + X(L, B) + X(B, BR) + X(BR, R) + X(R, T);
#undef X

			vertexNormals[y*ny + x] = Normal(Normalize(aggerateNormal));
			DebugerDrawLine(M, M + Normalize(aggerateNormal)/5, true);
		}
	}

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
	//DebugerDrawPoint(gridIntersect, 0, 0, 1, true);

	// Set up 3D DDA for ray
	float NextCrossingT[3], DeltaT[3];
	int Step[3], Out[3], Pos[3];
	for (int axis = 0; axis < 3; ++axis) {
		// Compute current voxel for axis
		Pos[axis] = posToVoxel(gridIntersect, axis);
		if (ray.d[axis] >= 0) {
			// Handle ray with positive direction for voxel stepping
			NextCrossingT[axis] = rayT + (voxelToPos(Pos[axis] + 1, axis) - gridIntersect[axis]) / ray.d[axis];
			DeltaT[axis] = 1.0 / (nVoxels[axis] * ray.d[axis]);
			Step[axis] = 1;
			Out[axis] = nVoxels[axis];
		} else {
			// Handle ray with negative direction for voxel stepping
			NextCrossingT[axis] = rayT + (voxelToPos(Pos[axis], axis) - gridIntersect[axis]) / ray.d[axis];
			DeltaT[axis] = -1.0 / (nVoxels[axis] * ray.d[axis]);
			Step[axis] = -1;
			Out[axis] = -1;
		}
	}

	// Walk ray through voxel grid
	Intersection intersection;
	float _tHit;
	bool hitSomething = false;
	for (;;) {
		int i = Pos[0], j = Pos[1];
		hitSomething = VoxelIntersector(r, i, j, &intersection, &_tHit);

		// no overlapping voxels in heightfield
		if (hitSomething) break;

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
	*tHit = _tHit;
	*rayEpsilon = intersection.rayEpsilon;
	*dg = intersection.dg;

	Point pp = ray(_tHit);
	//DebugerDrawPoint(pp, 1, 0, 0, true);
	//DebugerDrawPoint(intersection.dg.p, 0, 1, 0, false);
	
	return true;
}

bool HeightfieldImproved::IntersectP(const Ray &r) const {
	return false;
}

void HeightfieldImproved::GetShadingGeometry(const Transform &obj2world,
	const DifferentialGeometry &dg,	DifferentialGeometry *dgShading) const {
	dg.shape->GetShadingGeometry(obj2world, dg, dgShading);
}

/**
 * Debug using vdb
 */
void HeightfieldImproved::DebugerDrawHeightField(bool inWorldSpace) const {
	for (int j = 0; j < ny - 1; j++) {
		for (int i = 0; i < nx - 1; i++) {
			Point TL(voxelToPos(i, 0), voxelToPos(j, 1), getZ(i, j));
			Point TR(voxelToPos(i + 1, 0), voxelToPos(j, 1), getZ(i + 1, j));
			Point BR(voxelToPos(i + 1, 0), voxelToPos(j + 1, 1), getZ(i + 1, j + 1));
			Point BL(voxelToPos(i, 0), voxelToPos(j + 1, 1), getZ(i, j + 1));

			if (inWorldSpace) {
				TL = (*ObjectToWorld)(TL);
				TR = (*ObjectToWorld)(TR);
				BL = (*ObjectToWorld)(BL);
				BR = (*ObjectToWorld)(BR);
			}

			vdb_line(TL.x, TL.y, TL.z, TR.x, TR.y, TR.z);
			vdb_line(TR.x, TR.y, TR.z, BR.x, BR.y, BR.z);
			vdb_line(BR.x, BR.y, BR.z, BL.x, BL.y, BL.z);
			vdb_line(BL.x, BL.y, BL.z, TL.x, TL.y, TL.z);
		}
	}
}

void HeightfieldImproved::DebugerDrawPoint(Point pt, float r, float g, float b, bool inWorldSpace) const {
	if (inWorldSpace) {
		pt = (*ObjectToWorld)(pt);
	}
	vdb_color(r, g, b);
	vdb_point(pt.x, pt.y, pt.z);
	vdb_color(1, 1, 1);
}

void HeightfieldImproved::DebugerDrawTriangle(Point pt1, Point pt2, Point pt3, bool inWorldSpace) const {
	if (inWorldSpace) {
		pt1 = (*ObjectToWorld)(pt1);
		pt2 = (*ObjectToWorld)(pt2);
		pt3 = (*ObjectToWorld)(pt3);
	}

	vdb_line(pt1.x, pt1.y, pt1.z, pt2.x, pt2.y, pt2.z);
	vdb_line(pt2.x, pt2.y, pt2.z, pt3.x, pt3.y, pt3.z);
	vdb_line(pt3.x, pt3.y, pt3.z, pt1.x, pt1.y, pt1.z);
}

void HeightfieldImproved::DebugerDrawLine(Point pt1, Point pt2, bool inWorldSpace) const {
	if (inWorldSpace) {
		pt1 = (*ObjectToWorld)(pt1);
		pt2 = (*ObjectToWorld)(pt2);
	}
	vdb_color(0, 1, 0);
	vdb_line(pt1.x, pt1.y, pt1.z, pt2.x, pt2.y, pt2.z);
	vdb_color(1,1,1);
}
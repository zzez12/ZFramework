#pragma once

#include "GlobalDefs.h"
#include <vector>
#include "vectorN.h"

namespace ZMeshSpace
{
	enum ZPlaneDistanceType {
		PlaneDistanceRM2,
		PlaneDistanceNormD
	};
	enum ZPlaneType {
		PlaneType_Default,
		PlaneType_Polar
	};

	class ZPlane
	{
	private:
		Vec3f center_;
		Vec3f norm_;
		float d_;	// plane: nx+d=0
		Vec3f polarCoord_;	// (theta, phi, d)

		Vec3f corners[4];
		std::vector<Vec3f> initPoints_;
		std::vector<Vec3f> projectedPoints_;
		float averageEdgeLength_;
		ZPlaneType planeType_;
	public:
		ZPlane();
		ZPlane(const Vec3f& center, const Vec3f& norm);
		ZPlane(const Vec3f& v0, const Vec3f& v1, const Vec3f& v2);

		// return the rm2 error
		double buildPlane(std::vector<Vec3f>& points);
		double distance(const ZPlane& plane, ZPlaneDistanceType type = PlaneDistanceRM2) const;
		double distance(ZPlane* plane, ZPlaneDistanceType type=PlaneDistanceRM2);

		void buildPlaneByPolarCoord(const Vec3f& polarCoord, bool withCenter=false, const Vec3f& center=Vec3f(0,0,0));
		void buildPlaneByPolarCoord(const MATH::vector5f& polarCoordWithCenter);

		// setter and getters
		void setCenter(const Vec3f& c) {center_ = c;}
		void setNorm(const Vec3f& n) {norm_=n; norm_.unify();}
		Vec3f getCenter()const {return center_;}
		Vec3f getNorm()const {return norm_;}
		Vec3f getPolarCoord() const {return polarCoord_;}
		float getAverageEdgeLength() const {return averageEdgeLength_;}
		void setD(float d) {d_=d;}
		float getD() const {return d_;}
		ZPlaneType getPlaneType() const {return planeType_;}
		void setPlaneType(ZPlaneType type) {planeType_ = type;}

		void computeAverageEdgeLength();
		void computePolarCoordCorners();

		Vec3f projectToPlane(const Vec3f& input);
		const std::vector<Vec3f>& getProjectedPoints() const {return projectedPoints_;}
		const std::vector<Vec3f>& getInitPoints() const {return initPoints_;}

	private:
		double distance_rm2(ZPlane *plane);
		double distance_rm2(const ZPlane& plane) const;
		double distance_normD(ZPlane* plane);
		double distance_normD(const ZPlane& plane) const;

	private:
		void computePolar();
	public:
		float colorValue_;	// for rendering [0,1]
		Vec3f polarCoordCorners_[4];	// for rendering the plane (using polar coord.)
	};
}
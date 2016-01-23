#include "ZPlane.h"
#include "ZMathAlgorithms.h"
#include "ANNWarper.h"
#include "SphericalSystem.h"

namespace ZMeshSpace
{
	ZPlane::ZPlane()
	{
		setPlaneType(PlaneType_Default);
		averageEdgeLength_ = 0.01f;
		setCenter(Vec3f(0,0,0));
		setNorm(Vec3f(1,0,0));
	}

	ZPlane::ZPlane(const Vec3f& center, const Vec3f& norm)
	{
		averageEdgeLength_ = 0.01f;
		setCenter(center);
		setNorm(norm);
	}

	ZPlane::ZPlane(const Vec3f& v0, const Vec3f& v1, const Vec3f& v2)
	{
		averageEdgeLength_ = 0.01f;
		setCenter((v0+v1+v2)/3);
		setNorm((v1-v0)*(v2-v0));
	}

	void ZPlane::buildPlaneByPolarCoord(const Vec3f& polarCoord, bool withCenter, const Vec3f& center)
	{
		setPlaneType(PlaneType_Polar);
		polarCoord_ = polarCoord;
		setNorm(MATH::spherical2Cartesian(polarCoord));
		setD(polarCoord.z);
		if (withCenter)
			setCenter(center);
		else
			setCenter(getNorm()*(-getD()));
		computePolarCoordCorners();
	}

	void ZPlane::buildPlaneByPolarCoord(const MATH::vector5f& polarCoordWithCenter)
	{
		setPlaneType(PlaneType_Polar);
		setNorm(MATH::spherical2Cartesian(Vec3f(polarCoordWithCenter[0], polarCoordWithCenter[1], 0)));
		setCenter(Vec3f(polarCoordWithCenter[2], polarCoordWithCenter[3], polarCoordWithCenter[4]));
		setD(getNorm().DotProduct(getCenter())*(-1.f));
		computePolarCoordCorners();
	}

	void ZPlane::computePolar()
	{
		// n_ = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
		// polarCoord_ = (theta, phi, d_)

		if (!g_isSameValue(getNorm().length(), 1.f))
			return;
		if (g_isSameValue(norm_.z, 1))
		{
			polarCoord_ = Vec3f(0, 0, d_);
			return;
		}

		float theta, phi;
		theta = acos(norm_.z);
		float sinTheta = sqrt(1-norm_.z*norm_.z);//sin(theta);
// 		if (g_isSameValue(sinTheta, 0))
// 			std::cout << "zero divider!" << std::endl;
		float cosTheta = norm_.z;
		float sinPhi = norm_.y / sinTheta;
		float cosPhi = norm_.x / sinTheta;
		phi = acos(cosPhi);
		if (sinPhi<0)
			phi += Z_PI;
		polarCoord_ = Vec3f(theta, phi, d_);
// 		std::cout << "(" << norm_ << "," << d_ << ") --> ("
// 				  << polarCoord_ << ")" << std::endl;

		// prepare the polar-coordinate plane
		computePolarCoordCorners();

	}

	void ZPlane::computePolarCoordCorners()
	{
		float cosTheta = cos(polarCoord_.x);
		float sinTheta = sin(polarCoord_.x);
		float cosPhi = cos(polarCoord_.y);
		float sinPhi = sin(polarCoord_.y);
		Vec3f p0 = getCenter();//Vec3f(sinTheta*cosPhi, sinTheta*sinPhi, cosTheta) * (-d_);
		Vec3f m1(-cosPhi*cosTheta, -sinPhi*cosTheta, sinTheta);
		Vec3f m2 = norm_*m1;
		polarCoordCorners_[0] = p0 + m1 + m2;
		polarCoordCorners_[1] = p0 + m1 - m2;
		polarCoordCorners_[2] = p0 - m1 - m2;
		polarCoordCorners_[3] = p0 - m1 + m2;
	}

	double ZPlane::buildPlane(std::vector<Vec3f>& points)
	{
		double err=0;
		Vec3f ori(0,0,0);
		for (unsigned int i=0; i<points.size(); i++)
		{
			Vec3f& p = points[i];
			ori = ori + p;
		}
		ori = ori/points.size();
		std::vector<float> values;
		std::vector<Vec3f> axes;
		ZMathAlgorithms::calPCA(points, ori, values, axes);
		setCenter(ori);
		setNorm(axes[0]*axes[1]);
		float d = getNorm().DotProduct(getCenter())*(-1.f);
		setD(d);		
		
		// compute error
		for (unsigned int i=0; i<points.size(); i++)
		{
			Vec3f& p = points[i];
			err += (p-ori).DotProduct(getNorm());
		}
		err/=points.size();

		// compute projected points
		projectedPoints_.resize(points.size());
		for (int i=0; i<points.size(); i++)
		{
			projectedPoints_[i] = projectToPlane(points[i]);
		}
		computeAverageEdgeLength();

		//float d = getNorm().DotProduct(projectedPoints_[0])*(-1.f);
		// change to polar coordinate system
		//computePolar();
		polarCoord_ = MATH::cartesian2Spherical(norm_, d_);

		// test
// 		{
// 			Vec3f v[7];
// 			v[0] = getNorm();
// 			for (int i=0; i<6; i++)
// 			{
// 				if (i%2==0)
// 					v[i+1] = MATH::cartesian2Spherical(v[i], d_);
// 				else
// 					v[i+1] = MATH::spherical2Cartesian(v[i]);
// 			}
// 			Vec3f dif = v[6] - v[0];
// 			std::cout << ZConsoleTools::green << dif.length() << ZConsoleTools::white << "\n";
// 			if (dif.length()>10e-3)
// 				std::cout << v[0] << " --> " << v[6];
// 		}

		computePolarCoordCorners();
		// save init points
		initPoints_ = points;
		return err;
	}

	Vec3f ZPlane::projectToPlane(const Vec3f& input)
	{
		Vec3f p = input-getCenter();
		return input - getNorm() * p.DotProduct(getNorm());
	}

	void ZPlane::computeAverageEdgeLength()
	{
		float lengths = 0.f;
		for (int i=0; i<projectedPoints_.size()-1; i++)
		{
			lengths += (projectedPoints_[i]-projectedPoints_[i+1]).norm();
		}
		lengths /= (projectedPoints_.size()-1);
	}

	double ZPlane::distance_rm2(ZPlane *plane)
	{
		if (plane==NULL)
			return 0.0;

		AnnWarper annSearch;
		annSearch.init(initPoints_);
		double dis = 0;
		const std::vector<Vec3f>& points = plane->getInitPoints();
		for (int i=0; i<points.size(); i++)
		{
			std::vector<int> nnIdx;
			std::vector<float> nnDist;
			annSearch.queryPoint(points[0], 1, nnIdx, nnDist);
			dis += nnDist[0];
		}
		dis /= points.size();
		return dis;
	}

	double ZPlane::distance_rm2(const ZPlane& plane) const
	{
		if (initPoints_.empty()) return 0;
		AnnWarper annSearch;
		annSearch.init(initPoints_);
		double dis = 0;
		const std::vector<Vec3f>& points = plane.getInitPoints();
		for (int i=0; i<points.size(); i++)
		{
			std::vector<int> nnIdx;
			std::vector<float> nnDist;
			annSearch.queryPoint(points[0], 1, nnIdx, nnDist);
			dis += nnDist[0];
		}
		dis /= points.size();
		return dis;
	}

	double ZPlane::distance_normD(ZPlane* plane)
	{
		double dis=0;
		dis += acos(getNorm().DotProduct(plane->getNorm()));
		dis += abs(getD()-plane->getD());
		return dis;
	}

	double ZPlane::distance_normD(const ZPlane& plane) const
	{
		double dis=0;
		dis += acos(getNorm().DotProduct(plane.getNorm()));
		dis += abs(getD()-plane.getD());
		return dis;
	}

	double ZPlane::distance(ZPlane* plane, ZPlaneDistanceType type)
	{
		if (type==PlaneDistanceRM2)
			return distance_rm2(plane);
		else if (type==PlaneDistanceNormD)
			return distance_normD(plane);
		return 0;
	}

	double ZPlane::distance(const ZPlane& plane, ZPlaneDistanceType type) const 
	{
		if (type==PlaneDistanceRM2)
			return distance_rm2(plane);
		else if (type==PlaneDistanceNormD)
			return distance_normD(plane);

		return 0;
	}
}
#pragma once
#ifndef SPHERICAL_SYSTEM_H_
#define SPHERICAL_SYSTEM_H_

#include "../Math/vector3.h"
#include "../Eigen/Eigen/Dense"

namespace MATH
{
	template <typename T>
	static CVector3<T> cartesian2Spherical(const CVector3<T>& normal, float d)
	{
		CVector3<T> ret;
		float cosTheta = normal.z;
		if (g_isSameValue(cosTheta, 1.f))
			return CVector3<T>(0, 0, d);
		else if (g_isSameValue(cosTheta, -1.f))
			return CVector3<T>(Z_PI, 0, d);

		float sinTheta = sqrt(1-cosTheta*cosTheta);
		float cosPhi = normal.x / sinTheta;
		float sinPhi = normal.y / sinTheta;
		cosPhi = cosPhi<-1.f ? -1.f : (cosPhi > 1.f ? 1.f : cosPhi);
		sinPhi = sinPhi<-1.f ? -1.f : (sinPhi > 1.f ? 1.f : sinPhi);
		float theta = acos(cosTheta);
		float phi = acos(cosPhi);
		if (sinPhi<0)
			phi = Z_PI*2-phi;
		ret = CVector3<T>(theta, phi, d);
		return ret;
	}

	// NOTE: the d is same as spherical.z
	template <typename T>
	static CVector3<T> spherical2Cartesian(const CVector3<T>& spherical)
	{
		float cosTheta = cos(spherical.x);
		float sinTheta = sin(spherical.x);
		float cosPhi = cos(spherical.y);
		float sinPhi = sin(spherical.y);
		CVector3<T> ret(sinTheta*cosPhi, sinTheta*sinPhi, cosTheta);
		return ret;
	}

	static Eigen::VectorXf cartesian2Spherical(const Eigen::VectorXf& normal, float d)
	{
		Vec3f inV(normal(0), normal(1), normal(2));
		Vec3f ret = cartesian2Spherical(inV, d);
		Eigen::VectorXf retE(3);
		retE(0) = ret.x;
		retE(1) = ret.y;
		retE(2) = ret.z;
		return retE;
	}

	static Eigen::VectorXf spherical2Cartesian(const Eigen::VectorXf& spherical)
	{
		Vec3f sphV(spherical(0), spherical(1), spherical.size()==3 ? spherical(2) : 1.f);
		Vec3f retV = spherical2Cartesian(sphV);
		Eigen::VectorXf retE(3);
		retE(0) = retV.x;
		retE(1) = retV.y;
		retE(2) = retV.z;
		return retE;
	}
}

#endif//SPHERICAL_SYSTEM_H_
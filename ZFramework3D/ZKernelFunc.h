#pragma once
#ifndef ZKERNEL_FUNC_H_
#define ZKERNEL_FUNC_H_

#include "GlobalDefs.h"
#include "../Scene/Mesh3D.h"
#include "../Eigen/Eigen/Dense"

namespace ZMeshSpace
{
	class ZKernelFuncs
	{
	public:
		static float GaussianKernelFunc(float f, float sigma)
		{
			return exp(-0.5f*f*f/sigma/sigma);
		}

		static float GaussianKernelFuncA(const Eigen::VectorXf& v1, const Eigen::VectorXf& v2, float sigma)
		{
			float dist = (v1-v2).squaredNorm();
			return exp(-0.5f*dist/sigma/sigma);
		}

		static float GaussianKernelFuncA(const Eigen::Vector3f& v1, const Eigen::Vector3f& v2, float sigma)
		{
			float dist = (v1-v2).squaredNorm();
			return exp(-0.5f*dist/sigma/sigma);
		}

		static float GaussianKernelFuncA(const Vec3f& v1, const Vec3f& v2, float sigma)
		{
			float dist = (v1-v2).norm2();
			return exp(-0.5f*dist/sigma/sigma);
		}

		static float GaussianKernelFuncN(const Eigen::VectorXf& n1, const Eigen::VectorXf& n2, float sigma)
		{
			float dist = 1.f-n1.dot(n2)/n1.norm()/n2.norm();
			return exp(-0.5*dist*dist/sigma/sigma);
		}

		static float GaussianKernelFuncA3(const Eigen::VectorXf& v1, const Eigen::VectorXf& v2, float sigma)
		{
			Eigen::VectorXf v = v1-v2;
			float dist = v(0)*v(0)+v(1)*v(1)+v(2)*v(2);
			return exp(-0.5f*dist/sigma/sigma);
		}

		static float GaussianKernelFuncNormal3(const Eigen::VectorXf& n1, const Eigen::VectorXf& n2, float sigma)
		{
			Eigen::Vector3f n11(n1.data());
			Eigen::Vector3f n21(n2.data());
			n11.normalize();
			n21.normalize();
			float dist = 1.f-n11.dot(n21);
			return exp(-0.5*dist*dist/sigma/sigma);

		}

		static float GaussianKernelFunc(const Eigen::VectorXf& v, float sigma)
		{
			float length = v.squaredNorm();
			return exp(-0.5f*length/sigma/sigma);
		}

		static float GaussianKernelSphericalFunc(const Eigen::VectorXf& v1, const Eigen::VectorXf& v2, float sigma)
		{
			float cosTheta1 = cos(v1(0));
			float sinTheta1 = sin(v1(0));
			float cosPhi1 = cos(v1(1));
			float sinPhi1 = sin(v1(1));
			Eigen::Vector3f cv1(sinTheta1*cosPhi1, sinTheta1*sinPhi1, cosTheta1);
			float cosTheta2 = cos(v2(0));
			float sinTheta2 = sin(v2(0));
			float cosPhi2 = cos(v2(1));
			float sinPhi2 = sin(v2(1));
			Eigen::Vector3f cv2(sinTheta2*cosPhi2, sinTheta2*sinPhi2, cosTheta2);
			float length = 1.f-cv1.dot(cv2);
			return exp(-0.5f*length/sigma/sigma);
		}
	};
}

#endif//ZKERNEL_FUNC_H_
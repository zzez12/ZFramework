#pragma once
#ifndef ZMESH_JOINT_BILATERAL_FILTERING_H_
#define ZMESH_JOINT_BILATERAL_FILTERING_H_

#include "ZMeshFilter.h"
#include "AnnWarper_Eigen.h"

namespace ZMeshSpace
{
	class ZMeshJointBilateralFiltering : public ZMeshFilter
	{
	public:
		ZMeshJointBilateralFiltering(Mesh3D* mesh);
		virtual ~ZMeshJointBilateralFiltering();

		virtual bool apply();
		virtual bool apply(const Eigen::MatrixXf& input, const std::vector<bool>& tags);
		virtual bool apply(const Eigen::MatrixXf& input, const Eigen::VectorXf& weights, const std::vector<bool>& tags);
		virtual bool apply(const Eigen::MatrixXf& input, int spatialDim, int rangeDim, const Eigen::VectorXf& weights, const std::vector<bool>& tags);
		virtual bool apply(const Eigen::MatrixXf& input, const Eigen::MatrixXf& joint_input, const Eigen::VectorXf& weights, const std::vector<bool>& tags);

	private: 
		AnnWarper_Eigen *pAnnSpatialSearch_;
		KernelFuncA rangeKernelFunc_;

	private:
		void init();
	};
}


#endif//ZMESH_JOINT_BILATERAL_FILTERING_H_
#pragma once
#ifndef ZMESH_BALATERAL_FILTER_H_
#define ZMESH_BALATERAL_FILTER_H_

#include "ZMeshFilter.h"
#include "../Eigen/Eigen/Dense"
#include "AnnWarper_Eigen.h"

namespace ZMeshSpace
{
	class ZMeshBilateralFilter : public ZMeshFilter
	{
	public:
		ZMeshBilateralFilter(Mesh3D* mesh);
		virtual ~ZMeshBilateralFilter();

		virtual bool apply();
		virtual bool apply(const Eigen::MatrixXf& input, const std::vector<bool>& tags);
		virtual bool apply(const Eigen::MatrixXf& input, const Eigen::VectorXf& weights, const std::vector<bool>& tags);
		virtual bool apply(const Eigen::MatrixXf& input, int rangeDim, int spatialDim, const Eigen::VectorXf& weights, const std::vector<bool>& tags);

		void setAnnSearchHandle(AnnWarper_Eigen *handle){pAnnSearch_=handle;}
		void setRangeKernelFunc(KernelFuncA func) {rangeKernelFunc_=func;}

	private:
		AnnWarper_Eigen *pAnnSearch_;
		KernelFuncA rangeKernelFunc_;

	private:
		void init();
	};
}

#endif//ZMESH_BALATERAL_FILTER_H_
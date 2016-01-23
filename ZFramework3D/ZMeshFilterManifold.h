#pragma once
#ifndef ZMESHFILTER_MANIFOLD_H_
#define ZMESHFILTER_MANIFOLD_H_

#include "ZMeshFilter.h"
#include "AnnWarper_Eigen.h"

namespace ZMeshSpace
{
	class ZMeshFilterManifold : public ZMeshFilter
	{
	public:
		ZMeshFilterManifold(Mesh3D* mesh);
		virtual ~ZMeshFilterManifold();

		virtual bool apply();
		virtual bool apply(const Eigen::MatrixXf& input, const std::vector<bool>& tags);
		virtual void setPara(const std::string& name, float value);

		void setRangeWeight(int level, float weight);
		void updateRange();
		Eigen::VectorXf getPointClusterIds() {return verticeClusterIds;}

	private:
		void init(const Eigen::MatrixXf& input);
		void destroy();
		bool buildManifoldAndPerformFiltering(const Eigen::MatrixXf& input, 
			const Eigen::MatrixXf& etaK, const std::vector<bool>& clusterK,
			float sigma_s, float sigma_r,
			int currentTreeLevel);

		void computeTreeHeight();
		void initMinPixelDist2Manifold();
		Eigen::VectorXf computeMaxEigenVector(const Eigen::MatrixXf& inputMat, const std::vector<bool>& clusterK);

	private:
		Eigen::MatrixXf sum_w_ki_Psi_blur_;
		Eigen::VectorXf sum_w_ki_Psi_blur_0_;
		Eigen::VectorXf min_pixel_dist_to_manifold_squared_;

		std::vector<Eigen::MatrixXf> wki_Psi_blurs_;
		std::vector<Eigen::VectorXf> wki_Psi_blur_0s_;
		std::vector<float> rangeWeights_;
		Eigen::VectorXf verticeClusterIds;

		ZMeshFilter *pGaussianFilter_;
		AnnWarper_Eigen *pAnnSearch_;
	};
}

#endif//ZMESHFILTER_MANIFOLD_H_
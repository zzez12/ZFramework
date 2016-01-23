#pragma once
#ifndef ZMESHFILTER_WEIGHTED_GAUSSIAN_H_
#define ZMESHFILTER_WEIGHTED_GAUSSIAN_H_

#include "ZMeshFilterGaussian.h"

namespace ZMeshSpace
{
	class ZMeshFilterWeightedGaussian : public ZMeshFilterGaussian
	{
	public:
		ZMeshFilterWeightedGaussian(Mesh3D* mesh);
		virtual ~ZMeshFilterWeightedGaussian();

		bool apply(const Eigen::MatrixXf& input, const Eigen::VectorXf& weight, const std::vector<bool>& tags);
	};
}

#endif//ZMESHFILTER_WEIGHTED_GAUSSIAN_H_
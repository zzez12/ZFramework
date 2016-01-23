#pragma once
#ifndef ZMESHFILTER_GAUSSIAN_H_
#define ZMESHFILTER_GAUSSIAN_H_

#include "ZMeshFilter.h"

namespace ZMeshSpace
{
	class ZMeshFilterGaussian : public ZMeshFilter
	{
	public:
		ZMeshFilterGaussian(Mesh3D* mesh);
		virtual ~ZMeshFilterGaussian();

		virtual void setPara(const std::string& name, float value);
		virtual bool apply();
		virtual bool apply(const Eigen::MatrixXf& input, const std::vector<bool>& tags);
	};
}

#endif//ZMESHFILTER_GAUSSIAN_H_
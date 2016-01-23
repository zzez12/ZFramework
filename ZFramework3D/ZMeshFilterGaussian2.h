#pragma once
#ifndef ZMESHFILTER_GAUSSIAN_2_H_
#define ZMESHFILTER_GAUSSIAN_2_H_

#include "ZMeshFilter.h"
#include "AnnWarper_Eigen.h"

namespace ZMeshSpace
{
	class ZMeshFilterGaussian2 : public ZMeshFilter
	{
	public:
		ZMeshFilterGaussian2(Mesh3D* mesh);
		virtual ~ZMeshFilterGaussian2();

		virtual bool apply();
		virtual bool apply(const Eigen::MatrixXf& input, const std::vector<bool>& tags);
		virtual bool apply(const Eigen::MatrixXf& input, const Eigen::VectorXf& weights, const std::vector<bool>& tags);
		virtual void setPara(const std::string& name, float value);

		//void setAnnSearchHandle(AnnWarper3 *handle);
		void setAnnSearchHandle(AnnWarper_Eigen *handle){pAnnSearch_=handle;}

	private:
		//AnnWarper3 *pAnnSearch_;
		AnnWarper_Eigen *pAnnSearch_;

	private:
		void init();
	};
}

#endif//ZMESHFILTER_GAUSSIAN_2_H_
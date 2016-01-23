#include "ZMeshFilterGaussian2.h"

namespace ZMeshSpace
{
	ZMeshFilterGaussian2::ZMeshFilterGaussian2(Mesh3D* mesh)
		: ZMeshFilter(mesh, MESH_FILTER_GAUSSIAN2)
	{
		pAnnSearch_ = NULL;
		setKernelFunc(ZKernelFuncs::GaussianKernelFuncA);
		setPara(ZMeshFilterParaNames::SpatialSigma, 0);
	}

	ZMeshFilterGaussian2::~ZMeshFilterGaussian2()
	{

	}

	void ZMeshFilterGaussian2::setPara(const std::string& name, float value)
	{
		if (name.compare(ZMeshFilterParaNames::SpatialSigma)==0)
		{
			filterPara_.spatial_sigma = value;
		}
		else if (name.compare(ZMeshFilterParaNames::RangeSigma)==0)
		{
			filterPara_.range_sigma = value;
		}
	}

	bool ZMeshFilterGaussian2::apply()
	{
		return false;
	}

	bool ZMeshFilterGaussian2::apply(const Eigen::MatrixXf& input, const std::vector<bool>& tags)
	{
		if (pAnnSearch_==NULL)
			return false;

		int vSize = input.cols();
		output_ = input;
		pAnnSearch_->setFlags(tags);

		for (int i=0; i<vSize; i++)
		{
			if (!tags[i]) continue;

			Eigen::VectorXf v(input.row(i));

			Eigen::VectorXi nnIdx;
			Eigen::VectorXf nnDists;
			float searchRad = filterPara_.spatial_sigma*2;
			int queriedNum = pAnnSearch_->queryFRPoint(v, searchRad, 0, nnIdx, nnDists);

			// convolute them
			float sumWeight = 0;
			Eigen::VectorXf sumPos(input.cols());
			sumPos.fill(0);
			for (int j=1; j<queriedNum; j++) // NOT include the point itself
			{
				int idx = nnIdx(j);
				float distWeight = ZKernelFuncs::GaussianKernelFunc(sqrt(nnDists(j)), filterPara_.spatial_sigma);
				float rangeWeight = kernelFuncA_(v, input.row(idx), filterPara_.range_sigma);
				float weight = rangeWeight*distWeight;	// rangeWeight*spatialWeight
				sumWeight += weight;
				sumPos = sumPos + Eigen::VectorXf(input.row(idx))*weight;
			}
			if (!g_isZero(sumWeight))
				output_.row(i) = sumPos*(1.f/sumWeight);
			//else 
			//	std::cout << "Zero!" << i << "\n";
		}

		return true;
	}
	
	bool ZMeshFilterGaussian2::apply(const Eigen::MatrixXf& input, const Eigen::VectorXf& weights, const std::vector<bool>& tags)
	{
		if (pMesh_==NULL || !pMesh_->isvalid() || pAnnSearch_==NULL)
			return false;

		int vSize = pMesh_->get_num_of_vertex_list();
		output_ = input;
		pAnnSearch_->setFlags(tags);

		for (int i=0; i<vSize; i++)
		{
			if (!tags[i]) continue;

			Eigen::VectorXf v(input.row(i));

			Eigen::VectorXi nnIdx;
			Eigen::VectorXf nnDists;
			float searchRad = filterPara_.spatial_sigma*2;
			int queriedNum = pAnnSearch_->queryFRPoint(v, searchRad, 0, nnIdx, nnDists);

			// convolute them
			float sumWeight = 0;
			Eigen::VectorXf sumPos(input.cols());
			sumPos.fill(0);
			for (int j=1; j<queriedNum; j++) // NOT include the point itself
			{
				int idx = nnIdx(j);
				if (!tags[idx]) continue;
				float dist = ZKernelFuncs::GaussianKernelFunc(sqrt(nnDists(j)), filterPara_.spatial_sigma);
				float weight = kernelFuncA_(v, input.row(idx), filterPara_.range_sigma)*dist*weights(idx); // rangeWeight*spatialWeight*addtionalWeight
				sumWeight += weight;
				sumPos = sumPos + Eigen::VectorXf(input.row(idx))*weight;
			}
			if (!g_isZero(sumWeight))
				output_.row(i) = sumPos*(1.f/sumWeight);
// 			else 
// 				std::cout << "Zero!" << i << "\n";
		}

		return true;
	}
}
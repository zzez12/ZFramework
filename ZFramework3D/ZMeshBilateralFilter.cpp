#include "ZMeshBilateralFilter.h"
#include "ZKernelFunc.h"

namespace ZMeshSpace
{
	ZMeshBilateralFilter::ZMeshBilateralFilter(Mesh3D* mesh)
		: ZMeshFilter(mesh, MESH_FILTER_BILATERAL)
	{
		pAnnSearch_ = NULL;
		rangeDim_ = 3;
		spatialDim_ = 3;
		setKernelFunc(ZKernelFuncs::GaussianKernelFuncN);
	}

	ZMeshBilateralFilter::~ZMeshBilateralFilter()
	{

	}

	bool ZMeshBilateralFilter::apply()
	{
		return false;
	}

	bool ZMeshBilateralFilter::apply(const Eigen::MatrixXf& input, const std::vector<bool>& tags)
	{
		if (pAnnSearch_==NULL)
			return false;
		if (spatialDim_+rangeDim_!=input.cols())
			return false;

		int nSize = input.rows();
		output_ = input;
		pAnnSearch_->setFlags(tags);

		for (int i=0; i<nSize; i++)
		{
			if (!tags[i]) continue;

			Eigen::VectorXf v(input.row(i));
			Eigen::VectorXi nnIdx;
			Eigen::VectorXf nnDists;

			// query
			float searchRad = filterPara_.spatial_sigma*2;
			Eigen::VectorXf queryV(v.head(spatialDim_));
			Eigen::VectorXf rangeV(v.tail(rangeDim_));
			int queryNum = pAnnSearch_->queryFRPoint(queryV, searchRad, 0, nnIdx, nnDists);
			//int queryNum = queryMeshTool_->query(i, 20, nnIdx, nnDists);

			// convolute
			Eigen::VectorXf sumRange(rangeDim_);
			sumRange.fill(0);
			float sumWeight = 0;
			for (int j=1; j<queryNum; j++)
			{
				int idx = nnIdx[j];
				if (!tags[idx]) continue;
				Eigen::VectorXf rangeU(input.row(idx).tail(rangeDim_));
				float distWeight = ZKernelFuncs::GaussianKernelFunc(sqrt(nnDists(j)), filterPara_.spatial_sigma);
				// if kernelFuncA_==NULL, then only using spatial filter
				float rangeWeidht = kernelFuncA_ ? kernelFuncA_(rangeV, rangeU, filterPara_.range_sigma) : 1.f;	
				float weight = rangeWeidht*distWeight;
				sumWeight += weight;
				//sumRange += (rangeU-rangeV)*weight;
				sumRange += rangeU*weight;
			}
			if (!g_isZero(sumWeight))
				output_.row(i).tail(rangeDim_) = sumRange*(1.f/sumWeight);
		}

		return true;
	}

	bool ZMeshBilateralFilter::apply(const Eigen::MatrixXf& input, const Eigen::VectorXf& weights, const std::vector<bool>& tags)
	{
		if (pAnnSearch_==NULL)
			return false;
		if (getRangeDim()+getSpatialDim()!=input.cols())
			return false;

		int nSize = input.rows();
		output_ = input;
		pAnnSearch_->setFlags(tags);

		float searchRad = filterPara_.spatial_sigma*sqrt(3.0);
		for (int i=0; i<nSize; i++)
		{
			if (!tags[i]) continue;

			Eigen::VectorXf v(input.row(i));
			Eigen::VectorXi nnIdx;
			Eigen::VectorXf nnDists;

			// query
			Eigen::VectorXf queryV(v.head(spatialDim_));
			Eigen::VectorXf rangeV(v.tail(rangeDim_));
			int queryNum = pAnnSearch_->queryFRPoint(queryV, searchRad, 0, nnIdx, nnDists);
			//int queryNum = queryMeshTool_->query(i, 20, nnIdx, nnDists);

			// convolute
			Eigen::VectorXf sumRange(rangeDim_);
			sumRange.fill(0);
			float sumWeight = 0;
			for (int j=1; j<queryNum; j++)
			{
				int idx = nnIdx[j];
				if (!tags[idx]) continue;
				Eigen::VectorXf rangeU(input.row(idx).tail(rangeDim_));
				float distWeight = ZKernelFuncs::GaussianKernelFunc(sqrt(nnDists(j)), filterPara_.spatial_sigma);
				// if kernelFuncA_==NULL, then only using spatial filter
				float rangeWeidht = kernelFuncA_ ? kernelFuncA_(rangeV, rangeU, filterPara_.range_sigma) : 1.f;
				float weight = rangeWeidht*distWeight*weights(idx);
				//if (i==1)
				//	std::cout << rangeU << " * " << distWeight << "*" << rangeWeidht << "*" << weights(idx) << "\n";
				sumWeight += weight;
				sumRange += rangeU*weight;
			}
			if (!g_isZero(sumWeight))
				output_.row(i).tail(rangeDim_) = sumRange*(1.f/sumWeight);
		}

		return true;
	}

	bool ZMeshBilateralFilter::apply(const Eigen::MatrixXf& input, int spatialDim, int rangeDim, const Eigen::VectorXf& weights, const std::vector<bool>& tags)
	{
		setRangeDim(rangeDim);
		setSpatialDim(spatialDim);
		Eigen::VectorXf weights2;
		std::vector<bool> tags2;
		bool bW = (weights.size()==0);
		bool bT = (tags.empty());
		if (bW)
		{
			weights2.resize(input.rows());
			weights2.fill(1.f);
		}
		if (bT)
		{
			tags2.resize(input.rows(), true);
		}
		return apply(input, bW ? weights2 :weights, bT ? tags2 : tags);
	}
}
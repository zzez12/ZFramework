#include "ZMeshJointBilateralFiltering.h"

namespace ZMeshSpace
{
	ZMeshJointBilateralFiltering::ZMeshJointBilateralFiltering(Mesh3D* mesh)
		: ZMeshFilter(mesh, MESH_FILTER_JOINT_BILATERAL)
	{
		pAnnSpatialSearch_ = NULL;
		setRangeDim(3);
		setSpatialDim(3);
	}

	ZMeshJointBilateralFiltering::~ZMeshJointBilateralFiltering()
	{

	}

	void ZMeshJointBilateralFiltering::init()
	{

	}

	bool ZMeshJointBilateralFiltering::apply()
	{
		return false;
	}

	bool ZMeshJointBilateralFiltering::apply(const Eigen::MatrixXf& input, int rangeDim, int spatialDim, const Eigen::VectorXf& weights, const std::vector<bool>& tags)
	{
		Eigen::MatrixXf joint_input = input;
		return apply(input, joint_input, weights, tags);
	}

	bool ZMeshJointBilateralFiltering::apply(const Eigen::MatrixXf& input, const Eigen::VectorXf& weights, const std::vector<bool>& tags)
	{
		Eigen::MatrixXf joint_input = input;
		return apply(input, joint_input,  weights, tags);
	}

	bool ZMeshJointBilateralFiltering::apply(const Eigen::MatrixXf& input, const std::vector<bool>& tags)
	{
		Eigen::MatrixXf joint_input = input;
		Eigen::VectorXf weights(input.rows());
		weights.fill(1.f);
		return apply(input, joint_input, weights, tags);
	}

	bool ZMeshJointBilateralFiltering::apply(const Eigen::MatrixXf& input, const Eigen::MatrixXf& joint_input, 
		const Eigen::VectorXf& weights, const std::vector<bool>& tags)
	{
		setSpatialDim(input.cols());
		setRangeDim(joint_input.cols());

		if (pAnnSpatialSearch_==NULL)
			return false;

		int nSize = input.rows();
		output_ = input;
		pAnnSpatialSearch_->setFlags(tags);
		float searchRad = filterPara_.spatial_sigma*sqrt(3.0);
		for (int i=0; i<nSize; i++)
		{
			if (!tags[i]) continue;

			Eigen::VectorXf queryV(input.row(i));
			Eigen::VectorXi nnIdx;
			Eigen::VectorXf nnDists;

			// query
			Eigen::VectorXf rangeV(joint_input.row(i));
			int queryNum = pAnnSpatialSearch_->queryFRPoint(queryV, searchRad, 0, nnIdx, nnDists);

			// convolute
			Eigen::VectorXf sumRange(rangeDim_);
			sumRange.fill(0);
			float sumWeight = 0;
			for (int j=1; j<queryNum; j++)
			{
				int idx = nnIdx[j];
				Eigen::VectorXf rangeU(joint_input.row(idx));
				float distWeight = ZKernelFuncs::GaussianKernelFunc(nnDists(j), filterPara_.spatial_sigma);
				// if kernelFuncA_==NULL, then only using spatial filter
				float rangeWeidht = kernelFuncA_ ? kernelFuncA_(rangeV, rangeU, filterPara_.range_sigma) : 1.f;
				float weight = rangeWeidht*distWeight*weights(idx);
				sumWeight += weight;
				sumRange += rangeU*weight;
			}
			if (!g_isZero(sumWeight))
				output_.row(i) = sumRange*(1.f/sumWeight);
		}

		return true;
	}
}
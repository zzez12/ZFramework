#include "ZMeshFilterWeightedGaussian.h"

namespace ZMeshSpace
{
	ZMeshFilterWeightedGaussian::ZMeshFilterWeightedGaussian(Mesh3D* mesh)
		: ZMeshFilterGaussian(mesh)
	{
		type_ = MESH_FILTER_WEIGHTED_GAUSSIAN;
	}

	ZMeshFilterWeightedGaussian::~ZMeshFilterWeightedGaussian()
	{

	}

	bool ZMeshFilterWeightedGaussian::apply(const Eigen::MatrixXf& input, const Eigen::VectorXf& weights, const std::vector<bool>& tags)
	{
		int vSize = pMesh_->get_num_of_vertex_list();
		output_.resize(vSize, input.cols());
		output_.fill(0);
// 		output_ = input;

		// only average 1-ring neighborhood
		for (int i=0; i<vSize; i++)
		{
			if (!tags[i]) continue;

			HE_vert* v = pMesh_->get_vertex(i);
			Eigen::VectorXf sumPos(input.cols());//=input.col(v->m_id);
			sumPos.fill(0);
			float sumWeight = 0;
			HE_edge* e = v->m_pedge;
			do 
			{
				HE_vert* vNeibor = e->m_pvert;
				if (tags[vNeibor->m_id]) 
				{
					float weight = kernelFuncA_(input.row(v->m_id), input.row(vNeibor->m_id), filterPara_.spatial_sigma)*weights(vNeibor->m_id);
					sumWeight += weight;
					sumPos = sumPos + Eigen::VectorXf(input.row(vNeibor->m_id))*weight;
				}
				e = e->m_ppair->m_pnext;
			} while (e!=v->m_pedge && e!=NULL);
			if (!g_isZero(sumWeight))
				output_.row(i) = sumPos*(1.f/sumWeight);
			//else
			//	std::cout << "Zero!" << i << "\n";
		}
		return true;
	}
}
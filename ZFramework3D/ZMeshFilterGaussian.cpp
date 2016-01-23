#include "ZMeshFilterGaussian.h"

namespace ZMeshSpace
{
	ZMeshFilterGaussian::ZMeshFilterGaussian(Mesh3D* mesh)
		: ZMeshFilter(mesh, MESH_FILTER_GAUSSIAN)
	{
		setKernelFunc(ZKernelFuncs::GaussianKernelFuncA);
		setPara(ZMeshFilterParaNames::SpatialSigma, 0);
	}

	ZMeshFilterGaussian::~ZMeshFilterGaussian()
	{

	}

	void ZMeshFilterGaussian::setPara(const std::string& name, float value)
	{
		if (name.compare(ZMeshFilterParaNames::SpatialSigma)==0)
		{
			filterPara_.spatial_sigma = value;
		}
	}

	bool ZMeshFilterGaussian::apply()
	{
		
		return true;
	}

	bool ZMeshFilterGaussian::apply(const Eigen::MatrixXf& input, const std::vector<bool>& tags)
	{
		int vSize = pMesh_->get_num_of_vertex_list();
		output_.resize(vSize, input.cols());
		output_.fill(0);
//		output_ = input;

		// only average 1-ring neighborhood
		for (int i=0; i<vSize; i++)
		{
			if (!tags[i]) continue;
// 			if (i==531) 
// 				DebugBreak();

			HE_vert* v = pMesh_->get_vertex(i);
			Eigen::VectorXf sumPos(input.cols());//=input.row(v->m_id);
			sumPos.fill(0);
			float sumWeight = 0;
			HE_edge* e = v->m_pedge;
			do 
			{
				HE_vert* vNeibor = e->m_pvert;
				if (tags[vNeibor->m_id]) 
				{
					float weight = kernelFuncA_(input.row(v->m_id), input.row(vNeibor->m_id), filterPara_.spatial_sigma);
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
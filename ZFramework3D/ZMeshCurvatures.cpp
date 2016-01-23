#include "ZMeshCurvatures.h"
#include "MeshTools.h"

namespace ZMeshSpace
{
	ZMeshCurvatures::ZMeshCurvatures(Mesh3D* pMesh)
	{
		mesh_ = pMesh;
	}

	bool ZMeshCurvatures::apply()
	{
		if (mesh_==NULL || !mesh_->isvalid()) return false;

		int vSize = mesh_->get_num_of_vertex_list();
		int eSize = mesh_->get_num_of_edges_list();

		// weighted normals
		weightedNormals_.resize(vSize);
		for (int i=0; i<vSize; i++)
		{
			HE_vert* vert = mesh_->get_vertex(i);
			Eigen::Vector3f nor(0,0,0);
			HE_edge* edge = vert->m_pedge;
			do 
			{
				HE_face* f = edge->m_pface;
				if (f!=NULL)
				{
					Eigen::Vector3f fNor(f->m_vnormal);
					nor += fNor*MeshTools::getFaceArea(f);
				}
				edge = edge->m_ppair->m_pnext;
			} while (edge!=vert->m_pedge && edge!=NULL);
			nor.normalize();
			weightedNormals_[i] = nor;
		}

		// directional curvatures 
		directionalCurvatures_.resize(eSize);
		for (int i=0; i<eSize; i++)
		{
			HE_edge* edge = mesh_->get_edge(i);
			HE_vert* sVert = edge->m_ppair->m_pvert;
			HE_vert* tVert = edge->m_pvert;
			Eigen::Vector3f vji(tVert->m_vpos-sVert->m_vpos);
			float dc = weightedNormals_[sVert->m_id].dot(vji)*2.f/vji.squaredNorm();
			directionalCurvatures_[i] = dc;
		}

		// tensor of curvature
		curvatureTensor_.resize(vSize);
		for (int i=0; i<vSize; i++)
		{
			HE_vert* vert = mesh_->get_vertex(i);
			Eigen::Matrix3f mat;
			mat.fill(0);
			Eigen::Vector3f vi(vert->m_vpos);
			Eigen::Vector3f& norI = weightedNormals_[i];
			Eigen::Matrix3f nnt = Eigen::Matrix3f::Identity() - norI*norI.adjoint();
			float allAreas = 0;
			HE_edge* edge = vert->m_pedge;
			do 
			{
				HE_vert* vertJ = edge->m_pvert;
				Eigen::Vector3f vj(vertJ->m_vpos);
				Eigen::Vector3f Tij = nnt*(vi-vj);
				Tij.normalize();
				float wij = MeshTools::getFaceArea(edge) + MeshTools::getFaceArea(edge->m_ppair);
				allAreas += wij;
				float kij = directionalCurvatures_[edge->m_id];
				mat += Tij*Tij.adjoint()*wij*kij;
				edge = edge->m_ppair->m_pnext;
			} while (edge!=vert->m_pedge && edge!=NULL);
			mat *= (1.f/allAreas);
			curvatureTensor_[i] = mat;
		}

		// principal direction
		principalDirections1_.resize(vSize);
		principalDirections2_.resize(vSize);
		for (int i=0; i<vSize; i++)
		{
			HE_vert* vert = mesh_->get_vertex(i);
			Eigen::Matrix3f mat = curvatureTensor_[i];
			Eigen::Vector3f nor = weightedNormals_[i];
			Eigen::Vector3f WAdd, WMinus;
			WAdd = Eigen::Vector3f(1,0,0)+nor;
			WMinus = Eigen::Vector3f(1,0,0)-nor;
			Eigen::Vector3f Wi;
			if (WMinus.squaredNorm()>WAdd.squaredNorm())
			{
				Wi = WMinus/WMinus.norm();
			}
			else
			{
				Wi = WAdd/WAdd.norm();
			}
			Eigen::Matrix3f Qi = Eigen::Matrix3f::Identity() - Wi*Wi.adjoint()*2;
			Eigen::Vector3f Ti1 = Qi.col(1);
			Eigen::Vector3f Ti2 = Qi.col(2);
			Eigen::Matrix3f Si = Qi.adjoint()*mat*Qi;
	//		std::cout << i << ":\n" << Si << "\n";
			// compute the rotation angle
			Eigen::Matrix2f mi=Si.block<2,2>(1,1);//(Si(1,1), Si(1,2), Si(2,1), Si(2,2));
			Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> eigensolver(mi);
			if (eigensolver.info()!=Eigen::Success) 
			{
				std::cerr << "Wrong when computing eigens!\n";
				continue;
			}
			Eigen::Matrix2f eigenVecs = eigensolver.eigenvectors();
			float cosTheta = eigenVecs(0,0);
			float sinTheta = eigenVecs(1,0);
			principalDirections1_[i] = Ti1*cosTheta - Ti2*sinTheta;
			principalDirections2_[i] = Ti1*sinTheta + Ti2*cosTheta;
		}

		return true;
	}

	void ZMeshCurvatures::setToMesh()
	{
		int vSize = mesh_->get_num_of_vertex_list();
		float fMax = directionalCurvatures_.maxCoeff();
		float fMin = directionalCurvatures_.minCoeff();
		float fMaxx = fMax>abs(fMin) ? fMax : (-fMin);
		float fMinn = -fMaxx;
		for (int i=0; i<vSize; i++)
		{
			HE_vert* vert = mesh_->get_vertex(i);
			Eigen::Vector3f d1 = principalDirections1_[i];
			Eigen::Vector3f d2 = principalDirections2_[i];
			vert->m_vPrincipalDir1 = Vec3f(d1(0), d1(1), d1(2));
			vert->m_vPrincipalDir2 = Vec3f(d2(0), d2(1), d2(2));
			float fVal = directionalCurvatures_[i];
			vert->m_colorValue = (fVal-fMinn)*0.5/fMaxx;
		}
	}
}
#ifndef MESH_TOOLS_H_
#define MESH_TOOLS_H_

#include "../Scene/Mesh3D.h"
#include "../Eigen/Dense"

#include "SphericalSystem.h"

namespace ZMeshSpace
{
	class MeshTools
	{
	public:
		static HE_edge* getEdge(HE_vert* v0, HE_vert* v1)
		{
			if (v0==NULL || v1==NULL || v0==v1) return NULL;
			HE_edge* e = v0->m_pedge;
			do 
			{
				if (e->m_pvert==v1)
				{
					return e;
				}
				e = e->m_ppair->m_pnext;
			} while (e!=v0->m_pedge);
			return NULL;
		}

		static HE_edge* getEdge(Mesh3D* pMesh, int vId0, int vId1)
		{
			HE_vert* v0 = pMesh->get_vertex(vId0);
			HE_vert* v1 = pMesh->get_vertex(vId1);
			if (v0==NULL || v1==NULL)
				return NULL;
			else
				return getEdge(v0, v1);
		}

		static float getTriangleArea(const Vec3f& pos0, const Vec3f& pos1, const Vec3f& pos2)
		{
			Vec3f v10 = pos1-pos0;
			Vec3f v20 = pos2-pos0;
			return (v10*v20).length();
		}

		static float getFaceArea(HE_face* f)
		{
			if (f==NULL) return 0.f;

			float areas = 0;
			HE_edge* e = f->m_pedge;
			HE_vert* v0 = e->m_pvert;
			e = e->m_pnext;
			HE_vert* v1 = e->m_pvert;
			e = e->m_pnext;
			HE_vert* v2 =NULL;
			while (e!=f->m_pedge)
			{
				v2 = e->m_pvert;
				areas += getTriangleArea(v0->m_vpos, v1->m_vpos, v2->m_vpos);
				e = e->m_pnext;
				v0 = v1;
				v1 = v2;
			}
			return areas;
		}

		static float getFaceArea(HE_edge* e)
		{
			if (e==NULL || e->m_pface==NULL) return 0;
			else return getFaceArea(e->m_pface);
		}

		static float getLaplacianCotValue(HE_edge* e)
		{
			if (e->m_pnext==NULL) return 0;

			Vec3f v0 = e->m_ppair->m_pvert->m_vpos;
			Vec3f v1 = e->m_pvert->m_vpos;
			Vec3f v2 = e->m_pnext->m_pvert->m_vpos;
			Vec3f v02 = v0 - v2;
			Vec3f v12 = v1 - v2;
			float sinL = (v02*v12).length();
			if (g_isZero(sinL))
				return 0;
			else
				return v02.DotProduct(v12)/abs(sinL);
		}

		static void changeMeshVrtColors(Mesh3D* pMesh, const Eigen::VectorXf& colors)
		{
			if (pMesh==NULL || pMesh->get_num_of_vertex_list()!=colors.size())
				return;

			float dmax = colors.maxCoeff();
			float dmin = colors.minCoeff();
			std::cout << "fMax: " << dmax << "\tfMin: " << dmin << "\n";
			//std::cout << " colors: " << colors.head(50) << "\n";
			for (int i=0; i<pMesh->get_num_of_vertex_list(); i++)
			{
				HE_vert* vert = pMesh->get_vertex(i);
				float col = (colors(i)-dmin)/(dmax-dmin);
				vert->m_colorValue = col;
			}
		}

		static Eigen::MatrixXf getAllVerticePos(Mesh3D* pMesh)
		{
			if (pMesh==NULL || !pMesh->isvalid())
				return Eigen::MatrixXf();

			int vSize = pMesh->get_num_of_vertex_list();
			Eigen::MatrixXf ret(vSize, 3);
			for (int i=0; i<vSize; i++)
			{
				Vec3f pos = pMesh->get_vertex(i)->m_vpos;
				ret(i, 0) = pos.x;
				ret(i, 1) = pos.y;
				ret(i, 2) = pos.z;
			}
			return ret;
		}

		static Eigen::MatrixXf getAllVerticeNormals(Mesh3D* pMesh)
		{
			if (pMesh==NULL || !pMesh->isvalid())
				return Eigen::MatrixXf();

			int vSize = pMesh->get_num_of_vertex_list();
			Eigen::MatrixXf ret(vSize, 3);
			for (int i=0; i<vSize; i++)
			{
				Vec3f pos = pMesh->get_vertex(i)->m_vnormal;
				ret(i, 0) = pos.x;
				ret(i, 1) = pos.y;
				ret(i, 2) = pos.z;
			}
			return ret;
		}

		static Eigen::MatrixXf getAllFaceNormals(Mesh3D* pMesh)
		{
			if (pMesh==NULL || !pMesh->isvalid())
				return Eigen::MatrixXf();

			int fSize = pMesh->get_num_of_faces_list();
			Eigen::MatrixXf ret(fSize, 3);
			for (int i=0; i<fSize; i++)
			{
				Vec3f pos = pMesh->get_face(i)->m_vnormal;
				ret(i, 0) = pos.x;
				ret(i, 1) = pos.y;
				ret(i, 2) = pos.z;
			}
			return ret;
		}

		static Eigen::MatrixXf getAllFaceNormalsSpherical(Mesh3D* pMesh)
		{
			if (pMesh==NULL || !pMesh->isvalid())
				return Eigen::MatrixXf();

			int fSize = pMesh->get_num_of_faces_list();
			Eigen::MatrixXf ret(fSize, 2);
			for (int i=0; i<fSize; i++)
			{
				Vec3f pos = MATH::cartesian2Spherical(pMesh->get_face(i)->m_vnormal, 1);
				ret(i, 0) = pos.x;
				ret(i, 1) = pos.y;
			}
			return ret;
		}

		static Eigen::MatrixXf getAllFaceCenters(Mesh3D* pMesh)
		{
			if (pMesh==NULL || !pMesh->isvalid())
				return Eigen::MatrixXf();

			int fSize = pMesh->get_num_of_faces_list();
			Eigen::MatrixXf ret(fSize, 3);
			for (int i=0; i<fSize; i++)
			{
				HE_face* f = pMesh->get_face(i);
				Vec3f center(0,0,0);
				int count = 0;
				HE_edge* e = f->m_pedge;
				do 
				{
					center = center + e->m_pvert->m_vpos;
					count++;
					e = e->m_pnext;
				} while (e!=f->m_pedge);
				ret(i, 0) = center.x/count;
				ret(i, 1) = center.y/count;
				ret(i, 2) = center.z/count;
			}
			return ret;
		}

		static void setAllVerticePos(Mesh3D* pMesh, const Eigen::MatrixXf& newPos)
		{
			if (pMesh==NULL || pMesh->get_num_of_vertex_list()!=newPos.rows())
				return;

			int vSize = pMesh->get_num_of_vertex_list();
			for (int i=0; i<vSize; i++)
			{
				HE_vert* vert = pMesh->get_vertex(i);
				vert->m_vpos = Vec3f(newPos(i, 0), newPos(i, 1), newPos(i, 2));
			}
		}

		static void setAllFaceNormals(Mesh3D *pMesh, const Eigen::MatrixXf& newNormals)
		{
			if (pMesh==NULL || pMesh->get_num_of_faces_list()!=newNormals.rows())
				return;

			int fSize = pMesh->get_num_of_faces_list();
			for (int i=0; i<fSize; i++)
			{
				HE_face* face = pMesh->get_face(i);
				face->m_vnormal = Vec3f(newNormals(i, 0), newNormals(i, 1), newNormals(i, 2));
			}
		}

		static void setAllFaceNormals2(Mesh3D* pMesh, const Eigen::MatrixXf& newNormals)
		{
			if (pMesh==NULL || pMesh->get_num_of_faces_list()!=newNormals.rows())
				return;

			int fSize = pMesh->get_num_of_faces_list();
			for (int i=0; i<fSize; i++)
			{
				HE_face* face = pMesh->get_face(i);
				Vec3f newN = Vec3f(newNormals(i, 0), newNormals(i, 1), newNormals(i, 2));
				newN.unify();
				face->m_vnormal2 = newN;
			}
		}

		static void setAllFaceNormal2Spherical(Mesh3D* pMesh, const Eigen::MatrixXf& newNormals)
		{
			if (pMesh==NULL || pMesh->get_num_of_faces_list()!=newNormals.rows())
				return;

			int fSize = pMesh->get_num_of_faces_list();
			for (int i=0; i<fSize; i++)
			{
				HE_face* face = pMesh->get_face(i);
				Eigen::VectorXf newN = MATH::spherical2Cartesian(newNormals.row(i));
				face->m_vnormal2 = Vec3f(newN(0), newN(1), newN(2));
			}
		}

		static void setAllFaceColorValue(Mesh3D* pMesh, const Eigen::VectorXf& colors)
		{
			if (pMesh==NULL || pMesh->get_num_of_faces_list()!=colors.size())
				return;

			float dmax = colors.maxCoeff();
			float dmin = colors.minCoeff();

			for (int i=0; i<colors.size(); i++)
			{
				HE_face* f = pMesh->get_face(i);
				f->m_value = (colors(i)-dmin)/(dmax-dmin);;
			}
		}
	};
	
}

#endif//MESH_TOOLS_H_
#include "ZMeshProjector.h"
#include "GlobalDefs.h"

namespace ZMeshSpace
{
	ZMeshProjector::ZMeshProjector(Mesh3D *pMesh)
	{
		pAnnSearch_ = NULL;
		setMesh(pMesh);
	}

	void ZMeshProjector::setMesh(Mesh3D* pMesh)
	{
		pMesh_ = pMesh;
		initAnnSearch();
	}

	void ZMeshProjector::initAnnSearch()
	{
		if (!pMesh_)
		{
			SAFE_DELETE(pAnnSearch_);
			return;
		}

		SAFE_DELETE(pAnnSearch_);
		pAnnSearch_ = new AnnWarper;
		std::vector<Vec3f> points;
		for (int i=0; i<pMesh_->get_num_of_vertex_list(); i++)
		{
			points.push_back(pMesh_->get_vertex(i)->m_vpos);
		}
		pAnnSearch_->init(points);
	}

	int ZMeshProjector::project(const Vec3f& p, Vec3f& projectedPoint)
	{
		if (!pMesh_ || !pMesh_->isvalid())
		{
			std::cerr << "Warning: Invalid mesh3D set. MeshProjector::project() failed."<<"\n";
			return -1;
		}

		vId_ = fId_ = -1;

		// 1. find the nearest N points on the mesh
		int nearestNum = 8;
		std::vector<int> nnIdx;
		std::vector<float> nnDists;
		pAnnSearch_->queryPoint(p, nearestNum, nnIdx, nnDists);

		//test
		//projectedPoint = pMesh_->get_vertex(nnIdx[0])->position();
		//vId_ = nnIdx[0];
		//return pMesh_->get_vertex(nnIdx[0])->pedge_->pface_->id();

		// 2. get the corresponding faces
		std::vector<bool> faceTags(pMesh_->get_num_of_faces_list(), false);
		for (unsigned int i=0; i<nnIdx.size(); i++)
		{
			HE_vert *vert = pMesh_->get_vertex(nnIdx[i]);
			HE_edge *edge = vert->m_pedge;
			do 
			{
				HE_face* face = edge->m_pface;
				faceTags[face->m_id] = true;
				edge = edge->m_ppair->m_pnext;
			} while (edge!=vert->m_pedge && edge!=NULL);
		}

		// 3. find the exact project point
		float minDist = 10e10;
		Vec3f retProjP;
		//bool bFind = false;
		int retFaceId = -1;
		for (unsigned int i=0; i<pMesh_->get_num_of_faces_list(); i++)
		{
			if (!faceTags[i])
			{
				continue;
			}
			std::vector<Vec3f> facePoints;
			HE_face *face = pMesh_->get_face(i);
			HE_edge *edge = face->m_pedge;
			do 
			{
				facePoints.push_back(edge->m_pvert->m_vpos);
				edge = edge->m_pnext;
			} while (edge!=face->m_pedge);
			Vec3f projectPoint;
			float dist;
			if (projectToFace(p, facePoints, projectPoint, dist))
			{
				if (dist<minDist)
				{
					//bFind = true;
					retFaceId = i;
					minDist = dist;
					retProjP = projectPoint;
				}
			}
		}

		if (retFaceId!=-1)
		{
			projectedPoint = retProjP;
			fId_ = retFaceId;
			return retFaceId;
		}
		else
		{
			// project the point to the nearest vertex
			projectedPoint = pMesh_->get_vertex(nnIdx[0])->m_vpos;
			vId_ = nnIdx[0];
			return pMesh_->get_vertex(nnIdx[0])->m_pedge->m_pface->m_id;
		}

		return -1;
	}

	bool ZMeshProjector::projectToFace(const Vec3f& p, const std::vector<Vec3f>& facePoints, 
		Vec3f& projectedPoint, float& dist)
	{
		if (facePoints.size()!=3)
		{
			std::cerr << "Warning: Cannot project point to non-triangle faces. MeshProjector::projectToFace() failed." << "\n";
			return false;
		}

		const Vec3f& p1 = facePoints[0];
		const Vec3f& p2 = facePoints[1];
		const Vec3f& p3 = facePoints[2];

		Vec3f n = (p2-p1)*(p3-p1);
		n.unify();
		float signedDist = (p-p1).DotProduct(n);
		//dist = ((p-p1)^n).length();
		Vec3f h = n*signedDist;
		projectedPoint = p - h;
		dist = abs(signedDist);

		return isPointInFace(projectedPoint, facePoints);
		//return true;
	}

	bool ZMeshProjector::isPointInFace(const Vec3f& q, const std::vector<Vec3f>& facePoints)
	{
		if (facePoints.size()!=3)
		{
			std::cerr << "Waring: Cannot compute for non-triangle faces. MeshProjector::isPointInFace() failed." << "\n";
			return false;
		}

		const Vec3f& p1 = facePoints[0];
		const Vec3f& p2 = facePoints[1];
		const Vec3f& p3 = facePoints[2];

		Vec3f qp1 = p1-q;
		Vec3f qp2 = p2-q;
		Vec3f qp3 = p3-q;

		Vec3f n12 = qp1 * qp2;
		Vec3f n23 = qp2 * qp3;
		Vec3f n31 = qp3 * qp1;
		n12.unify();
		n23.unify();
		n31.unify();

		if (g_isSameValue(n12.length(), 0))
		{
			std::cerr << "Length=0!" << "\n";
			return g_isSameValue(abs(n23.DotProduct(n31)), 1.f);
		}
		if (g_isSameValue(n23.length(), 0))
		{
			std::cerr << "Length=0!" << "\n";
			return g_isSameValue(abs(n31.DotProduct(n12)), 1.f);
		}
		if (g_isSameValue(n31.length(), 0))
		{
			std::cerr << "Length=0!" << "\n";
			return g_isSameValue(abs(n12.DotProduct(n23)), 1.f);
		}
		return g_isSameValue(n12.DotProduct(n23), 1.f) && g_isSameValue(n23.DotProduct(n31), 1.f);
	}

	ZMeshProjector::~ZMeshProjector()
	{
		SAFE_DELETE(pAnnSearch_);
	}
}
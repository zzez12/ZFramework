#include "ZQueryMesh.h"

namespace ZMeshSpace
{
	ZQueryMesh::ZQueryMesh(Mesh3D* pMesh, QueryMeshType type/* =QueryMesh_Face_Neighbor */)
		: pMesh_(pMesh), type_(type)
	{
		pAnnSearch_ = NULL;
	}


	ZQueryMesh::~ZQueryMesh()
	{

	}

	int ZQueryMesh::query(int queryId, int k, Eigen::VectorXi& outputIds, Eigen::VectorXf& dists)
	{
		if (type_==QueryMesh_Face_Neighbor)
		{
			return queryNeighborFaces(queryId, k, outputIds, dists);
		}
		return -1;
	}

	int ZQueryMesh::query(const Eigen::VectorXf& queryPos, float rad, Eigen::VectorXi& outputIds, Eigen::VectorXf& dists)
	{
		if (type_==QueryMesh_Face_Neighbor)
		{
			return queryNeighborFaces(queryPos, rad, outputIds, dists);
		}
		return -1;
	}

	int ZQueryMesh::queryNeighborFaces(int queryId, int k, Eigen::VectorXi& outputFaceIds, Eigen::VectorXf& dists)
	{
		HE_face* queryFace = pMesh_->get_face(queryId);
		if (queryFace==NULL) return false;

		HE_edge* edge = queryFace->m_pedge;
		std::vector<int> faceIds;
		faceIds.push_back(queryId);
		getNeighborFaces(queryFace, k, faceIds);
		int nSize = faceIds.size();
		outputFaceIds.resize(nSize);
		dists.resize(nSize);
		for (int i=0; i<faceIds.size(); i++)
		{
			outputFaceIds(i) = faceIds[i];
			dists(i) = (pMesh_->get_face(outputFaceIds(i))->m_vCenter - queryFace->m_vCenter).norm();
		}
		return -1;
	}

	int ZQueryMesh::queryNeighborFaces(const Eigen::VectorXf& queryPos, float rad, Eigen::VectorXi& outputIds, Eigen::VectorXf& dists)
	{
		if (pAnnSearch_==NULL)
			return -1;

		if(!pAnnSearch_->queryFRPoint(queryPos, rad*rad, 0, outputIds, dists))
			return -1;
		return outputIds.rows();
	}

	bool ZQueryMesh::getNeighborFaces(HE_face* startFace, int stop, std::vector<int>& faceIds)
	{
		if (faceIds.size()>=stop)
			return true;

		HE_edge* edge = startFace->m_pedge;
		std::vector<HE_face*> neighFaces;
		do 
		{
			HE_face* innerFace = edge->m_ppair->m_pface;
			if (innerFace)
				neighFaces.push_back(innerFace);

			edge = edge->m_pnext;
		} while (edge!=startFace->m_pedge);

		// add this neighbors first
		for (int i=0; i<neighFaces.size(); i++)
		{
			faceIds.push_back(neighFaces[i]->m_id);
		}

		// then check the new added faces neighbors
		for (int i=0; i<neighFaces.size(); i++)
		{
			if (getNeighborFaces(neighFaces[i], stop, faceIds))
				return true;
		}
		return false;
	}
}
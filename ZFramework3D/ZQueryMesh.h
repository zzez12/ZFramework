#pragma once
#ifndef ZQUERY_MESH_H_
#define ZQUERY_MESH_H_

#include "GlobalDefs.h"
#include "../Scene/Mesh3D.h"
#include "../Eigen/Eigen/Dense"
#include "../AnnWarper_Eigen.h"

namespace ZMeshSpace
{
	enum QueryMeshType
	{
		QueryMesh_Face_Neighbor = 0,
		QueryMesh_Vertice_Neighbor
	};

	class ZQueryMesh
	{
	public:
		ZQueryMesh(Mesh3D* pMesh, QueryMeshType type=QueryMesh_Face_Neighbor);
		virtual ~ZQueryMesh();

		void setQueryType(QueryMeshType type){type_=type;}
		void setSearchTool(AnnWarper_Eigen *pSearch) {pAnnSearch_=pSearch;}

		int query(int queryId, int k, Eigen::VectorXi& outputIds, Eigen::VectorXf& dists);
		int query(const Eigen::VectorXf& queryPos, float rad, Eigen::VectorXi& outputIds, Eigen::VectorXf& dists);

	private:
		int queryNeighborFaces(int queryId, int k, Eigen::VectorXi& outputFaceIds, Eigen::VectorXf& dists);
		int queryNeighborFaces(const Eigen::VectorXf& queryPos, float rad, Eigen::VectorXi& outputIds, Eigen::VectorXf& dists);

		bool getNeighborFaces(HE_face* startFace, int stop, std::vector<int>& faceIds);

	private:
		QueryMeshType type_;
		Mesh3D* pMesh_;

		AnnWarper_Eigen *pAnnSearch_;
	};
}

#endif//ZQUERY_MESH_H_
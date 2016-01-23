#pragma once

#include "GlobalDefs.h"
#include "../Scene/Mesh3D.h"
#include "ANNWarper.h"


namespace ZMeshSpace
{
	/************************************************************************/
	/* Project 3D points to 3D mesh                                                                     */
	/************************************************************************/
	class ZMeshProjector
	{
	public:
		ZMeshProjector(Mesh3D* pMesh=NULL);
		virtual ~ZMeshProjector();

		void setMesh(Mesh3D* pMesh);

		/**
		* Calculates the object coordinates (projectedPoint) of
		* coordinate (p) when projected onto the given mesh. 
		* Returns faceId on successful projection. Else return -1.
		*/
		virtual int project(const Vec3f& p, Vec3f& projectedPoint);
		void getResult(int& vId, int& fId) {vId = vId_; fId = fId_;}

		AnnWarper* getAnnSearch() {return pAnnSearch_;}


	protected:
		void initAnnSearch();
		bool projectToFace(const Vec3f& p, const std::vector<Vec3f>& facePoints, 
			Vec3f& projectedPoint, float& dist);
		bool isPointInFace(const Vec3f& p, const std::vector<Vec3f>& facePoints);

	protected:
		Mesh3D *pMesh_;
		AnnWarper *pAnnSearch_;
		int vId_;
		int fId_;
	};
}
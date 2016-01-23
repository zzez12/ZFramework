#pragma once

#ifndef ZEXPMAPGENERATEWARP_H_
#define ZEXPMAPGENERATEWARP_H_

#include "GlobalDefs.h"
#include "../Scene/Mesh3D.h"
#include "../ExLib/expmap/ExpMapGenerator.h"
#include "../ExLib/expmap/VFTriangleMesh.h"
#include "ZMeshParser.h"

namespace ZMeshSpace
{
	class ZExpMapGenerateWarp
	{
	public:
		ZExpMapGenerateWarp(Mesh3D* pMesh=NULL);
		~ZExpMapGenerateWarp();

		void destroy();

		void setMesh(Mesh3D* pMesh);
		void setHitVertex(HE_vert* vert);
		void setHitVertex(const Vec3f& vPos, const Vec3f& vNormal);
		//void computeExpMap();
		void validateExpMap();
		Vec2f getUV(int vId);

	private:
		Mesh3D* mesh_;
		rms::VFTriangleMesh* vfMesh_;
		rms::ExpMapGenerator expMapGen_;
		rms::IMeshBVTree bvTree_;

		bool bExpMapInitialized_;
		bool bExpMapValid_;

		bool bSmoothNormals_;
		bool bUpwindAverage_;
		bool bUseConstantNormal_;
		bool bPreserveProjectedLengths_;

		float fDecalRadius_;
		float fBoundaryWidth_;
		float fMaxEdgeLength_;

		rms::Frame3f vSeedFrame_;
	};
}

#endif //ZEXPMAPGENERATEWARP_H_
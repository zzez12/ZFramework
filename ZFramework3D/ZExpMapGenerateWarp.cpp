#include "ZExpMapGenerateWarp.h"

namespace ZMeshSpace
{
	ZExpMapGenerateWarp::ZExpMapGenerateWarp(Mesh3D* pMesh/* =NULL */)
	{
		mesh_ = pMesh; 
		vfMesh_ = NULL;
		bExpMapInitialized_ = false;
		bExpMapValid_ = false;
		fDecalRadius_ = 0.3f;
		bSmoothNormals_ = true;
		bUpwindAverage_ = true;
		bUseConstantNormal_ = false;
		bPreserveProjectedLengths_ = true;
	}

	ZExpMapGenerateWarp::~ZExpMapGenerateWarp(){}

	void ZExpMapGenerateWarp::destroy()
	{
		if (vfMesh_)
			delete vfMesh_;
		vfMesh_ = NULL;
		expMapGen_.Reset();
		bvTree_.Clear();
		bExpMapInitialized_ = false;
		bExpMapValid_ = false;
		fDecalRadius_ = 0.3f;
		bSmoothNormals_ = true;
		bUpwindAverage_ = true;
		bUseConstantNormal_ = false;
		bPreserveProjectedLengths_ = true;
	}

	void ZExpMapGenerateWarp::setMesh(Mesh3D* pMesh)
	{
		destroy();

		mesh_ = pMesh;
		zMesh3D2VFMesh(mesh_, vfMesh_);

		bvTree_.SetMesh(vfMesh_);
		expMapGen_.SetSurface(vfMesh_, &bvTree_);

		fMaxEdgeLength_ = vfMesh_->GetMaxEdgeLength();
		if (!vfMesh_->HasUVSet(0))
			vfMesh_->AppendUVSet();
		vfMesh_->InitializeUVSet(0);

		fBoundaryWidth_ = vfMesh_->GetMaxEdgeLength() * 1.0f;

		//if (bExpMapInitialized_)
		//	bExpMapValid_ = false;
		bExpMapInitialized_ = true;
	}

	void ZExpMapGenerateWarp::setHitVertex(HE_vert* vert)
	{
		setHitVertex(vert->m_vpos, vert->m_vnormal);
	}

	void ZExpMapGenerateWarp::setHitVertex(const Vec3f& vPos, const Vec3f& vNormal)
	{
		Wml::Vector3f vHit(vPos);
		Wml::Vector3f vNor(vNormal);
		vSeedFrame_.Origin() = vHit;
		vSeedFrame_.AlignZAxis(vNor);
	}

	void ZExpMapGenerateWarp::validateExpMap()
	{
		if (!bExpMapInitialized_)
			return;
// 		if (bExpMapValid_)
// 			return;
		expMapGen_.SetUseNeighbourNormalSmoothing(bSmoothNormals_);
		expMapGen_.SetUseUpwindAveraging(bUpwindAverage_);
		expMapGen_.SetUseConstantNormal(bUseConstantNormal_);
		expMapGen_.SetPreserveProjectedLengths(bPreserveProjectedLengths_);

		float fParaRadius = fDecalRadius_ + fBoundaryWidth_;
		expMapGen_.SetSurfaceDistances(vSeedFrame_.Origin(), 1.1f*fMaxEdgeLength_,
			fParaRadius, &vSeedFrame_, &vSeedFrame_.Z());
		expMapGen_.CopyVertexUVs(vfMesh_, 0);

		bExpMapValid_ = true;
	}

	Vec2f ZExpMapGenerateWarp::getUV(int vId)
	{
		Wml::Vector2f uv;
		vfMesh_->GetUV(vId, 0, uv);
		return Vec2f(uv.X(), uv.Y());
	}
}
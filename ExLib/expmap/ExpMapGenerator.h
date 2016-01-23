// Ryan Schmidt
// Copyright (c) 2006-2010
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt

#ifndef __RMS_EXPMAP_GENERATOR_H
#define __RMS_EXPMAP_GENERATOR_H

#include <list>
#include <set>
#include <map>
#include <limits>
#include <queue>

#include "Frame.h"
#include "ParticleGrid.h"
#include "IMeshBVTree.h"
#include "WmlExtPlane3.h"


namespace rms {


class ExpMapParticle
{
public:
	struct ListEntry {
		ExpMapParticle * pParticle;
		ListEntry * pNext;
	};

	enum PropagationState {
		Frozen,
		Active,
		Inactive
	};

	ExpMapParticle(
		const Wml::Vector3f & vPosition = Wml::Vector3f::ZERO,
		const Wml::Vector3f & vNormal = Wml::Vector3f::ZERO);
	~ExpMapParticle();

	Wml::Vector3f & Position() { return m_vPosition; }
	Wml::Vector3f & Normal() { return m_vNormal; }

	float & SurfaceDistance() { return m_fSurfaceDistance; }
	float SurfaceDistance() const { return m_fSurfaceDistance; }

	Frame3f & WorldFrame() { return m_WorldFrame; }
	Wml::Vector2f & SurfaceVector() { return m_vSurfaceVector; }

	unsigned int & VertexID() { return m_nVertexID; }

	ListEntry * GetNeighbourList() { return m_pNeighbourList; }
	void AddNeighbour( ListEntry * pNeighbour );
	void ClearNeighbourList();

	PropagationState & State() { return m_eState; }
	ExpMapParticle * NearestParticle() { return m_pNearest; }
	void SetNearestParticle( ExpMapParticle * pNearest ) { m_pNearest = pNearest; }

	void Clear() {
		SurfaceDistance() = std::numeric_limits<float>::max();
		State() = ExpMapParticle::Inactive;
		SetNearestParticle(NULL);
		SurfaceVector() = Wml::Vector2f::ZERO;
	}

	bool m_bNbrFlag;
	bool m_bNormalSmoothed;

protected:
	Wml::Vector3f m_vPosition;
	Wml::Vector3f m_vNormal;

	float m_fSurfaceDistance;
	Wml::Vector2f m_vSurfaceVector;

	Frame3f m_WorldFrame;

	ListEntry * m_pNeighbourList;
	ExpMapParticle * m_pNearest;

	unsigned int m_nVertexID;

	PropagationState m_eState;

public:
	static Wml::Vector2f INVALID_PARAM;


	bool operator<( const ExpMapParticle & pParticle2 ) const
		{ return SurfaceDistance() < pParticle2.SurfaceDistance();  }


	// used to clear internal members - this is important sometimes...
	virtual void Reset() {
		m_fSurfaceDistance = 0.0f;
		m_pNeighbourList = NULL;
		m_nVertexID = IMesh::InvalidID;
		m_eState = ExpMapParticle::Inactive;
		m_pNearest = NULL;
		m_bNbrFlag = false;
	}
};


/*
 * This class is just a wrapper for an ExpMapParticle * that we
 * can put into STL classes (like std::multiset)
 */
class ParticleQueueWrapper
{
public:
	ParticleQueueWrapper( ExpMapParticle * pParticle )
		{ m_pParticle = pParticle; }
	ParticleQueueWrapper( float fSearchDistance )
		{ m_pParticle = NULL; m_fSearchDistance = fSearchDistance; }
	ExpMapParticle * Particle() { return m_pParticle; }
	float QueueValue() const 
		{ if ( m_pParticle == NULL ) return m_fSearchDistance; else return m_pParticle->SurfaceDistance(); }
	bool operator<( const ParticleQueueWrapper & pParticle2 ) const
		{ float f1 = QueueValue(); float f2 = pParticle2.QueueValue(); return f1 < f2; }
protected:
	ExpMapParticle * m_pParticle;
	float m_fSearchDistance;
};



class ExpMapGenerator
{
public:
	ExpMapGenerator();
	~ExpMapGenerator();

	void SetSurface( IMesh * pMesh, IMeshBVTree * pMeshBVTree );

	void SetUseNeighbourNormalSmoothing( bool bEnable );
	bool GetUseNeighbourNormalSmoothing() { return m_bUseNeighbourNormalSmoothing; }

	void SetUseUpwindAveraging( bool bEnable ) { m_bUseUpwindAveraging = bEnable; }
	bool GetUseUpwindAveraging() { return m_bUseUpwindAveraging; }

	void SetUseConstantNormal( bool bEnable ) { m_bUseSeedFrameAtParticles = bEnable; }
	bool GetUseConstantNormal() { return m_bUseSeedFrameAtParticles; }

	void SetPreserveProjectedLengths( bool bEnable ) { m_bPreserveProjectedLengths = bEnable; }
	bool GetPreserveProjectedLengths() { return m_bPreserveProjectedLengths; }

	void Reset();

	void SetSurfaceDistances( const Wml::Vector3f & vSeedPoint, float fNeighbourThreshold,
		float fStopDistance,
		const Frame3f * pLastFrame,
		const Wml::Vector3f * pSeedNormal );

	void CopyVertexUVs( IMesh * pMesh, IMesh::UVSetID nSetID );


protected:
	IMesh * m_pMesh;
	IMeshBVTree * m_pMeshBVTree;

	float m_fMaxEdgeLength;
	float m_fAvgEdgeLength;

	bool m_bUseMeshNeighbours;
	bool m_bPrecomputeNeighbours;
	bool m_bUseUpwindAveraging;
	bool m_bUseNeighbourNormalSmoothing;
	bool m_bUseSeedFrameAtParticles;
	bool m_bPreserveProjectedLengths;


	// vertex to particle map
	std::vector< unsigned int > m_vVtxToParticleMap;

	// neighbour list setup
	bool m_bNeighbourListsValid;
	float m_fMaxNbrDist;
	std::vector<ExpMapParticle *> m_vNeighbourBuf;
	void InitializeNeighbourLists(float fRadiusThreshold);
	void ClearNeighbourLists();
	void FindNeighbours( ExpMapParticle * pParticle,
						 std::vector<ExpMapParticle *> & vParticles);
	void FindNeighbours( const Wml::Vector3f & vPoint, const Wml::Vector3f & vNormal,
						 std::vector<ExpMapParticle *> & vParticles, 
						 float fRadiusThreshold, ExpMapParticle * pSkip = NULL );
	MemoryPool<ExpMapParticle::ListEntry> m_NeighbourMemoryPool;

	// neighbor list init
	ExpMapParticle::ListEntry * GetNeighbourList( ExpMapParticle * pParticle );

	void SetSmoothNormal( ExpMapParticle * pParticle, ExpMapParticle::ListEntry * pNbrs, bool bEnable );


	void ComputeExpMap( float fStopDistance );

	void RemoveNeighbours( ExpMapParticle * pParticle, std::multiset< ParticleQueueWrapper > & pq );
	void UpdateNeighbours( ExpMapParticle * pParticle, std::multiset< ParticleQueueWrapper > & pq );

	void PropagateFrameFromNearest( ExpMapParticle * pParticle );
	void PropagateFrameFromNearest_Average( ExpMapParticle * pParticle );

	void PrecomputePropagationData( ExpMapParticle * pCenterParticle, Wml::ExtPlane3f & vTangentPlane, 
										Frame3f & vCenterWorldFrame, Wml::Matrix2f & matFrameRotate );
	Wml::Vector2f ComputeSurfaceVector( ExpMapParticle * pCenterParticle, 
											ExpMapParticle * pNbrParticle,
											Wml::ExtPlane3f & vTangentPlane, 
											Frame3f & vCenterWorldFrame, 
											Wml::Matrix2f & matFrameRotate );



	//! seed particle
	ExpMapParticle * m_pSeedParticle;

	ExpMapParticle * InitializeSeedParticle( const Wml::Vector3f & vPosition,
		const Wml::Vector3f * pSeedNormal, const Frame3f * pLastSeedPointFrame);

	// particle grid (for point-based neighbour finding)
	ParticleGrid<ExpMapParticle *> m_particleGrid;
	bool m_bParticleGridValid;
	void InitializeParticleGrid(float fCellSize);

	// last-touched cache to speed up clears, etc
	std::vector<ExpMapParticle *> m_vLastParticles;
	void InitializeParticles( ExpMapParticle * pSkip );


	// Particle management
public:
	virtual unsigned int AddParticle( ExpMapParticle * pParticle ) {
		m_vParticles.push_back( pParticle );
		return (unsigned int)(m_vParticles.size()-1);
	}

	virtual ExpMapParticle * AllocateParticle() { 
		ExpMapParticle * pParticle = m_vParticleMemoryPool.Allocate();
		pParticle->Reset();
		return pParticle;
	}

	virtual void ClearParticles() { 
		m_vParticles.clear();
		m_vParticleMemoryPool.ClearAll();
	}

	class Iterator {
	public:
		Iterator( DynamicVector<ExpMapParticle *> * pVector, unsigned int nCur ) 
			{ m_pVector = pVector; m_nCount = (unsigned int)m_pVector->size(); m_nCur = nCur; }
		~Iterator() {}

		ExpMapParticle * operator *() 
			{ return (*m_pVector)[m_nCur]; }
		void operator++() 
			{ if ( m_nCur < m_nCount ) m_nCur++; }
		void operator++(int nPostfix) 
			{ if ( m_nCur < m_nCount ) m_nCur++; }
		bool operator==( const Iterator & p2 ) const 
			{ return m_nCur == p2.m_nCur; }
		bool operator!=( const Iterator & p2 ) const
			{ return m_nCur != p2.m_nCur; }

	protected:
		DynamicVector<ExpMapParticle *> * m_pVector;
		unsigned int m_nCount;
		unsigned int m_nCur;
	};

	Iterator Begin() { return Iterator( &m_vParticles, 0 ); }
	Iterator End() { return Iterator( &m_vParticles, (unsigned int)m_vParticles.size() ); }

protected:
	DynamicVector<ExpMapParticle *> m_vParticles;
	MemoryPool<ExpMapParticle> m_vParticleMemoryPool;


};


} // end namespace rms




#endif // __RMS_EXPMAP_GENERATOR_H
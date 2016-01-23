// Ryan Schmidt
// Copyright (c) 2006-2010
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt

//#include <stdafx.h>
#include "ExpMapGenerator.h"
//#include "rmsprofile.h"

#include <limits>
#include <list>
#include <set>
#include <algorithm>

#include "VectorUtil.h"
#include "MeshUtils.h"
#include "VFTriangleMesh.h"

using namespace rms;


Wml::Vector2f ExpMapParticle::INVALID_PARAM = Wml::Vector2f( std::numeric_limits<float>::max(), std::numeric_limits<float>::max() );

ExpMapParticle::ExpMapParticle(const Wml::Vector3f & vPosition,
							   const Wml::Vector3f & vNormal)
							   : m_vPosition(vPosition), m_vNormal(vNormal)
{
	m_fSurfaceDistance = 0.0f;
	m_vSurfaceVector = Wml::Vector2f::ZERO;
	ClearNeighbourList();
	m_nVertexID = std::numeric_limits<unsigned int>::max();

	m_eState = ExpMapParticle::Inactive;
	m_pNearest = NULL;
	m_bNbrFlag = false;
	m_bNormalSmoothed = false;
}

ExpMapParticle::~ExpMapParticle()
{
	ClearNeighbourList();
}

void ExpMapParticle::AddNeighbour( ListEntry * pNeighbour )
{
	if ( m_pNeighbourList == NULL ) {
		m_pNeighbourList = pNeighbour;
		pNeighbour->pNext = NULL;
	} else {
		pNeighbour->pNext = m_pNeighbourList;
		m_pNeighbourList = pNeighbour;
	}
}

void ExpMapParticle::ClearNeighbourList()
{
	// cleared using memory pool in particle system!
	m_pNeighbourList = NULL;
}




/*
 * ExpMapGenerator
 */


ExpMapGenerator::ExpMapGenerator()
{
	m_bNeighbourListsValid = false;
	m_bParticleGridValid = false;
	m_fMaxNbrDist = 0.0f;
	m_pSeedParticle = NULL;

	m_pMesh = NULL;
	m_pMeshBVTree = NULL;

	m_bUseMeshNeighbours = false;
	m_bPrecomputeNeighbours = false;
	if ( ! m_bUseMeshNeighbours )
		m_bPrecomputeNeighbours = true;
	m_bUseNeighbourNormalSmoothing = false;
	m_bUseUpwindAveraging = false;
	m_bUseSeedFrameAtParticles = false;
	m_bPreserveProjectedLengths = true;

}

ExpMapGenerator::~ExpMapGenerator()
{
}

void ExpMapGenerator::SetSurface( IMesh * pMesh, IMeshBVTree * pMeshBVTree )
{
	m_pMesh = pMesh;
	m_pMeshBVTree = pMeshBVTree;

	float fMin = 0;
	MeshUtils::GetEdgeLengthStats(pMesh, fMin, m_fMaxEdgeLength, m_fAvgEdgeLength);
	if ( m_fMaxEdgeLength == 0 )
	{
//		DebugBreak();
	}

	Reset();

	unsigned int nMaxVert = pMesh->GetVertexCount();
	m_vVtxToParticleMap.resize(0);
	m_vVtxToParticleMap.resize( nMaxVert, IMesh::InvalidID );

	IMesh::IVtxIterator curv(pMesh->BeginIVertices()), endv(pMesh->EndIVertices());
	while ( curv != endv ) {
		IMesh::VertexID nID = *curv;	curv++;

		// create particle
		ExpMapParticle * p = AllocateParticle();
		pMesh->GetVertex( nID, p->Position(), & p->Normal() );

		// associate w/ vertex
		p->VertexID() = nID;
		m_vVtxToParticleMap[nID] = AddParticle(p);

		// set arbitrary surface frame
		p->WorldFrame() = Frame3f( p->Position(), p->Normal() );
	}

	if ( m_bUseNeighbourNormalSmoothing )
		SetUseNeighbourNormalSmoothing(true);
}


void ExpMapGenerator::Reset()
{
	m_bParticleGridValid = false;
	ClearNeighbourLists();
	ClearParticles();
	m_vLastParticles.resize(0);
}


void ExpMapGenerator::SetUseNeighbourNormalSmoothing( bool bEnable )
{
	m_bUseNeighbourNormalSmoothing = bEnable;
}

void ExpMapGenerator::SetSmoothNormal( ExpMapParticle * pParticle, ExpMapParticle::ListEntry * pNbrs, bool bEnable )
{
	if ( bEnable && pParticle->m_bNormalSmoothed )
		return;
	else if ( !bEnable && !pParticle->m_bNormalSmoothed)
		return;

	Wml::Vector3f vPos, vNorm, vNbrPos, vNbrNorm;
	if ( bEnable ) {
		float fWeightSum = 0.0f;
		Wml::Vector3f vAverage = Wml::Vector3f::ZERO;
		m_pMesh->GetVertex( pParticle->VertexID(), vPos, &vNorm);
		ExpMapParticle::ListEntry * pCur = pNbrs;
		while ( pCur != NULL ) {
			IMesh::VertexID nID = pCur->pParticle->VertexID();
			Wml::Vector3f vNbrPos, vNbrNorm;
			m_pMesh->GetVertex( nID, vNbrPos, &vNbrNorm);
			float fWeight = 1.0f / ( (vPos - vNbrPos).Length()  + (0.0001f*m_fMaxEdgeLength) );
			vAverage += fWeight * vNbrNorm;
			fWeightSum += fWeight;
			pCur = pCur->pNext;
		}
		vAverage /= fWeightSum;
		vAverage.Normalize();
		pParticle->Normal() = vAverage;
	} else {
		m_pMesh->GetNormal( pParticle->VertexID(), pParticle->Normal() );
	}

	pParticle->WorldFrame() = Frame3f( pParticle->Position(), pParticle->Normal() );
	pParticle->m_bNormalSmoothed = bEnable;
}




void ExpMapGenerator::CopyVertexUVs( IMesh * pMesh, IMesh::UVSetID nSetID )
{
	pMesh->ClearUVSet( nSetID );

	Wml::Vector2f vTexCoord;
	size_t nLastCount = m_vLastParticles.size();
	for ( unsigned int i = 0; i < nLastCount; ++i ) {
		ExpMapParticle * pParticle = m_vLastParticles[i];
		if ( pParticle->VertexID() != IMesh::InvalidID )
			pMesh->SetUV( pParticle->VertexID(), nSetID, pParticle->SurfaceVector() );
	}

}








// neighbour list setup code


ExpMapParticle::ListEntry * ExpMapGenerator::GetNeighbourList( ExpMapParticle * pParticle )
{
	if ( pParticle->GetNeighbourList() == NULL ) {
		m_vNeighbourBuf.resize(0);
		if ( m_bUseMeshNeighbours ) 
			FindNeighbours(pParticle, m_vNeighbourBuf);
		else
			FindNeighbours( pParticle->Position(), pParticle->Normal(), m_vNeighbourBuf, m_fMaxNbrDist, pParticle );	

		size_t nCount = m_vNeighbourBuf.size();
		for ( unsigned int i = 0; i < nCount; ++i ) {
			ExpMapParticle::ListEntry * pEntry = m_NeighbourMemoryPool.Allocate();
			pEntry->pParticle = m_vNeighbourBuf[i];
			pParticle->AddNeighbour(pEntry);
		}
	}
	
	if ( pParticle != m_pSeedParticle )
		SetSmoothNormal(pParticle, pParticle->GetNeighbourList(), m_bUseNeighbourNormalSmoothing);	

	return pParticle->GetNeighbourList();
}


class  NeighborTriBuffer : public IMesh::NeighborTriCallback
{
public:
	NeighborTriBuffer() {}

	const std::vector<IMesh::TriangleID> & Triangles() { return m_vTriangles; }

	virtual void BeginTriangles() 
		{ m_vTriangles.resize(0); }
	virtual void NextTriangle( IMesh::TriangleID tID )
		{ m_vTriangles.push_back(tID); }

protected:
	std::vector<IMesh::TriangleID> m_vTriangles;
};


void ExpMapGenerator::FindNeighbours( ExpMapParticle * pParticle,
									  std::vector<ExpMapParticle *> & vParticles)
{
	vParticles.resize(0);

	IMesh::VertexID vID = pParticle->VertexID();

	NeighborTriBuffer vBuffer;
	m_pMesh->NeighbourIteration(vID, &vBuffer);
	const std::vector<IMesh::TriangleID> & vTris = vBuffer.Triangles();
	size_t nCount = vTris.size();
	for ( unsigned int i = 0; i < nCount; ++i ) {
		IMesh::VertexID nTri[3];
		m_pMesh->GetTriangle( vTris[i], nTri);

		// check each vertex
		for ( int j = 0; j < 3; ++j ) {
			if ( nTri[j] == vID )
				continue;
			ExpMapParticle * pParticle = m_vParticles[ m_vVtxToParticleMap[nTri[j]] ];
			bool bFound = false;
			for ( unsigned int k = 0; k < vParticles.size() && ! bFound; ++k ) {
				if ( vParticles[k] == pParticle )
					bFound = true;
			}
			if ( ! bFound )
				vParticles.push_back( pParticle );
		}
	}
}




#define USE_KNN

bool ParticleDistanceSort( ExpMapParticle * p1, ExpMapParticle * p2 )
{
	return p1->SurfaceDistance() < p2->SurfaceDistance();
}


void ExpMapGenerator::FindNeighbours( const Wml::Vector3f & vPoint, const Wml::Vector3f & vNormal,
									  std::vector<ExpMapParticle *> & vParticles, float fRadiusThreshold,
									  ExpMapParticle * pSkip )
{
	//float fDistSqrThreshold = fRadiusThreshold*fRadiusThreshold;

	std::vector<ExpMapParticle *> vNeighbours;
	const unsigned int sMaxNbrs = 15;
//	const unsigned int sMaxNbrs = 8;

rinse_and_repeat:

	vNeighbours.resize(0);
	ParticleGrid<ExpMapParticle *>::BoxIterator itr( &m_particleGrid, vPoint, fRadiusThreshold );
	if ( itr.Done() ) {
		fRadiusThreshold *= 1.5f;
		goto rinse_and_repeat;
	}
	while ( ! itr.Done() ) {
		ExpMapParticle * pTmp = *itr;
		++itr;

		if ( pTmp == pSkip )
			continue;

		float fDistSqr = ( pTmp->Position() - vPoint ).SquaredLength();
		//float fNormalDir = pTmp->Normal().Dot( vNormal );
		//if ( fDistSqr < fDistSqrThreshold && fNormalDir > -0.1f) {
		//if ( fDistSqr < fDistSqrThreshold ) {
		if ( true ) {
			pTmp->SurfaceDistance() = fDistSqr;
			vNeighbours.push_back( pTmp );
		} 
	}

	size_t nFoundNbrs = vNeighbours.size();
	if ( nFoundNbrs < sMaxNbrs/2 ) {
		for ( unsigned int k = 0; k < nFoundNbrs; ++k )
			vNeighbours[k]->SurfaceDistance() = std::numeric_limits<float>::max();
		fRadiusThreshold *= 1.5f;
		goto rinse_and_repeat;
	}

	std::sort( vNeighbours.begin(), vNeighbours.end(), ParticleDistanceSort );
	unsigned int nMaxNbrs = (unsigned int)std::min(sMaxNbrs, (unsigned int)nFoundNbrs );
	for ( unsigned int k = 0; k < nFoundNbrs; ++k )
		vNeighbours[k]->SurfaceDistance() = std::numeric_limits<float>::max();

	vParticles.resize(nMaxNbrs);
	for ( unsigned int i = 0; i < nMaxNbrs; ++i ) {
		vParticles[i] = vNeighbours[i];
	}
}


void ExpMapGenerator::InitializeNeighbourLists( float fRadiusThreshold )
{
	if (! m_bUseMeshNeighbours ) {
		// find max nbr dist
		m_fMaxNbrDist = m_fAvgEdgeLength;
		InitializeParticleGrid(m_fMaxNbrDist);
	}

	if ( m_bPrecomputeNeighbours && ! m_bNeighbourListsValid ) {
		// force nbr generation for all particles
		size_t nParticleCount = m_vParticles.size();	
		for ( unsigned int i = 0; i < nParticleCount; ++i ) 
			GetNeighbourList(m_vParticles[i]);

		m_bNeighbourListsValid = true;
	}
}



void ExpMapGenerator::ClearNeighbourLists()
{
	size_t nParticleCount = m_vParticles.size();
	for ( unsigned int i = 0; i < nParticleCount; ++i ) {
		m_vParticles[i]->ClearNeighbourList();
	}
	m_NeighbourMemoryPool.ClearAll();

	m_bNeighbourListsValid = false;
}









// particle grid setup code

void ExpMapGenerator::InitializeParticleGrid(float fCellSize)
{
	if (m_bParticleGridValid == true )
		return;

	Wml::AxisAlignedBox3f partBounds;
	size_t nParticleCount = m_vParticles.size();	
	for ( unsigned int i = 0; i < nParticleCount; ++i ) {
		ExpMapParticle * pCur = m_vParticles[i];
		if ( i == 0 )
			partBounds = Wml::AxisAlignedBox3f( 
			pCur->Position().X(), pCur->Position().X(),
			pCur->Position().Y(), pCur->Position().Y(),
			pCur->Position().Z(), pCur->Position().Z() );
		else
			rms::Union( partBounds, pCur->Position() );
	}
	// dilate by one cell
	for ( int k = 0; k < 3; ++k ) {
		partBounds.Min[k] -= fCellSize;
		partBounds.Max[k] += fCellSize;
	}

	m_particleGrid.Initialize( rms::Center(partBounds), fCellSize );

	// add each particle
	for ( unsigned int i = 0; i < nParticleCount; ++i ) {
		ExpMapParticle * pCur = m_vParticles[i];

		m_particleGrid.AddParticle( pCur, pCur->Position() );
	}

	m_bParticleGridValid = true;
}	


void ExpMapGenerator::InitializeParticles( ExpMapParticle * pSkip )
{
	if ( m_vLastParticles.empty() ) {
		size_t nParticleCount = m_vParticles.size();
		for ( unsigned int i = 0; i < nParticleCount; ++i ) {
			if ( m_vParticles[i] == pSkip )
				continue;
			m_vParticles[i]->Clear();
		}
	} else {
		size_t nCount = m_vLastParticles.size();
		for ( unsigned int i = 0; i < nCount; ++i ) {
			if ( m_vLastParticles[i] == pSkip )
				continue;
			m_vLastParticles[i]->Clear();
		}
		m_vLastParticles.resize(0);
	}
}


















void ExpMapGenerator::SetSurfaceDistances( const Wml::Vector3f & vSeedPoint, float fNeighbourThreshold,
													   float fStopDistance,
													   const Frame3f * pLastFrame,
													   const Wml::Vector3f * pSeedNormal )
{
	//std::cout << "SetSurfaceDistances\n";
	//_RMSTUNE_start(4);

	// create seed particle
	InitializeNeighbourLists(fNeighbourThreshold);
	InitializeSeedParticle( vSeedPoint, pSeedNormal, pLastFrame );

	// compute approximate geodesic distances to seed particle
	ComputeExpMap( fStopDistance );

	//_RMSTUNE_end(4);
	//_RMSInfo("Total time was %f\n", _RMSTUNE_time(4) );
}




void ExpMapGenerator::ComputeExpMap( float fStopDistance )
{
	// set all particle distances to max and state to inactive
	size_t nParticleCount = m_vParticles.size();
	InitializeParticles( m_pSeedParticle );

	// now initialize pq
	std::multiset< ParticleQueueWrapper > pq;
	UpdateNeighbours( m_pSeedParticle, pq );

	// run until pq is empty
	float fStopDistSquare = fStopDistance / (float)sqrt(2.0f);
	unsigned int nTouched = 0;
	while ( ! pq.empty() ) {

		// pop front	
		//ExpMapParticle * pFront = pq.begin()->Particle();
		ExpMapParticle * pFront = (const_cast<ParticleQueueWrapper&>(*pq.begin())).Particle();
		pq.erase( pq.begin() );

		// freeze particle
		pFront->State() = ExpMapParticle::Frozen;
		m_vLastParticles.push_back( pFront );

		// set frame for pFront
		if ( m_bUseUpwindAveraging ) 
			PropagateFrameFromNearest_Average( pFront );
		else
			PropagateFrameFromNearest( pFront );

		if ( pFront->SurfaceDistance() > fStopDistance )
			continue;

		// Square-culling. This gives a significant speed-up...
		if ( (float)abs(pFront->SurfaceVector().X()) > fStopDistSquare || 
			 (float)abs(pFront->SurfaceVector().Y()) > fStopDistSquare )
			 continue;

		// update neighbours
		RemoveNeighbours( pFront, pq );
		UpdateNeighbours( pFront, pq );
		++nTouched;
	}
	//_RMSInfo("Touched %d particles while updating\n", nTouched);

	// mark any left-over particles for cleanup
	std::multiset<ParticleQueueWrapper>::iterator cur(pq.begin()), end(pq.end());
	while ( cur != end ) {
		//m_vLastParticles.push_back( (*cur).Particle() );
		m_vLastParticles.push_back( (const_cast<ParticleQueueWrapper&>(*pq.begin())).Particle() );
		++cur;
	}
}





void ExpMapGenerator::RemoveNeighbours( ExpMapParticle * pParticle, std::multiset< ParticleQueueWrapper > & pq )
{
	ExpMapParticle::ListEntry * pCur = GetNeighbourList( pParticle );
	//if ( pCur == NULL ) DebugBreak();
	while ( pCur != NULL ) {

		ExpMapParticle * pCurParticle = pCur->pParticle;
		pCur = pCur->pNext;

		if ( pCurParticle->State() != ExpMapParticle::Active )
			continue;

		// find entry in pq
#if 1
		std::multiset<ParticleQueueWrapper>::iterator found( 
			pq.find( ParticleQueueWrapper( pCurParticle->SurfaceDistance() ) ) );
		if ( found != pq.end() ) {

// 			while ( (*found).Particle() != pCurParticle &&
// 						(*found).QueueValue() == pCurParticle->SurfaceDistance() )
			while( (const_cast<ParticleQueueWrapper&>(*found)).Particle() != pCurParticle
				&& (const_cast<ParticleQueueWrapper&>(*found)).QueueValue()==pCurParticle->SurfaceDistance() )
				++found;

			// [RMS: this should always happen...]
			//ASSERT( (*found).Particle() == pCurParticle );
			ASSERT((const_cast<ParticleQueueWrapper&>(*found)).Particle() == pCurParticle);
			if ( (const_cast<ParticleQueueWrapper&>(*found)).Particle() == pCurParticle )
			//if ( (*found).Particle() == pCurParticle )
				pq.erase( found );
		} else {
			ASSERT( found != pq.end() );
		}

#else 
		std::multiset<ParticleQueueWrapper>::iterator cur( pq.begin() );
		bool bFound = false;
		while ( !bFound && cur != pq.end() ) {
			if ( (*cur).Particle() == pCurParticle ) {
				pq.erase( cur );
				bFound = true;
			} else
				++cur;
		}
		ASSERT( bFound );
#endif

	}
}


void ExpMapGenerator::UpdateNeighbours( ExpMapParticle * pParticle, std::multiset< ParticleQueueWrapper > & pq )
{
	// iterate through neighbours, updating particle distances and pushing onto pq
	ExpMapParticle::ListEntry * pCur = GetNeighbourList( pParticle );
	//if ( pCur == NULL ) DebugBreak();
	while ( pCur != NULL ) {

		ExpMapParticle * pCurParticle = pCur->pParticle;
		pCur = pCur->pNext;

		// skip inactive particles
		if ( pCurParticle->State() == ExpMapParticle::Frozen )
			continue;

		// set active state
		pCurParticle->State() = ExpMapParticle::Active;

		// compute new distance
		float fDistToPoint = (pParticle->Position() - pCurParticle->Position()).Length();
		float fSurfDist = fDistToPoint + pParticle->SurfaceDistance();

		// update particle distance and/or nearest particle
		bool bUpdated = false;
		if ( fSurfDist < pCurParticle->SurfaceDistance() ) {
			pCurParticle->SetNearestParticle( pParticle );
			pCurParticle->SurfaceDistance() = fSurfDist;
			bUpdated = true;
		}

		// re-insert particle into priority queue
		pq.insert( ParticleQueueWrapper(pCurParticle) );
	}	

}


ExpMapParticle * ExpMapGenerator::InitializeSeedParticle( const Wml::Vector3f & vPosition,
														  const Wml::Vector3f * pSeedNormal,
														  const Frame3f * pLastSeedPointFrame )
{
	if ( m_pSeedParticle == NULL )
		m_pSeedParticle = new ExpMapParticle();

	ExpMapParticle * pSeedParticle = m_pSeedParticle;
	pSeedParticle->Reset();
	pSeedParticle->Position() = vPosition;
	pSeedParticle->Normal() = *pSeedNormal;
	pSeedParticle->State() = ExpMapParticle::Frozen;


	if ( m_bUseMeshNeighbours ) {
		// just use 3 nearest mesh neighbours...

		Wml::Vector3f vNearest;
		IMesh::TriangleID tID;
		if ( ! m_pMeshBVTree->FindNearest( vPosition, vNearest, tID ) )
		{
			//DebugBreak();
		}
		IMesh::VertexID nTri[3], nNbrTri[3];
		m_pMesh->GetTriangle( tID, nTri );

		// add direct nbrs
		m_vNeighbourBuf.resize(0);
		for ( int j = 0; j < 3; ++j ) {
			m_vNeighbourBuf.push_back( m_vParticles[ m_vVtxToParticleMap[nTri[j]] ] );
			m_vNeighbourBuf.back()->m_bNbrFlag = true;
		}

		// add each of their one-rings
		NeighborTriBuffer vBuffer;
		for ( int j = 0; j < 3; ++j ) {
			m_pMesh->NeighbourIteration( nTri[j], &vBuffer );
			const std::vector<IMesh::TriangleID> & vTriangles = vBuffer.Triangles();
			size_t nCount = vTriangles.size();
			for ( unsigned int i = 0; i < nCount; ++i ) {
				m_pMesh->GetTriangle( vTriangles[i], nNbrTri );
				for ( int k = 0; k < 3; ++k ) {
					ExpMapParticle * pParticle = m_vParticles[ m_vVtxToParticleMap[nNbrTri[j]] ];
					if ( pParticle->m_bNbrFlag == false ) {
						m_vNeighbourBuf.push_back(pParticle);
						m_vNeighbourBuf.back()->m_bNbrFlag = true;
					}
				}
			}
		}

	} else {
		// multiplying distance by 2 here is a hack, to try and fix
		// some problems where we don't get enough neighbours around the seed point 
		// if it is in the middle of a triangle or something like that...
		m_vNeighbourBuf.resize(0);
		FindNeighbours( pSeedParticle->Position(), pSeedParticle->Normal(), m_vNeighbourBuf, m_fMaxNbrDist*2.0f );	
	}

	//// just use 3 nearest mesh neighbours...
	//Wml::Vector3f vNearest;
	//IMesh::TriangleID tID;
	//if ( ! m_pMeshBVTree->FindNearest( vPosition, vNearest, tID ) )
	//	DebugBreak();
	//IMesh::VertexID nTri[3], nNbrTri[3];
	//m_pMesh->GetTriangle( tID, nTri );
	//
	//m_vNeighbourBuf.resize(0);

	//// add direct nbrs
	//for ( int j = 0; j < 3; ++j ) {
	//	m_vNeighbourBuf.push_back( m_vParticles[ m_vVtxToParticleMap[nTri[j]] ] );
	//	m_vNeighbourBuf.back()->m_bNbrFlag = true;
	//}

	//// add each of their one-rings
	//NeighborTriBuffer vBuffer;
	//for ( int j = 0; j < 3; ++j ) {
	//	m_pMesh->NeighbourIteration( nTri[j], &vBuffer );
	//	const std::vector<IMesh::TriangleID> & vTriangles = vBuffer.Triangles();
	//	size_t nCount = vTriangles.size();
	//	for ( unsigned int i = 0; i < nCount; ++i ) {
	//		m_pMesh->GetTriangle( vTriangles[i], nNbrTri );
	//		for ( int k = 0; k < 3; ++k ) {
	//			ExpMapParticle * pParticle = m_vParticles[ m_vVtxToParticleMap[nNbrTri[j]] ];
	//			if ( pParticle->m_bNbrFlag == false ) {
	//				m_vNeighbourBuf.push_back(pParticle);
	//				m_vNeighbourBuf.back()->m_bNbrFlag = true;
	//			}
	//		}
	//	}
	//}

	// ok now add them all as neighbours, and clear neighbour flags
	size_t nCount = m_vNeighbourBuf.size();
	for ( unsigned int i = 0; i < nCount; ++i ) {
		ExpMapParticle::ListEntry * pEntry = m_NeighbourMemoryPool.Allocate();
		pEntry->pParticle = m_vNeighbourBuf[i];
		pSeedParticle->AddNeighbour( pEntry );
		m_vNeighbourBuf[i]->m_bNbrFlag = false;
	}

	// estimate smooth normal  (only works if IMesh is VFTriangleMesh...)
	VFTriangleMesh * pVFMesh = dynamic_cast<VFTriangleMesh *>(m_pMesh);
	if ( pVFMesh && m_bUseNeighbourNormalSmoothing ) {
		// get average normal...
		std::vector<IMesh::VertexID> vNbrs;
		for ( unsigned int k = 0; k < m_vNeighbourBuf.size(); ++k )
			vNbrs.push_back(m_vNeighbourBuf[k]->VertexID());

		Wml::Vector3f vNewNormal(0,0,0);
		float fWeightSum = 0.0f;
		for ( unsigned int i = 0; i < m_vNeighbourBuf.size(); ++i) {
			ExpMapParticle * pTriParticle = m_vNeighbourBuf[i];
			float fWeight = 1.0f / ((pTriParticle->Position() - vPosition).Length() + (0.0001f*m_fMaxEdgeLength) );
			fWeightSum += fWeight;
			vNewNormal += fWeight * MeshUtils::GetAverageNormal( *pVFMesh, pTriParticle->VertexID() );
		}
		vNewNormal /= fWeightSum;
		vNewNormal.Normalize();
		pSeedParticle->Normal() = vNewNormal;
	}


	// compute 3D frame at seed point
	Frame3f startWorldFrame( pSeedParticle->Position(), pSeedParticle->Normal() );

	// if this argument was passed non-null, try to minimize the tangent-frame rotation of the
	// current seed point wrt the last one
	if ( pLastSeedPointFrame != NULL ) {
		Frame3f vLastFrame( *pLastSeedPointFrame );
		vLastFrame.AlignZAxis( startWorldFrame );
		vLastFrame.Origin() = startWorldFrame.Origin();
		startWorldFrame = vLastFrame;
	}
	pSeedParticle->WorldFrame() = startWorldFrame;

	// [TEST]
	//for ( unsigned int k = 0; k < m_vParticles.size(); ++k ) {
	//	m_vParticles[k]->Normal() = pSeedParticle->Normal();
	//	m_vParticles[k]->WorldFrame() = Frame3f( m_vParticles[k]->Position(), m_vParticles[k]->Normal() );
	//}
	// [TEST]

	return pSeedParticle;
}







void ExpMapGenerator::PropagateFrameFromNearest( ExpMapParticle * pParticle )
{
	ExpMapParticle * pCenterParticle = pParticle->NearestParticle();

	Wml::ExtPlane3f vTangentPlane;
	Frame3f vCenterWorldFrame;
	Wml::Matrix2f matFrameRotate;
	PrecomputePropagationData( pCenterParticle, vTangentPlane, vCenterWorldFrame, matFrameRotate );

	pParticle->SurfaceVector() = ComputeSurfaceVector( pCenterParticle, pParticle, 
		vTangentPlane, vCenterWorldFrame, matFrameRotate );

	// check error and re-set distance if it is too high
#if 1
	float fVecLengthSqr = pParticle->SurfaceVector().SquaredLength();
	float fError = fabs(
		fVecLengthSqr / (pParticle->SurfaceDistance()*pParticle->SurfaceDistance()) - 1.0f );
	static const float s_fMaxErr = 2*0.5f - 0.5f*0.5f;		// 0.5f == 50% - arbitrary threshold...
	if ( fError > s_fMaxErr ) {
		if ( true ) {
			pParticle->SurfaceVector().Normalize();
			pParticle->SurfaceVector() *= pParticle->SurfaceDistance();
			//_RMSInfo("Fix!\n");
		}
	}
#endif
}




struct NbrInfo
{
	float fNbrWeight;
	Wml::Vector2f vNbrUV;
};

void ExpMapGenerator::PropagateFrameFromNearest_Average( ExpMapParticle * pParticle )
{
	ExpMapParticle * pNearest = pParticle->NearestParticle();

	if ( pParticle->SurfaceDistance() == 0.0f )  {
		// pathological case where seed point == input point (or other duplicate points)
		pParticle->SurfaceVector() = pNearest->SurfaceVector();
		return;
	}

	static std::vector<NbrInfo> vNbrs;
	Wml::ExtPlane3f vTangentPlane;
	Frame3f vCenterWorldFrame;
	Wml::Matrix2f matFrameRotate;

	float fWeightSum = 0.0f;
	vNbrs.resize(0);

	//! need to look for nearest particle, because if it is seed
	//! point, it's not in nbr lists...
	if ( pNearest->State() != ExpMapParticle::Frozen )
	{
	//	DebugBreak();
	}

	bool bSawNearest = false;
	ExpMapParticle::ListEntry * pCurNbr = GetNeighbourList(pParticle);
	while ( pCurNbr != NULL ) {

		ExpMapParticle * pCenterParticle = pCurNbr->pParticle;
		if ( pCenterParticle == pNearest )
			bSawNearest = true;

		if ( pCenterParticle->State() == ExpMapParticle::Frozen ) {
			PrecomputePropagationData( pCenterParticle, vTangentPlane, vCenterWorldFrame, matFrameRotate );

			NbrInfo info;
			info.vNbrUV = ComputeSurfaceVector( pCenterParticle, pParticle, 
				vTangentPlane, vCenterWorldFrame, matFrameRotate );

			bool bSkip = false;
			if ( ! _finite(info.vNbrUV.Length()) )
				bSkip = true;

			info.fNbrWeight = 1.0f / ( ( pCenterParticle->Position() - pParticle->Position() ).Length() + (0.00001f*m_fMaxEdgeLength) );
			if ( ! _finite(info.fNbrWeight) )
				bSkip = true;

			// weight by geodesic delta - tries to avoid bias from neigbours w/ same "time-of-arrival" (ie at same front-timestep)
			// [RMS: not sure if this has a significant effect - note that need to enable below as well....
			// [RMS] also this is broken for particles at the same 3D position, which seems to happen w/ surface tree...
			//float fGeoDelta = fabs(pCenterParticle->SurfaceVector().Length() - info.vNbrUV.Length());
			//info.fNbrWeight *= fGeoDelta;

			if ( !bSkip ) {
				fWeightSum += info.fNbrWeight;
				vNbrs.push_back(info);
			}
		}

		pCurNbr = pCurNbr->pNext;
	}

	if ( ! _finite(fWeightSum) )
	{
		//DebugBreak();
	}

	// add un-seen "nearest" particle
	if ( ! bSawNearest ) {
		PrecomputePropagationData( pNearest, vTangentPlane, vCenterWorldFrame, matFrameRotate );
		NbrInfo info;
		info.vNbrUV =  ComputeSurfaceVector( pNearest, pParticle, 
				vTangentPlane, vCenterWorldFrame, matFrameRotate );

		info.fNbrWeight = 1.0f / ( ( pNearest->Position() - pParticle->Position() ).Length() + (0.00001f*m_fMaxEdgeLength) );

		// weight by geo delta
		//float fGeoDelta = fabs(pNearest->SurfaceVector().Length() - info.vNbrUV.Length());
		//info.fNbrWeight *= fGeoDelta;

		fWeightSum += info.fNbrWeight;
		vNbrs.push_back(info);
	}

	if ( ! _finite(fWeightSum) )
	{
		//DebugBreak();
	}

	// take weighted average to get real UV
	Wml::Vector2f vUV = Wml::Vector2f::ZERO;
	size_t nCount = vNbrs.size();
	for ( unsigned int i = 0; i < nCount; ++i ) {
		float fWeight = vNbrs[i].fNbrWeight / fWeightSum;
		vUV += fWeight * vNbrs[i].vNbrUV;
	}

	if ( ! _finite(vUV.Length()) )
	{
		//DebugBreak();
	}

	pParticle->SurfaceVector() = vUV;


	// check error and re-set distance if it is too high
#if 0
	float fVecLengthSqr = pParticle->SurfaceVector().SquaredLength();
	float fError = fabs(
		fVecLengthSqr / (pParticle->SurfaceDistance()*pParticle->SurfaceDistance()) - 1.0f );
	static const float s_fMaxErr = 2*0.5f - 0.5f*0.5f;		// 0.5f == 50% - arbitrary threshold...
	if ( fError > s_fMaxErr ) {
		if ( true ) {
			pParticle->SurfaceVector().Normalize();
			pParticle->SurfaceVector() *= pParticle->SurfaceDistance();
//			_RMSInfo("Fix!\n");
		}
	}
#endif

}




void ExpMapGenerator::PrecomputePropagationData( ExpMapParticle * pCenterParticle,
															 Wml::ExtPlane3f & vTangentPlane, 
															 Frame3f & vCenterWorldFrame,
															 Wml::Matrix2f & matFrameRotate )
{
	// compute surface-tangent plane at particle
	Wml::Vector3f vNormal( m_bUseSeedFrameAtParticles ?
		 m_pSeedParticle->Normal() : pCenterParticle->Normal() );

	// compute plane at this point
	vTangentPlane = Wml::ExtPlane3f( vNormal, pCenterParticle->Position() );

	// compute 3D frame at center point
	//vCenterWorldFrame = Frame3f( pCenterParticle->Position(), vNormal );
	vCenterWorldFrame = m_bUseSeedFrameAtParticles ?
		m_pSeedParticle->WorldFrame() : pCenterParticle->WorldFrame();

	// rotate seed world frame Z into this particle's Z
	Frame3f seedWorldFrame( m_pSeedParticle->WorldFrame() );
	seedWorldFrame.AlignZAxis( vCenterWorldFrame );

	// compute cos(angle) between axis
	Wml::Vector3f vCenterAxisX( vCenterWorldFrame.Axis( Frame3f::AxisX ) );
	Wml::Vector3f vSeedAxisX( seedWorldFrame.Axis( Frame3f::AxisX ) );
	float fCosTheta = vCenterAxisX.Dot(vSeedAxisX);

	// compute rotated min-dist vector for this particle
	float fTmp = 1 - fCosTheta*fCosTheta;
	if ( fTmp < 0 ) fTmp = 0;		// need to clamp so that sqrt works...
	float fSinTheta = (float)sqrt(fTmp);
	Wml::Vector3f vCross = vCenterAxisX.Cross(vSeedAxisX);
	if ( vCross.Dot( vNormal ) < 0 )	// get the right sign...
		fSinTheta = -fSinTheta;

	// create transposed 2D rotation matrix
	matFrameRotate = Wml::Matrix2f( fCosTheta, fSinTheta, -fSinTheta, fCosTheta );
}

Wml::Vector2f ExpMapGenerator::ComputeSurfaceVector( ExpMapParticle * pCenterParticle, 
																 ExpMapParticle * pNbrParticle,
																 Wml::ExtPlane3f & vTangentPlane, 
																 Frame3f & vCenterWorldFrame, 
																 Wml::Matrix2f & matFrameRotate )
{
	// project point into plane
	Wml::Vector3f vPlanePoint = m_bPreserveProjectedLengths ? 
		vTangentPlane.RotatePointIntoPlane( pNbrParticle->Position() ) : 
		vTangentPlane.ProjectPointToPlane( pNbrParticle->Position() );

	// project point into coord system of frame
	vPlanePoint -= pCenterParticle->Position();
	vCenterWorldFrame.ToFrameLocal(vPlanePoint);

	// now we can project into surface frame simply by dropping z (which should be 0 anyway,
	//   since the vector lies in the plane!)
	Wml::Vector2f vSurfaceFrame( vPlanePoint.X(), vPlanePoint.Y() );

	// reverse vector so it points back to current particle
	vSurfaceFrame *= -1.0f;

	// transform local vector into coord system of initial surface reference frame
	//  and add accumulated surface vector
	return pCenterParticle->SurfaceVector() + (matFrameRotate * vSurfaceFrame);
}

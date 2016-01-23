// Ryan Schmidt
// Copyright (c) 2006-2010
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt

//#include "StdAfx.h"
#include "VFTriangleMesh.h"
#include "VectorUtil.h"

using namespace rms;

VFTriangleMesh::VFTriangleMesh(void)
{
}

VFTriangleMesh::~VFTriangleMesh(void)
{
}

VFTriangleMesh::VFTriangleMesh( IMesh & copy, std::vector<VertexID> * vVertMap )
{
	Copy(copy, vVertMap);
}



IMesh::TriangleID VFTriangleMesh::AppendTriangle( IMesh::VertexID v1, IMesh::VertexID v2, IMesh::VertexID v3 )
{ 
	ASSERT( m_vVertices.isValid(v1) && m_vVertices.isValid(v2) && m_vVertices.isValid(v3) );

	// insert new triangle
	TriangleID tID = (TriangleID)m_vTriangles.insert( Triangle(v1,v2,v3) ); 

	// increment reference counts
	m_vVertices.increment( v1 );
	m_vVertices.increment( v2 );
	m_vVertices.increment( v3 );

	// add to triangle lists
	AddTriEntry( tID, v1 );
	AddTriEntry( tID, v2 );
	AddTriEntry( tID, v3 );

	return tID;
}



void VFTriangleMesh::SetTriangle( IMesh::TriangleID tID, IMesh::VertexID v1, IMesh::VertexID v2, IMesh::VertexID v3 )
{
	if ( m_vTriangles.isValid(tID) ) {
		Triangle & t = m_vTriangles[tID];

		// decrement existing reference counts
		ASSERT( m_vVertices.isValid(t.nVertices[0]) && m_vVertices.isValid(t.nVertices[1]) && m_vVertices.isValid(t.nVertices[2]) );
		RemoveTriEntry( tID, t.nVertices[0] );
		RemoveTriEntry( tID, t.nVertices[1] );
		RemoveTriEntry( tID, t.nVertices[2] );
		m_vVertices.decrement(t.nVertices[0]);
		m_vVertices.decrement(t.nVertices[1]);
		m_vVertices.decrement(t.nVertices[2]);

		// set new IDs
		t.nVertices[0] = v1;
		t.nVertices[1] = v2;
		t.nVertices[2] = v3;

		// increment new reference counts
		ASSERT( m_vVertices.isValid(v1) && m_vVertices.isValid(v2) && m_vVertices.isValid(v3) );
		m_vVertices.increment(v1);
		m_vVertices.increment(v2);
		m_vVertices.increment(v3);
		AddTriEntry( tID, v1 );
		AddTriEntry( tID, v2 );
		AddTriEntry( tID, v3 );
	}

}


void VFTriangleMesh::Clear( bool bFreeMem )
{
	IMesh::Clear(bFreeMem);
	m_vVertices.clear( bFreeMem );
	m_vTriangles.clear( bFreeMem );
	m_VertDataMemPool.ClearAll();
	m_VertListPool.Clear(bFreeMem);
}



void VFTriangleMesh::RemoveVertex( VertexID vID )
{
	if ( ! m_vVertices.isValid( vID ) )
		return;

	VFTriangleMesh::Vertex & v = m_vVertices[vID];

	// [RMS] HACK! This case shouldn't happen, but it does in ::Weld() because
	//  SetTriangle() doesn't remove un-referenced vertices (which it really
	//   shouldn't, since we might be performing mesh surgery stuff....
	if ( v.pData->vTriangles.pFirst == NULL ) {
		if ( m_vVertices.refCount( vID ) == 1 )
			m_vVertices.remove( vID );
//		else
//			DebugBreak();
	} else { 
		// remove each attached face
		while ( m_vVertices.isValid( vID ) && v.pData->vTriangles.pFirst != NULL )
			RemoveTriangle(v.pData->vTriangles.pFirst->data );
	}

	// vertex should be gone now because no attached faces remain! 
	ASSERT( ! m_vVertices.isValid( vID ) );
}


void VFTriangleMesh::RemoveTriangle( TriangleID tID )
{
	Triangle & t = m_vTriangles[tID];

	// decrement existing reference counts
	ASSERT( m_vVertices.isValid(t.nVertices[0]) && m_vVertices.isValid(t.nVertices[1]) && m_vVertices.isValid(t.nVertices[2]) );
	for ( int i = 0; i < 3; ++i ) {
		RemoveTriEntry( tID, t.nVertices[i] );
		m_vVertices.decrement(t.nVertices[i]);

		// remove vertex if refcount == 1  (means that it is only referenced by self, so is safe to delete)
		if ( m_vVertices.refCount( t.nVertices[i] ) == 1 )
			m_vVertices.remove( t.nVertices[i]  );
	}

	// remove triangle
	m_vTriangles.remove( tID );
}

void VFTriangleMesh::RemoveUnreferencedGeometry()
{
	
}


//! initialize vertex neighbour iteration
void VFTriangleMesh::BeginVtxTriangles( VtxNbrItr & v )
{
	ASSERT( m_vVertices.isValid(v.vID) );
	VFTriangleMesh::Vertex & vtx = m_vVertices[v.vID];
	if ( vtx.pData == NULL )
		v.nData[0] = IMesh::InvalidID;
	else
		v.nData[0] = (unsigned long long)(vtx.pData->vTriangles.pFirst);
	v.nData[1] = 1234567890;
}

//! (possibly) un-ordered iteration around one-ring of a vertex. Returns InvalidID when done
IMesh::TriangleID VFTriangleMesh::GetNextVtxTriangle( VtxNbrItr & v )
{
	if ( v.nData[0] == IMesh::InvalidID )
		return IMesh::InvalidID;
	
	TriListEntry * pCur = (TriListEntry *)v.nData[0];
	IMesh::TriangleID tID = pCur->data;

	if ( pCur->pNext == NULL )
		v.nData[0] = IMesh::InvalidID;
	else
		v.nData[0] = (unsigned long long)(pCur->pNext);
	
	return tID;
}


bool VFTriangleMesh::IsBoundaryVertex( VertexID vID )
{
	// RMS TODO: it should be possible to do this just using
	//  a 3-vertex buffer (ie determine direction by looking
	//  at last and second-last vertices...)
	//  (Maybe even 2-vertex buffer?)

	Vertex & v = m_vVertices[vID];

	// count triangles and make a list of them
	int nCount = 0;
	TriListEntry * pCur = v.pData->vTriangles.pFirst;
	while ( pCur != NULL ) {
		pCur = pCur->pNext;
		nCount++;
	}
	if ( nCount == 0 )
		return true;

	std::vector< TriangleID > vTris;
	vTris.resize( nCount-1 );
	pCur = v.pData->vTriangles.pFirst->pNext;
	int i = 0;
	while ( pCur != NULL ) {
		vTris[ i++ ] = pCur->data;
		pCur = pCur->pNext;
	}

	// pick first edge
	VertexID vCurID = InvalidID;
	VertexID vStopID = InvalidID;
	pCur = v.pData->vTriangles.pFirst;
	VertexID * pTri = m_vTriangles[ pCur->data ].nVertices;
	for ( int i = 0; i < 3; ++i ) {
		if ( pTri[i] == vID ) {
			vCurID = pTri[ (i+1) % 3 ];
			vStopID = pTri[ (i+2) % 3];
			break;
		} else if ( pTri [ (i+1) % 3 ] == vID ) {
			vCurID = pTri[i];
			vStopID = pTri[ (i+2) % 3 ];
			break;
		}
	}
	nCount--;   // used up first tri

	// loop until we hit vStopID
	while ( vCurID != vStopID ) {
		
		// find unused tri w/ [vID, vCurID, X]
		int nCurTri = InvalidID;
		for ( int i = 0; i < nCount; ++i ) {
			if ( vTris[i] == InvalidID )
				continue;
			VertexID * pTri = m_vTriangles[ vTris[i] ].nVertices;
			if ( pTri[0] == vCurID || pTri[1] == vCurID || pTri[2] == vCurID ) {
				nCurTri = i;
				break;
			}
		}
		if ( nCurTri == InvalidID )
			return true;			// 1-ring is not connected - must be invalid!

		// mark tri as done
		VertexID * pTri = m_vTriangles[ vTris[nCurTri] ].nVertices;
		vTris[ nCurTri ] = InvalidID;

		// go to next tri in one-ring
		if ( pTri[0] == vID ) 
			vCurID = ( pTri[1] == vCurID ) ? pTri[2] : pTri[1];
		else if ( pTri[1] == vID )
			vCurID = ( pTri[0] == vCurID ) ? pTri[2] : pTri[0];
		else
			vCurID = ( pTri[0] == vCurID ) ? pTri[1] : pTri[0];
	}

	return false;
}


void VFTriangleMesh::FindNeighbours( TriangleID tID, TriangleID vNbrs[3] )
{
	Triangle & t = m_vTriangles[tID];
	for ( int i = 0; i < 3; ++i ) {

		VertexID v1 = t.nVertices[i];
		VertexID v2 = t.nVertices[ (i+1) % 3];

		// iterate over triangles of v1, looking for another tri with edge [v1,v2]
		Vertex & v = m_vVertices[v1];
		TriListEntry * pCur = v.pData->vTriangles.pFirst;
		TriListEntry * pLast = NULL;
		bool bFound = false;
		while ( pCur != NULL && ! bFound ) {
			pLast = pCur;
			TriangleID tCurID = pCur->data; pCur = pCur->pNext;
			if ( tCurID == tID )
				continue;
			VertexID * vTri2 = m_vTriangles[tCurID].nVertices;
			if ( vTri2[0] == v1 && ( vTri2[1] == v2 || vTri2[2] == v2 ) )
				bFound = true;
			else if ( vTri2[1] == v1 && ( vTri2[0] == v2 || vTri2[2] == v2 ) )
				bFound = true;
			else if ( vTri2[2] == v1 && ( vTri2[0] == v2 || vTri2[1] == v2 ) )
				bFound = true;
		}
		if ( bFound )
			vNbrs[i] = pLast->data;
		else
			vNbrs[i] = IMesh::InvalidID;
	}
}


void VFTriangleMesh::Copy( VFTriangleMesh & mesh )
{
	Clear(false);
	Wml::Vector3f vVertex, vNormal;
	vertex_iterator curv(mesh.BeginVertices()), endv(mesh.EndVertices());
	while ( curv != endv ) {
		VertexID vID = *curv;  ++curv;
		mesh.GetVertex(vID, vVertex, &vNormal);
		VertexID vNew = AppendVertex(vVertex, &vNormal);
//		if ( vNew != vID )
//			DebugBreak();
	}

	VertexID nTri[3];
	triangle_iterator curt(mesh.BeginTriangles()), endt(mesh.EndTriangles());
	while ( curt != endt ) { 
		TriangleID tID = *curt; ++curt;
		mesh.GetTriangle(tID, nTri);
		TriangleID tNew = AppendTriangle(nTri[0], nTri[1], nTri[2]);
//		if ( tNew != tID )
//			DebugBreak();
	}
}

void VFTriangleMesh::CopyVertInfo( VFTriangleMesh & mesh )
{
//	if ( mesh.GetVertexCount() != GetVertexCount() )
//		DebugBreak();
	Wml::Vector3f vVertex, vNormal;
	vertex_iterator curv(mesh.BeginVertices()), endv(mesh.EndVertices());
	while ( curv != endv ) {
		VertexID vID = *curv;  ++curv;
		mesh.GetVertex(vID, vVertex, &vNormal);
		SetVertex(vID, vVertex, &vNormal);
	}
}


void VFTriangleMesh::Copy( IMesh & copy, std::vector<VertexID> * vVertMap )
{
	Clear(false);

	std::vector< IMesh::VertexID > vInternalVertMap;
	std::vector< IMesh::VertexID > & vMapV = (vVertMap) ? *vVertMap : vInternalVertMap;
	vMapV.resize( copy.GetVertexCount() );

	Wml::Vector3f vVertex, vNormal;
	IMesh::IVtxIterator curv(copy.BeginIVertices()), endv(copy.EndIVertices());
	while ( curv != endv ) {
		VertexID vID = *curv;  ++curv;
		copy.GetVertex(vID, vVertex, &vNormal );
		vMapV[vID] = AppendVertex(vVertex, &vNormal);
	}

	TriangleID nTri[3];
	IMesh::ITriIterator curt(copy.BeginITriangles()), endt(copy.EndITriangles());
	while ( curt != endt ) {
		TriangleID tID = *curt;  ++curt;
		copy.GetTriangle(tID, nTri);
		AppendTriangle( vMapV[nTri[0]], vMapV[nTri[1]], vMapV[nTri[2]] );
	}
}

bool VFTriangleMesh::Append( VFTriangleMesh & mAppend, std::vector<IMesh::VertexID> & vVertexMap )
{
	unsigned int nMaxID = mAppend.GetMaxVertexID();
	std::vector<IMesh::VertexID> vNewIDs;
	bool bUseInternal = (vVertexMap.size() < nMaxID);
//	if ( bUseInternal )
//		DebugBreak();
	std::vector<IMesh::VertexID> & vMapV = 
		bUseInternal ? vNewIDs : vVertexMap;
	vMapV.resize( nMaxID );

	Wml::Vector3f vVertex, vNormal;
	vertex_iterator curv( mAppend.BeginVertices() ), endv( mAppend.EndVertices() );
	while ( curv != endv ) {
		IMesh::VertexID vID =  *curv;  ++curv;
		mAppend.GetVertex( vID, vVertex, &vNormal );
		IMesh::VertexID vNewID = AppendVertex( vVertex, &vNormal );
		vMapV[vID] = vNewID;
	}

	IMesh::VertexID nTri[3];
	triangle_iterator curt( mAppend.BeginTriangles() ), endt( mAppend.EndTriangles() );
	while ( curt != endt ) {
		IMesh::TriangleID tID =  *curt;  ++curt;
		mAppend.GetTriangle( tID, nTri );
		AppendTriangle( 
			vMapV[ nTri[0] ], vMapV[ nTri[1] ], vMapV[ nTri[2] ] );
	}

	return ! bUseInternal;
}



void VFTriangleMesh::Weld( VertexID vKeep, VertexID vDiscard )
{
	ASSERT( IsVertex(vKeep) && IsVertex(vDiscard) );
	
	//DebugBreak();		// This function is broken! One problem is that the TriListEntry iteration
						// will break once the triangle is changed...but that still doesn't explain
						// all the broken-ness...
	//_RMSInfo("Start weld\n");

	Vertex & v = m_vVertices[vDiscard];
	TriListEntry * pCur = v.pData->vTriangles.pFirst;
	while ( pCur != NULL && IsVertex(vDiscard) ) {
		TriangleID tID = pCur->data;
		pCur = pCur->pNext;
		VertexID nTri[3];
		GetTriangle(tID, nTri);
		bool bModified = false;
		for ( int j = 0; j < 3; ++j ) {
			if ( nTri[j] == vDiscard ) {
				nTri[j] = vKeep;
				bModified = true;
			}
		}
		if ( bModified ) {
			SetTriangle(tID, nTri[0], nTri[1], nTri[2]);
			//_RMSInfo(" Removed Tri - new count %d %d\n", GetTriangleCount(vDiscard), GetTriangleCount(vKeep));
			pCur = v.pData->vTriangles.pFirst;
		}
	}

	// should be unreferenced now...
	//_RMSInfo("  Vtx has %d %d triangles\n",GetTriangleCount(vDiscard), GetTriangleCount(vKeep));
	if ( IsVertex(vDiscard) && GetTriangleCount(vDiscard) == 0 )
		RemoveVertex(vDiscard);
}


void VFTriangleMesh::SplitEdge( VertexID e1, VertexID e2 )
{
	TriangleID t[2] = { InvalidID, InvalidID };
	int ti = 0;

	Vertex * pV1 = &m_vVertices[e1];
	Vertex * pV2 = &m_vVertices[e2];

	// find two triangles with edge e1e2
	TriListEntry * pCur = pV1->pData->vTriangles.pFirst;
	while ( pCur != NULL ) {
		Triangle * pTri = &m_vTriangles[pCur->data];
		if ( (pTri->nVertices[0] == e1 || pTri->nVertices[1] == e1 || pTri->nVertices[2] == e1)
			&& (pTri->nVertices[0] == e2 || pTri->nVertices[1] == e2 || pTri->nVertices[2] == e2) )
			t[ti++] = pCur->data;
		pCur = pCur->pNext;
	}
	if ( ti != 2 )
		return;

	// append new vertex
	Wml::Vector3f vInterp = 0.5f * (pV1->vVertex + pV2->vVertex);
	Wml::Vector3f nInterp = 0.5f * (pV1->vNormal + pV2->vNormal);
	VertexID vNew = AppendVertex(vInterp, &nInterp);
	
	// update triangles
	for ( int j = 0; j < 2; ++j ) {

		// figure out which index is which
		VertexID nTri[3];
		GetTriangle(t[j],nTri);
		VertexID nE1, nE2, nOther;
		for ( int k = 0; k < 3; ++k ) {
			if ( nTri[k] == e1 ) nE1 = k;
			else if ( nTri[k] == e2 ) nE2 = k;
			else nOther = k;
		}
		
		// set triangles
		nTri[nE1] = vNew;
		AppendTriangle(nTri[0], nTri[1], nTri[2]);
		nTri[nE1] = e1;  nTri[nE2] = vNew;
		SetTriangle(t[j], nTri[0], nTri[1], nTri[2]);
	}
}



void VFTriangleMesh::ReverseOrientation()
{
	triangle_iterator curt(BeginTriangles()), endt(EndTriangles());
	while ( curt != endt ) {
		TriangleID tID = *curt;  ++curt;
		Triangle & tri = m_vTriangles[tID];
		VertexID tmp = tri.nVertices[2];
		tri.nVertices[2] = tri.nVertices[1];
		tri.nVertices[1] = tmp;
	}
}



void VFTriangleMesh::GetBoundingBox( Wml::AxisAlignedBox3f & bounds )
{
	bounds.Min[0] = std::numeric_limits<float>::max();

	vertex_iterator curv( BeginVertices() ), endv( EndVertices() );
	while ( curv != endv ) {
		VertexID vID = *curv;  ++curv;
		Wml::Vector3f vVertex;
		GetVertex(vID, vVertex);
		if ( bounds.Min[0] == std::numeric_limits<float>::max() )
			bounds = Wml::AxisAlignedBox3f( vVertex.X(), vVertex.X(), vVertex.Y(), vVertex.Y(), vVertex.Z(), vVertex.Z() );
		else
			rms::Union( bounds, vVertex);
	}
}


float VFTriangleMesh::GetMaxEdgeLength()
{
	float fMaxEdgeLength = 0.0f;

	triangle_iterator curt( BeginTriangles()), endt(EndTriangles());
	while ( curt != endt ) {
		TriangleID nID = *curt;
		curt++;

		Wml::Vector3f vVertices[3];
		GetTriangle(nID, vVertices);

		for ( int i = 0; i < 3; ++i ) {
			float fLen = (vVertices[i] - vVertices[(i+1)%3]).SquaredLength();
			if ( fLen > fMaxEdgeLength )
				fMaxEdgeLength = fLen;
		}
	}
	return (float)sqrt(fMaxEdgeLength);	
}





void VFTriangleMesh::ClearBit( unsigned int nBit )
{
	RefCountedVector<Vertex>::item_iterator 
		curv(m_vVertices.begin_items()), endv(m_vVertices.end_items());
	while ( curv != endv ) {
		(*curv).vecData.nBits &= ~(1<<nBit);
		++curv;
	}
}


bool VFTriangleMesh::ReadOBJ( const char * pFilename, std::string & errString )
{
	Clear(false);

	ifstream in(pFilename);
	if(!in){
		errString = string("Cannot open file ") + string(pFilename);
		cerr << errString << endl;
		return false;
	}

	string command;
	char c1,c2;
	unsigned int tv1, tv2, tv3, tn1, tn2, tn3, tt1, tt2, tt3;
	char linebuf[1024];

	Wml::Vector3f fvec = Wml::Vector3f::ZERO;

	bool bHasTextures = false;

	// need to save normals separately and then match to vertices (maya
	//  "optimizes" the mesh...argh!)
	std::vector<Wml::Vector3f> vNormals;

	unsigned int ivec[3];
	while(in){
		ostrstream s;
		in >> command;
		if(!in)
			continue;
		switch(command.c_str()[0]){

		case 'v':    
			in >> fvec[0] >> fvec[1];
			if(!in)
				continue;
			switch(command.c_str()[1]){
			case '\0':  // vertex
				in >> fvec[2];
				AppendVertex(fvec);
				break;
			case 'n': // vertex normal
				in >> fvec[2];
				fvec.Normalize();
				vNormals.push_back( fvec );
				break;
			case 't':
//				AppendVertexData(NULL, NULL, NULL, fVec);
				bHasTextures = true;
				break;
			default:
				string err("Got unknown OBJ command ");
				err += command;
				cerr << err << endl;
			}
		break;

		case 'f':
			if ( bHasTextures ) {
				in >> tv1 >> c1 >> tt1 >> c2 >> tn1;
				in >> tv2 >> c1 >> tt2 >> c2 >> tn2;
				in >> tv3 >> c1 >> tt3 >> c2 >> tn3;
				
			} else {
				in >> tv1 >> c1 >> c2 >> tn1;
				in >> tv2 >> c1 >> c2 >> tn2;
				in >> tv3 >> c1 >> c2 >> tn3;
			}
			ivec[0] = tv1-1; ivec[1] = tv2-1; ivec[2] = tv3-1;
			AppendTriangle(ivec[0], ivec[1], ivec[2]);

			// set proper normals
			SetNormal( ivec[0], vNormals[ tn1-1 ] );
			SetNormal( ivec[1], vNormals[ tn2-1 ] );
			SetNormal( ivec[2], vNormals[ tn3-1 ] );
			break;

		default:
			in.getline(linebuf, 1023, '\n');
		}
	}

	return true;
}


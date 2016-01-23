// Ryan Schmidt
// Copyright (c) 2006-2010
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt

#ifndef __RMS_VF_TRIANGLE_MESH_H__
#define __RMS_VF_TRIANGLE_MESH_H__

#include <Wm4Vector3.h>
#include <Wm4AxisAlignedBox3.h>
#include <Wm4Vector2.h>
#include <Wm4ColorRGBA.h>
#include <set>
#include <map>
#include <vector>
#include <hash_map>
#include <limits>

#include "IMesh.h"
#include "MemoryPool.h"
#include "RefCountedVector.h"
#include "SparseArray.h"

using namespace std;

namespace rms {

class VFTriangleMesh : public IMesh
{
public:
	VFTriangleMesh(void);
	VFTriangleMesh( IMesh & copy, std::vector<VertexID> * vVertMap = NULL );
	~VFTriangleMesh(void);

	bool ReadOBJ( const char * pFilename, std::string & errString );


/*
 * IMesh read interface (mandatory)
 */
	virtual void GetVertex( VertexID vID, Wml::Vector3f & vVertex, Wml::Vector3f * pNormal = NULL ) const;
	virtual void GetNormal( VertexID vID, Wml::Vector3f & vNormal ) const;
	virtual unsigned int GetVertexCount() const;

	virtual void GetTriangle( TriangleID tID, VertexID vTriangle[3]  ) const;
	virtual void GetTriangle( TriangleID tID, Wml::Vector3f vTriangle[3], Wml::Vector3f * pNormals = NULL ) const;
	virtual unsigned int GetTriangleCount() const;

/*
 *  IMesh write interface (optional)
 */
	virtual VertexID AppendVertex( const Wml::Vector3f & vVertex, const Wml::Vector3f * pNormal = NULL );

	virtual TriangleID AppendTriangle( VertexID v1, VertexID v2, VertexID v3 );
	virtual void SetTriangle( TriangleID tID, VertexID v1, VertexID v2, VertexID v3 );

	virtual void Clear( bool bFreeMem );

/*
 * IMesh mesh info interface - has default implementation
 */
	//! initialize vertex neighbour iteration
	virtual void BeginVtxTriangles( VtxNbrItr & v );

	//! (possibly) un-ordered iteration around one-ring of a vertex. Returns InvalidID when done
	virtual TriangleID GetNextVtxTriangle( VtxNbrItr & v );

	//! determine if a vertex is on the mesh boundary (if there is one)
	virtual bool IsBoundaryVertex( VertexID vID );

/*
 *  functions that should probably go into IMesh
 */
	inline virtual void SetVertex( VertexID vID, const Wml::Vector3f & vVertex, const Wml::Vector3f * pNormal = NULL );
	virtual void SetNormal( VertexID vID, const Wml::Vector3f & vNormal );
	virtual VertexID GetTriVertex( TriangleID tID, int nIndex ) const;
	virtual void GetTriVertex( TriangleID tID, int nIndex, Wml::Vector3f & vVertex ) const;
	virtual unsigned int GetTriangleCount( IMesh::VertexID vID ) const;
	inline void GetTriangleNormal( TriangleID tID, Wml::Vector3f & vNormal );

	// stuff that is sort of specific to this class...
	inline unsigned int GetMaxVertexID() const;
	inline unsigned int GetMaxTriangleID() const;
	inline bool IsVertex( VertexID v ) const;
	inline bool IsTriangle( TriangleID t ) const;

	virtual void Copy( VFTriangleMesh & pMesh );
	virtual void CopyVertInfo( VFTriangleMesh & pMesh );
	virtual void Copy( IMesh & pCopy, std::vector<VertexID> * vVertMap = NULL );

	// if vVertexMap is not large enough ( GetMaxVertexID() ), an internal vector is used, and false returns
	virtual bool Append( VFTriangleMesh & mAppend, std::vector<IMesh::VertexID> & vVertexMap );

	virtual void GetColor( VertexID vID, Wml::ColorRGBA & cColor ) const;
	virtual void SetColor( VertexID vID, const Wml::ColorRGBA & cColor );

	// mesh removal functions
	void RemoveUnreferencedGeometry();
	void RemoveVertex( VertexID vID );
	void RemoveTriangle( TriangleID tID );

	//! neighbour tris are listed in order for edges [0,1],  [1,2],  [2,0]. Maybe be invalid, if no nbr
	void FindNeighbours( TriangleID tID, TriangleID vNbrs[3] );

	// mesh editing operations
	void Weld( VertexID vKeep, VertexID vDiscard );
	void SplitEdge( VertexID e1, VertexID e2 );
	void ReverseOrientation();

	// mesh info
	void GetBoundingBox( Wml::AxisAlignedBox3f & bounds );
	float GetMaxEdgeLength();

	// bitmask functions. Bits [0:15] are reserved for internal use of VFTriangleMesh.
	//  Bits [16:] are available for callers, but be careful...
	void ClearBit( unsigned int nBit );
	void ClearBit( VertexID vID, unsigned int nBit );
	void SetBit( VertexID vID, unsigned int nBit );
	bool GetBit( VertexID vID, unsigned int nBit );


protected:

	//! memory pool for fast allocation of per-vertex triangle list elements
	ListPool<TriangleID> m_VertListPool;

	// per-vertex triangle list management
	inline void ClearTriList( VertexID nVertID );
	inline void AddTriEntry( TriangleID nTriID, VertexID nVertID );
	inline void RemoveTriEntry( TriangleID nTriID, VertexID nVertID );

	typedef ListPool<TriangleID> TriListPool;
	typedef ListPool<TriangleID>::Entry TriListEntry;
	struct VertexData {
		ListPool<TriangleID>::List vTriangles;
		Wml::ColorRGBA vColor;
	};

	MemoryPool<VertexData> m_VertDataMemPool;

	struct Vertex {
		Wml::Vector3f vVertex;		//! vertex location
		Wml::Vector3f vNormal;		//! vertex normal 
		VertexData * pData;			//! other per-vertex data
		unsigned int nBits;			//! insanely-useful per-vertex bitmask

		Vertex() : pData(NULL) { }
		Vertex( const Wml::Vector3f & v, const Wml::Vector3f & n )
			: vVertex(v), vNormal(n), pData(NULL) {}
	};

	struct Triangle {
		VertexID nVertices[3];		//! vertices of this triangle

		Triangle() {};
		Triangle( VertexID v1, VertexID v2, VertexID v3 )
			{ nVertices[0] = v1; nVertices[1] = v2; nVertices[2] = v3; }
	};

	RefCountedVector<Vertex> m_vVertices;
	RefCountedVector<Triangle> m_vTriangles;

public:
	typedef RefCountedVector<Vertex>::index_iterator vertex_iterator;
	inline vertex_iterator BeginVertices()
		{ return m_vVertices.begin_indexes(); }
	inline vertex_iterator EndVertices()
		{ return m_vVertices.end_indexes(); }

	typedef RefCountedVector<Triangle>::index_iterator triangle_iterator;
	inline triangle_iterator BeginTriangles()
		{ return m_vTriangles.begin_indexes(); }
	inline triangle_iterator EndTriangles()
		{ return m_vTriangles.end_indexes(); }


protected:
/*
 * IMesh iterator interface
 */
	inline virtual void * ivtx_make_iterator(bool bStart)
		{ return new vertex_iterator( (bStart) ? BeginVertices() : EndVertices() ); }
	inline virtual void * ivtx_make_iterator( void * pFromItr )
		{ return new vertex_iterator( * ((vertex_iterator *)pFromItr) ); }
	inline virtual void ivtx_free_iterator( void * pItr ) 
		{ delete (vertex_iterator *)pItr; }
	inline virtual void ivtx_set( void * pItr, void * pTo )
		{ *((vertex_iterator *)pItr) = *((vertex_iterator *)pTo); }
	inline virtual void ivtx_goto_next( void * pItr )
		{ ++(*((vertex_iterator *)pItr)); }
	inline virtual bool ivtx_equal( void * pItr1, void * pItr2 )
		{ return *((vertex_iterator *)pItr1) == *((vertex_iterator *)pItr2); }
	inline virtual VertexID ivtx_value( void * pItr )
		{ return **((vertex_iterator *)pItr); }

	inline virtual void * itri_make_iterator(bool bStart)
		{ return new triangle_iterator( (bStart) ? BeginTriangles() : EndTriangles() ); }
	inline virtual void * itri_make_iterator( void * pFromItr )
		{ return new triangle_iterator( * ((triangle_iterator *)pFromItr) ); }
	inline virtual void itri_free_iterator( void * pItr ) 
		{ delete (triangle_iterator *)pItr; }
	inline virtual void itri_set( void * pItr, void * pTo )
		{ *((triangle_iterator *)pItr) = *((triangle_iterator *)pTo); }
	inline virtual void itri_goto_next( void * pItr )
		{ ++(*((triangle_iterator *)pItr)); }
	inline virtual bool itri_equal( void * pItr1, void * pItr2 )
		{ return *((triangle_iterator *)pItr1) == *((triangle_iterator *)pItr2); }
	inline virtual TriangleID itri_value( void * pItr )
		{ return **((triangle_iterator *)pItr); }
};



inline IMesh::VertexID VFTriangleMesh::AppendVertex( const Wml::Vector3f & vVertex, const Wml::Vector3f * pNormal )
{ 
	VertexID vNewID = 
		(VertexID)m_vVertices.insert( Vertex(vVertex, (pNormal) ? *pNormal : Wml::Vector3f::UNIT_Z ) );
	if ( m_vVertices[vNewID].pData == NULL ) {
		m_vVertices[vNewID].pData = m_VertDataMemPool.Allocate();
		m_vVertices[vNewID].pData->vTriangles = m_VertListPool.GetList();
	}

	// clear lists...
	VertexData & v = * m_vVertices[vNewID].pData;
	m_VertListPool.Clear( v.vTriangles );

	return vNewID;
}

inline void VFTriangleMesh::SetVertex( VertexID vID, const Wml::Vector3f & vVertex, const Wml::Vector3f * pNormal )
{
	Vertex & v = m_vVertices[vID]; 
	v.vVertex = vVertex;
	if ( pNormal )
		v.vNormal = *pNormal;
}

inline void VFTriangleMesh::SetNormal( VertexID vID, const Wml::Vector3f & vNormal )
{
	Vertex & v = m_vVertices[vID]; 
	v.vNormal = vNormal;
}

inline void VFTriangleMesh::GetVertex( IMesh::VertexID vID, Wml::Vector3f & vVertex, Wml::Vector3f * pNormal ) const 
{ 
	const Vertex & v = m_vVertices[vID]; 
	vVertex = v.vVertex; 
	if ( pNormal ) 
		*pNormal = v.vNormal; 
}

inline void VFTriangleMesh::GetNormal( IMesh::VertexID vID, Wml::Vector3f & vNormal ) const
{ 
	const Vertex & v = m_vVertices[vID]; 
	vNormal = v.vNormal; 
}

inline unsigned int VFTriangleMesh::GetVertexCount() const
{ 
	return m_vVertices.size(); 
}


inline void VFTriangleMesh::GetTriangle( IMesh::TriangleID tID, IMesh::VertexID vTriangle[3]  ) const
{
	const Triangle & t = m_vTriangles[tID];
	memcpy( vTriangle, t.nVertices, sizeof(TriangleID)*3 );
}

inline void VFTriangleMesh::GetTriangle( IMesh::TriangleID tID, Wml::Vector3f vTriangle[3], Wml::Vector3f * pNormals ) const
{
	const Triangle & t = m_vTriangles[tID];
	if ( pNormals ) {
		GetVertex( t.nVertices[0], vTriangle[0], & pNormals[0] );
		GetVertex( t.nVertices[1], vTriangle[1], & pNormals[1] );
		GetVertex( t.nVertices[2], vTriangle[2], & pNormals[2] );
	} else {
		GetVertex( t.nVertices[0], vTriangle[0], NULL );
		GetVertex( t.nVertices[1], vTriangle[1], NULL );
		GetVertex( t.nVertices[2], vTriangle[2], NULL );
	}
}

inline IMesh::VertexID VFTriangleMesh::GetTriVertex( IMesh::TriangleID tID, int nIndex ) const
{
	const Triangle & t = m_vTriangles[tID];
	return t.nVertices[ nIndex ];
}

inline void VFTriangleMesh::GetTriVertex( IMesh::TriangleID tID, int nIndex, Wml::Vector3f & vVertex ) const
{
	const Triangle & t = m_vTriangles[tID];
	GetVertex(tID, vVertex);
}

inline unsigned int VFTriangleMesh::GetTriangleCount() const
{
	return m_vTriangles.size();
}

inline unsigned int VFTriangleMesh::GetTriangleCount( IMesh::VertexID vID ) const
{ 
//	VertTriEntry *pCur = m_vVertices[vID].pTriList;		// RMS VD
	TriListEntry * pCur = m_vVertices[vID].pData->vTriangles.pFirst;

	int nCount = 0;
	while ( pCur != NULL ) {
		++nCount;
		pCur = pCur->pNext;
	}
	return nCount;
}

inline unsigned int VFTriangleMesh::GetMaxVertexID() const
{
	return m_vVertices.max_index();
}

inline unsigned int VFTriangleMesh::GetMaxTriangleID() const
{
	return m_vTriangles.max_index();
}

inline bool VFTriangleMesh::IsVertex( VertexID v ) const
{
	return v != InvalidID && m_vVertices.isValid(v);
}

inline bool VFTriangleMesh::IsTriangle( TriangleID t ) const
{
	return t != InvalidID && m_vTriangles.isValid(t);
}


inline void VFTriangleMesh::GetTriangleNormal( TriangleID tID, Wml::Vector3f & vNormal )
{
	VertexID * verts = m_vTriangles[tID].nVertices;
	Wml::Vector3f e1( m_vVertices[verts[1]].vVertex - m_vVertices[verts[0]].vVertex );
	Wml::Vector3f e2( m_vVertices[verts[2]].vVertex - m_vVertices[verts[0]].vVertex );
	e1.Normalize();
	e2.Normalize();
	vNormal = e1.Cross(e2);
	vNormal.Normalize();
}


inline void VFTriangleMesh::AddTriEntry( TriangleID nTriID, VertexID nVertID )
{
	ASSERT( m_vVertices.isValid(nVertID) );
	VertexData & v = * m_vVertices[nVertID].pData;
	m_VertListPool.Insert( v.vTriangles, nTriID );
}

inline void VFTriangleMesh::RemoveTriEntry( TriangleID nTriID, VertexID nVertID )
{
	ASSERT( m_vVertices.isValid(nVertID) );
	VertexData & v = * m_vVertices[nVertID].pData;
	m_VertListPool.Remove( v.vTriangles, nTriID );
}

inline void VFTriangleMesh::ClearTriList( VertexID nVertID )
{
	ASSERT( m_vVertices.isValid(nVertID) );
	VertexData & v = * m_vVertices[nVertID].pData;
	m_VertListPool.Clear( v.vTriangles );
}

inline void VFTriangleMesh::GetColor( VertexID vID, Wml::ColorRGBA & cColor ) const
{
	ASSERT( m_vVertices.isValid(vID) );
	VertexData & v = * m_vVertices[vID].pData;
	cColor = v.vColor;
}
inline void VFTriangleMesh::SetColor( VertexID vID, const Wml::ColorRGBA & cColor )
{
	ASSERT( m_vVertices.isValid(vID) );
	VertexData & v = * m_vVertices[vID].pData;
	v.vColor = cColor;
}



inline void VFTriangleMesh::ClearBit( VertexID vID, unsigned int nBit )
{
	ASSERT( m_vVertices.isValid(vID) );
	m_vVertices[vID].nBits &= ~(1<<nBit);
}

inline void VFTriangleMesh::SetBit( VertexID vID, unsigned int nBit )
{
	ASSERT( m_vVertices.isValid(vID) );
	m_vVertices[vID].nBits |= (1<<nBit);
}

inline bool VFTriangleMesh::GetBit( VertexID vID, unsigned int nBit )
{
	return ( m_vVertices[vID].nBits & (1<<nBit) ) != 0;
}


} // namespace rmsmesh


#endif  // __RMS_VF_TRIANGLE_MESH_H__
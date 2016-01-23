// Ryan Schmidt
// Copyright (c) 2006-2010
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt

#ifndef __RMS_ITRIANGLEMESH_3D__
#define __RMS_ITRIANGLEMESH_3D__

#include <vector>
#include <Wm4Vector3.h>
#include <Wm4Vector2.h>
#include "SparseArray.h"


namespace rms
{


class IMesh
{
public:
	typedef unsigned int VertexID;
	typedef unsigned int TriangleID;
	typedef unsigned int UVSetID;

	static unsigned int InvalidID;

/*
 * IMesh read interface (mandatory)
 */
	virtual void GetVertex( VertexID vID, Wml::Vector3f & vVertex, Wml::Vector3f * pNormal = NULL ) const = 0;
	virtual void GetNormal( VertexID vID, Wml::Vector3f & vNormal ) const = 0;
	virtual unsigned int GetVertexCount() const = 0;

	virtual void GetTriangle( TriangleID tID, VertexID vTriangle[3]  ) const = 0;
	virtual void GetTriangle( TriangleID tID, Wml::Vector3f vTriangle[3], Wml::Vector3f * pNormals = NULL ) const = 0;
	virtual unsigned int GetTriangleCount() const = 0;

/*
 *  IMesh write interface (optional)
 */
	virtual VertexID AppendVertex( const Wml::Vector3f & vVertex, const Wml::Vector3f * pNormal = NULL ) { return InvalidID; }

	virtual TriangleID AppendTriangle( VertexID v1, VertexID v2, VertexID v3 ) { return InvalidID; }
	virtual void SetTriangle( TriangleID tID, VertexID v1, VertexID v2, VertexID v3 ) {}

	virtual void Clear( bool bFreeMem ) {
		for ( unsigned int i = 0; i < m_vUVSets.size(); ++i )
			m_vUVSets[i]->Clear();
		//	if ( m_vUVSets[i] )
		//		delete m_vUVSets[i];
		//m_vUVSets.clear();
	}


	struct VtxNbrItr {
		VertexID vID;
		unsigned long long nData[2];
		VtxNbrItr() {}
		VtxNbrItr( VertexID id ) { vID = id; }
	};

/*
 * IMesh mesh info interface - has default implementation
 */
	//! initialize vertex neighbour iteration
	virtual void BeginVtxTriangles( VtxNbrItr & v )  {}

	//! (possibly) un-ordered iteration around one-ring of a vertex. Returns InvalidID when done
	virtual TriangleID GetNextVtxTriangle( VtxNbrItr & v ) { return InvalidID; }

	//! determine if a vertex is on the mesh boundary (if there is one)  (depends on GetNextVtxTriangle)
	virtual bool IsBoundaryVertex( VertexID vID );


/*
 * callback versions that actually work...
 */
	class NeighborTriCallback {
	public:
		virtual void BeginTriangles() {}
		virtual void NextTriangle( TriangleID tID ) = 0;
		virtual void EndTriangles() {}
	};

	virtual void NeighbourIteration( VertexID vID, NeighborTriCallback * pCallback );

/*
 * IMesh UV interface - has default implementation
 */
	inline virtual bool HasUVSet( UVSetID nSetID ) 
		{ return nSetID < m_vUVSets.size(); }
	inline virtual UVSetID AppendUVSet( ) 
		{ m_vUVSets.push_back( new UVSet() ); return (unsigned int)m_vUVSets.size()-1; }
	inline virtual void InitializeUVSet( UVSetID nSetID )
		{ m_vUVSets[nSetID]->Initialize(this); }
	inline virtual bool GetUV( VertexID vID, UVSetID nSetID, Wml::Vector2f & vUV ) const
		{ return m_vUVSets[nSetID]->GetUV(vID, vUV); }
	inline virtual bool GetTriangleUV( TriangleID tID, UVSetID nSetID, Wml::Vector2f vUV[3] ) const
		{ VertexID nTri[3]; GetTriangle(tID, nTri);
		  return m_vUVSets[nSetID]->GetUV(nTri, vUV); }
	inline virtual void SetUV( VertexID vID, UVSetID nSetID, const Wml::Vector2f & vUV )
		{ m_vUVSets[nSetID]->SetUV(vID, vUV); }
	inline virtual void ClearUVSet( UVSetID nSetID )
		{ m_vUVSets[nSetID]->Clear(); }

/*
 * IMesh iterator interface
 */
public:

	/*
	 * iterators
	 */
	class IVtxIterator {
	public:
		inline IVtxIterator( const IVtxIterator & i2 ) { 
			m_pMesh = i2.m_pMesh;
			m_itr = m_pMesh->ivtx_make_iterator( i2.m_itr ); }
		inline virtual ~IVtxIterator() 
			{ if ( m_itr ) m_pMesh->ivtx_free_iterator( m_itr ); }
		inline const IVtxIterator & operator=( const IVtxIterator & i2 )
			{ m_pMesh->ivtx_set(m_itr, i2.m_itr); return *this; }
		
		inline VertexID operator*() 
			{ return m_pMesh->ivtx_value( m_itr ); }

		//! postfix operator doesn't return value, avoids memory allocation for copy
		inline void operator++(int)  // postfix
			{ m_pMesh->ivtx_goto_next(m_itr); }
		inline IVtxIterator & operator++() // prefix
			{ m_pMesh->ivtx_goto_next(m_itr); return *this; }

		inline bool operator==( const IVtxIterator & i2 ) 
			{ return m_pMesh->ivtx_equal( m_itr, i2.m_itr ); }
		inline bool operator!=( const IVtxIterator & i2 ) 
			{ return ! m_pMesh->ivtx_equal( m_itr, i2.m_itr ); }

	protected:
		IMesh * m_pMesh;
		void * m_itr;
	
		//! if bStart, return start iterator, else return end iterator
		inline IVtxIterator(IMesh * pMesh, bool bStart) : m_pMesh(pMesh)
			{ m_itr = m_pMesh->ivtx_make_iterator(bStart); }
		friend class IMesh;
	};


	class ITriIterator {
	public:
		inline ITriIterator( const ITriIterator & i2 ) {
			m_pMesh = i2.m_pMesh;
			m_itr = m_pMesh->itri_make_iterator( i2.m_itr ); }
		inline virtual ~ITriIterator() { 
			if ( m_itr ) m_pMesh->itri_free_iterator( m_itr ); }
		inline const ITriIterator & operator=( const ITriIterator & i2 ) { 
			m_pMesh->itri_set(m_itr, i2.m_itr); return *this; }
		
		inline TriangleID operator*() 
			{ return m_pMesh->itri_value( m_itr ); }

		//! postfix operator doesn't return value, avoids memory allocation for copy
		inline void operator++(int)  // postfix
			{ m_pMesh->itri_goto_next(m_itr); }
		inline ITriIterator & operator++() // prefix
			{ m_pMesh->itri_goto_next(m_itr); return *this; }

		inline bool operator==( const ITriIterator & i2 ) 
			{ return m_pMesh->itri_equal( m_itr, i2.m_itr ); }
		inline bool operator!=( const ITriIterator & i2 ) 
			{ return ! m_pMesh->itri_equal( m_itr, i2.m_itr ); }

	protected:
		IMesh * m_pMesh;
		void * m_itr;
	
		//! if bStart, return start iterator, else return end iterator
		inline ITriIterator(IMesh * pMesh, bool bStart) : m_pMesh(pMesh)
			{ m_itr = m_pMesh->itri_make_iterator(bStart); }
		friend class IMesh;
	};


	inline IVtxIterator BeginIVertices()   { return IVtxIterator(this, true); }
	inline IVtxIterator EndIVertices()     { return IVtxIterator(this, false); }

	inline ITriIterator BeginITriangles()  { return ITriIterator(this, true); }
	inline ITriIterator EndITriangles()    { return ITriIterator(this, false); }

protected:
/*
 * IMesh iterator interface (mandatory)
 */
	virtual void * ivtx_make_iterator(bool bStart) = 0;
	virtual void * ivtx_make_iterator( void * pFromIterator ) = 0;
	virtual void ivtx_set( void * pIterator, void * pToIteratorVal ) = 0;
	virtual void ivtx_free_iterator( void * pIterator ) = 0;
	virtual void ivtx_goto_next( void * pIterator ) = 0;
	virtual bool ivtx_equal( void * pIterator1, void * pIterator2 ) = 0;
	virtual VertexID ivtx_value( void * pIterator ) = 0;

	virtual void * itri_make_iterator(bool bStart) = 0;
	virtual void * itri_make_iterator( void * pFromIterator ) = 0;
	virtual void itri_set( void * pIterator, void * pToIteratorVal ) = 0;
	virtual void itri_free_iterator( void * pIterator ) = 0;
	virtual void itri_goto_next( void * pIterator ) = 0;
	virtual bool itri_equal( void * pIterator1, void * pIterator2 ) = 0;
	virtual TriangleID itri_value( void * pIterator ) = 0;




public:
/*
 * UVSet class
 */
	class UVSet {
	public:
		UVSet() {}
		UVSet( IMesh * pMesh )
			{ Initialize(pMesh); }

		inline void Initialize( IMesh * pMesh )
			{ Clear(); m_vUV.resize(pMesh->GetVertexCount()); }

		inline void Clear()
			{ m_vUV.clear(false); }

		inline void SetUV( VertexID vID, const Wml::Vector2f & vUV ) 
			{ m_vUV.set(vID, vUV); }
		
		inline bool HasUV( VertexID vID ) const
			{ return m_vUV.has(vID); }

		inline bool GetUV( VertexID vID, Wml::Vector2f & vUV ) const
			{ if ( m_vUV.has(vID) )
				{ vUV = m_vUV[vID]; return true; } 
			  return false; 
			}
		
		bool GetUV( VertexID vID[3], Wml::Vector2f vUV[3] ) const
			{ return GetUV(vID[0],vUV[0]) && GetUV(vID[1],vUV[1]) && GetUV(vID[2],vUV[2]); }

	protected:
		SparseArray<Wml::Vector2f> m_vUV;
	};

	UVSet & GetUVSet( UVSetID setID ) { return * m_vUVSets[setID]; }
	
protected:
	std::vector< UVSet * > m_vUVSets;



};


} // end namespace rms


#endif  // __RMS_ITRIANGLEMESH_3D__
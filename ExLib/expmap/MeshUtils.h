// Ryan Schmidt
// Copyright (c) 2006-2010
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt

#pragma once

#include "VFTriangleMesh.h"
#include "IMeshBVTree.h"
#include "Wm4AxisAlignedBox2.h"
#include "Wm4Matrix3.h"

namespace rms {

class MeshUtils
{
public:


	//! bClosed flag is only set if bOrdered is true
	static bool VertexOneRing( VFTriangleMesh & mesh, IMesh::VertexID vID, std::vector<IMesh::VertexID> & vOneRing );

	//! bClosed flag is only set if bOrdered is true
	static bool TriangleOneRing( VFTriangleMesh & mesh, IMesh::VertexID vID, std::vector<IMesh::TriangleID> & vOneRing );

	static bool GetTriVerts( IMesh::VertexID nTri[3], IMesh::VertexID vID, IMesh::VertexID & other1, IMesh::VertexID & other2 );

	static Wml::Vector3f MeshLaplacian( VFTriangleMesh & mesh, IMesh::VertexID vID, std::vector<IMesh::VertexID> & vOneRing, std::vector<float> & vWeights );

	static void GetEdgeLengthStats( IMesh * pMesh, float & fMin, float & fMax, float & fAvg );

	// currently only does one-ring...
	static Wml::Vector3f GetAverageNormal( VFTriangleMesh & mesh, IMesh::VertexID vID, unsigned int nKRings = 1 );
	static void SmoothNormals( VFTriangleMesh & mesh );

	enum NormalEstMode {
		UniformFaceAvg,
		AreaWeightedFaceAvg
	};
	static Wml::Vector3f EstimateNormal( VFTriangleMesh & mesh, IMesh::VertexID vID, NormalEstMode eMode = AreaWeightedFaceAvg );
	static void EstimateNormals( VFTriangleMesh & mesh, NormalEstMode eMode = AreaWeightedFaceAvg, bool bSkipBoundary = false, Wml::Vector3f * pBuffer = NULL );

	static Wml::Vector3f FaceNormal( VFTriangleMesh & mesh, IMesh::TriangleID tID );
	static Wml::Vector3f FaceInterpNormal( IMesh * pMesh, IMesh::TriangleID tID, const Wml::Vector3f & vTriPoint );

	//! note: scale is currently not applied to normals...
	static void ScaleMesh( VFTriangleMesh & mesh, const Wml::Vector3f & vScale, const Wml::Vector3f & vCenter = Wml::Vector3f::ZERO );
	static void TranslateMesh( VFTriangleMesh & mesh, const Wml::Vector3f & vTranslate );
	static void RotateMesh( VFTriangleMesh & mesh, const Wml::Matrix3f & mRotate, const Wml::Vector3f & vCenter = Wml::Vector3f::ZERO );

	static float LoopLength( const VFTriangleMesh & mesh, const std::vector<IMesh::VertexID> & vLoop );
	static float LoopLengthUV( const VFTriangleMesh & mesh, const std::vector<IMesh::VertexID> & vLoop, IMesh::UVSetID nSetID );
	static void ScaleMeshByLoop( VFTriangleMesh & mesh, const std::vector<IMesh::VertexID> & vLoop, float fTargetPerimeter );

	static float MeshArea( VFTriangleMesh & mesh );
	static float MeshUVArea( VFTriangleMesh & mesh, IMesh::UVSetID nUVSet = 0 );


	static void ScaleMeshUV( VFTriangleMesh & mesh, IMesh::UVSetID nSetID, const Wml::Vector2f & vScale, const Wml::Vector2f & vCenter = Wml::Vector2f::ZERO );


	//! appends faces of vID to vSelection
	static void SelectFaces( VFTriangleMesh & mesh, IMesh::VertexID & vID, std::set<IMesh::TriangleID> & vSelection );

	//! appends all faces connected to vVerts to vSelection
	static void SelectFaces( VFTriangleMesh & mesh, const std::vector<IMesh::VertexID> & vVerts, std::set<IMesh::TriangleID> & vSelection );

	// TODO: this would be more efficient if we could find the boundary verts of a selection...
	//! appends all one-ring neighbours of current selection. pNewVerts returns added verts, if requested
	static void GrowSelection( VFTriangleMesh & mesh, 
							   std::set<IMesh::TriangleID> & vSelection, 
							   std::set<IMesh::VertexID> * pNewVerts );


	//! returns the two vertices from nTri which != vOpp  (returns false if nTri does not contain vOpp)
	static bool PickTriVerts( IMesh::VertexID nTri[3], IMesh::VertexID vOpp, IMesh::VertexID & vOther1, IMesh::VertexID & vOther2 );


private:
	MeshUtils(void);
	~MeshUtils(void);
};


/*
 * This one-ring neighbour iteration constructs a list [parents] of all
 * the neighbouring vertices of [vertid] which are in [knownverts]
 */
class MakeOneRingListCallback : public IMesh::NeighborTriCallback
{
public:
	MakeOneRingListCallback( IMesh::VertexID vertid, VFTriangleMesh * mesh) {
		vID = vertid;
		pMesh = mesh;
	}
	IMesh::VertexID vID;
	std::set<IMesh::VertexID> vOneRing;
	VFTriangleMesh * pMesh;

	virtual void NextTriangle( IMesh::TriangleID tID ) {
		IMesh::VertexID nTri[3];
		pMesh->GetTriangle(tID, nTri);
		for ( int k = 0; k < 3; ++k ) {
			if ( nTri[k] != vID )
				vOneRing.insert( nTri[k] );
		}
	}
};




/*
 * simple classes for one-ring triangle iterations
 */

class NeighborTriPassThruCallback : public IMesh::NeighborTriCallback
{
public:
	NeighborTriPassThruCallback( IMesh::NeighborTriCallback * pPassThru )
		{ m_pPassThru = pPassThru; }

	virtual void NextTriangle( IMesh::TriangleID tID )
		{ m_pPassThru->NextTriangle(tID); }

protected:
	IMesh::NeighborTriCallback * m_pPassThru;
};


class NeighborTriBuffer : public IMesh::NeighborTriCallback
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







}  // end namespace rms


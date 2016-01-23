// Ryan Schmidt
// Copyright (c) 2006-2010
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt


//#include <stdafx.h>
#include "MeshUtils.h"

#include "VectorUtil.h"

using namespace rms;

MeshUtils::MeshUtils(void)
{
}

MeshUtils::~MeshUtils(void)
{
}



bool MeshUtils::VertexOneRing( VFTriangleMesh & mesh, IMesh::VertexID vID, std::vector<IMesh::VertexID> & vOneRing )
{
	std::set<IMesh::VertexID> vNbrs;

	VFTriangleMesh::VtxNbrItr itr(vID);
	mesh.BeginVtxTriangles(itr);
	IMesh::TriangleID tID = mesh.GetNextVtxTriangle(itr);
	while ( tID != IMesh::InvalidID ) {
		IMesh::VertexID nTri[3];
		mesh.GetTriangle( tID, nTri );
		for ( int k = 0; k < 3; ++k )
			if ( nTri[k] != vID )
				vNbrs.insert(nTri[k]);
		tID = mesh.GetNextVtxTriangle(itr);
	}
	vOneRing.insert(vOneRing.end(), vNbrs.begin(), vNbrs.end());
	return true;
}



//! bClosed flag is only set if bOrdered is true
bool MeshUtils::TriangleOneRing( VFTriangleMesh & mesh, IMesh::VertexID vID, std::vector<IMesh::TriangleID> & vOneRing )
{
	VFTriangleMesh::VtxNbrItr itr(vID);
	mesh.BeginVtxTriangles(itr);
	IMesh::TriangleID tID = mesh.GetNextVtxTriangle(itr);
	while ( tID != IMesh::InvalidID ) {
		vOneRing.push_back(tID);
		tID = mesh.GetNextVtxTriangle(itr);
	}
	return true;	
}


bool MeshUtils::GetTriVerts( IMesh::VertexID nTri[3], IMesh::VertexID vID, IMesh::VertexID & other1, IMesh::VertexID & other2 )
{
	if ( nTri[0] == vID ) {
		other1 = nTri[1];
		other2 = nTri[2];
	} else if ( nTri[1] == vID ) {
		other1 = nTri[0];
		other2 = nTri[2];
	} else if ( nTri[2] == vID ) {
		other1 = nTri[0];
		other2 = nTri[1];
	} else
		return false;
	return true;
}








Wml::Vector3f MeshUtils::MeshLaplacian( VFTriangleMesh & mesh, IMesh::VertexID vID, std::vector<IMesh::VertexID> & vOneRing, std::vector<float> & vWeights )
{
	Wml::Vector3f vCenter, vNbr;
	mesh.GetVertex(vID, vCenter);
	size_t nNbrs = vOneRing.size();
	Wml::Vector3f vLaplacian(Wml::Vector3f::ZERO);
	for ( unsigned int i = 0; i < nNbrs; ++i ) {
		mesh.GetVertex(vOneRing[i], vNbr);
		vLaplacian += vWeights[i] * ( vNbr - vCenter );
	}
	return vLaplacian;
}


void MeshUtils::GetEdgeLengthStats( IMesh * pMesh, float & fMin, float & fMax, float & fAvg )
{
	fMax = std::numeric_limits<float>::min();
	fMin = std::numeric_limits<float>::max();
	fAvg = 0;
	
	int nCount = 0;
	IMesh::ITriIterator curt( pMesh->BeginITriangles()), endt( pMesh->EndITriangles());
	while ( curt != endt ) {
		IMesh::TriangleID nID = *curt;
		curt++;
		++nCount;

		Wml::Vector3f vVertices[3];
		pMesh->GetTriangle(nID, vVertices);

		for ( int i = 0; i < 3; ++i ) {
			float fLen = (vVertices[i] - vVertices[(i+1)%3]).Length();
			if ( fLen > fMax )
				fMax = fLen;
			if ( fLen < fMin )
				fMin = fLen;
			fAvg += fLen / 3.0f;
		}
	}
	fAvg /= (float)nCount;
}





bool MeshUtils::PickTriVerts( IMesh::VertexID nTri[3], IMesh::VertexID vOpp, IMesh::VertexID & vOther1, IMesh::VertexID & vOther2 )
{
	if ( nTri[0] == vOpp ) {
		vOther1 = nTri[1];
		vOther2 = nTri[2];
		return true;
	} else if ( nTri[1] == vOpp ) {
		vOther1 = nTri[2];
		vOther2 = nTri[0];
		return true;
	} else if ( nTri[2] == vOpp ) {
		vOther1 = nTri[0]; 
		vOther2 = nTri[1];
		return true;
	} else
		return false;
}

Wml::Vector3f MeshUtils::GetAverageNormal( VFTriangleMesh & mesh, IMesh::VertexID vID, unsigned int nKRings )
{
	Wml::Vector3f vVertex, vNormal, vAvgNormal(Wml::Vector3f::ZERO);
	float fWeightSum = 0.0f;

	std::vector<IMesh::VertexID> vOneRing;
	VertexOneRing(mesh, vID, vOneRing);
	size_t nCount = vOneRing.size();
	for ( unsigned int k = 0; k < nCount; ++k ) {
		mesh.GetVertex(vOneRing[k], vVertex, &vNormal);
		vAvgNormal += vNormal;
		fWeightSum += 1.0f;
	}	

	vAvgNormal *= 1.0f / fWeightSum;
	vAvgNormal.Normalize();

	return vAvgNormal;
}


void MeshUtils::SmoothNormals( VFTriangleMesh & mesh )
{
	VFTriangleMesh::vertex_iterator curv(mesh.BeginVertices()), endv(mesh.EndVertices());
	while ( curv != endv ) {
		IMesh::VertexID vID = *curv++;
		Wml::Vector3f vSmoothNorm = GetAverageNormal(mesh, vID);
		mesh.SetNormal(vID, vSmoothNorm);
	}
}



Wml::Vector3f MeshUtils::EstimateNormal( VFTriangleMesh & mesh, IMesh::VertexID vID, NormalEstMode eMode )
{
	Wml::Vector3f vEstimate(Wml::Vector3f::ZERO), vTri[3], vNormal;

	IMesh::VtxNbrItr itr(vID);
	mesh.BeginVtxTriangles(itr);
	IMesh::TriangleID tID = mesh.GetNextVtxTriangle(itr);
	while ( tID != IMesh::InvalidID ) {
		mesh.GetTriangle(tID, vTri);

		if ( eMode == AreaWeightedFaceAvg ) {
			float fWeight;
			vNormal = Normal(vTri[0], vTri[1], vTri[2], &fWeight );
			vEstimate += fWeight * vNormal;
		} else if ( eMode == UniformFaceAvg ) {
			vNormal = Normal(vTri[0], vTri[1], vTri[2]);
			vEstimate += vNormal;
		}

		tID = mesh.GetNextVtxTriangle(itr);
	}

	vEstimate.Normalize();
	return vEstimate;
}


void MeshUtils::EstimateNormals( VFTriangleMesh & mesh, NormalEstMode eMode, bool bSkipBoundary, Wml::Vector3f * pBuffer )
{
	VFTriangleMesh::vertex_iterator curv(mesh.BeginVertices()), endv(mesh.EndVertices());
	while ( curv != endv ) {
		IMesh::VertexID vID = *curv++;
		if ( bSkipBoundary && mesh.IsBoundaryVertex(vID) )
			continue;
		Wml::Vector3f vEstimate = EstimateNormal(mesh, vID, eMode);
		if ( pBuffer )
			pBuffer[vID] = vEstimate;
		else
			mesh.SetNormal(vID, vEstimate);
	}	
}



Wml::Vector3f MeshUtils::FaceNormal( VFTriangleMesh & mesh, IMesh::TriangleID tID )
{
	// [TODO] should pick most robust set of edges here...

	Wml::Vector3f vTri[3];
	mesh.GetTriangle(tID, vTri);
	Wml::Vector3f vEdge1(vTri[1] - vTri[0]);					vEdge1.Normalize();
	Wml::Vector3f vEdge2(vTri[2] - vTri[0]);					vEdge2.Normalize();
	Wml::Vector3f vFaceNormal( vEdge1.Cross(vEdge2) );
	vFaceNormal.Normalize();
	return vFaceNormal;
}


Wml::Vector3f MeshUtils::FaceInterpNormal( IMesh * pMesh, IMesh::TriangleID tID, const Wml::Vector3f & vTriPoint )
{
	Wml::Vector3f vTri[3], vNorm[3];
	pMesh->GetTriangle( tID, vTri, vNorm );
	return rms::InterpNormal(vTri[0],vTri[1],vTri[2], vNorm[0], vNorm[1], vNorm[2], vTriPoint);
}






void MeshUtils::ScaleMesh( VFTriangleMesh & mesh, const Wml::Vector3f & vScale, const Wml::Vector3f & vCenter )
{
	// TODO: handle normals properly...

	Wml::Vector3f vVertex, vNormal;
	VFTriangleMesh::vertex_iterator curv(mesh.BeginVertices()), endv(mesh.EndVertices());
	while ( curv != endv ) {
		IMesh::VertexID vID = *curv++;

		mesh.GetVertex( vID, vVertex, &vNormal );
		vVertex -= vCenter;
		vVertex.X() *= vScale.X();  vVertex.Y() *= vScale.Y();   vVertex.Z() *= vScale.Z();
		vVertex += vCenter;

		mesh.SetVertex( vID, vVertex, &vNormal );
	}
}


void MeshUtils::TranslateMesh( VFTriangleMesh & mesh, const Wml::Vector3f & vTranslate )
{
	Wml::Vector3f vVertex;
	VFTriangleMesh::vertex_iterator curv(mesh.BeginVertices()), endv(mesh.EndVertices());
	while ( curv != endv ) {
		IMesh::VertexID vID = *curv++;
		mesh.GetVertex( vID, vVertex );
		vVertex += vTranslate;
		mesh.SetVertex( vID, vVertex );
	}
}

void MeshUtils::RotateMesh( VFTriangleMesh & mesh, const Wml::Matrix3f & mRotate, const Wml::Vector3f & vCenter )
{
	Wml::Vector3f vVertex, vNormal;
	VFTriangleMesh::vertex_iterator curv(mesh.BeginVertices()), endv(mesh.EndVertices());
	while ( curv != endv ) {
		IMesh::VertexID vID = *curv++;

		mesh.GetVertex( vID, vVertex, &vNormal );
		vVertex -= vCenter;
		vVertex = mRotate * vVertex;
		vVertex += vCenter;
		vNormal = mRotate * vNormal;

		mesh.SetVertex( vID, vVertex, &vNormal );
	}
}




float MeshUtils::LoopLength( const VFTriangleMesh & mesh, const std::vector<IMesh::VertexID> & vLoop )
{
	float fLength = 0.0f;
	size_t nCount = vLoop.size();
	Wml::Vector3f v1,v2;
	mesh.GetVertex(vLoop[0], v1);
	for ( unsigned int i = 1; i < nCount; ++i ) {
		mesh.GetVertex(vLoop[(i+1)%nCount], v2);
		fLength += (v1-v2).Length();
		v1 = v2;
	}
	return fLength;
}
float MeshUtils::LoopLengthUV( const VFTriangleMesh & mesh, const std::vector<IMesh::VertexID> & vLoop, IMesh::UVSetID nSetID )
{
	float fLength = 0.0f;
	size_t nCount = vLoop.size();
	Wml::Vector2f v1,v2;
	if (! mesh.GetUV(vLoop[0], nSetID, v1) )
		return -1.0f;
	for ( unsigned int i = 1; i < nCount; ++i ) {
		if ( ! mesh.GetUV(vLoop[(i+1)%nCount], nSetID, v2) )
			return -1.0f;
		fLength += (v1-v2).Length();
		v1 = v2;
	}
	return fLength;
}



void MeshUtils::ScaleMeshByLoop( VFTriangleMesh & mesh, const std::vector<IMesh::VertexID> & vLoop, float fTargetPerimeter )
{
	float fLength = LoopLength(mesh, vLoop);
	float fScale = fTargetPerimeter / fLength;
	ScaleMesh(mesh, Wml::Vector3f(fScale,fScale,fScale));
}



float MeshUtils::MeshArea( VFTriangleMesh & mesh )
{
	float fSum = 0.0f;
	VFTriangleMesh::triangle_iterator curt(mesh.BeginTriangles()), endt(mesh.EndTriangles());
	while ( curt != endt ) {
		IMesh::TriangleID tID = *curt++;
		Wml::Vector3f vTri[3];
		mesh.GetTriangle(tID, vTri);
		fSum += Area( vTri[0], vTri[1], vTri[2] );
	}
	return fSum;
}

float MeshUtils::MeshUVArea( VFTriangleMesh & mesh, IMesh::UVSetID nUVSet )
{
	float fSum = 0.0f;
	VFTriangleMesh::triangle_iterator curt(mesh.BeginTriangles()), endt(mesh.EndTriangles());
	while ( curt != endt ) {
		IMesh::TriangleID tID = *curt++;
		Wml::Vector2f vTri[3];
		if ( mesh.GetTriangleUV(tID, nUVSet, vTri) )
			fSum += Area( vTri[0], vTri[1], vTri[2] );
	}
	return fSum;
}



void MeshUtils::ScaleMeshUV( VFTriangleMesh & mesh, IMesh::UVSetID nSetID, const Wml::Vector2f & vScale, const Wml::Vector2f & vCenter )
{
	if ( ! mesh.HasUVSet(nSetID) )
		return;

	IMesh::UVSet & uvset = mesh.GetUVSet(nSetID);

	Wml::Vector2f vUV;
	VFTriangleMesh::vertex_iterator curv(mesh.BeginVertices()), endv(mesh.EndVertices());
	while ( curv != endv ) {
		IMesh::VertexID vID = *curv++;

		if ( uvset.GetUV(vID, vUV) ) {
			vUV -= vCenter;
			vUV.X() *= vScale.X();  vUV.Y() *= vScale.Y();
			vUV += vCenter;
			uvset.SetUV(vID, vUV);
		}
	}
}









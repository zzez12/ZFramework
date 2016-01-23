// Ryan Schmidt
// Copyright (c) 2006-2010
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt

//#include "StdAfx.h"
#include ".\vfmeshrenderer.h"

#include <gl/GL.h>
#include <gl/GLU.h>

#include "VectorUtil.h"

using namespace rms;

VFMeshRenderer::VFMeshRenderer( VFTriangleMesh * pMesh )
{
	m_pMesh = pMesh;

	m_eNormalMode = FaceNormals;
	m_eColorMode = NoColors;
	m_eTexture2DMode = NoTexture2D;
	
	m_bDrawNormals = false;

}

VFMeshRenderer::~VFMeshRenderer(void)
{
}

void VFMeshRenderer::Render()
{
	IMesh::VertexID vTri[3];
	Wml::Vector3f vVtx[3], vNorm[3];
	Wml::Vector2f vTriUV[3];
	Wml::Vector2f vInvalidUV[3];
	vInvalidUV[0] = vInvalidUV[1] = vInvalidUV[2] = Wml::Vector2f(99999.0f, 99999.0f);
	Wml::Vector2f * pCurUV;
	Wml::ColorRGBA cColor;

	VFTriangleMesh::triangle_iterator 
		cur( m_pMesh->BeginTriangles() ),
		end( m_pMesh->EndTriangles() );

	glBegin(GL_TRIANGLES);
	while ( cur != end ) {
		IMesh::TriangleID nTriID = *cur;	cur++;

		m_pMesh->GetTriangle( nTriID, vTri );
		m_pMesh->GetTriangle( nTriID, vVtx, vNorm );

		if ( m_eTexture2DMode == VertexTexture2D || m_eTexture2DMode == VertexTexture2D_Required ) {
			pCurUV = ( m_pMesh->GetTriangleUV( nTriID, 0, vTriUV )) ? vTriUV : vInvalidUV;
		}
		if ( m_eTexture2DMode == VertexTexture2D_Required && pCurUV == vInvalidUV )
			continue;

		if ( m_eNormalMode == FaceNormals ) 
			glNormal3fv( rms::Normal(vVtx[0], vVtx[1], vVtx[2]) );

		for ( int j = 0 ; j < 3; ++j ) {
			if ( m_eNormalMode == VertexNormals )
				glNormal3fv( vNorm[j] );
			if ( m_eColorMode == VertexColors ) {
				m_pMesh->GetColor( vTri[j], cColor );
				glColor4fv( cColor );
			}
			if ( m_eTexture2DMode == VertexTexture2D || m_eTexture2DMode == VertexTexture2D_Required )
				glTexCoord2fv( pCurUV[j] );
			glVertex3fv( vVtx[j] );
		}

	}
	glEnd();


	if ( m_bDrawNormals )
		RenderNormals();
}

void VFMeshRenderer::RenderNormals()
{
	float fScaleFactor = 0.025f;

	Wml::Vector3f vVtx[3], vNorm[3];

	glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT);
	glDisable(GL_LIGHTING);
	glLineWidth(2.0f);

	glColor4f(0.5f, 0.5f, 0.5f, 1.0f);

	glBegin(GL_LINES);
	VFTriangleMesh::triangle_iterator 
		cur( m_pMesh->BeginTriangles() ),
		end( m_pMesh->EndTriangles() );
	while ( cur != end ) {
		IMesh::TriangleID nTriID = *cur;	cur++;
		m_pMesh->GetTriangle( nTriID, vVtx, vNorm );

		for ( int j = 0 ; j < 3; ++j ) {
			glVertex3fv( vVtx[j] );
			glVertex3fv( vVtx[j] + fScaleFactor * vNorm[j] );
		}

	}
	glEnd();

	glPopAttrib();
}





void VFMeshRenderer::Render_UV()
{
	Wml::Vector2f vTriUV[3];
	Wml::Vector2f vInvalidUV[3];
	vInvalidUV[0] = vInvalidUV[1] = vInvalidUV[2] = Wml::Vector2f(99999.0f, 99999.0f);
	Wml::Vector2f * pCurUV;

	VFTriangleMesh::triangle_iterator 
		cur( m_pMesh->BeginTriangles() ),
		end( m_pMesh->EndTriangles() );

	glBegin(GL_TRIANGLES);
	while ( cur != end ) {
		IMesh::TriangleID nTriID = *cur;	cur++;

		pCurUV = ( m_pMesh->GetTriangleUV( nTriID, 0, vTriUV )) ? vTriUV : vInvalidUV;
		if ( pCurUV == vInvalidUV )
			continue;

		for ( int j = 0 ; j < 3; ++j ) {
			glTexCoord2fv( pCurUV[j] );
			glVertex2fv( pCurUV[j] );
		}

	}
	glEnd();
}

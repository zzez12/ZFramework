// Ryan Schmidt
// Copyright (c) 2006-2010
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt

#ifndef __RMS_VFMESH_RENDERER_H__
#define __RMS_VFMESH_RENDERER_H__

#include "VFTriangleMesh.h"

namespace rms {

class VFMeshRenderer
{
public:
	VFMeshRenderer( VFTriangleMesh * pMesh );
	~VFMeshRenderer(void);

	void Render();
	void Render_UV();

	enum NormalMode {
		NoNormals,
		VertexNormals,
		FaceNormals
	};
	NormalMode GetNormalMode() { return m_eNormalMode; }
	void SetNormalMode( NormalMode eMode ) { m_eNormalMode = eMode; }

	enum ColorMode {
		NoColors,
		VertexColors
	};
	ColorMode GetColorMode() { return m_eColorMode; }
	void SetColorMode( ColorMode eMode ) { m_eColorMode = eMode; }


	enum Texture2DMode {
		NoTexture2D,
		VertexTexture2D,
		FaceTexture2D,
		VertexTexture2D_Required
	};
	Texture2DMode GetTexture2DMode() { return m_eTexture2DMode; }
	void SetTexture2DMode( Texture2DMode eMode ) { m_eTexture2DMode = eMode; }

	void SetModes( NormalMode eNormalMode, Texture2DMode eTex2DMode ) 
		{ m_eNormalMode = eNormalMode; m_eTexture2DMode = eTex2DMode; }

	void SetDrawNormals( bool bDrawNormals ) { m_bDrawNormals = bDrawNormals; }
	bool GetDrawNormals() { return m_bDrawNormals; }

protected:
	VFTriangleMesh * m_pMesh;

	NormalMode m_eNormalMode;
	ColorMode m_eColorMode;
	Texture2DMode m_eTexture2DMode;

	bool m_bDrawNormals;
	void RenderNormals();
};

} // end namespace rmsmesh


#endif // __RMS_VFMESH_RENDERER_H__
#pragma once

#include "../ZFramework3D/GlobalDefs.h"
#include "ZRenderer.h"
#include "ZMeshAlgorithms.h"
#include "ZPlane.h"
#include "RenderOption.h"

namespace ZMeshSpace
{
	class Mesh3D;
}

class MeshRenderer : public ZRenderer
{
public:

	MeshRenderer();

	void draw();
	void setRenderType(int type);
	void setRenderType(const std::string& type);

	void setRenderWireFrame(bool b) {renderOpt_.show_wire_frames=b;}
	RenderOption& getRenderOption() {return renderOpt_;}

private:
	//MeshRenderType drawType_;
	RenderOption renderOpt_;

private:
	void beginToDraw();
	void endDraw();

	void drawMesh(ZMeshSpace::Mesh3D *mesh);
	void drawIsoPlanes(const ZMeshSpace::ZIsoLineContainer& isoLineContainer, int drawWhich);

	void drawSolidFlat(ZMeshSpace::Mesh3D *mesh);
	void drawSolidSmooth(ZMeshSpace::Mesh3D *mesh);
	void drawPoints(ZMeshSpace::Mesh3D *mesh);
	void drawWireFrame(ZMeshSpace::Mesh3D *mesh);
	void drawMeshByVerticeColors(ZMeshSpace::Mesh3D *mesh);
	void drawMeshByFaceColors(ZMeshSpace::Mesh3D *mesh);

	void drawIsoLines(ZMeshSpace::ZMeshAlgorithms *algo);
	void drawIsoLine(const ZMeshSpace::ZEigenIsoLine& isoLine);

	void drawPlanePoints(const ZMeshSpace::ZIsoLineContainer& isoLineContainer, int drawWhich);
	void drawPlanes(const ZMeshSpace::ZIsoLineContainer& isoLineContainer, int drawWhich);
	void drawPlaneLines(const ZMeshSpace::ZIsoLineContainer& isoLineContainer, int drawWhich);
	void drawPlaneLines_ByClusterId(const ZMeshSpace::ZIsoLineContainer& isoLineContainer, int drawWhich);
	void drawPlanes_PlaneLists(const ZMeshSpace::ZIsoLineContainer& isoLineContainer, int which);
	void drawPlanes_PlaneList(const ZMeshSpace::ZIsoLineContainer& isoLineContainer, const std::vector<int>& ids, float * color);

 	void drawPlanePoints(const ZMeshSpace::ZPlane& plane);
	void drawPlane(const ZMeshSpace::ZPlane& plane, float* color);
	void drawPlaneLines(const ZMeshSpace::ZPlane& plane, float* color);

	void drawPlaneByPolarCoord(const ZMeshSpace::ZPlane& plane, float* color);

	void drawMeshPointByType(ZMeshSpace::Mesh3D *mesh);
	void drawMSEdges(ZMeshSpace::ZQuadrangulation* quadAlgo);

	void drawFaceNormals(ZMeshSpace::Mesh3D *mesh);
};
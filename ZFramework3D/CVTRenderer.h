#pragma once
#ifndef ZCVT_RENDERER_H_
#define ZCVT_RENDERER_H_

#include "ZCVTAlgorithms.h"
#include "RenderOption.h"
#include "ZRenderer.h"
#include "CGAL_Algorithms.h"

class CVTRenderer : public ZRenderer
{
public:
	CVTRenderer();

	void draw();
	void setRenderType(int type);
	void setRenderType(const std::string& type);

	RenderOption& getRenderOption() {return renderOpt_;}

private:
	RenderOption renderOpt_;

private:
	void beginToDraw();
	void endDraw();

	void drawLines(ZCVT::ZCVTData* data);
	void drawCVT(CGAL_Algorithms::Delaunay_triangulation_2& dt);
	void drawCVT(ZCVT::ZCVTData* data);
};

#endif//ZCVT_RENDERER_H_
#include "CVTRenderer.h"
#include "ZDataManager.h"
#include "Color.h"
#include "qgl.h"
#include <gl/glut.h>

#ifndef _RTYPE_
#define _RTYPE_(a) (std::string)(render_type(a))
#endif

CVTRenderer::CVTRenderer()
{
	setRenderType(RENDER_TYPE_CVT_POINTS);
}

void CVTRenderer::setRenderType(int type)
{
	renderOpt_.strRenderMethod = render_type(type);
}

void CVTRenderer::setRenderType(const std::string& type)
{
	// add a temp code for testing pull-request and gerrit
	renderOpt_.strRenderMethod = type;
}

void CVTRenderer::beginToDraw()
{
	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);
}

void CVTRenderer::endDraw()
{
	glEnable(GL_DEPTH_TEST);
}

void CVTRenderer::draw()
{
	ZCVT::ZCVTData* data = ZDataManager::getDataManager()->getCVTData();
	ZCVT::ZCVTAlgorithms* algo = ZDataManager::getDataManager()->getCVTAlgorithmHandler();
	if (data==NULL) return;
	
	if (renderOpt_.show_cvt_init_samples)
	{
		drawLines(data);
	}

	if (renderOpt_.show_cvt_delaunay || renderOpt_.show_cvt_voronoi)
	{
		//drawCVT(algo->getDelaunay());
		//drawCVT(algo->getVoronoi());
		drawCVT(data);
	}

	// draw constrian bounding Box
	if (renderOpt_.show_cvt_clippingRegion)
	{
		glColor4f(1.f, 1.f, 1.f, 0.8f);
		glBegin(GL_LINE_LOOP);
		glVertex2f(0, 0);
		glVertex2f(1, 0);
		glVertex2f(1, 1);
		glVertex2f(0, 1);
		glEnd();
	}
}

void CVTRenderer::drawLines(ZCVT::ZCVTData* data)
{
	std::vector<ZCVT::ZLine>& lines = data->getLines();
	glEnable(GL_LINE_STIPPLE);
	glLineStipple(1, 0x2222);
	glColor3f(1.f, 1.f, 1.f);
	glBegin(GL_LINES);
	for (int i=0; i<lines.size(); i++)
	{
		ZCVT::ZLine& line = lines[i];
		glVertex2fv(line.p0);
		glVertex2fv(line.p1);
	}
	glEnd();
	glDisable(GL_LINE_STIPPLE);

	glPointSize(4.f);
	glBegin(GL_POINTS);
	int subSize = data->getSubSize();
	for (int i=0; i<lines.size(); i++)
	{
		ZCVT::ZLine& line = lines[i];
		Vec2f step = (line.p1-line.p0)/subSize;
		for (int j=0; j<=subSize; j++)
		{
			glVertex2fv(line.p0 + step*j);
		}
	//	glVertex2fv(line.p0);
	//	glVertex2fv(line.p1);
	}
	glEnd();
	glPointSize(1.f);
}

void CVTRenderer::drawCVT(CGAL_Algorithms::Delaunay_triangulation_2& dt)
{
	if (renderOpt_.show_cvt_delaunay)
	{
		glColor3f(0.f, 0.f, 1.f);
		CGAL_Algorithms::Delaunay_triangulation_2::Finite_faces_iterator fit;
		for (fit=dt.finite_faces_begin(); fit!=dt.finite_faces_end(); fit++)
		{
			glBegin(GL_LINE_LOOP);
			glVertex2f(fit->vertex(0)->point().hx(), fit->vertex(0)->point().hy());
			glVertex2f(fit->vertex(1)->point().hx(), fit->vertex(1)->point().hy());
			glVertex2f(fit->vertex(2)->point().hx(), fit->vertex(2)->point().hy());
			glEnd();
		}
	}

	if (renderOpt_.show_cvt_voronoi)
	{
		glColor3f(0.f, 1.f, 0.f);
		glEnable(GL_LINE_STIPPLE);
		glLineStipple(1, 0x3333);
		CGAL_Algorithms::Delaunay_triangulation_2::Finite_edges_iterator eit;
		for (eit=dt.finite_edges_begin(); eit!=dt.finite_edges_end(); eit++)
		{
			CGAL::Object o = dt.dual(eit);
			if (CGAL::object_cast<CGAL_Algorithms::Segment_2>(&o))	// Line
			{
				glBegin(GL_LINES);
				glVertex2f( CGAL::object_cast<CGAL_Algorithms::Segment_2>(&o)->source().hx(),
					CGAL::object_cast<CGAL_Algorithms::Segment_2>(&o)->source().hy());
				glVertex2f( CGAL::object_cast<CGAL_Algorithms::Segment_2>(&o)->target().hx(),
					CGAL::object_cast<CGAL_Algorithms::Segment_2>(&o)->target().hy());
				glEnd();
			}
			else if (CGAL::object_cast<CGAL_Algorithms::Ray_2>(&o))	// ray
			{
				glBegin(GL_LINES);
				glVertex2f( CGAL::object_cast<CGAL_Algorithms::Ray_2>(&o)->source().hx(),
					CGAL::object_cast<CGAL_Algorithms::Ray_2>(&o)->source().hy());
				glVertex2f( CGAL::object_cast<CGAL_Algorithms::Ray_2>(&o)->point(1).hx(),
					CGAL::object_cast<CGAL_Algorithms::Ray_2>(&o)->point(1).hy());
				glEnd();
			}
		}
		glDisable(GL_LINE_STIPPLE);
	}

}

void CVTRenderer::drawCVT(ZCVT::ZCVTData* data)
{
	ZCVT::ZVoronoi* vor = data->getVoronoi();
	if (renderOpt_.show_cvt_vertices)
	{
		glColor3f(0.6f, 0.6f, 0.6f);
		glPointSize(4.f);
		glBegin(GL_POINTS);
		for (int i=0; i<vor->voronoiVerticeSize(); i++)
		{
			glVertex2fv(vor->getVertex(i));
		}
		glEnd();
		glPointSize(1.f);
	}

	if (renderOpt_.show_cvt_delaunay)
	{
	}

	if (renderOpt_.show_cvt_voronoi)
	//if (renderOpt_.show_cvt_delaunay)
	{
		glColor3f(1.f, 0.f, 0.f);
		glEnable(GL_LINE_STIPPLE);
		glLineStipple(1, 0x3333);
		for (int i=0; i<vor->voronoiCellSize(); i++)
		{
			ZCVT::ZVoronoiCell& cell = vor->cell(i);
			glBegin(GL_LINE_LOOP);
			for (int j=0; j<cell.vertices_.size(); j++)
			{
				int vId = cell.vertices_[j];
				glVertex2fv(vor->getVertex(vId));
			}
			glEnd();
		}
		glDisable(GL_LINE_STIPPLE);
	}

	// draw voronoi centers
	if (renderOpt_.show_cvt_voronoi_vertices)
	{
		glColor3f(0.6f, 0.6f, 0.f);
		glPointSize(2.f);
		glBegin(GL_POINTS);
		for (int i=0; i<vor->voronoiCellSize(); i++)
		{
			glVertex2fv(vor->cell(i).center_);
		}
		glEnd();
		glPointSize(1.f);
	}
}

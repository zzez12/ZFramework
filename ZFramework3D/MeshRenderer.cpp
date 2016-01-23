#include "MeshRenderer.h"
#include "ZDataManager.h"
#include "qgl.h"
#include <gl/glut.h>

using namespace ZMeshSpace;

#ifndef _RTYPE_
#define _RTYPE_(a) (std::string)(render_type(a))
#endif

MeshRenderer::MeshRenderer()
{
	setRenderType(RENDER_TYPE_SOLID_SMOOTH);
}

void MeshRenderer::setRenderType(int type)
{
	renderOpt_.strRenderMethod = render_type(type);
}

void MeshRenderer::setRenderType(const std::string& type)
{
	renderOpt_.strRenderMethod = type;
}

void MeshRenderer::beginToDraw()
{
	if (renderOpt_.translucent)
	{
		glEnable(GL_BLEND);
		//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glBlendFunc(GL_ONE, GL_ONE);
		glDisable(GL_DEPTH_TEST);

	}
}

void MeshRenderer::endDraw()
{
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
}

void MeshRenderer::draw()
{
	ZMeshSpace::Mesh3D* mesh = ZDataManager::getDataManager()->getMesh();
	ZMeshSpace::ZMeshAlgorithms* algo = ZDataManager::getDataManager()->getAlgorithmHandler();
	if (!mesh) return;

	beginToDraw();
	drawMesh(mesh);
	if (renderOpt_.show_iso_plane) 
	{
		int which = -1;
		if (renderOpt_.show_iso_line_by_id)
		{
			which = algo->getCurrentIsoLineIdx();
		}
		drawIsoPlanes(algo->getIsoLineContainer(), which);

	}
	if (renderOpt_.show_iso_plane_plane_list)
	{
		int which = -1;
		if (renderOpt_.show_iso_line_by_id)
		{
			which = algo->getCurrentIsoLineIdx();
		}
		drawPlanes_PlaneLists(algo->getIsoLineContainer(), which);
	}

	if (renderOpt_.show_wire_frames)
	{
		drawWireFrame(mesh);
	}
	if (renderOpt_.show_iso_lines)
	{
		drawIsoLines(algo);
	}
	if (renderOpt_.show_vertex_tags)
	{
		drawMeshPointByType(mesh);
		drawMSEdges(algo->getQuadrangulation());
	}
	if (renderOpt_.show_face_normals)
	{
		drawFaceNormals(mesh);
	}
	//if (renderOpt_.show_iso_plane)
	//{
	//	drawPlane(algo);
	//}
	//if (renderOpt_.show_iso_plane_points)
	//{
	//	drawPlanePoints(algo);
	//}
	endDraw();
}

void MeshRenderer::drawMesh(ZMeshSpace::Mesh3D *mesh)
{
	if (renderOpt_.strRenderMethod == _RTYPE_(RENDER_TYPE_NONE))
	{
		;
	}
	else if (renderOpt_.strRenderMethod == _RTYPE_(RENDER_TYPE_POINTS))
	{
		drawPoints(mesh);
	}
	else if (renderOpt_.strRenderMethod == _RTYPE_(RENDER_TYPE_SOLID_FLAT))
	{
		drawSolidFlat(mesh);
	}
	else if (renderOpt_.strRenderMethod == _RTYPE_(RENDER_TYPE_SOLID_SMOOTH))
	{
		drawSolidSmooth(mesh);
	}
	else if (renderOpt_.strRenderMethod == _RTYPE_(RENDER_TYPE_VERTEX_COLOR))
	{
		drawMeshByVerticeColors(mesh);
	}
	else if (renderOpt_.strRenderMethod == _RTYPE_(RENDER_TYPE_FACE_COLOR))
	{
		drawMeshByFaceColors(mesh);
	}
}

void MeshRenderer::drawIsoPlanes(const ZMeshSpace::ZIsoLineContainer& isoLineContainer, int which)
{
	if (renderOpt_.strRenderIsoPlaneMethod == _RTYPE_(RENDER_TYPE_IP_POINTS))
	{
		drawPlanePoints(isoLineContainer, which);
	}
	else if (renderOpt_.strRenderIsoPlaneMethod == _RTYPE_(RENDER_TYPE_IP_LINES))
	{
		drawPlaneLines(isoLineContainer, which);
	}
	else if (renderOpt_.strRenderIsoPlaneMethod == _RTYPE_(RENDER_TYPE_IP_PLANES))
	{
		drawPlanes(isoLineContainer, which);
	}
	else if (renderOpt_.strRenderIsoPlaneMethod == _RTYPE_(RENDER_TYPE_IP_LINES_BY_CLUSTERID))
	{
		drawPlaneLines_ByClusterId(isoLineContainer, which);
	}
}

void MeshRenderer::drawPoints(ZMeshSpace::Mesh3D *mesh)
{

}

void MeshRenderer::drawSolidFlat(ZMeshSpace::Mesh3D *mesh)
{
	glEnable(GL_LIGHTING);
	glShadeModel(GL_FLAT);
	for (int i=0; i<mesh->get_num_of_faces_list(); i++)
	{
		HE_face* face = mesh->get_face(i);
		HE_edge* edge = face->m_pedge;
		glBegin(GL_POLYGON);
		do 
		{
			HE_vert* v = edge->m_pvert;
			glColor4fv(v->m_color);
			glNormal3fv(v->m_vnormal);
			glVertex3fv(v->m_vpos);
			edge = edge->m_pnext;
		} while (edge!=face->m_pedge);
		glEnd();
	}
}

void MeshRenderer::drawSolidSmooth(ZMeshSpace::Mesh3D *mesh)
{
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	for (int i=0; i<mesh->get_num_of_faces_list(); i++)
	{
		HE_face* face = mesh->get_face(i);
		HE_edge* edge = face->m_pedge;
		glBegin(GL_POLYGON);
		do 
		{
			HE_vert* v = edge->m_pvert;
			glColor4fv(v->m_color);
			glNormal3fv(v->m_vnormal);
			glVertex3fv(v->m_vpos);
			edge = edge->m_pnext;
		} while (edge!=face->m_pedge);
		glEnd();
	}
}

void MeshRenderer::drawMeshByVerticeColors(ZMeshSpace::Mesh3D *mesh)
{
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	//glShadeModel(GL_FLAT);
	for (int i=0; i<mesh->get_num_of_faces_list(); i++)
	{
		HE_face* face = mesh->get_face(i);
		HE_edge* edge = face->m_pedge;
		glBegin(GL_POLYGON);
		do 
		{
			HE_vert* v = edge->m_pvert;
			float rgb[4];
			colorCoding(v->m_colorValue, rgb);
			rgb[3] = 0.7f;
			glColor3fv(rgb);
			glNormal3fv(v->m_vnormal);
			glVertex3fv(v->m_vpos);
			edge = edge->m_pnext;
		} while (edge!=face->m_pedge);
		glEnd();
	}

}

void MeshRenderer::drawMeshByFaceColors(ZMeshSpace::Mesh3D *mesh)
{
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	glShadeModel(GL_FLAT);
	for (int i=0; i<mesh->get_num_of_faces_list(); i++)
	{
		HE_face* face = mesh->get_face(i);
		HE_edge* edge = face->m_pedge;
		glBegin(GL_POLYGON);
		float rgb[4];
		colorCoding(face->m_value, rgb);
		do 
		{
			HE_vert* v = edge->m_pvert;
			rgb[3] = 0.7f;
			glColor3fv(rgb);
			glNormal3fv(v->m_vnormal);
			glVertex3fv(v->m_vpos);
			edge = edge->m_pnext;
		} while (edge!=face->m_pedge);
		glEnd();
	}
}

void MeshRenderer::drawIsoLine(const ZMeshSpace::ZEigenIsoLine& isoLine)
{
	glLineWidth(3.f);
	//glColor3f(1.f, 0.f, 0.f);
	glBegin(GL_LINE_STRIP);
	int size = isoLine.linePoints.size();
	for (int j=0; j<isoLine.linePoints.size(); j++)
	{
		const ZMeshSpace::ZMeshSurfacePoint& p = isoLine.linePoints[j];
		float color[4];
		colorCoding(1.f*j/size, color);
		color[3] = isoLine.isGood_ ? 1.f : 0.3f;
		glColor4fv(color);
		glVertex3fv(p.pos);
	}
	glEnd();
}

void MeshRenderer::drawIsoLines(ZMeshSpace::ZMeshAlgorithms *algo)
{
	if (algo->getIsoLineSize()==0)
		return;
	if (algo->getCurrentIsoLineIdx()==-1)
	{
		for (int i=0; i<algo->getIsoLineSize(); i++)
		{
			const ZMeshSpace::ZEigenIsoLine& line = algo->getIsoLine(i);
			drawIsoLine(line);
		}
	}
	else
	{
		const ZMeshSpace::ZEigenIsoLine& line = algo->getIsoLine(algo->getCurrentIsoLineIdx());
		drawIsoLine(line);
	}
}

void MeshRenderer::drawWireFrame(ZMeshSpace::Mesh3D *mesh)
{
	glDisable(GL_LIGHTING);
	glLineWidth(1.f);
	for (int i=0; i<mesh->get_num_of_faces_list(); i++)
	{
		HE_face* f = mesh->get_face(i);
		glBegin(GL_LINE_STRIP);
		HE_edge* e = f->m_pedge;
		do 
		{
			//glColor4fv(e->m_pvert->m_color);
			glColor4fv(CColor::gray);
			glVertex3fv(e->m_pvert->m_vpos);
			e = e->m_pnext;
		} while (e!=f->m_pedge);
		glEnd();
	}
}

void MeshRenderer::drawPlanePoints(const ZMeshSpace::ZIsoLineContainer& isoLineContainer, int which)
{
	if (which==-1)
	{
		for (int i=0; i<isoLineContainer.getIsoLineSize(); i++)
		{
			drawPlanePoints(isoLineContainer.getIsoLine(i).plane);
		}
	}
	else
	{
		drawPlanePoints(isoLineContainer.getIsoLine(which).plane);
	}

}

void MeshRenderer::drawPlaneLines(const ZMeshSpace::ZIsoLineContainer& isoLineContainer, int which)
{
	int size = isoLineContainer.getIsoLineSize();
	glDisable(GL_LIGHTING);
	if (which==-1)
	{
		for (int i=0; i<size; i++)
		{
			const ZMeshSpace::IsoLineData& isoLineData = isoLineContainer.getIsoLine(i);
			float rgb[3];
			colorCoding(1.f*i/size, rgb);
			drawPlaneLines(isoLineData.plane, rgb);
		}
	}
	else
	{
		const ZMeshSpace::IsoLineData& isoLineData = isoLineContainer.getIsoLine(which);
		float rgb[3];
		colorCoding(1.f*which/size, rgb);
		drawPlaneLines(isoLineData.plane, rgb);
		drawPlaneByPolarCoord(isoLineData.plane, rgb);
		// draw plane norm
		Vec3f c = isoLineData.plane.getCenter();
		Vec3f n = isoLineData.plane.getNorm();
		glBegin(GL_LINES);
		glColor4f(1.f, 0.f, 0.f, 0.8f);
		glVertex3fv(c);
		glVertex3fv((c+n*0.3));
		glEnd();
	}
	glEnable(GL_LIGHTING);
}

void MeshRenderer::drawPlaneLines_ByClusterId(const ZMeshSpace::ZIsoLineContainer& isoLineContainer, int clusterId)
{
	int size = isoLineContainer.getIsoLineSize();
	glDisable(GL_LIGHTING);
	if (clusterId==-1)
	{
		for (int i=0; i<size; i++)
		{
			const ZMeshSpace::IsoLineData& isoLineData = isoLineContainer.getIsoLine(i);
			float rgb[3];
			colorCoding(1.f*isoLineData.clusterId/isoLineContainer.getKmeansK(), rgb);
			drawPlaneLines(isoLineData.plane, rgb);
			//drawPlaneByPolarCoord(isoLineData.plane, rgb);
		}
		if (renderOpt_.show_kmeans_plane)
		{
			for (int i=0; i<isoLineContainer.getKmeansK(); i++)
			{
				float rgb[3];
				colorCoding(1.f*i/isoLineContainer.getKmeansK(), rgb);
				if (renderOpt_.show_kmeans_center_plane)
					drawPlaneByPolarCoord(isoLineContainer.getKmeansPlanes()[i], rgb);
			}
		}

	}
	else
	{
		for (int i=0; i<size; i++)
		{
			const ZMeshSpace::IsoLineData& isoLineData = isoLineContainer.getIsoLine(i);
			if (isoLineData.clusterId!=clusterId) continue;

			float rgb[3];
			colorCoding(1.f*isoLineData.clusterId/isoLineContainer.getKmeansK(), rgb);
			drawPlaneLines(isoLineData.plane, NULL);
		}	
		if (renderOpt_.show_kmeans_plane)
		{
			float rgb[3];
			colorCoding(1.f*clusterId/isoLineContainer.getKmeansK(), rgb);
			if (renderOpt_.show_kmeans_center_plane)
				drawPlaneByPolarCoord(isoLineContainer.getKmeansPlanes()[clusterId], rgb);
		}

	}
	glEnable(GL_LIGHTING);
}

void MeshRenderer::drawPlanes(const ZMeshSpace::ZIsoLineContainer& isoLineContainer, int which)
{
	int size = isoLineContainer.getIsoLineSize();
	glDisable(GL_LIGHTING);
	if (which==-1)
	{
		for (int i=0; i<size; i++)
		{
			float rgb[3];
			colorCoding(1.f*i/size, rgb);
			drawPlane(isoLineContainer.getIsoLine(i).plane, rgb);
		}
	}
	else
	{
		float rgb[3];
		colorCoding(1.f*which/size, rgb);
		drawPlane(isoLineContainer.getIsoLine(which).plane, rgb);
	}

	glEnable(GL_LIGHTING);
}

void MeshRenderer::drawPlanePoints(const ZMeshSpace::ZPlane& plane)
{
	bool bUseSphere = false;
	if (bUseSphere)
	{
		glEnable(GL_LIGHTING);
		float ballRad = plane.getAverageEdgeLength()*0.6;
		GLUquadricObj *qobj = gluNewQuadric();
		gluQuadricDrawStyle(qobj, GLU_FILL);
		gluQuadricNormals(qobj, GLU_SMOOTH);
		const std::vector<Vec3f>& points = plane.getProjectedPoints();
		for (int i=0; i<points.size(); i++)
		{
			const Vec3f& p = points[i];
			float rgb[3];
			colorCoding(1.f*i/points.size(), rgb);
			glColor3fv(rgb);
			glTranslatef(p.x, p.y, p.z);
			gluSphere(qobj, ballRad, 5, 5);
			glTranslatef(-p.x, -p.y, -p.z);
		}
		gluDeleteQuadric(qobj);	
	}
	else
	{
		glDisable(GL_LIGHTING);
		glPointSize(4.0);
		const std::vector<Vec3f>& points = plane.getProjectedPoints();
		glBegin(GL_POINTS);
		for (int i=0; i<points.size(); i++)
		{
			const Vec3f& p = points[i];
			float rgb[3];
			colorCoding(1.f*i/points.size(), rgb);
			glColor3fv(rgb);
			glVertex3fv(p);
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}																								
}

void MeshRenderer::drawPlaneLines(const ZMeshSpace::ZPlane& plane, float *color)
{
	const std::vector<Vec3f>& points = plane.getProjectedPoints();
	if (color!=NULL)
	{
		glLineWidth(renderOpt_.line_width);
		glBegin(GL_LINE_STRIP);
		for (int i=0; i<points.size(); i++)
		{
			const Vec3f& p = points[i];
			//float rgb[3];
			//colorCoding(1.f*i/points.size(), rgb);
			glColor3fv(color);
			glVertex3fv(p);
		}
		glEnd();
	}
	else
	{
		float rgb[3];
		colorCoding(plane.colorValue_, rgb);
		glLineWidth(renderOpt_.line_width);
		glBegin(GL_LINE_STRIP);
		for (int i=0; i<points.size(); i++)
		{
			const Vec3f& p = points[i];
			glColor3fv(rgb);
			glVertex3fv(p);
		}
		glEnd();
	}

}

void MeshRenderer::drawPlane(const ZMeshSpace::ZPlane& plane, float *color)
{
	const std::vector<Vec3f>& points = plane.getProjectedPoints();
	glBegin(GL_POLYGON);
	for (int i=0; i<points.size(); i++)
	{
		const Vec3f& p = points[i];
		glColor3fv(color);
		glVertex3fv(p);
	}
	glEnd();
}

void MeshRenderer::drawPlaneByPolarCoord(const ZMeshSpace::ZPlane& plane, float* color)
{
// 	if (plane.getPlaneType()!=ZMeshSpace::PlaneType_Polar)
// 		return;
	float col[4];
	col[0] = color[0]; col[1] = color[1]; col[2] = color[2];
	col[3] = 0.7f;
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_ONE, GL_ONE);
	//glDisable(GL_DEPTH_TEST);
	// 	std::cout << std::endl;
// 	std::cout << ZConsoleTools::blue;
// 	std::cout << "P:" << std::endl;
// 	std::cout << "n:" << plane.getNorm() << " --l:" << plane.getNorm().length() << std::endl;
// 	std::cout << "np: " << plane.getPolarCoord() << std::endl;
	glBegin(GL_POLYGON);
	for (int i=0; i<4; i++)
	{
		glColor4fv(col);
		glVertex3fv(plane.polarCoordCorners_[i]);
//		std::cout << plane.polarCoordCorners_[i] << std::endl;
	}
//	std::cout << ZConsoleTools::white;
	glEnd();
	//glDisable(GL_BLEND);
}

void MeshRenderer::drawPlanes_PlaneLists(const ZMeshSpace::ZIsoLineContainer& isoLineContainer, int which)
{
	glDisable(GL_LIGHTING);
	const std::vector<std::vector<int>>& idxx = isoLineContainer.getPlaneLists();
	float rgb[3];
	if (which==-1)
	{
		for (int i=0; i<idxx.size(); i++)
		{
			colorCoding(1.f*i/idxx.size(), rgb);
			drawPlanes_PlaneList(isoLineContainer, idxx[i], rgb);
		}		
	}
	else if (which<idxx.size())
	{
		colorCoding(1.f*which/idxx.size(), rgb);
		drawPlanes_PlaneList(isoLineContainer, idxx[which], NULL);
	}
	glEnable(GL_LIGHTING);
}

void MeshRenderer::drawPlanes_PlaneList(const ZMeshSpace::ZIsoLineContainer& isoLineContainer, const std::vector<int>& idx, float *rgb)
{
	for (int i=0; i<idx.size(); i++)
	{
		const ZMeshSpace::IsoLineData& isoLineData = isoLineContainer.getIsoLine(idx[i]);
		float col[3];
		if (rgb==NULL) colorCoding(1.f*i/idx.size(), col);
		drawPlaneLines(isoLineData.plane, rgb==NULL ? col : rgb);
	}
}

void MeshRenderer::drawMeshPointByType(ZMeshSpace::Mesh3D *mesh)
{
	glDisable(GL_LIGHTING);
	glPointSize(4.0);
	glBegin(GL_POINTS);
	for (int i=0; i<mesh->get_num_of_vertex_list(); i++)
	{
		HE_vert* v = mesh->get_vertex(i);
		if (v->m_flag==4) continue;
		const Vec3f& p = v->m_vpos;
		float rgb[3];
		colorCodingInt(v->m_flag, rgb);
		//std::cout << rgb[0] << " " << rgb[1] << " " << rgb[2] << "\n";
		glColor3fv(rgb);
		glVertex3fv(p);
	}
	glEnd();
	glEnable(GL_LIGHTING);
}

void MeshRenderer::drawMSEdges(ZMeshSpace::ZQuadrangulation* quadAlgo)
{
	glDisable(GL_LIGHTING);

	ZMeshSpace::Mesh3D* mesh = quadAlgo->getMesh();
	int eSize = quadAlgo->getMSEdgeSize();
	float rgb[3];
	CColor(0.8, 0.5, 0.8).fillColor(rgb);
	glColor3fv(rgb);
	glLineWidth(2.0);
	for (int i=0; i<eSize; i++)
	{
		ZMeshSpace::ZQuadrangulation::MeshMSEdge* edge = quadAlgo->getMSEdge(i);
		Vec3f v1 = mesh->get_vertex(edge->vert1_->id_)->m_vpos;
		Vec3f v2 = mesh->get_vertex(edge->vert2_->id_)->m_vpos;
		glBegin(GL_LINES);
		glVertex3fv(v1);
		glVertex3fv(v2);
		glEnd();
	}

	glEnable(GL_LIGHTING);
}

void MeshRenderer::drawFaceNormals(ZMeshSpace::Mesh3D *mesh)
{
	float length = mesh->m_averageEdgeLength*0.3f;

	glDisable(GL_LIGHTING);
	int fSize = mesh->get_num_of_faces_list();
	glBegin(GL_LINES);
	for (int i=0; i<fSize; i++)
	{
		HE_face* f = mesh->get_face(i);
		Vec3f c = f->m_vCenter;
		Vec3f p = c + f->m_vnormal*length;
		Vec3f p2 = c + f->m_vnormal2*length;
		glColor3fv(CColor::green);
		glVertex3fv(c);
		glVertex3fv(p);
		glColor3fv(CColor::red);
		glVertex3fv(c);
		glVertex3fv(p2);
	}	
	glEnd();
	glEnable(GL_LIGHTING);

}
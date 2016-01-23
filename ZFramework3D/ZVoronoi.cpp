#include "ZVoronoi.h"
#include <map>
#include <fstream>
#include "ZLine.h"

namespace ZCVT
{
	ZVoronoi::ZVoronoi(){}

	void ZVoronoi::clear()
	{
		voronoiCells_.clear();
		voronoiEdges_.clear();
		voronoiVertice_.clear();
	}

	bool ZVoronoi::isValid()
	{
		return !voronoiCells_.empty() && !voronoiVertice_.empty();
	}

	void ZVoronoi::setClippingRect(float xMin, float xMax, float yMin, float yMax)
	{
		clipRect_[0] = xMin;
		clipRect_[1] = xMax;
		clipRect_[2] = yMin;
		clipRect_[3] = yMax;
	}

	void ZVoronoi::getVoronoiCells(CGAL_Algorithms::Delaunay_triangulation_2& dt)
	{
		clear();

		// delaunay faces <--> voronoi vertices
		int dtFaceNum = dt.number_of_faces();
		voronoiVertice_.resize(dtFaceNum);
		std::map<CGAL_Algorithms::Triangulation_2::Face_handle, int> faceMap;
		CGAL_Algorithms::Delaunay_triangulation_2::Face_iterator fit = dt.faces_begin();

		for (int k=0; fit!=dt.faces_end(); fit++, k++)
		{
			CGAL_Algorithms::Point_2 p = dt.dual(fit);
			voronoiVertice_[k] = Vec2f(p.hx(), p.hy());
			faceMap[fit] = k;
		}
		// test infinite
		//voronoiVertice_.push_back(Vec2f(-1,-1));
		//std::map<CGAL_Algorithms::Triangulation_2::Vertex_handle, int> vrtMap;
		//CGAL_Algorithms::Triangulation_2::Vertex_iterator vit = dt.vertices_begin();
		//for (int k=0; vit!=dt.vertices_end(); vit++, k++) vrtMap[vit] = k;

		// delaynay vertices <--> voronoi cells
		// delaunay edges <--> voronoi edges
		//CGAL_Algorithms::Iso_rectangle_2 rect(clipRect_[0], clipRect_[2], clipRect_[1], clipRect_[3]);
		int dtVertexNum = dt.number_of_vertices();
		CGAL_Algorithms::Triangulation_2::Vertex_iterator vit = dt.vertices_begin();
		voronoiCells_.resize(dtVertexNum);
		//vit=dt.vertices_begin();
		CGAL_Algorithms::Triangulation_2::Face_circulator fcIter, temp, tmp2;
		//for (int k=0; vit!=dt.vertices_end(); vit++, k++)
		for (; vit!=dt.vertices_end(); vit++)
		{
			unsigned int k = vit->info();
			ZVoronoiCell& cell = voronoiCells_[k];
			cell.id_ = k;
			fcIter = vit->incident_faces();
			temp = fcIter;
			CGAL_For_all(fcIter, temp)
			{
				if (dt.is_infinite(fcIter))
				{
					tmp2 = fcIter;
					tmp2++;
					if (dt.is_infinite(tmp2)) {
						tmp2 = fcIter; tmp2--;
						if (dt.is_infinite(tmp2)) {
							std::cerr << "Both sides are infinite!\n";
						}
					}
					CGAL_Algorithms::Edge_2 edge = CGAL_Algorithms::getEdge(dt, fcIter, tmp2);
					CGAL::Object o = dt.dual(edge);
					if (CGAL::object_cast<CGAL_Algorithms::Ray_2>(&o))
					{
						const CGAL_Algorithms::Ray_2* ray2 = CGAL::object_cast<CGAL_Algorithms::Ray_2>(&o);
						voronoiVertice_.push_back(Vec2f(ray2->point(1).hx(), ray2->point(1).hy()));
						cell.vertices_.push_back(voronoiVertice_.size()-1);
					}
					else
					{
						std::cout << " Error!\n";
					}
				}
				else
				{
					cell.vertices_.push_back(faceMap[fcIter]);
				}
			}
			cell.edgeIdx_.resize(cell.vertices_.size(), -1);
		}

		buildVoronoiEdges(dt);
		clipped(0, 1, 0, 1);
		computeCenters();
	}

	void ZVoronoi::computeCenters()
	{
		for (int i=0; i<voronoiCells_.size(); i++)
		{
			ZVoronoiCell& cell = voronoiCells_[i];
			Vec2f p(0,0);
			int count = cell.vertices_.size();
			for (int j=0; j<count; j++)
			{
				p += getVertex(cell.vertices_[j]);
			}
			if (count==0) {
				std::cerr << " count=0 in " << i << " cell\n";
				cell.center_ = p;
				continue;
			}
			p = p*(1.f/count);
			cell.center_ = p;
		}
	}

	void ZVoronoi::buildVoronoiEdges(CGAL_Algorithms::Delaunay_triangulation_2& dt)
	{
		if (voronoiCells_.empty() || voronoiVertice_.empty()) return;

		for (int i=0; i<voronoiCells_.size(); i++)
		{
			ZVoronoiCell& cell = voronoiCells_[i];
			for (int j=0; j<cell.vertices_.size(); j++)
			{
				int v0 = cell.vertices_[j];
				int v1 = cell.vertices_[(j+1)%cell.vertices_.size()];
				ZVoronoiEdge edge;
				edge.v0 = v0;
				edge.v1 = v1;
				edge.f0 = edge.f1 = -1;
				std::vector<ZVoronoiEdge>::iterator findIter = std::find_if(voronoiEdges_.begin(), voronoiEdges_.end(), findSameEdge(edge));
				if (findIter== voronoiEdges_.end())
				{
					edge.id_ = voronoiEdges_.size();
					edge.f0 = i;
					voronoiEdges_.push_back(edge);
				}
				else
				{
					findIter->f1 = i;
				}
			}
		}

		for (int i=0; i<voronoiEdges_.size(); i++)
		{
			ZVoronoiEdge& edge = voronoiEdges_[i];
			if (edge.f0<voronoiCells_.size())
			{
				ZVoronoiCell& cell0 = voronoiCells_[edge.f0];
				for (int j=0; j<cell0.vertices_.size(); j++)
				{
					int v0 = cell0.vertices_[j];
					int v1 = cell0.vertices_[(j+1)%cell0.vertices_.size()];
					if ( (v0==edge.v0 && v1==edge.v1) || (v0==edge.v1 && v1==edge.v0) )
					{
						cell0.edgeIdx_[j] = edge.id_;
						break;
					}
				}
			}
			if (edge.f1<voronoiCells_.size() && edge.f1>=0)
			{
				ZVoronoiCell& cell1 = voronoiCells_[edge.f1];
				for (int j=0; j<cell1.vertices_.size(); j++)
				{
					int v0 = cell1.vertices_[j];
					int v1 = cell1.vertices_[(j+1)%cell1.vertices_.size()];
					if ( (v0==edge.v0 && v1==edge.v1) || (v0==edge.v1 && v1==edge.v0) )
					{
						cell1.edgeIdx_[j] =edge.id_;
						break;
					}
				}
			}
		}
	}

	bool ZVoronoi::isPointInRect(const Vec2f& p)
	{
		return p.x<=clipRect_[1] && p.x>=clipRect_[0] && p.y <=clipRect_[3] && p.y>=clipRect_[2];
	}

	void ZVoronoi::clipped(float xMin, float xMax, float yMin, float yMax)
	{
		setClippingRect(xMin, xMax, yMin, yMax);
		//std::vector<int> tagVoronoiCells(voronoiCells_.size(), 0);
		std::vector<int> bVerticesInside(voronoiVertice_.size(), 0);
		for (int i=0; i<voronoiVertice_.size(); i++)
			bVerticesInside[i] = isPointInRect(voronoiVertice_[i]) ? 1 : 0;

		struct tmp_edge {
			int id;
			int v0, v1;
			int u0, u1;
			bool bBothOutside;
		};

		std::vector<tmp_edge> edgeStructs(voronoiEdges_.size());
		int nBeforeClippingVerticeSize = voronoiVertice_.size();
		for (int i=0; i<voronoiEdges_.size(); i++)
		{
			ZVoronoiEdge& edge = voronoiEdges_[i];
			tmp_edge& tmpEdge = edgeStructs[i];
			// initial tmp_edge
			tmpEdge.id = edge.id_;
			tmpEdge.v0 = edge.v0;
			tmpEdge.v1 = edge.v1;
			tmpEdge.u0 = tmpEdge.u1 = -1;
			tmpEdge.bBothOutside = false;

			Vec2f p0 = getVertex(edge.v0);
			Vec2f p1 = getVertex(edge.v1);
			bool b0 = bVerticesInside[edge.v0]==1;
			bool b1 = bVerticesInside[edge.v1]==1;
			if ( b0 && b1 ) continue;
			if ( b0^b1 )
			{
				Vec2f intersection;
				if (intersect(intersection, p0, p1))
				{
					voronoiVertice_.push_back(intersection);
					if (b0) {
						//edge.v1 = voronoiVertice_.size()-1;
						tmpEdge.u1 = voronoiVertice_.size()-1;
					}
					else {
						//edge.v0 = voronoiVertice_.size()-1;
						tmpEdge.u0 = voronoiVertice_.size()-1;
					}
				}
			}
			if (!b0 && !b1) tmpEdge.bBothOutside = true;
		}

		// add corner points
		voronoiVertice_.push_back(Vec2f(clipRect_[0], clipRect_[2]));
		voronoiVertice_.push_back(Vec2f(clipRect_[1], clipRect_[2]));
		voronoiVertice_.push_back(Vec2f(clipRect_[1], clipRect_[3]));
		voronoiVertice_.push_back(Vec2f(clipRect_[0], clipRect_[3]));

		// modify the voronoi cells
		for (int i=0; i<voronoiCells_.size(); i++)
		{
			ZVoronoiCell& cell = voronoiCells_[i];
			std::vector<int> newVerticesLoop;
			for (int j=0; j<cell.vertices_.size(); j++)
			{
				int eIdx = cell.edgeIdx_[j];
				tmp_edge& edge = edgeStructs[eIdx];
				if (edge.bBothOutside) continue;
				bool bReverse = (edge.v1==cell.vertices_[j]);
				if (!bReverse && edge.u0==-1 && edge.u1!=-1) {
					newVerticesLoop.push_back(edge.v0);
					newVerticesLoop.push_back(edge.u1);
				} 
				else if (bReverse && edge.u0!=-1 && edge.u1==-1) {
					newVerticesLoop.push_back(edge.v1);
					newVerticesLoop.push_back(edge.u0);
				}
				else if (!bReverse && edge.u0!=-1 && edge.u1==-1) {
					newVerticesLoop.push_back(edge.u0);
				}
				else if (bReverse && edge.u0==-1 && edge.u1!=-1) {
					newVerticesLoop.push_back(edge.u1);
				}
				else {
					newVerticesLoop.push_back(cell.vertices_[j]);
				}
				//if (edge.u0!=-1) newVerticesLoop.push_back(edge.u0);
				//else if (edge.u1!=-1) newVerticesLoop.push_back(edge.u1);
				//else newVerticesLoop.push_back(edge.v0);
			}
			// add corner vertices
			std::vector<int> finalVoronoiCellV = newVerticesLoop;
			for (int j=0; j<newVerticesLoop.size(); j++)
			{
				int v0 = newVerticesLoop[j];
				int v1 = newVerticesLoop[(j+1)%newVerticesLoop.size()];
				if (v0>=nBeforeClippingVerticeSize && v1>=nBeforeClippingVerticeSize)
				{
					Vec2f p0 = getVertex(v0);
					Vec2f p1 = getVertex(v1);
					int i0 = onLine(p0);
					int i1 = onLine(p1);
					int addPIdx = -1;
					if ((i0==0 && i1==2) || (i0==2 && i1==0)) {
						addPIdx = voronoiVertice_.size()-4;
					}
					if ((i0==1 && i1==2) || (i0==2 && i1==1)) {
						addPIdx = voronoiVertice_.size()-3;
					}
					if ((i0==1 && i1==3) || (i0==3 && i1==1)) {
						addPIdx = voronoiVertice_.size()-2;
					}
					if ((i0==0 && i1==3) || (i0==3 && i1==0)) {
						addPIdx = voronoiVertice_.size()-1;
					}
					if (addPIdx!=-1) {
						finalVoronoiCellV.insert(finalVoronoiCellV.begin()+j+1, addPIdx);
						break;
					}
				}
			}
			cell.vertices_ = finalVoronoiCellV;
		}
	}

	/**
	       3
		 ---------
		 |       |
		0|       |1
		 ---------
		   2
	*/
	int ZVoronoi::onLine(const Vec2f& p)
	{
		if (p.x==clipRect_[0]) return 0;
		if (p.x==clipRect_[1]) return 1;
		if (p.y==clipRect_[2]) return 2;
		if (p.y==clipRect_[3]) return 3;
		return -1;
	}

	bool ZVoronoi::intersect(Vec2f& intersection, const Vec2f& p, const Vec2f& q)
	{
		float xMin = clipRect_[0];
		float xMax = clipRect_[1];
		float yMin = clipRect_[2];
		float yMax = clipRect_[3];
		static ZLine lines[4] = {ZLine(xMin, yMin, xMax, yMin),
								ZLine(xMax, yMin, xMax, yMax),
								ZLine(xMax, yMax, xMin, yMax),
								ZLine(xMin, yMax, xMin, yMin)};

		ZLine testLine(p, q);
		for (int i=0; i<4; i++)
		{
			if (lines[i].intersect(testLine, intersection))
				return true;
		}
		return false;
	}

	int ZVoronoi::cellNeighborIdx(int cellIdx, int neiIdx)
	{
		ZVoronoiCell& vorCell = cell(cellIdx);
		ZVoronoiEdge& vorEdge = voronoiEdges_[vorCell.edgeIdx_[neiIdx]];
		if (vorEdge.f0==cellIdx) return vorEdge.f1;
		else if (vorEdge.f1==cellIdx) return vorEdge.f0;
		else if (vorEdge.f0==-1 || vorEdge.f1==-1) return -1;
		else {
			std::cerr << "Error with voronoi edges! " << cellIdx << "->" << vorCell.edgeIdx_[neiIdx] << "\n";
			return -2;
		}
	}

	void ZVoronoi::saveVoronoiToFile(const char* fileName)
	{
		std::ofstream fout(fileName);

		for (int i=0; i<voronoiVertice_.size(); i++)
		{
			fout << "v " << voronoiVertice_[i] << " 0\n";
		}
		for (int i=0; i<voronoiCells_.size(); i++)
		{
			ZVoronoiCell& cell = voronoiCells_[i];
			fout << "f ";
			for (int j=0; j<cell.vertices_.size(); j++)
			{
				fout << cell.vertices_[j]+1 << " ";
			}
			fout << "\n";
		}

		fout.close();
	}
}
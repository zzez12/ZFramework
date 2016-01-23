#pragma once
#ifndef ZVORONOI_H_
#define ZVORONOI_H_

#include "GlobalDefs.h"
#include <vector>
#include "CGAL_Algorithms.h"

namespace ZCVT
{
	class ZVoronoiCell
	{
	public:
		int id_;
		std::vector<int> vertices_;
		std::vector<int> edgeIdx_;
		std::vector<int> neighborCellIdx_;
		Vec2f center_;
	};

	class ZVoronoiEdge
	{
	public:
		int id_;
		int v0, v1;	// the two end point indices
		int f0, f1;	// the two voronoi cell indices	
	};

	class findSameEdge
	{
	public:
		findSameEdge(const ZVoronoiEdge& edge) {
			v0 = edge.v0;
			v1 = edge.v1;
		}

		bool operator()(const ZVoronoiEdge& edge) {
			return (v0==edge.v0&&v1==edge.v1) || (v0==edge.v1&&v1==edge.v0);
		}

	private:
		int v0, v1;
	};


	class ZVoronoi
	{
	//public:
	private:
		std::vector<ZVoronoiCell> voronoiCells_;
		std::vector<Vec2f> voronoiVertice_;

	private:
		std::vector<ZVoronoiEdge> voronoiEdges_;
		float clipRect_[4];	// [clipXMin_, clipXMax_, clipYMin_, clipYMax_]

	public:
		ZVoronoi();

		bool isValid();
		void clear();
		void setClippingRect(float xMin, float xMax, float yMin, float yMax);
		void getVoronoiCells(CGAL_Algorithms::Delaunay_triangulation_2& dt);
		void clipped(float xMin, float xMax, float yMin, float yMax);
		void saveVoronoiToFile(const char* fileName);

		int voronoiVerticeSize() {return voronoiVertice_.size();}
		Vec2f getVertex(int idx) {return voronoiVertice_[idx];}

		void computeCenters();

		int voronoiCellSize() {return voronoiCells_.size();}
		ZVoronoiCell& cell(int idx) {return voronoiCells_[idx];}
		int cellNeighborSize(int idx) {return cell(idx).edgeIdx_.size();}
		int cellNeighborIdx(int cellIdx, int neiIdx);

	private:
		void buildVoronoiEdges(CGAL_Algorithms::Delaunay_triangulation_2& dt);

	private:
		bool isPointInRect(const Vec2f& p);
		bool intersect(Vec2f& intersection, const Vec2f& p, const Vec2f& q);
		int onLine(const Vec2f& p);
	};

// /	struct tmp_OutsideVertex
// 	{
// 		int id;
// 		int 
// 	};
}

#endif//ZVORONOI_H_
#pragma once
#ifndef ZDELAUNAY_H_
#define ZDELAUNAY_H_

#include "CGAL_Algorithms.h"

namespace ZCVT
{
	class ZDelaunayTriangle
	{
	public:
		int id_;
		int vertices_[3];
	};

	class ZDelaunay
	{
	public:
		ZDelaunay();

		void getDelaunay(CGAL_Algorithms::Delaunay_triangulation_2& dt);

	private:
		std::vector<Vec2f> delaunayVertices_;
		std::vector<ZDelaunayTriangle> delaunayTriangles_;
	};
}

#endif//ZDELAUNAY_H_
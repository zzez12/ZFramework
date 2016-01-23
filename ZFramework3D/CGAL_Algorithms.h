#pragma once
#ifndef CGAL_ALGORITHMS_H_
#define CGAL_ALGORITHMS_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <iterator>

#include "GlobalDefs.h"

namespace CGAL_Algorithms
{
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef K::Point_2 Point_2;
	typedef K::Iso_rectangle_2 Iso_rectangle_2;
	typedef K::Segment_2 Segment_2;
	typedef K::Ray_2 Ray_2;
	typedef K::Line_2 Line_2;
	typedef CGAL::Delaunay_triangulation_2<K>  Delaunay_triangulation_2;
	typedef Delaunay_triangulation_2::Triangulation Triangulation_2;
	typedef Triangulation_2::Edge Edge_2;

	struct Cropped_voronoi_from_delaunay{
		std::list<Segment_2> m_cropped_vd;
		Iso_rectangle_2 m_bbox;

		Cropped_voronoi_from_delaunay(const Iso_rectangle_2& bbox):m_bbox(bbox){}

		template <class RSL>
		void crop_and_extract_segment(const RSL& rsl){
			CGAL::Object obj = CGAL::intersection(rsl,m_bbox);
			const Segment_2* s=CGAL::object_cast<Segment_2>(&obj);
			if (s) m_cropped_vd.push_back(*s);
		}

		void operator<<(const Ray_2& ray)    { crop_and_extract_segment(ray); }
		void operator<<(const Line_2& line)  { crop_and_extract_segment(line); }
		void operator<<(const Segment_2& seg){ crop_and_extract_segment(seg); }
	};

	static void delaunay(const std::vector<Vec2f>& input, Delaunay_triangulation_2& dt2)
	{
		std::vector<Point_2> points;
		for (int i=0; i<input.size(); i++)
		{
			points.push_back(Point_2(input[i].x, input[i].y));
		}

		dt2.clear();
		//insert points into the triangulation
		dt2.insert(points.begin(), points.end());

		//construct a rectangle
//		Iso_rectangle_2 bbox(0, 0, 1, 1);
//		Cropped_voronoi_from_delaunay vor(bbox);

		//extract the cropped Voronoi diagram
//		dt2.draw_dual(vor);
		//print the cropped Voronoi diagram as segments
// 		std::copy(vor.m_cropped_vd.begin(),vor.m_cropped_vd.end(),
// 			std::ostream_iterator<Segment_2>(std::cout,"\n"));
// 		std::ofstream fout("voronoiOut.txt");
// 		std::ofstream fout1("input.txt");
// 		std::copy(input.begin(), input.end(), std::ostream_iterator<Vec2f>(fout1, "\n"));
// 		std::copy(vor.m_cropped_vd.begin(),vor.m_cropped_vd.end(),
// 			std::ostream_iterator<Segment_2>(fout,"\n"));
	}

	static Edge_2 getEdge(Triangulation_2& t, Triangulation_2::Face_handle f0, Triangulation_2::Face_handle f1)
	{
		int which = 0;
		for (;which<3; which++)
		{
			if (f0->neighbor(which)==f1)
				break;
		}
		if (which==3)
			return Edge_2();
		return Edge_2(f0, which);
	}
}

#endif//CGAL_ALGORITHMS_H_
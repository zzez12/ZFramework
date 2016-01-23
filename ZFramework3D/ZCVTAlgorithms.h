#pragma once
#ifndef ZCVT_ALGORITHMS_H_
#define ZCVT_ALGORITHMS_H_

#include "ZCVTData.h"
#include "./LpCVT/combinatorics/mesh.h"
#include "CGAL_Algorithms.h"
#include "ZVoronoi.h"

namespace ZCVT
{
	class ZCVTAlgorithms
	{
	private:
		ZCVTData *pData_;
	public:
		ZCVTAlgorithms();
		~ZCVTAlgorithms();

		void init(ZCVTData* data);
		void destroy();

		ZCVTData* getData() {return pData_;}
		CGAL_Algorithms::Delaunay_triangulation_2& getDelaunay() {return dt_;}
		//ZVoronoi* getVoronoi() {return &voronoi_;}

		bool lineCVT_init(int nSamples);
		bool lineCVT_run();
		bool lineCVT_iterate();

	private:
		void get_combinatorics(Geex::Mesh* M, const std::vector<Geex::vec3>& pts, 
			std::vector<int>& I, std::vector<Geex::vec3>& C, std::vector<int>& F, bool volume);
		void compute_F_g(Geex::Mesh* M, const std::vector<Geex::vec3>& pts, unsigned int p, bool volume);

		Vec2f computeMassCenter(int idx);

		//float getVoronoiCellL2Energy(ZVoronoi* vor, int cellIdx, const Vec2f& c0);
		//float getCVTLineEnergy();
		//float getRegularizationEnergy();

		//float getSimplexEnergy(const Vec2f& p0, const Vec2f& p1, const Vec2f& p2);

		void buildCVTLineCellNeighbors();

		void optimizeByHLBFGS(int sampleN, double* initX, int numIter, int M, int T, ZCVTData* pData);

	private:
		CGAL_Algorithms::Delaunay_triangulation_2 dt_;
		//ZVoronoi voronoi_;
	};
}

#endif//ZCVT_ALGORITHMS_H_
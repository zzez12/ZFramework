#include "ZCVTAlgorithms.h"
#include "CGAL_Algorithms.h"

#include <LpCVT/combinatorics/delaunay.h>
#include <LpCVT/combinatorics/RVD.h>
#include <LpCVT/combinatorics/clipped_VD.h>
#include <LpCVT/algebra/F_Lp.h>
#include <LpCVT/common/line_stream.h>
#include "ZCVTEnergy.h"

#include "hlbfgs/HLBFGS.h"
#include "mathUtils.h"

namespace ZCVT
{

	/**
     * Used by get_combinatorics() in volume mode
     */
    class MemorizeIndices {
    public:
        MemorizeIndices(
            std::vector<int>& I_in,
            std::vector<Geex::vec3>& C_in
        ) : I(I_in), C(C_in) {
            I.resize(0) ;
            C.resize(0) ;
        }

        void operator() (
            unsigned int i, 
            int j,
            const Geex::VertexEdge& v1, 
            const Geex::VertexEdge& v2, 
            const Geex::VertexEdge& v3
        ) const {
            I.push_back(i) ;
            I.push_back(v1.sym[2]) ;
            I.push_back(v1.sym[1]) ;
            I.push_back(v1.sym[0]) ;
            I.push_back(v2.sym[2]) ;
            I.push_back(v2.sym[1]) ;
            I.push_back(v2.sym[0]) ;
            I.push_back(v3.sym[2]) ;
            I.push_back(v3.sym[1]) ;
            I.push_back(v3.sym[0]) ;
            C.push_back(v1) ;
            C.push_back(v2) ;
            C.push_back(v3) ;
        }
    private:
        mutable std::vector<int>& I ;
        mutable std::vector<Geex::vec3>& C ;
    } ;

    /**
     * Used by get_combinatorics() in surface mode
     */
    class MemorizeIndicesAndFacets{
    public:
        MemorizeIndicesAndFacets(
            const Geex::RestrictedVoronoiDiagram& RVD_in,
            std::vector<int>& I_in,
            std::vector<Geex::vec3>& C_in,
            std::vector<int>& F_in
        ) : RVD(RVD_in), I(I_in), C(C_in), F(F_in) {
            I.resize(0) ;
            C.resize(0) ;
            F.resize(0) ;
        }

        void operator() (
            unsigned int i, 
            const Geex::VertexEdge& v1, 
            const Geex::VertexEdge& v2, 
            const Geex::VertexEdge& v3
        ) const {
            I.push_back(i) ;
            I.push_back(v1.sym[2]) ;
            I.push_back(v1.sym[1]) ;
            I.push_back(v1.sym[0]) ;
            I.push_back(v2.sym[2]) ;
            I.push_back(v2.sym[1]) ;
            I.push_back(v2.sym[0]) ;
            I.push_back(v3.sym[2]) ;
            I.push_back(v3.sym[1]) ;
            I.push_back(v3.sym[0]) ;
            F.push_back(RVD.current_facet()) ;
            C.push_back(v1) ;
            C.push_back(v2) ;
            C.push_back(v3) ;
        }
    private:
        const Geex::RestrictedVoronoiDiagram& RVD ;
        mutable std::vector<int>& I ;
        mutable std::vector<Geex::vec3>& C ;
        mutable std::vector<int>& F ;
    } ;

	ZCVTAlgorithms::ZCVTAlgorithms()
	{
		init(NULL);
	}

	ZCVTAlgorithms::~ZCVTAlgorithms()
	{
		destroy();
	}

	void ZCVTAlgorithms::init(ZCVTData* data)
	{
		pData_ = data;
	}

	void ZCVTAlgorithms::destroy()
	{
	}

	bool ZCVTAlgorithms::lineCVT_init(int nSamples)
	{
		if (pData_==NULL) return false;
		return pData_->init(nSamples);
	}

	bool ZCVTAlgorithms::lineCVT_run()
	{
		if (pData_==NULL) return false;
		// delaunay
		std::vector<ZLine> lines = pData_->getLines();
		std::vector<Vec2f> points;
		int subSize = pData_->getSubSize();
		for (int i=0; i<lines.size(); i++)
		{
			//points.push_back(lines[i].p0);
			//points.push_back(lines[i].p1);
			std::vector<int> sampleIdx;
			for (int j=0; j<=subSize; j++)
			{
				points.push_back(pData_->getSamplePoint(i, j));
				sampleIdx.push_back(points.size()-1);
			}
			pData_->setCorresponding(i, sampleIdx);
		}
		CGAL_Algorithms::delaunay(points, dt_);
		//voronoi_.getVoronoiCells(dt_);
		pData_->getVoronoi()->getVoronoiCells(dt_);
		buildCVTLineCellNeighbors();

		return false;
	}

	bool ZCVTAlgorithms::lineCVT_iterate()
	{
		// compute the centers of each voronoi cell
		//std::vector<Vec2f> points;
		//for (int i=0; i<pData_->getVoronoi()->voronoiCellSize(); i++)
		//{
		//	ZVoronoiCell& cell = pData_->getVoronoi()->cell(i);
		//	points.push_back(cell.center_);
		//}

		//CGAL_Algorithms::delaunay(points, dt_);
		//pData_->getVoronoi()->getVoronoiCells(dt_);
		int sampleN = pData_->numLines()*2*2;
		std::vector<double> x(sampleN);
		int numIter = 1000;
		int M = 7;
		int T = 0;

		// init x
		for (int i=0; i<pData_->numLines(); i++)
		{
			pData_->getLineEndPoints(i, &x[i*4]);
		}

		optimizeByHLBFGS(sampleN, &x[0], numIter, M, T, pData_);

		// copy the result to the voronoi cell
		for (int i=0; i<pData_->numLines(); i++)
		{
			pData_->setLineEndPoints(i, Vec2f(x[i*4], x[i*4+1]), Vec2f(x[i*4+2], x[i*4+3]));
		}

		return true;
	}

	void ZCVTAlgorithms::optimizeByHLBFGS(int sampleN, double* initX, int numIter, int M, int T, ZCVTData* pData)
	{
		double parameter[20];
		int info[20];
		INIT_HLBFGS(parameter, info);

		//info[3] = 1;
		info[4] = numIter;
		info[6] = T;
		info[7] =  0;	// without hessian
		info[10] = 0;
		info[11] = 1;

		HLBFGS(sampleN, M, initX, ZCVTEnergy::evalFunc, 0, HLBFGS_UPDATE_Hessian, ZCVTEnergy::newIteration, parameter, info, (void*)pData);

	}

	void ZCVTAlgorithms::buildCVTLineCellNeighbors()
	{
		if (pData_==NULL || pData_->getVoronoi()==NULL || !pData_->getVoronoi()->isValid())
			return;

		int nSampleSize = pData_->numSampleSize();
		for (int i=0; i<nSampleSize; i++)
		{
			ZVoronoiCell& vorCell = pData_->getVoronoi()->cell(i);
			std::vector<int> neighbors;
			for (int j=0; j<vorCell.edgeIdx_.size(); j++)
			{
				int neighb = pData_->getVoronoi()->cellNeighborIdx(i, j);
				if (neighb>=0)
					neighbors.push_back(neighb);
			}
			vorCell.neighborCellIdx_ = neighbors;
		}
		//std::vector<int> vrt2lineCorres(nSampleSize, 0);
		//for (int i=0; i<nSampleSize; i++)
		//{
		//	vrt2lineCorres[i] = pData_->sampleOnLine(i);
		//}

		//int nLineSize = pData_->numLines();
		//int nSubSize = pData_->getSubSize();
		//for (int i=0; i<nLineSize; i++)
		//{
		//	ZCVTLineCell* lineCell = pData_->getLineCell(i);
		//	std::vector<int> inLineIdx;
		//	std::vector<int> neiLineIdx;
		//	for (int j=0; j<=nSubSize; j++)
		//	{
		//		int vorCellIdx = pData_->getVoronoiCellIdx(i, j);
		//		inLineIdx.push_back(vorCellIdx);
		//		int neiSize = pData_->getVoronoi()->cellNeighborSize(vorCellIdx);
		//		for (int k=0; k<neiSize; k++)
		//		{
		//			neiLineIdx.push_back(pData_->getVoronoi()->cellNeighborIdx(vorCellIdx, k));
		//		}
		//	}
		//	// remove duplicates
		//	std::sort(neiLineIdx.begin(), neiLineIdx.end());
		//	std::unique(neiLineIdx.begin(), neiLineIdx.end());

		//	// remove the inLineIdx's
		//	std::vector<int> lineNeighbors;
		//	MATH_Utils::sub(lineNeighbors, neiLineIdx, inLineIdx);
		//	lineCell->neighborLineIds_ = lineNeighbors;
		//}
	}

	Vec2f ZCVTAlgorithms::computeMassCenter(int idx)
	{
		ZVoronoiCell& cell = pData_->getVoronoi()->cell(idx);
		float area = 0;
		for (int i=0; i<cell.vertices_.size(); i++)
		{
			Vec2f p0 = pData_->getVoronoi()->getVertex(cell.vertices_[i]);
			Vec2f p1 = pData_->getVoronoi()->getVertex(cell.vertices_[(i+1)%cell.vertices_.size()]);
			area += p0.x*p1.y - p0.y*p1.x;
		}
		if (area<0) area=-area;

		// unfinished!!
		return Vec2f();
	}

	/*float ZCVTAlgorithms::getVoronoiCellL2Energy(ZVoronoi* vor, int cellIdx, const Vec2f& c0)
	{
		ZVoronoiCell& cell = vor->cell(cellIdx);
		if (cell.vertices_.size()<=2) return 0;
		float ret=0;

		float area = 0;
		for (int i=0; i<cell.vertices_.size(); i++)
		{
			Vec2f p0 = vor->getVertex(cell.vertices_[i]);
			Vec2f p1 = vor->getVertex(cell.vertices_[(i+1)%cell.vertices_.size()]);
			area += p0.x*p1.y - p0.y*p1.x;
		}
		if (area<0) area=-area;

		Vec2f v0 = vor->getVertex(cell.vertices_[0])-c0;
		for (int i=1; i<=cell.vertices_.size()-2; i++)
		{
			Vec2f v1 = vor->getVertex(cell.vertices_[i])-c0;
			Vec2f v2 = vor->getVertex(cell.vertices_[i+1])-c0;
			ret += getSimplexEnergy(v0, v1, v2);
		}

		return ret*area/6;
	}

	float ZCVTAlgorithms::getSimplexEnergy(const Vec2f& p0, const Vec2f& p1, const Vec2f& p2)
	{
		return p0.dot(p0) + p1.dot(p1) + p2.dot(p2) + p0.dot(p1) + p1.dot(p2) + p2.dot(p0);
	}

	float ZCVTAlgorithms::getCVTLineEnergy()
	{
		float ret = 0;
		for (int i=0; i<pData_->numLines(); i++)
		{
			for (int j=0; j<=pData_->getSubSize(); j++)
			{
				int cellIdx = pData_->getVoronoiCellIdx(i, j);
				Vec2f xk = pData_->getSamplePoint(i, j);
				ret += getVoronoiCellL2Energy(pData_->getVoronoi(), cellIdx, xk);
			}
		}

		return ret;
	}*/

	/**
     * Gets the combinatorics of the integration simplices,
     * i.e. 10 integers per integration simplex.
     * (see Section 3.1 in the paper)
     * Returns also the array of C vertices (three per integration simplex).
     * Since they are easy to get during the combinatorial phase, they are
     * computed here and kept for the algebraic phase.
     *
     * In 2D mode (volume = false), returns also the array F.
     *   F[i] indicates the facet that contains the i-th integration simplex.
     *
     */
    void ZCVTAlgorithms::get_combinatorics(
		Geex::Mesh* M, const std::vector<Geex::vec3>& pts, 
		std::vector<int>& I, std::vector<Geex::vec3>& C, std::vector<int>& F, bool volume
    ) {
        Geex::Delaunay* delaunay = Geex::Delaunay::create("CGAL") ;
        delaunay->set_vertices(pts) ;
        if(volume) {
            Geex::ClippedVoronoiDiagram CVD(delaunay,M) ;
            CVD.for_each_triangle(MemorizeIndices(I,C)) ;
        } else {
            Geex::RestrictedVoronoiDiagram RVD(delaunay,M) ;
            RVD.set_symbolic(true) ;
            RVD.for_each_triangle(MemorizeIndicesAndFacets(RVD,I,C,F)) ;
        }
        delete delaunay ;
    }

	/**
     * Computes F_{L_p} and its gradient.
     */
    void ZCVTAlgorithms::compute_F_g(Geex::Mesh* m, const std::vector<Geex::vec3>& pts, unsigned int p, bool volume) {
        std::cerr << "nb pts = " << pts.size() << "   nb facets = " << m->nb_facets() << std::endl ;
        std::vector<int> I ;
        std::vector<Geex::vec3> C ;
        std::vector<int> F ;
        get_combinatorics(m, pts, I, C, F, volume) ;
        unsigned int nb_integration_simplices = (unsigned int)I.size() / 10 ;
        std::vector<Geex::mat3> M(nb_integration_simplices) ;
        for(unsigned int i=0; i<M.size(); i++) {
            M[i].load_identity() ; 
                // or replace with anisotropy field
                //   In 2D: use F[i] to retreive the index of the facet that contains
                //      the current integration simplex (and access an array of per-facet anisotropy).
                //   In 3D: use geometric search from the centroid of the current
                //      integration simplex.
        }
        std::vector<Geex::plane3> Q(m->nb_facets()) ;
        for(unsigned int i=0; i<m->nb_facets(); i++) {
            Q[i] = m->facet_plane(i) ;
        }
        std::vector<double> g(pts.size() * 3) ;
        double f = compute_F_Lp(volume, p, m, I, C, pts, Q, M, g) ;
        double gnorm = 0.0 ;
        for(unsigned int i=0; i<g.size(); i++) {
            gnorm += g[i]*g[i] ;
        }
        gnorm = ::sqrt(gnorm) ;
        std::cerr.precision(16);
        std::cerr << (volume ? "volume " : "surface ") 
                  << "F_L" << p << ":" 
                  << "f=" << std::scientific << f << "  g=" << gnorm << std::endl ;
    }
}

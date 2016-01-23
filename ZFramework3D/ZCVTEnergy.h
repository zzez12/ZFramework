#pragma once
#ifndef ZCVTENERGY_H_
#define ZCVTENERGY_H_

#include "ZCVTData.h"
#include "ZVoronoi.h"

namespace ZCVT
{
	class ZCVTEnergy
	{
	public:
		static float getSimplexEnergy(const Vec2f& p0, const Vec2f& p1, const Vec2f& p2)
		{
			return p0.dot(p0) + p1.dot(p1) + p2.dot(p2) + p0.dot(p1) + p1.dot(p2) + p2.dot(p0);
		}

		static Vec2f getCellSampleCurrentPos(ZCVTData* pData, int cellIdx, double* x)
		{
			int lineIdx = pData->sampleOnLine(cellIdx);
			int lineInternalIdx= pData->sampleOnLineInterval(cellIdx);
			Vec2f p0 = Vec2f(x[lineIdx*4+0], x[lineIdx*4+1]);
			Vec2f p1 = Vec2f(x[lineIdx*4+2], x[lineIdx*4+3]);
			return p0 + (p1-p0)*(1.f*lineInternalIdx/pData->getSubSize());
		}

		static float getVoronoiCellL2Energy(ZCVTData* pData, int cellIdx, const Vec2f& c0, double *x)
		{
			ZVoronoiCell& cell = pData->getVoronoi()->cell(cellIdx);
			if (cell.vertices_.size()<=2) return 0;
			float ret=0;

			float area = 0;
			for (int i=0; i<cell.vertices_.size(); i++)
			{
				int idx0 = cell.vertices_[i];
				int idx1 = cell.vertices_[(i+1)%cell.vertices_.size()];
				Vec2f p0 = pData->getVoronoi()->getVertex(cell.vertices_[i]);//getCellSampleCurrentPos(pData, idx0, x);//
				Vec2f p1 = pData->getVoronoi()->getVertex(cell.vertices_[(i+1)%cell.vertices_.size()]);//getCellSampleCurrentPos(pData, idx1, x);//
				area += p0.x*p1.y - p0.y*p1.x;
			}
			if (area<0) area=-area;

			Vec2f v0 = pData->getVoronoi()->getVertex(cell.vertices_[0])-c0;//getCellSampleCurrentPos(pData, cell.vertices_[0], x)-c0;//
			for (int i=1; i<=cell.vertices_.size()-2; i++)
			{
				Vec2f v1 = pData->getVoronoi()->getVertex(cell.vertices_[i])-c0;//getCellSampleCurrentPos(pData, cell.vertices_[i], x)-c0;//
				Vec2f v2 = pData->getVoronoi()->getVertex(cell.vertices_[i+1])-c0;//getCellSampleCurrentPos(pData, cell.vertices_[i+1], x)-c0;//
				ret += getSimplexEnergy(v0, v1, v2);
			}

			return ret*area/6;
		}

		static void evalFunc(int sampleN, double *x, double *prev_x, double *f, double *g, void *pData)
		{
			static double para_gamma = 1.0;
			ZCVTData* cvtData = (ZCVTData*)pData;
			*f = 0;
			int subSize = cvtData->getSubSize();
			int nLines = cvtData->numLines();
			int nSampleSize = cvtData->numSampleSize();

			// cell centers, precomputed for using them many times
			std::vector<Vec2f> cellCenters(nSampleSize);
			for (int i=0; i<nLines; i++)
			{
				int id0 = i*2;
				int id1 = i*2+1;
				//Vec2f x0 = Vec2f(x[id0*2], x[id0*2+1]);
				//Vec2f x1 = Vec2f(x[id1*2], x[id1*2+1]);
				//Vec2f step = (x1-x0)*(1.f/subSize);
				for (int j=0; j<=subSize; j++)
				{
					//Vec2f xk = x0+step*j;
					Vec2f sumNeigh(0,0);
					int vorId = cvtData->getVoronoiCellIdx(i, j);
					ZVoronoiCell& cell = cvtData->getVoronoi()->cell(vorId);
					for (int k=0; k<cell.neighborCellIdx_.size(); k++)
					{
						int neiIdx = cell.neighborCellIdx_[k];
						sumNeigh += getCellSampleCurrentPos(cvtData, neiIdx, x);//Vec2f(x[neiIdx*2], x[neiIdx*2+1]);
					}
					if (cell.neighborCellIdx_.size()>0)
						sumNeigh = sumNeigh*(1.f/cell.neighborCellIdx_.size());
					cellCenters[vorId] = sumNeigh;
				}
			}

			// mk, ck
			std::vector<double> mk(nSampleSize, 0);
			std::vector<Vec2f> ck(nSampleSize);
			for (int i=0; i<nLines; i++)
			{
				int id0 = i*2;
				int id1 = i*2+1;
				//Vec2f x0 = Vec2f(x[id0*2], x[id0*2+1]);
				//Vec2f x1 = Vec2f(x[id1*2], x[id1*2+1]);
				//Vec2f step = (x1-x0)*(1.f/subSize);
				for (int j=0; j<=subSize; j++)
				{
					//Vec2f xk = x0+step*j;
					Vec2f sumNeigh(0,0);
					int vorId = cvtData->getVoronoiCellIdx(i, j);
					ZVoronoiCell& cell = cvtData->getVoronoi()->cell(vorId);
					double area = 0;
					for (int k=0; k<cell.vertices_.size(); k++)
					{
						int neiIdx0 = cell.vertices_[k];
						int neiIdx1 = cell.vertices_[(k+1)%cell.vertices_.size()];
						//sumNeigh += Vec2f(x[neiIdx*2], x[neiIdx*2+1]);
						Vec2f p0 = cvtData->getVoronoi()->getVertex(neiIdx0);//getCellSampleCurrentPos(cvtData, neiIdx0, x);//Vec2f(x[neiIdx0*2], x[neiIdx0*2+1]);
						Vec2f p1 = cvtData->getVoronoi()->getVertex(neiIdx1);//getCellSampleCurrentPos(cvtData, neiIdx1, x);//Vec2f(x[neiIdx1*2], x[neiIdx1*2+1]);
						area += p0.x*p1.y - p0.y*p1.x;
					}
					mk[vorId] = area>0 ? area : -area;
					ck[vorId] = cell.center_;
				}
			}

			// cvt 
			for (int i=0; i<nLines; i++)
			{
				int id0 = i*2;
				int id1 = i*2+1;
				Vec2f x0 = Vec2f(x[id0*2], x[id0*2+1]);
				Vec2f x1 = Vec2f(x[id1*2], x[id1*2+1]);
				Vec2f step = (x1-x0)*(1.f/subSize);
				for (int j=0; j<=subSize; j++)
				{
					Vec2f xk = x0 + step*j;
					int vorId = cvtData->getVoronoiCellIdx(i, j);
					//*f += getVoronoiCellL2Energy(cvtData->getVoronoi(), vorId, xk);	// cvt energy
					*f += getVoronoiCellL2Energy(cvtData, vorId, xk, x);
				}
			}


			// regularization
			for (int i=0; i<nLines; i++)
			{
				int id0 = i*2;
				int id1 = i*2+1;
				Vec2f x0 = Vec2f(x[id0*2], x[id0*2+1]);
				Vec2f x1 = Vec2f(x[id1*2], x[id1*2+1]);
				Vec2f step = (x1-x0)*(1.f/subSize);
				for (int j=0; j<=subSize; j++)
				{
					Vec2f xk = x0+step*j;
					Vec2f deltaK = xk-ck[cvtData->getVoronoiCellIdx(i, j)];
					*f += deltaK.dot(deltaK)*para_gamma;
				}
			}

			// deviation
			for (int i=0; i<nLines; i+=2)
			{
				int idxI = i*2;
				int idxJ = i*2+1;
				g[idxI*2] = g[idxI*2+1] = g[idxJ*2] = g[idxJ*2+1] = 0;
				Vec2f xi = Vec2f(x[idxI*2], x[idxI*2+1]);
				Vec2f xj = Vec2f(x[idxJ*2], x[idxJ*2+1]);
				Vec2f step = (xj-xi)*(1.f/subSize);
				for(int j=0; j<=subSize; j++)
				{
					int idxK = cvtData->getVoronoiCellIdx(i, j);
					Vec2f pk = xi + step*j;
					// cvt part
					double lamdaK = 1.f*j/subSize;
					Vec2f deltaK = (pk-ck[idxK])*(2*mk[idxK]);
					g[idxI*2]	+= deltaK.x*(lamdaK);
					g[idxI*2+1] += deltaK.y*(lamdaK);
					g[idxJ*2]	+= deltaK.x*(1-lamdaK);
					g[idxJ*2+1] += deltaK.y*(1-lamdaK);
				}

				// regularization part
				Vec2f regPI = xi - ck[idxI];
				ZVoronoiCell& cellI = cvtData->getVoronoi()->cell(idxI);
				Vec2f neiI(0,0);
				for (int k=0; k<cellI.neighborCellIdx_.size(); k++)
				{
					int neiJIdx = cellI.neighborCellIdx_[k];
					int neineiJSize = cvtData->getVoronoi()->cell(neiJIdx).neighborCellIdx_.size();
					if (neineiJSize>1)
						//neiI = neiI + (Vec2f(x[neiJIdx*2], x[neiJIdx*2+1])-cellCenters[neiJIdx])*(1.f/neineiJSize);
						neiI += (getCellSampleCurrentPos(cvtData, neiJIdx, x) - ck[neiJIdx])*(1.f/neineiJSize);
				}
				regPI = (regPI - neiI)*2;
				g[idxI*2]	+= regPI.x*para_gamma;
				g[idxI*2+1] += regPI.y*para_gamma;

				Vec2f regPJ = xj - ck[idxJ];
				ZVoronoiCell& cellJ = cvtData->getVoronoi()->cell(idxJ);
				Vec2f neiJ(0,0);
				for (int k=0; k<cellJ.neighborCellIdx_.size(); k++)
				{
					int neiJIdx = cellJ.neighborCellIdx_[k];
					int neineiJSize = cvtData->getVoronoi()->cell(neiJIdx).neighborCellIdx_.size();
					if (neineiJSize>1)
						//neiJ = neiJ + (Vec2f(x[neiJIdx*2], x[neiJIdx*2+1])-cellCenters[neiJIdx])*(1.f/neineiJSize);
						neiJ += (getCellSampleCurrentPos(cvtData, neiJIdx, x) - ck[neiJIdx])*(1.f/neineiJSize);
				}
				regPJ = (regPJ - neiJ)*2;
				g[idxJ*2]	+= regPJ.x*para_gamma;
				g[idxJ*2+1] += regPJ.y*para_gamma;
			}
		}

		static void newIteration(int iter, int call_iter, double *x, double *f, double *g, double *gnorm, void *pData)
		{
			std::cout << " Iter=" << iter << ", Call_iter=" << call_iter << ", f=" << *f << ", gnorm=" << *gnorm << "\n";
		}

	};
}

#endif//ZCVTENERGY_H_
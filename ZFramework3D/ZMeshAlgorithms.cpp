#include "ZMeshAlgorithms.h"
#include <vector>
#include <map>
#include <algorithm>
#include <list>
#include "../Eigen/spmatrixNew.h"
#include "../Eigen/speigenNew.h"
#include "ZQuadrangulation.h"
#include "MeshTools.h"

#include "ZMeshFilterManifold.h"
#include "ZMeshFilterGaussian.h"
#include "ZMeshFilterGaussian2.h"
#include "ZMeshBilateralFilter.h"
#include "AnnWarper_Eigen.h"

#include "randomGenerator.h"
#include "ZFileHelper.h"

using namespace std;

namespace ZMeshSpace
{
	//int ZMeshAlgorithms::maxEigens = 20;


	ZMeshAlgorithms::ZMeshAlgorithms(Mesh3D* mesh/* =NULL */)
	{
		maxEigens = 20;
		setMesh(mesh);
		laplacianType_ = COTANGENT;
		currentEigenIdx_ = -1;
		currentIsoLineIdx_ = -1;
		pMeshFilterManifold_ = NULL;
	}

	void ZMeshAlgorithms::reset()
	{
		eigenDatas_.clear();
		//isoLines_.clear();
		isoLineContainer_.clearData();
		currentEigenIdx_ = -1;
		currentIsoLineIdx_ = -1;
		setMesh(NULL);
		SAFE_DELETE(pMeshFilterManifold_);
	}

	void ZMeshAlgorithms::setMesh(Mesh3D* mesh)
	{
		mesh_ = mesh;
		meshProjector_.setMesh(mesh);
		isoLineContainer_.setMesh(mesh);
		if (mesh_)
			meshInitPositions_ = MeshTools::getAllVerticePos(mesh);
	}

	float ZMeshAlgorithms::laplacianWeight(HE_edge* he)
	{
		HE_vert* vi = he->m_pvert;
		HE_vert* vj = he->m_ppair->m_pvert;
		float w1 = 0, w2 = 0;
		int count = 0;
		HE_vert* vi1 = NULL;
		HE_vert* vj1 = NULL;
		switch (laplacianType_) 
		{
		case UNIFORM:
			return -1;
		case COTANGENT:
			if (he->m_pface!=NULL) 
			{
				vi1 = he->m_pnext->m_pvert;
				w1 = 1*computeCotangent(vi1->m_vpos, vi->m_vpos, vj->m_vpos); // cot laplacian
				count ++;
			}			
			if (he->m_ppair->m_pface!=NULL) 
			{
				vj1 = he->m_ppair->m_pnext->m_pvert;
				w2 = 1*computeCotangent(vj1->m_vpos, vi->m_vpos, vj->m_vpos); // cot laplacian
				count ++;
			}			
			return (w1+w2)/count;
		}
		return 0;
	}

	float ZMeshAlgorithms::laplacianWeight(HE_vert* v)
	{
		float w=0;
		HE_edge* e = v->m_pedge;
		//if (e==NULL) return 0.f;
		do 
		{
			w += e->m_laplacianWeight;
// 			if (e->m_ppair==NULL) {
// 				std::cerr << "NULL pair!" << std::endl;
// 				break;
// 			}
			e = e->m_ppair->m_pnext;
		} while (e!=NULL && e!=v->m_pedge);
		//v->m_laplacianWeight = -1.f*w;
		return -1.f*w;
	}

	void ZMeshAlgorithms::computeLaplacian()
	{
		if (mesh_==NULL) return;

		// edge laplacian
		for (int i=0; i<mesh_->get_num_of_edges_list(); i++)
		{
			HE_edge* e = mesh_->get_edge(i);
			e->m_laplacianWeight = laplacianWeight(e);
		}
		// vertices laplacian
		for (int i=0; i<mesh_->get_num_of_vertex_list(); i++)
		{
			HE_vert* v = mesh_->get_vertex(i);
			v->m_laplacianWeight = laplacianWeight(v);
		}

		// test: show the weights
// 		if (eigenDatas_.empty()) this->eigenDatas_.push_back(EigenData());
// 		EigenData& data = this->eigenDatas_[0];
// 		data.eigVec.resize(mesh_->get_num_of_vertex_list());
// 		for (int j=0; j<mesh_->get_num_of_vertex_list(); j++)
// 		{
// 			HE_vert* v = mesh_->get_vertex(j);
// 			data.eigVec[j] = v->m_laplacianWeight;
// 		}
// 		data.maxVal = *max_element(data.eigVec.begin(), data.eigVec.end());
// 		data.minVal = *min_element(data.eigVec.begin(), data.eigVec.end());
// 
// 		// normalize
// 		for (int j=0; j<mesh_->get_num_of_vertex_list(); j++)
// 		{
// 			data.eigVec[j] = (data.eigVec[j] - data.minVal)/(data.maxVal-data.minVal);
// 		}
// 		setMeshVrtColors(0);
	}

	bool ZMeshAlgorithms::saveLaplacian(const char* fileName)
	{
		if (!mesh_) return false;
		std::ofstream fout(fileName);
		fout << mesh_->get_num_of_vertex_list() << "\n";
		for (int i=0; i<mesh_->get_num_of_edges_list(); i++)
		{
			HE_edge* edge = mesh_->get_edge(i);
			fout << edge->m_pvert->m_id << " " << edge->m_ppair->m_pvert->m_id
				<< " " << edge->m_laplacianWeight << "\n";
		}
		for (int i=0; i<mesh_->get_num_of_vertex_list(); i++)
		{
			HE_vert* v = mesh_->get_vertex(i);
			fout << v->m_id << " " << v->m_id << v->m_laplacianWeight << "\n";
		}
		fout.close();
		return true;
	}

	bool ZMeshAlgorithms::computeEigens()
	{
		
		// build the laplacian matrix
		int nVrt = mesh_->get_num_of_vertex_list();
		numcNew::RowMatSym<double> rowMap(nVrt, nVrt);
		for (int i=0; i<nVrt; i++)
		{
			//map<int, double> oneRow;
			HE_vert* vert = mesh_->get_vertex(i);
			HE_edge* e = vert->m_pedge;
			int count = 0;
			do 
			{
				rowMap(i, e->m_pvert->m_id) = e->m_laplacianWeight;
				//rowMap(i, e->m_pvert->m_id) = -1.0;
				count++;
				e = e->m_ppair->m_pnext;
			} while (e!=vert->m_pedge);
			//rowMap(i, vert->m_id) = vert->m_laplacianWeight;
			rowMap(i, i) = vert->m_laplacianWeight;
			//rowMap(i, vert->m_id) = count;
		}

		// build CSRMatrix
		numcNew::CSRMatrix<numcNew::mklReal> mat;
		numcNew::CreateCSRMatrixFromRowMap(mat, rowMap);
		mat.mMtype = numcNew::CSRMatrix<numcNew::mklReal>::RealSymmIndef;
		//numcNew::RowMatSym<numcNew::mklReal> symMat(rowMap);

		// solve eigen problem
		numcNew::EigenProblem<numcNew::mklReal> solver(mat);
		numcNew::mklReal *eigenVals = new numcNew::mklReal[maxEigens];
		numcNew::mklReal *eigenVecs = new numcNew::mklReal[maxEigens*nVrt];
		//solver.setParam(numcNew::EigenProblem<numcNew::mklReal>::LM);
		solver.setParam(numcNew::EigenProblem<numcNew::mklReal>::SM);
		solver.solve(maxEigens, eigenVals, eigenVecs);

		saveLaplacian("lap.txt");
		saveEigens("eigens.txt", maxEigens, nVrt, eigenVals, eigenVecs);
		cout << "Eigens computed!\n";

		// store eigens
		eigenDatas_.resize(maxEigens);

		for (int i=0; i<maxEigens; i++)
		{
			ZEigenData& data = eigenDatas_[maxEigens-i-1];
			//EigenData& data = eigenDatas_[i];
			data.eigValue = eigenVals[i];
			data.eigVec.resize(nVrt);
			double maxValue = 0;
			for (int j=0; j<nVrt; j++)
			{
				data.eigVec[j] = eigenVecs[j+i*nVrt];
			}
			data.normalize();
		}
		cout << "Eigens saved.\n";

		currentEigenIdx_ = 0;
		setMeshVrtColors(currentEigenIdx_);

		delete [] eigenVals;
		delete [] eigenVecs;
		return true;
	}

	void ZMeshAlgorithms::setMeshVrtColors()
	{
		setMeshVrtColors(currentEigenIdx_);
	}

	void ZMeshAlgorithms::setMeshVrtColors(int eigIdx)
	{
		if (mesh_==NULL || eigIdx>=eigenDatas_.size())
			return;
		cout << eigIdx << " --> " << eigenDatas_[eigIdx].eigValue << "--maxVec: " << eigenDatas_[eigIdx].maxVal
			 << " minVec: " << eigenDatas_[eigIdx].minVal << "\n";
		int nVrt = mesh_->get_num_of_vertex_list();
		// change colors of mesh vertices
		for (int i=0; i<nVrt; i++)
		{
			HE_vert* v = mesh_->get_vertex(i);
			v->m_colorValue = eigenDatas_[eigIdx].normalizeValues[i];
		}
	}


	bool ZMeshAlgorithms::saveEigens(const char* fileName, int eSize, int nCols, double* eigenValues, double* eigenVecs)
	{
		ofstream fout(fileName);
		fout << eSize << " " << nCols << "\n";
		for (int i=0; i<eSize; i++)
		{
			fout << eigenValues[i] << "\n";
			double* vals = eigenVecs + i*nCols;
			for (int j=0; j<nCols; j++)
			{
				fout << vals[j] << " ";
			}
			fout << "\n";
		}
		fout.close();
		return true;
	}

	bool ZMeshAlgorithms::loadEigenVector(const char* fileName)
	{
		if (mesh_==NULL) return false;
		ifstream fin(fileName);
		int nVrt = mesh_->get_num_of_vertex_list();
		if (eigenDatas_.empty())
			eigenDatas_.push_back(ZEigenData());
		ZEigenData& data = eigenDatas_[0];
		data.eigVec.resize(nVrt);
		for (int i=0; i<nVrt; i++)
		{
			fin >> data.eigVec[i];
		}
		data.maxVal = *max_element(data.eigVec.begin(), data.eigVec.end());
		data.minVal = *min_element(data.eigVec.begin(), data.eigVec.end());
		for (int i=0; i<nVrt; i++)
		{
			data.eigVec[i] = (data.eigVec[i] - data.minVal)/(data.maxVal-data.minVal);
		}
		fin.close();
		setMeshVrtColors(0);
		return true;
	}
		
	bool ZMeshAlgorithms::saveEigenVectors(const char* fileName)
	{
		if (mesh_==NULL || !mesh_->isvalid()) return false;
		ofstream fout(fileName);
		int vSize = mesh_->get_num_of_vertex_list();
		fout << maxEigens << " " << vSize << "\n";
		for (int i=0; i<maxEigens; i++)
		{
			ZEigenData& data = eigenDatas_[i];
			fout << i << " " << data.eigValue << "\n";
			for (int j=0; j<vSize; j++)
			{
				fout << data.eigVec[j] << " ";
			}
			fout << "\n";
		}
		fout.close();
		return true;
	}

	bool ZMeshAlgorithms::loadEigenVectors(const char* fileName)
	{
		if (mesh_==NULL || !mesh_->isvalid()) return false;
		ifstream fin(fileName);
		int vSize;
		fin >> maxEigens >> vSize;
		if (vSize!=mesh_->get_num_of_vertex_list()) return false;

		eigenDatas_.resize(maxEigens);
		for (int i=0; i<maxEigens; i++)
		{
			ZEigenData& data = eigenDatas_[i];
			data.eigVec.resize(vSize, 0);
			int dummy;
			fin >> dummy >> data.eigValue;
			for (int j=0; j<vSize; j++)
			{
				fin >> data.eigVec[j];
			}
			data.normalize();
		}
		fin.close();
		setMeshVrtColors(0);
		return true;
	}
// 
// 	void ZMeshAlgorithms::EigenData::normalize()
// 	{
// 		maxVal = *max_element(eigVec.begin(), eigVec.end());
// 		minVal = *min_element(eigVec.begin(), eigVec.end());
// 		normalizeValues.resize(eigVec.size());
// 		for (int i=0; i<eigVec.size(); i++)
// 		{
// 			normalizeValues[i] = (eigVec[i]-minVal)/(maxVal-minVal);
// 		}
// 	}

// 	void ZMeshAlgorithms::EigenData::extractIsoBands(Mesh3D* mesh)
// 	{
// 		if (mesh==NULL) return;
// 		if (eigVec.empty() || eigVec.size()!= mesh->get_num_of_vertex_list()) return;
// 
// 		int bandSize = 10;
// 		double bandWidth = (maxVal - minVal)/bandSize;
// 		int nVrt = mesh->get_num_of_vertex_list();
// 
// 		isoBands.resize(bandSize+1);
// 		normalizeValues.resize(nVrt);
// 		for (int i=0; i<bandSize+1; i++)
// 		{
// 			isoBands[i].bandMin = bandWidth*i + minVal;
// 			isoBands[i].bandMax = bandWidth*(i+1) + minVal;
// 		}
// 		for (int i=0; i<nVrt; i++)
// 		{
// 			double s = eigVec[i];
// 			int bandId = (s-minVal)/bandWidth;
// 			isoBands[bandId].vrtIdxs.push_back(i);
// 			normalizeValues[i] = 1.0*bandId/bandSize;
// 			//mesh->get_vertex(i)->m_colorValue = normalizeValues[i];
// 		}
// 	}

// 	void ZMeshAlgorithms::computeIsoLine(HE_vert* startVert, double value, ZEigenIsoLine& outLine)
// 	{
// 		EigenData& data = eigenDatas_[getCurrentEigenIdx()];
// 
// 		std::vector<bool> traversedF(mesh_->get_num_of_faces_list(), false);
// 		std::vector<bool> traversedV(mesh_->get_num_of_vertex_list(), false);
// 		ZMeshSurfacePoint startP;
// 		startP.pos = startVert->m_vpos;
// 		startP.edge = startVert->m_pedge;
// 		startP.face = startVert->m_pedge->m_pface;
// 		startP.vert = startVert;
// 		startP.type = VERTEX_POINT;
// 		//traversed[startVert->m_id] = true;
// 
// 		outLine.linePoints.push_back(startP);
// 		ZMeshSurfacePoint nextP(startP);
// 
// 		bool stop = false;
// 		do 
// 		{
// 			// find next surface point with value 
// 			if (nextP.type==VERTEX_POINT)
// 			{
// 				HE_edge* e = nextP.vert->m_pedge;
// 				HE_edge* prevE = e->m_pprev->m_ppair;
// 				if (prevE==NULL)
// 				{
// 					prevE = e;
// 					e = e->m_ppair->m_pnext;
// 				}
// 				//e = e->m_ppair->m_pnext;
// 				HE_face* f = prevE->m_pface;
// 				bool bFind = false;
// 				do 
// 				{
// 					double dPrev = data.eigVec[prevE->m_pvert->m_id];
// 					double dV = data.eigVec[e->m_pvert->m_id];
// 					double d = (dPrev-value)*(dV-value);
// 					ZMeshSurfacePoint p;
// 					if (g_isZero(dPrev-value) && !traversedV[prevE->m_pvert->m_id])
// 					{
// 						p.type = VERTEX_POINT;
// 						p.vert = prevE->m_pvert;
// 						p.pos = prevE->m_pvert->m_vpos;
// 						traversedV[prevE->m_pvert->m_id] = true;
// 						bFind = true;
// 					}
// 					else if (g_isZero(dV-value) && !traversedV[e->m_pvert->m_id])
// 					{
// 						p.type = VERTEX_POINT;
// 						p.vert = e->m_pvert;
// 						p.pos = e->m_pvert->m_vpos;
// 						traversedV[e->m_pvert->m_id] = true;
// 						bFind = true;
// 					}
// 					else if (d<0 && !traversedF[e->m_pface->m_id])
// 					{
// 						p.type = EDGE_POINT;
// 						double ddPrev = abs(dPrev-value);
// 						double ddV = abs(dV-value);
// 						p.pos = e->m_pvert->m_vpos * (ddPrev/(ddPrev+ddV))
// 							+ prevE->m_pvert->m_vpos * (ddV/(ddPrev+ddV));
// 						p.edge = findEdge(prevE->m_pvert, e->m_pvert);
// 						traversedF[e->m_pface->m_id] = true;
// 						bFind = true;
// 					}
// 					if (bFind)
// 					{
// 						outLine.linePoints.push_back(p);
// 						nextP = p;
// 						break;
// 					}
// 					prevE = e;
// 					e = e->m_ppair->m_pnext;
// 				} while (e!=nextP.vert->m_pedge);
// 				if (!bFind) stop = true;
// 			}
// 			else if (nextP.type==EDGE_POINT)
// 			{
// 				HE_edge* e = nextP.edge;
// 				if (e->m_pface && !traversedF[e->m_pface->m_id])
// 				{
// 					ZMeshSurfacePoint p;
// 					findIsoLineInTriangle(e, data.eigVec, value, nextP, p);
// 					outLine.linePoints.push_back(p);
// 					nextP = p;
// 					if (p.vert!=NULL)
// 						traversedV[p.vert->m_id] = true;
// 					traversedF[e->m_pface->m_id] = true;
// 				}
// 				else if (e->m_ppair->m_pface && !traversedF[e->m_ppair->m_pface->m_id])
// 				{
// 					ZMeshSurfacePoint p;
// 					findIsoLineInTriangle(e->m_ppair, data.eigVec, value, nextP, p);
// 					outLine.linePoints.push_back(p);
// 					nextP = p;
// 					if (p.vert!=NULL)
// 						traversedV[p.vert->m_id] = true;
// 					traversedF[e->m_ppair->m_pface->m_id] = true;
// 				}
// 				else
// 				{
// 					stop = true;
// 				}
// 			}
// 			if (nextP.vert==startVert) stop = true;
// 			//else if (outLine.linePoints.size()>200) stop = true;
// 		} while (!stop);
// 
// 		//std::cout << "IsoLine, Point Count: " << outLine.linePoints.size() << "\n";
// 	}
// 
// 	void ZMeshAlgorithms::computeIsoLine(ZMeshSurfacePoint& startP, double value, ZEigenIsoLine& outLine)
// 	{
// 		//std::cout << "--start from edge..\n";
// 		EigenData& data = eigenDatas_[getCurrentEigenIdx()];
// 
// 		std::vector<bool> traversedF(mesh_->get_num_of_faces_list(), false);
// 		std::vector<bool> traversedV(mesh_->get_num_of_vertex_list(), false);
// 
// 		outLine.linePoints.push_back(startP);
// 		ZMeshSurfacePoint nextP(startP);
// 
// 		bool stop = false;
// 		do 
// 		{
// 			// find next surface point with value 
// 			if (nextP.type==VERTEX_POINT)
// 			{
// 				HE_edge* e = nextP.vert->m_pedge;
// 				HE_edge* prevE = e->m_pprev->m_ppair;
// 				if (prevE==NULL)
// 				{
// 					prevE = e;
// 					e = e->m_ppair->m_pnext;
// 				}
// 				//e = e->m_ppair->m_pnext;
// 				HE_face* f = prevE->m_pface;
// 				bool bFind = false;
// 				do 
// 				{
// 					double dPrev = data.eigVec[prevE->m_pvert->m_id];
// 					double dV = data.eigVec[e->m_pvert->m_id];
// 					double d = (dPrev-value)*(dV-value);
// 					ZMeshSurfacePoint p;
// 					if (g_isZero(dPrev-value) && !traversedV[prevE->m_pvert->m_id])
// 					{
// 						p.type = VERTEX_POINT;
// 						p.vert = prevE->m_pvert;
// 						p.pos = prevE->m_pvert->m_vpos;
// 						traversedV[prevE->m_pvert->m_id] = true;
// 						bFind = true;
// 					}
// 					else if (g_isZero(dV-value) && !traversedV[e->m_pvert->m_id])
// 					{
// 						p.type = VERTEX_POINT;
// 						p.vert = e->m_pvert;
// 						p.pos = e->m_pvert->m_vpos;
// 						traversedV[e->m_pvert->m_id] = true;
// 						bFind = true;
// 					}
// 					else if (d<0 && !traversedF[e->m_pface->m_id])
// 					{
// 						p.type = EDGE_POINT;
// 						double ddPrev = abs(dPrev-value);
// 						double ddV = abs(dV-value);
// 						p.pos = e->m_pvert->m_vpos * (ddPrev/(ddPrev+ddV))
// 							+ prevE->m_pvert->m_vpos * (ddV/(ddPrev+ddV));
// 						p.edge = findEdge(prevE->m_pvert, e->m_pvert);
// 						traversedF[e->m_pface->m_id] = true;
// 						bFind = true;
// 					}
// 					if (bFind)
// 					{
// 						outLine.linePoints.push_back(p);
// 						nextP = p;
// 						break;
// 					}
// 					prevE = e;
// 					e = e->m_ppair->m_pnext;
// 				} while (e!=nextP.vert->m_pedge);
// 				if (!bFind) stop = true;
// 			}
// 			else if (nextP.type==EDGE_POINT)
// 			{
// 				HE_edge* e = nextP.edge;
// 				if (e->m_pface && !traversedF[e->m_pface->m_id])
// 				{
// 					ZMeshSurfacePoint p;
// 					findIsoLineInTriangle(e, data.eigVec, value, nextP, p);
// 					outLine.linePoints.push_back(p);
// 					nextP = p;
// 					if (p.vert!=NULL)
// 						traversedV[p.vert->m_id] = true;
// 					traversedF[e->m_pface->m_id] = true;
// 				}
// 				else if (e->m_ppair->m_pface && !traversedF[e->m_ppair->m_pface->m_id])
// 				{
// 					ZMeshSurfacePoint p;
// 					findIsoLineInTriangle(e->m_ppair, data.eigVec, value, nextP, p);
// 					outLine.linePoints.push_back(p);
// 					nextP = p;
// 					if (p.vert!=NULL)
// 						traversedV[p.vert->m_id] = true;
// 					traversedF[e->m_ppair->m_pface->m_id] = true;
// 				}
// 				else
// 				{
// 					stop = true;
// 				}
// 			}
// 			if (nextP.edge==startP.edge || (nextP.vert!=NULL && nextP.vert==startP.vert)) stop = true;
// 			//else if (outLine.linePoints.size()>200) stop = true;
// 		} while (!stop);
// 	}

	void ZMeshAlgorithms::computeIsoLines(int idx, double value)
	{

		ZEigenData& data = eigenDatas_[idx];
		if (abs(data.eigValue)<0.00001) 
		{
			cerr << "Eigen value is 0: " << data.eigValue << "\n";
			return;
		}
		isoLineContainer_.computeIsoLines(data.eigVec, value);
	}

	void ZMeshAlgorithms::computeProjectionPlanes()
	{
		isoLineContainer_.computeProjectionPlanes();
	}

	const ZEigenIsoLine& ZMeshAlgorithms::getIsoLine(int idx)
	{
		return isoLineContainer_.getIsoLine(idx);
	}

	int ZMeshAlgorithms::getIsoLineSize()
	{
		return isoLineContainer_.getIsoLineSize();
	}

	void ZMeshAlgorithms::clearIsoLineData()
	{
		isoLineContainer_.clearData();
	}

	void ZMeshAlgorithms::intersectMeshByPlane(ZPlane* plane)
	{

	}

	bool ZMeshAlgorithms::initQuadrangulation()
	{
		std::vector<Vec3f> localU(mesh_->get_num_of_vertex_list(), Vec3f(1,0,0));
		std::vector<Vec3f> localV(mesh_->get_num_of_vertex_list(), Vec3f(0,1,0));
		std::vector<double> initScaleFiled;
		if(!eigenDatas_.empty())
		{
			for (int i=0; i<eigenDatas_[4].eigVec.size(); i++)
			{
				//initScaleFiled.push_back(eigenDatas_[4].eigVec[i]);
				initScaleFiled.push_back(1);
			}
		}
		meshQuadrangulation.init(mesh_,localU, localV);
		meshQuadrangulation.setVrtScaleFiled(initScaleFiled);

		return true;
	}

	bool ZMeshAlgorithms::runQuadrangulation(float faceScale, float lambda)
	{
		//ZQuadrangulation algo;
// 		algo.testingEigen();
		std::cout << "EigenValue: " << eigenDatas_[4].eigValue << "\n";
		std::cout << "Lambda: " << lambda << "\n";
		float lambda0 = 1.0/sqrt(abs(eigenDatas_[4].eigValue))*lambda;
		meshQuadrangulation.setLambda(lambda0);
		std::vector<float> faceFactors(mesh_->get_num_of_faces_list(), faceScale);
		meshQuadrangulation.setFaceFactors(faceFactors);
		meshQuadrangulation.go();
		//meshQuadrangulation.smoothScaleField();
		meshQuadrangulation.tagVertices();
		meshQuadrangulation.generateMSEdges();
		MeshTools::changeMeshVrtColors(mesh_, meshQuadrangulation.getVrtScaleField());

		return true;
	}

	bool ZMeshAlgorithms::runBilateralFiltering(float sigma_r, float sigma_s)
	{

		int nSize = mesh_->get_num_of_faces_list();
		int spatialDim = 3;
		int rangeDim = 3;

		Eigen::MatrixXf faceNormals = MeshTools::getAllFaceNormals(mesh_);
		//Eigen::MatrixXf faceNormals = MeshTools::getAllFaceNormalsSpherical(mesh_);
		Eigen::MatrixXf faceCenters = MeshTools::getAllFaceCenters(mesh_);
		Eigen::MatrixXf verticePos = MeshTools::getAllVerticePos(mesh_);
		Eigen::MatrixXf faceAttributes(nSize, spatialDim+rangeDim);
		faceAttributes.block(0, 0, nSize, spatialDim) = faceCenters;
		faceAttributes.block(0, spatialDim, nSize, rangeDim) = faceNormals;
		AnnWarper_Eigen annSearch;
		annSearch.init(faceCenters);
		ZMeshBilateralFilter filter(mesh_);
		filter.setAnnSearchHandle(&annSearch);
		float sigma_spatial = mesh_->m_minEdgeLength*sigma_s;
		float sigma_range = cos(sigma_r*Z_PI/180.f);
		//filter.setKernelFunc(NULL);
		filter.setPara(ZMeshFilterParaNames::SpatialSigma, sigma_spatial);
		filter.setPara(ZMeshFilterParaNames::RangeSigma, sigma_range);
		filter.apply(faceAttributes, std::vector<bool>(nSize, true));
		Eigen::MatrixXf output = filter.getResult();

		MeshTools::setAllFaceNormals2(mesh_, output.block(0, 3, nSize, 3));
		//MeshTools::setAllFaceNormal2Spherical(mesh_, output.block(0, spatialDim, nSize, rangeDim));

		return true;
	}

	bool ZMeshAlgorithms::runManifoldSmooth(float sigma_r, float sigma_s)
	{
		SAFE_DELETE(pMeshFilterManifold_);
		int nSize = mesh_->get_num_of_faces_list();
		//int nSize = mesh_->get_num_of_vertex_list();
		int spatialDim = 3;
		int rangeDim = 3;
		Eigen::MatrixXf faceNormals = MeshTools::getAllFaceNormals(mesh_);
		//Eigen::MatrixXf faceNormals = MeshTools::getAllFaceNormalsSpherical(mesh_);
		Eigen::MatrixXf faceCenters = MeshTools::getAllFaceCenters(mesh_);
		Eigen::MatrixXf verticePos = MeshTools::getAllVerticePos(mesh_);
		Eigen::MatrixXf faceAttributes(nSize, spatialDim+rangeDim);
		faceAttributes.block(0, 0, nSize, spatialDim) = faceCenters;
		faceAttributes.block(0, 3, nSize, rangeDim) = faceNormals;
// 		faceAttributes.block(0, 0, nSize, spatialDim) = verticePos;
// 		faceAttributes.block(0, 3, nSize, rangeDim) = verticePos;
		//ZMeshFilterManifold filter(mesh_);
		pMeshFilterManifold_ = new ZMeshFilterManifold(mesh_);
		float sigma_spatial = mesh_->m_minEdgeLength*sigma_s;
		float sigma_range = sin(sigma_r*Z_PI/180.f);
		pMeshFilterManifold_->setPara(ZMeshFilterParaNames::SpatialSigma, sigma_spatial);
		pMeshFilterManifold_->setPara(ZMeshFilterParaNames::RangeSigma, sigma_range);
		pMeshFilterManifold_->setRangeDim(rangeDim);
		pMeshFilterManifold_->setSpatialDim(spatialDim);
		CHECK_FALSE_AND_RETURN(pMeshFilterManifold_->apply(faceAttributes, std::vector<bool>(nSize, true)));
		Eigen::MatrixXf output = pMeshFilterManifold_->getResult();
		meshNewNormals_ = output.block(0, spatialDim, nSize, rangeDim);
		MeshTools::setAllFaceNormals2(mesh_, meshNewNormals_);
		MeshTools::setAllFaceColorValue(mesh_, pMeshFilterManifold_->getPointClusterIds());
		//MeshTools::setAllVerticePos(mesh_, output.block(0, spatialDim, nSize, rangeDim));
		//MeshTools::setAllFaceNormal2Spherical(mesh_, output.block(0, spatialDim, nSize, rangeDim));

		ZFileHelper::saveEigenMatrixToFile("oldNormal.txt", faceNormals);
		ZFileHelper::saveEigenMatrixToFile("newNormal.txt", output.block(0,spatialDim, nSize, rangeDim));

		return true;
	}

	void ZMeshAlgorithms::resetMeshPositions()
	{
		MeshTools::setAllVerticePos(mesh_, meshInitPositions_);
		mesh_->update_mesh();
	}

	void ZMeshAlgorithms::setRangeLevelWeight(int level, float weight)
	{
		if (pMeshFilterManifold_)
		{
			pMeshFilterManifold_->setRangeWeight(level, weight);
		}
	}

	void ZMeshAlgorithms::updateRange()
	{
		pMeshFilterManifold_->updateRange();
		Eigen::MatrixXf output = pMeshFilterManifold_->getResult();
		meshNewNormals_ = output.block(0, pMeshFilterManifold_->getSpatialDim(), pMeshFilterManifold_->getPointSize(), pMeshFilterManifold_->getRangeDim());
		MeshTools::setAllFaceNormals2(mesh_, meshNewNormals_);
		//MeshTools::setAllFaceNormal2Spherical(mesh_, meshNewNormals_);
	}

	void ZMeshAlgorithms::updateVerticesPos()
	{
		meshNewPositions_ = computeNewPositions(meshNewNormals_);
		MeshTools::setAllVerticePos(mesh_, meshNewPositions_);
		mesh_->update_mesh();
	}

	Eigen::MatrixXf ZMeshAlgorithms::computeNewPositions(const Eigen::MatrixXf& newNormals)
	{
		int vSize = mesh_->get_num_of_vertex_list();
		Eigen::MatrixXf newPositions(vSize, 3);
		newPositions.fill(0);
		const int maxIterTime = 1;
		float moveDist = 0;
		for (int iterTime=0; iterTime<maxIterTime; iterTime++)
		{
			for (int i=0; i<vSize; i++)
			{
				HE_vert* v = mesh_->get_vertex(i);
				HE_edge* e = v->m_pedge;
				if (e==NULL) continue;
				Eigen::VectorXf oldPos(3);
				oldPos(0) = v->m_vpos.x; oldPos(1) = v->m_vpos.y; oldPos(2) = v->m_vpos.z;
				Eigen::VectorXf newPos(3);
				newPos.fill(0);
				int count = 0;
				do 
				{
					HE_face* f = e->m_pface;
					if (f!=NULL) 
					{
						//Eigen::VectorXf newN(newNormals.row(f->m_id));
						Eigen::VectorXf newN(3);
						newN(0) = f->m_vnormal2.x; newN(1) = f->m_vnormal2.y; newN(2) = f->m_vnormal2.z;
						float weight = f->m_vnormal2.DotProduct(f->m_vCenter-v->m_vpos);//newN.dot(Eigen::VectorXf(faceCenters.row(f->m_id))-oldPos);
						newPos += newN*weight;
						count ++;
					}
					e = e->m_ppair->m_pnext;
				} while (e!=v->m_pedge && e!=NULL);
				//count = 6;
				//if (i<50)
				//	std::cout << i << " " << newPos.norm() << "\n";
				newPos = oldPos + newPos*(1.f/count);
				moveDist += newPos.norm()/count;
				newPositions.row(i) = newPos;
			}
		}
		
		moveDist /= vSize;
		//std::ofstream fout("positions.txt");
		//fout<<newPositions;
		//fout.close();
		std::cout << "Moved: " << moveDist << "\n";
		return newPositions;
	}

	void ZMeshAlgorithms::addNoiseToMesh(float sigma)
	{
		if (mesh_==NULL)
			return;
	
		int vSize = mesh_->get_num_of_vertex_list();
		float edgeLength = mesh_->m_averageEdgeLength;
		float sigmaF = edgeLength*sigma;
		Random_Generator::Start();
		for (int i=0; i<vSize; i++)
		{
			HE_vert* v = mesh_->get_vertex(i);
			Vec3f n = v->m_vnormal;
			float randomF = Random_Generator::NormalRandom(0, sigmaF, -edgeLength, edgeLength);//(1.f*(rand()%10000)/10000-0.5)*edgeLength*sigma;
			v->m_vpos = v->m_vpos + v->m_vnormal*randomF;
		}
		mesh_->update_mesh();
	}

	void ZMeshAlgorithms::updateVerticePos(const Eigen::MatrixXf& newNormals, const Eigen::MatrixXf& faceCenters)
	{
		MeshTools::setAllVerticePos(mesh_, computeNewPositions(newNormals));
		mesh_->update_mesh();
	}

	void ZMeshAlgorithms::test0()
	{
		int idx = getCurrentEigenIdx();
		ZEigenData& data = eigenDatas_[idx];
		//isoLines_.clear();
// 		for (int i=0; i<10; i++)
// 		{
// 			//double dValue = (data.maxVal+data.minVal)/2;
// 			double dValue = data.minVal + (data.maxVal-data.minVal)/9*i;
// 			computeIsoLines(idx, dValue);		
// 		}
		isoLineContainer_.clear();
		isoLineContainer_.extract(data);

	}

	void ZMeshAlgorithms::test1()
	{
		isoLineContainer_.clear();
		isoLineContainer_.extract2(eigenDatas_);
	}
}


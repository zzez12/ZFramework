#include "ZQuadrangulation.h"
#include "MeshTools.h"
#include <map>
#include <fstream>

namespace ZMeshSpace
{
	// inner class
	void ZQuadrangulation::MeshExtData::clear()
	{
		vrtData_.clear();
		edgeData_.clear();
	}

	//////////////////////////////////////////////////////////////////////////
	// class ZQuadrangulation
	ZQuadrangulation::ZQuadrangulation(Mesh3D* mesh)
	{
		pMesh_ = mesh;
		paras_.omega1_ = 0.5;
		paras_.omega2_ = 1.f;
		paras_.lambda_ = 1.f;	// test value
		paras_.threshold_stopIteration_ = 1.f*10e-7;

	}

	ZQuadrangulation::~ZQuadrangulation()
	{
		destroy();
	}

	void ZQuadrangulation::destroy()
	{
		pMesh_ = NULL;
		solver_.clear();
		meshExtData_.clear();
		expMapGen_.destroy();
	}

	void ZQuadrangulation::buildVertex2RingNeighbor()
	{
		if (pMesh_==NULL || !pMesh_->isvalid())
			return;

		int vSize = pMesh_->get_num_of_vertex_list();
		meshExtData_.resetNeighbors(vSize);
		for (int i=0; i<vSize; i++)
		{
			MeshVrt2RingNeighbor* vrtData = meshExtData_.getVrt2RingNeighbor(i);
			HE_vert* vert = pMesh_->get_vertex(i);
			HE_edge* edge = vert->m_pedge;
			do 
			{
				HE_vert* innerVrt = edge->m_pvert;
				vrtData->addNeighborVrt(innerVrt);
				HE_edge* innerEdge = innerVrt->m_pedge;
				do 
				{
					if (innerEdge->m_pvert!=vert)
						vrtData->addNeighborVrt(innerEdge->m_pvert);
					innerEdge = innerEdge->m_ppair->m_pnext;
				} while (innerEdge!=innerVrt->m_pedge && innerVrt!=NULL);
				edge = edge->m_ppair->m_pnext;
			} while (edge!=vert->m_pedge && edge!=NULL);
		}
	}

	void ZQuadrangulation::init(Mesh3D* mesh, const std::vector<Vec3f>& localU, const std::vector<Vec3f>& localV, const std::vector<float>& faceFactors)
	{
		destroy();
		pMesh_ = mesh;
		expMapGen_.setMesh(pMesh_);
		buildVertex2RingNeighbor();

		initExtData(localU, localV, faceFactors);
	}

	void ZQuadrangulation::setLocalUV(const std::vector<Vec3f>& localU, const std::vector<Vec3f>& localV)
	{
		//initExtData(localU, localV);
		resetLocalCoordinates(localU, localV);
	}

	void ZQuadrangulation::initExtData(const std::vector<Vec3f>& localU, const std::vector<Vec3f>& localV, const std::vector<float>& faceFactors)
	{
		if (!pMesh_ || !pMesh_->isvalid())
			return;

		if (localU.size()!=pMesh_->get_num_of_vertex_list() || localV.size()!=pMesh_->get_num_of_vertex_list())
			return;

		meshExtData_.clear();
		// alloc spaces
		int vSize = pMesh_->get_num_of_vertex_list();
		int eSize = pMesh_->get_num_of_edges_list();
		int fSize = pMesh_->get_num_of_faces_list();
		meshExtData_.resetData(vSize, eSize, fSize);
		setFaceFactors(faceFactors);
		resetVrtAreas();
		resetLocalCoordinates(localU, localV);
		resetLaplacianCoord();
		std::cout << "Local coordinate done!\n";
	}

	void ZQuadrangulation::resetLaplacianCoord()
	{
		int eSize = pMesh_->get_num_of_edges_list();
		for (int i=0; i<eSize; i++)
		{
			MeshEdgeExtData* edgeData = meshExtData_.getEdgeAttri(i);
			edgeData->id_ = i;
			edgeData->cotAlpha_ = MeshTools::getLaplacianCotValue(pMesh_->get_edge(i));
		}
		for (int i=0; i<eSize; i++)
		{
			MeshEdgeExtData* edgeData = meshExtData_.getEdgeAttri(i);
			MeshEdgeExtData* edgePairData = meshExtData_.getEdgeAttri(pMesh_->get_edge(i)->m_ppair->m_id);
			edgeData->laplacianW_ = (edgeData->cotAlpha_+edgePairData->cotAlpha_);
		}
	}

	void ZQuadrangulation::setFaceFactors(const std::vector<float>& faceFactors)
	{
		if (!pMesh_) return;
		int fSize = pMesh_->get_num_of_faces_list();
		std::vector<float> fFacetors = faceFactors;
		if (faceFactors.empty()) fFacetors.resize(fSize, 1.f);
		if (fSize!=faceFactors.size()) return;
		if (meshExtData_.getFaceAttriSize()!=fSize) meshExtData_.resizeFaceAttris(fSize);
		for (int i=0; i<fSize; i++)
		{
			MeshFaceExtData* face = meshExtData_.getFaceAttri(i);
			face->id_ = i;
			face->scaleFactor_ = fFacetors[i];
		}
		resetVrtAreas();
		resetLocalTransformation();
	}

	void ZQuadrangulation::resetVrtAreas()
	{
		int vSize = pMesh_->get_num_of_vertex_list();
		for (int i=0; i<vSize; i++)
		{
			HE_vert* vert = pMesh_->get_vertex(i);
			HE_edge* edge = vert->m_pedge;
			MeshVrtLocalCoord* vrtData = (MeshVrtLocalCoord*)meshExtData_.getVrtAttri(i);

			float areas = 0.f;
			do
			{
				HE_face* f = edge->m_pface;
				areas += MeshTools::getFaceArea(f)*meshExtData_.getFaceAttri(f->m_id)->scaleFactor_;
				edge = edge->m_ppair->m_pnext;
			} while (edge!=vert->m_pedge && edge!=NULL);
			vrtData->area_ = areas/3;
			vrtData->area_sqrt_ = sqrt(vrtData->area_);	
		}


	}

	void ZQuadrangulation::resetLocalCoordinates(const std::vector<Vec3f>& localU, const std::vector<Vec3f>& localV)
	{
		// alloc spaces
		int vSize = pMesh_->get_num_of_vertex_list();
		int eSize = pMesh_->get_num_of_edges_list();
		int fSize = pMesh_->get_num_of_faces_list();

		// assign local frame to the vertices
		int vrtAddSize = 0;
		for (int i=0; i<vSize; i++) 
		{
			HE_vert* v = pMesh_->get_vertex(i);
			Vec3f vPos = v->m_vpos;
			MeshVrtLocalCoord* vrt = (MeshVrtLocalCoord*)meshExtData_.getVrtAttri(i);

			vrt->type_ = VERTEX_DEFAULT;
			vrt->id_ = i;
			// project the directions to the tangent plane
			//vrt->directionU_ = localU[i];
			//vrt->directionV_ = localV[i];
			Vec3f vNormal = v->m_vnormal;
			Vec3f U = localU[i];
			Vec3f V = localV[i];
			Vec3f vProjU = U - vNormal*(U.DotProduct(vNormal));
			Vec3f vProjV = V - vNormal*(V.DotProduct(vNormal));
			if (!g_isZero(vProjU.length()))
			{
				vProjU.unify();
				vrt->directionU_ = vProjU;
				vrt->directionV_ = vNormal*vProjU;
			}
			else if (!g_isZero(vProjV.length()))
			{
				vProjV.unify();
				vrt->directionU_ = vProjV*vNormal;
				vrt->directionV_ = vProjV;
			}
			vrt->directionU_.unify();
			vrt->directionV_.unify();
 		}
		std::cout << "Begin to build local parameterization..\n";

#define NOT_USE_EXPMAP 0
#if NOT_USE_EXPMAP
		for (int i=0; i<vSize; i++) 
		{
			HE_vert* v = pMesh_->get_vertex(i);
			Vec3f vPos = v->m_vpos;
			MeshVrtLocalCoord* vrt = (MeshVrtLocalCoord*)meshExtData_.getVrtAttri(i);
			HE_edge* edge = v->m_pedge;


			// add current vertices
			vrt->addNeighborVrt(i, Vec3f(0,0,0));

			// get 1-ring vertices
			edge = v->m_pedge;
			do 
			{
				int vId = edge->m_pvert->m_id;
				Vec3f vDir = edge->m_pvert->m_vpos - vPos;
				vrt->addNeighborVrt(vId, vDir);

				edge = edge->m_ppair->m_pnext;
			} while (edge!=v->m_pedge && edge!=NULL);
			// make sure the neighbors have at least 6 vertices
			if (vrt->neighborIds_.size()<8)
			{
				bool bFinish = false;
				edge = v->m_pedge;
				do 
				{
					HE_vert* innerVrt = edge->m_pvert;
					HE_edge* innerEdge = innerVrt->m_pedge;
					do 
					{
						int nId = innerEdge->m_pvert->m_id;
						if (std::find(vrt->neighborIds_.begin(), vrt->neighborIds_.end(), nId)==vrt->neighborIds_.end())
						{
							Vec3f vDir = innerEdge->m_pvert->m_vpos - vPos;
							vrt->addNeighborVrt(nId, vDir);
						}
						innerEdge = innerEdge->m_ppair->m_pnext;
						if (vrt->neighborIds_.size()>=8)
						{
							bFinish = true;
							break;
						}
					} while (innerEdge!=innerVrt->m_pedge && innerEdge!=NULL);
					edge = edge->m_ppair->m_pnext;
					if (bFinish) break;
				} while (edge!=v->m_pedge && edge!=NULL);

				if (vrt->neighborIds_.size()<6)
				{
					std::cerr << "Error: Still less than 6 neighbors!" << std::endl;
				}
			}
		}
#else
		for (int i=0; i<vSize; i++) 
		{
			if (i%40==0)
				std::cout << (char)13 << i*100/vSize << "%..." << i << "/" << vSize;
			HE_vert* v = pMesh_->get_vertex(i);
			Vec3f vPos = v->m_vpos;
			MeshVrtLocalCoord* vrt = (MeshVrtLocalCoord*)meshExtData_.getVrtAttri(i);
			HE_edge* edge = v->m_pedge;

			//std::cout << "update uv..\n";
			expMapGen_.setHitVertex(v);
			expMapGen_.validateExpMap();
			//std::cout << " done!\n";
			// add current vertices
			//vrt->addNeighborVrt(i, expMapGen_.getUV(i));
			vrt->addNeighborVrt(i, Vec2f(0,0));

			// get 1-ring vertices
			edge = v->m_pedge;
			do 
			{
				int vId = edge->m_pvert->m_id;
				Vec3f vDir = edge->m_pvert->m_vpos - vPos;
				vrt->addNeighborVrt(vId, expMapGen_.getUV(vId));

				edge = edge->m_ppair->m_pnext;
			} while (edge!=v->m_pedge && edge!=NULL);
			// make sure the neighbors have at least 6 vertices
			if (vrt->neighborIds_.size()<8)
			{
				bool bFinish = false;
				edge = v->m_pedge;
				do 
				{
					HE_vert* innerVrt = edge->m_pvert;
					HE_edge* innerEdge = innerVrt->m_pedge;
					do 
					{
						int nId = innerEdge->m_pvert->m_id;
						if (std::find(vrt->neighborIds_.begin(), vrt->neighborIds_.end(), nId)==vrt->neighborIds_.end())
						{
							Vec3f vDir = innerEdge->m_pvert->m_vpos - vPos;
							vrt->addNeighborVrt(nId, expMapGen_.getUV(nId));
						}
						innerEdge = innerEdge->m_ppair->m_pnext;
						if (vrt->neighborIds_.size()>=8)
						{
							bFinish = true;
							break;
						}
					} while (innerEdge!=innerVrt->m_pedge && innerEdge!=NULL);
					edge = edge->m_ppair->m_pnext;
					if (bFinish) break;
				} while (edge!=v->m_pedge && edge!=NULL);

				if (vrt->neighborIds_.size()<6)
				{
					std::cerr << "Error: Still less than 6 neighbors!" << std::endl;
				}
			}
		}
		std::cout << vSize << " neighbor done!\n";
#endif

		
	}

	void ZQuadrangulation::resetLocalTransformation()
	{
		int vSize = pMesh_->get_num_of_vertex_list();
		//std::ofstream fout("neighborTrans.txt");
		// compute the operators (transformation matrix)
		for (int i=0; i<vSize; i++)
		{
			MeshVrtLocalCoord* vrt = (MeshVrtLocalCoord*)meshExtData_.getVrtAttri(i);
			int neighborSize = vrt->neighborIds_.size();
			// matrix U (dim:6*neighSize)
			//std::vector<std::vector<float>> matU(6, std::vector<float>(neighborSize, 0));
			Eigen::MatrixXf matU(6, neighborSize);
			Eigen::MatrixXf matD(neighborSize, neighborSize);
			matD.fill(0);
			for (int y=0; y<neighborSize; y++)
			{
				int neighborId = vrt->neighborIds_[y];
				//float sqrtArea = sqrt(meshExtData_.getVrtAttri(neighborId)->area_);
				float u = vrt->neighborCoords_[y][0];
				float v = vrt->neighborCoords_[y][1];
				matU(0,y) = 1;
				matU(1,y) = u;
				matU(2,y) = v;
				matU(3,y) = u*u*0.5f;
				matU(4,y) = u*v;
				matU(5,y) = v*v*0.5f;
				matD(y,y) = meshExtData_.getVrtAttri(neighborId)->area_;
			}
			Eigen::Matrix<float, 6, 6> m6 = matU*matD*(matU.transpose());
			vrt->neighborTransMat_ = m6.inverse()*matU*matD;
			//fout << i << " U=\n" << matU << "\n";
			//fout << i << " D=\n" << matD << "\n";
			//fout << i << " M=\n" << vrt->neighborTransMat_ << "\n";
		}
		//fout.close();
	}

	void ZQuadrangulation::go()
	{
		int maxIter = 20;
		//generateInitField();
		//initLinearSystem();
		initLinearSystem2();
		vrtScaleField_.fill(1);
		for (int i=0; i<maxIter; i++)
		{
			float err = computeOnce();
			if (err<paras_.threshold_stopIteration_)
				break;
			std::cout << " ..Iteration: " << i << " Error: " << err << "\n";
		}
		bool bSaveRetToFile = true;
		if (bSaveRetToFile)
		{
			std::ofstream fout("field.txt");
			fout << vrtScaleField_ ;
			fout.close();
		}
	}

	void ZQuadrangulation::generateInitField()
	{
		// using setVrtScaleField() for initialization now (testing)
		recomputeQuadraticCoordinate();
	}

	void ZQuadrangulation::initLinearSystem()
	{
		// prepare matrix
		int nSize = pMesh_->get_num_of_vertex_list();
		numcNew::CSRMatrix<numcNew::mklReal> mat;
		numcNew::RowMat<numcNew::mklReal> rowC(nSize, nSize);
		for (int i=0; i<nSize; i++)
		{
			//std::map<int, double> oneRow;
			//MeshVrtLocalCoord* vrtData = (MeshVrtLocalCoord*)meshExtData_.getVrtAttri(i);
// 			HE_vert* vert = pMesh_->get_vertex(i);
// 			HE_edge* edge = vert->m_pedge;
// 			do 
// 			{
// 				int nId = edge->m_pvert->m_id;
// 				// laplacian term
// 				float fLij = getLaplacianTerm(i, nId);
// 
// 				// orientation term
// 				float fQij = getOrientationTerm(i, nId);
// 
// 				// alignment term (undo)
// 				float fAij = 0;
// 
// 				//oneRow[nId] = fLij + paras_.omega1_*fQij + paras_.omega2_*fAij;
// 				rowC(i, nId) = fLij + paras_.omega1_*fQij + paras_.omega2_*fAij;
// 				edge = edge->m_ppair->m_pnext;
// 			} while (edge!=vert->m_pedge && edge!=NULL);
			int neiSize = meshExtData_.getNeighborSize(i);
			for (int j=0; j<neiSize; j++)
			{
				HE_vert* vNeighbor = meshExtData_.getNeighbor(i, j);
				int nId = vNeighbor->m_id;
				// laplacian term
				float fLij = getLaplacianTerm(i, nId);
				 
				// orientation term
				float fQij = getOrientationTerm(i, nId);
				 
				// alignment term (undo)
				float fAij = 0;
				rowC(i, nId) = fLij + paras_.omega1_*fQij + paras_.omega2_*fAij;
			}
			float fLii = getLaplacianTerm(i, i);
			float fQii = getOrientationTerm(i, i);
			float fAii = 0;
			//oneRow[i] = fLii + paras_.omega1_*fQii + paras_.omega2_*fAii;
			//rowC.push_back(oneRow);
			rowC(i, i) = fLii + paras_.omega1_*fQii + paras_.omega2_*fAii;
		}

		solver_.clear();
		numcNew::CreateCSRMatrixFromRowMap(mat, rowC);
		solver_.init(&mat);
	}

	void ZQuadrangulation::initLinearSystem2()
	{
		// prepare matrix
		// cholmod solver
		int nSize = pMesh_->get_num_of_vertex_list();
		SpMatrix matL(nSize, nSize);
		SpMatrix matDinverse(nSize, nSize);
		SpMatrix matD(nSize, nSize);
		SpMatrix matOrient(nSize, nSize);
		SpMatrix matAlign(nSize, nSize);
		std::vector<SpMatrixTriplet> dataL, dataOri, dataAlign, dataDinverse, dataD;
		dataL.reserve(nSize*8);
		dataOri.reserve(nSize*8);
		dataDinverse.reserve(nSize);
		dataD.reserve(nSize);
		// laplacian term
		std::cout << " ..Laplacian term...";
		for (int i=0; i<nSize; i++)
		{
			MeshVrtLocalCoord* vrtData = (MeshVrtLocalCoord*)meshExtData_.getVrtAttri(i);
			HE_vert* vert = pMesh_->get_vertex(i);
			HE_edge* edge = vert->m_pedge;
			do 
			{
			 	int nId = edge->m_pvert->m_id;
				float fLij = getLaplacian(i, nId);
				dataL.push_back(SpMatrixTriplet(i, nId, fLij));
			 
			 	edge = edge->m_ppair->m_pnext;
			} while (edge!=vert->m_pedge && edge!=NULL);
			dataL.push_back(SpMatrixTriplet(i, i, getLaplacian(i, i)));
			dataDinverse.push_back(SpMatrixTriplet(i, i, 1.0/vrtData->area_));
			dataD.push_back(SpMatrixTriplet(i, i, vrtData->area_));
		}
		std::cout << "done!\n";

		// Orientation term
		std::cout << " .. Orientation term..";
		for (int i=0; i<nSize; i++)
		{
			MeshVrtLocalCoord* vrtData = (MeshVrtLocalCoord*)meshExtData_.getVrtAttri(i);
			for (int j=0; j<vrtData->neighborIds_.size(); j++)
			{
				int nId = vrtData->neighborIds_[j];
				float fOij = getOrientationTerm(i, nId);
				dataOri.push_back(SpMatrixTriplet(i, nId, fOij));
				//std::cout << (char)13 << fOij;
			}
		}
		std::cout << "done!\n";

		matL.setFromTriplets(dataL.begin(), dataL.end());
		matOrient.setFromTriplets(dataOri.begin(), dataOri.end());
		matDinverse.setFromTriplets(dataDinverse.begin(), dataDinverse.end());
		matD.setFromTriplets(dataD.begin(), dataD.end());

		SpMatrix matA = matL.adjoint()*matDinverse*matL + matL*(2*paras_.lambda_) + matD * (paras_.lambda_*paras_.lambda_);
		SpMatrix matO = matOrient.adjoint()*matOrient*paras_.omega1_;

// 		std::ofstream fout("matrix.txt");
// 		fout << matA << "\n";
// 		fout.close();
// 
// 		std::ofstream fout2("matrixO.txt");
// 		fout2 << matOrient << "\n";
// 		fout2.close();

		cholmodSolver_.compute(matA+matO);
		//cholmodSolver_.compute(matA);
		if (cholmodSolver_.info()!=Eigen::Success)
		{
			std::cerr << "Error in cholmod init\n";
		}
		std::cout << "Linear System has been built!\n";
	}

	void ZQuadrangulation::recomputeQuadraticCoordinate()
	{
		int vSize = pMesh_->get_num_of_vertex_list();
		for (int i=0; i<vSize; i++)
		{
			MeshVrtLocalCoord* vrt = (MeshVrtLocalCoord*)meshExtData_.getVrtAttri(i);
			Eigen::MatrixXf fMat(vrt->neighborIds_.size(), 1);
			for (int j=0; j<vrt->neighborIds_.size(); j++)
			{
				int id = vrt->neighborIds_[j];
				MeshVrtAttri* neiVrt = meshExtData_.getVrtAttri(id);
				fMat(j, 0) = vrtScaleField_(id)*meshExtData_.getVrtAttri(id)->area_sqrt_;
// 				if(neiVrt->type_==VERTEX_DEFAULT)
// 				{
// 					fMat(j, 0) = vrtScaleField_(id)*meshExtData_.getVrtAttri(id)->area_sqrt_;
// 				}
// 				else	// VERTEX_ADDED
// 				{
// 					MeshVrtOnEdge* neiVrt2 = (MeshVrtOnEdge*)neiVrt;
// 					float scale = neiVrt2->scale_;
// 					fMat(j, 0) = (vrtScaleField_(neiVrt2->vId0_)*scale + vrtScaleField_(neiVrt2->vId1_)*(1-scale))*meshExtData_.getVrtAttri(id)->area_sqrt_;
// 				}
			}
			vrt->operatorMat_ = vrt->neighborTransMat_*fMat;
		}
		// testing
		//saveTransformations("trans.txt");
	}

	double ZQuadrangulation::computeOnce()
	{
		int nSize = pMesh_->get_num_of_vertex_list();

		// 0. recompute the local quadratic coordinates
		//recomputeQuadraticCoordinate();


		// 1. Ay = J_k^T --> y
		//std::vector<double> y(nSize);
		//solver_.solve(&vrtScaleField_[0], &y[0], 1);
		//lsSolver_.solve(&vrtScaleField_[0], &y[0], 1);
		Eigen::VectorXf y;
		y = cholmodSolver_.solve(vrtScaleField_);


		// 2. mu = -1/(J_k*y) --> mu
		double miu = 0;
		for (int i=0; i<nSize; i++)
		{
			miu += vrtScaleField_(i)*y(i);
		}
		//miu *= -1;

		// 3. A*(f_k+1) = -mu*J_k^T
// 		std::vector<double> f(nSize);
// 		std::vector<double> b(nSize);
// 		for (int i=0; i<nSize; i++)
// 		{
// 			b[i] = miu*vrtScaleField_[i];
// 		}
		Eigen::VectorXf b = vrtScaleField_*miu;
		Eigen::VectorXf f;
		//solver_.solve(&b[0], &f[0], 1);
		f = cholmodSolver_.solve(b);

		// error
		double err = 0;
		for (int i=0; i<nSize; i++)
		{
			err += meshExtData_.getVrtAttri(i)->area_ * (vrtScaleField_(i) - f(i));
		}
		err /= nSize;

		// 4. f_(k+1)
//		std::cout << err << " --done!\n";
		vrtScaleField_ = f;
// 		std::ofstream fout("field.txt");
// 		fout << vrtScaleField_;
// 		fout.close();
		return err;
	}

	float ZQuadrangulation::getLaplacian(int vId0, int vId1)
	{
		if (vId0==vId1)
		{
			HE_vert* vert = pMesh_->get_vertex(vId0);
			float ret = 0;
			HE_edge* edge = vert->m_pedge;
			do 
			{
				//ret += (meshExtData_.getEdgeAttri(edge->m_id)->cotAlpha_ + meshExtData_.getEdgeAttri(edge->m_ppair->m_id)->cotAlpha_)*0.5;
				ret += meshExtData_.getEdgeAttri(edge->m_id)->laplacianW_;
				edge = edge->m_ppair->m_pnext;
			} while (edge!=vert->m_pedge && edge!=NULL);
			return ret*(-1.f);
		}
		else
		{
			HE_edge* e = MeshTools::getEdge(pMesh_, vId0, vId1);
			if (e==NULL) return 0;
			//float ret = (meshExtData_.getEdgeAttri(e->m_id)->cotAlpha_ + meshExtData_.getEdgeAttri(e->m_ppair->m_id)->cotAlpha_)*0.5;
			return meshExtData_.getEdgeAttri(e->m_id)->laplacianW_;
		}
	}

	float ZQuadrangulation::getLaplacian(HE_edge* e)
	{
		if (e==NULL) return 0;
		//float ret = (meshExtData_.getEdgeAttri(e->m_id)->cotAlpha_ + meshExtData_.getEdgeAttri(e->m_ppair->m_id)->cotAlpha_)*0.5;
		return meshExtData_.getEdgeAttri(e->m_id)->laplacianW_;
	}

	float ZQuadrangulation::getLaplacianTerm(int vId0, int vId1)
	{
		if (vId0==vId1)
		{
			HE_vert* vert = pMesh_->get_vertex(vId0);
			HE_edge* edge = vert->m_pedge;
			MeshVrtLocalCoord* vrtData = (MeshVrtLocalCoord*)meshExtData_.getVrtAttri(vId0);
			float ret = 0;
			do 
			{
				int vId = edge->m_pvert->m_id;
				float lap = getLaplacian(vId0, vId);
				ret += lap*lap/meshExtData_.getVrtAttri(vId)->area_;
				edge = edge->m_ppair->m_pnext;
			} while (edge!=vert->m_pedge && edge!=NULL);
			ret += getLaplacian(vId0, vId1)*paras_.lambda_*2;
			ret += vrtData->area_*paras_.lambda_*paras_.lambda_;
			return ret;
		}
		else
		{
			HE_vert* v0 = pMesh_->get_vertex(vId0);
			HE_vert* v1 = pMesh_->get_vertex(vId1);
			MeshVrtLocalCoord* vrtData0 = (MeshVrtLocalCoord*)meshExtData_.getVrtAttri(vId0);
			float ret = 0;

			HE_edge* edge = v0->m_pedge;
			do 
			{
				HE_vert* v = edge->m_pvert;
				HE_edge* e2 = MeshTools::getEdge(v, v1);
				if (e2!=NULL)	// v0 and v1 connect at v
				{
					ret += getLaplacian(edge)*getLaplacian(e2)/meshExtData_.getVrtAttri(v->m_id)->area_;
				}
				edge = edge->m_ppair->m_pnext;
			} while (edge!=v0->m_pedge&& edge!=NULL);

			ret += getLaplacian(vId0, vId1)*paras_.lambda_*2;
			return ret;
		}
		
	}

	float ZQuadrangulation::getOrientationTerm(int vId0, int vId1)
	{
// 		float ret = meshExtData_.getVrtAttri(vId0)->area_sqrt_ * meshExtData_.getVrtAttri(vId1)->area_sqrt_;
// 		ret *= ((MeshVrtLocalCoord*)meshExtData_.getVrtAttri(vId0))->operatorMat_[4];
// 		ret	*= ((MeshVrtLocalCoord*)meshExtData_.getVrtAttri(vId1))->operatorMat_[4];
// 		return ret;

		float ret = 0;
		MeshVrtLocalCoord* vrt = (MeshVrtLocalCoord*)meshExtData_.getVrtAttri(vId0);
		int neiIdx = -1;
		for (int i=0; i<vrt->neighborIds_.size(); i++)
		{
			if (vrt->neighborIds_[i]==vId1)
			{
				neiIdx = i;
				break;
			}
		}
		if (neiIdx!=-1)
		{
			ret += vrt->neighborTransMat_(4, neiIdx) * vrt->area_sqrt_;
		}
		else
		{
			std::cout << "Cannot find term (" << vId0 << "," << vId1 << ")\n";
		}
		return ret;
	}

	float ZQuadrangulation::getLaplacianTerm2(int vId0, int vId1)
	{
		float Lij = getLaplacian(vId0, vId1);
		float ret = 0;
		float weight = meshExtData_.getVrtAttri(vId0)->area_sqrt_;
		ret += Lij/weight;
		if (vId0==vId1)
			ret += weight;
		return ret;
	}

	float ZQuadrangulation::getOrientationTerm2(int vId0, int vId1)
	{
		return 0;
	}

	void ZQuadrangulation::saveTransformations(const char* fileName)
	{
		std::ofstream fout(fileName);
		for (int i=0; i<meshExtData_.vrtSize_; i++)
		{
			fout << "/***\n";
			MeshVrtLocalCoord* vrt = (MeshVrtLocalCoord*)meshExtData_.getVrtAttri(i);
			//fout << vrt->neighborTransMat_ << "\n";
			fout << vrt->operatorMat_ << "\n";
			fout << "**/\n";
		}
		fout.close();
	}

	void ZQuadrangulation::testingEigen()
	{
// 		Eigen::MatrixXf mat(6, 6);
// 		mat.fill(0);
// 		for (int i=0; i<6; i++)
// 		{
// 			for (int j=0; j<6; j++)
// 			{
// 				mat(i,j) = rand()%1000;
// 			}
// 		}
// 		std::cout << mat << "\n";
// 		Eigen::MatrixXf mat2 = mat.inverse();
// 		std::cout << "Inverse:\n" << mat2 << "\n";
// 		Eigen::MatrixXf mat3 = mat*mat2;
// 		std::cout << "M*M2=\n" << mat3 << "\n";

		typedef Eigen::Triplet<double> T;
		std::vector<T> tripletList;
		tripletList.reserve(11);
		for (int i=0; i<9; i++)
		{
			tripletList.push_back(T(i, i, i+1));
		}
		tripletList.push_back(T(0,1,3));
		tripletList.push_back(T(1,0,3));
		//Eigen::SparseMatrixType mat4(9, 9);
		SpMatrix mat4(9,9);
		//mat4.
		mat4.setFromTriplets(tripletList.begin(), tripletList.end());
		std::cout << "Sparse Matrix:\n " << mat4 << "\n";
		SpMatrix mat5 = mat4*mat4.adjoint();
		std::cout << " m*m^T=\n" << mat5 << "\n";

		Eigen::CholmodSimplicialLDLT<SpMatrix> cholmod_solver;
		cholmod_solver.compute(mat5);
		if (cholmod_solver.info()!=Eigen::Success)
		{
			std::cout << "Wrong!\n";
		}

		Eigen::VectorXf x,b(9);
		for (int i=0; i<9; i++)
		{
			b(i) = i+1;
		}
		x = cholmod_solver.solve(b);
		std::cout << "After solve: x=\n" << x << "\n";
	}

	void ZQuadrangulation::smoothScaleField()
	{
		int vSize = pMesh_->get_num_of_vertex_list();
		Eigen::VectorXf newVF(vSize);
		for (int i=0; i<vSize; i++)
		{
			HE_vert* vert = pMesh_->get_vertex(i);
			HE_edge* edge = vert->m_pedge;
			double dSum = 0;
			do 
			{
				dSum += vrtScaleField_[edge->m_pvert->m_id];
				edge = edge->m_ppair->m_pnext;
			} while (edge!=vert->m_pedge && edge!=NULL);
			newVF(i) = vrtScaleField_(i) + (dSum/vert->m_degree - vrtScaleField_(i))*0.1;
		}
		vrtScaleField_ = newVF;
	}

	void ZQuadrangulation::tagVertices()
	{
		if (pMesh_==NULL) return;

		int vSize = vrtScaleField_.size();
		for (int i=0; i<vSize; i++)
		{
			HE_vert* vert = pMesh_->get_vertex(i);
			MeshVrtLocalCoord* vertData = (MeshVrtLocalCoord*)meshExtData_.getVrtAttri(i);
			float fCenter = vrtScaleField_[i];
			
			// check whether maximum or minimum
			int tag = 0;	// 0->default, 1->bigger, 2->smaller, 3->stop
			std::vector<HE_vert*> biggerSet, smallerSet;
			HE_edge* edge = vert->m_pedge;
			do 
			{
				HE_vert* v = edge->m_pvert;
				bool bBig = vrtScaleField_[v->m_id] > fCenter;
				if (bBig) biggerSet.push_back(v);
				else smallerSet.push_back(v);
				if (tag==0)
				{
					tag = bBig ? 1 : 2;
				}
				else if (tag==1)
				{
					tag = bBig ? 1 : 3;
				}
				else if (tag==2)
				{
					tag = bBig ? 3 : 2;
				}
				else if (tag==3)
				{
					;
				}
				edge = edge->m_ppair->m_pnext;
			} while (edge!=NULL && edge!=vert->m_pedge);
			if (tag==1)
			{
				vertData->type_ = VERTEX_MAXIMUM;
				continue;
			}
			else if (tag==2)
			{
				vertData->type_ = VERTEX_MINIMUM;
				continue;
			}

			// check regular or saddle
			int notConnected = 0;
			int smallSize = smallerSet.size();
			for (int j=0; j<smallSize; j++)
			{
				if (MeshTools::getEdge(smallerSet[j], smallerSet[(j+1)%smallSize])==NULL)
				{
					notConnected ++;
				}
			}
			if (notConnected<=1)
				vertData->type_ = VERTEX_REGULAR;
			else 
				vertData->type_ = VERTEX_SADDLE;
		}

		// save to mesh
		for (int i=0; i<vSize; i++)
		{
			pMesh_->get_vertex(i)->m_flag = meshExtData_.getVrtAttri(i)->type_;
		}
	}

	void ZQuadrangulation::generateMSEdges()
	{
		meshExtData_.clearMSEdges();

		std::vector<MeshMSEdge>& msEdges = meshExtData_.getMSEdges();
		int vSize = pMesh_->get_num_of_vertex_list();
		for (int i=0; i<vSize; i++)
		{
			MeshVrtLocalCoord* vertData = (MeshVrtLocalCoord*)meshExtData_.getVrtAttri(i);
			if (vertData->type_!=VERTEX_SADDLE) continue;

			//std::cout << i << "/vSize:" << vSize << "\n";
			// start from saddle points
			findMSNeighbors(i, true);
		}

		// add the edges to the vertices
		for (int i=0; i<msEdges.size(); i++)
		{
			MeshMSEdge& msEdge = msEdges[i];
			msEdge.vert1_->msEdges_.push_back(&msEdge);
			msEdge.vert2_->msEdges_.push_back(&msEdge);
		}
	}

	void ZQuadrangulation::findMSNeighbors(int vId, bool bAscending)
	{
		MeshVrtLocalCoord* vertData = (MeshVrtLocalCoord*)meshExtData_.getVrtAttri(vId);
		HE_vert* vert = pMesh_->get_vertex(vId);
		float vertValue = vrtScaleField_[vId];

		if (vertData->type_==VERTEX_SADDLE)
		{
			std::vector<HE_vert*> lowerNeighbors, higherNeighbors;
			HE_edge* edge = vert->m_pedge;
			do 
			{
				HE_vert* v = edge->m_pvert;
				if (vrtScaleField_[v->m_id]>vertValue)
					higherNeighbors.push_back(v);
				else
					lowerNeighbors.push_back(v);
				edge = edge->m_ppair->m_pnext;
			} while (edge!=NULL && edge!=vert->m_pedge);	

			std::vector<std::vector<HE_vert*>> highChains;
			findConnectedChains(higherNeighbors, highChains);
			for (int i=0; i<highChains.size(); i++)
			{
				std::vector<HE_vert*> oneChain = highChains[i];
				HE_vert *smallVert, *bigVert;
				getBigSmallVerts(oneChain, bigVert, smallVert);
				// add to the vertex
				if (bigVert!=NULL)
				{
					MeshVrtAttri* bigVertData = meshExtData_.getVrtAttri(bigVert->m_id);
					meshExtData_.getMSEdges().push_back(MeshMSEdge(vertData, bigVertData));
					//vertData->msEdges_.push_back(meshExtData_.getLastMSEdge());
					if (bigVertData->type_==VERTEX_REGULAR)
						findMSNeighbors(bigVert->m_id, true);
				}
				if (smallVert!=NULL)
				{
					MeshVrtAttri* smallVertData = meshExtData_.getVrtAttri(smallVert->m_id);
					meshExtData_.getMSEdges().push_back(MeshMSEdge(vertData, smallVertData));
					//vertData->msEdges_.push_back(meshExtData_.getLastMSEdge());
					if (smallVertData->type_==VERTEX_REGULAR)
						findMSNeighbors(smallVert->m_id, false);
				}
			}
			std::vector<std::vector<HE_vert*>> lowChains;
			findConnectedChains(lowerNeighbors, lowChains);
			for (int i=0; i<lowChains.size(); i++)
			{
				std::vector<HE_vert*> oneChain = lowChains[i];
				HE_vert *smallVert, *bigVert;
				getBigSmallVerts(oneChain, bigVert, smallVert);
				// add to the vertex
				if (bigVert!=NULL)
				{
					MeshVrtAttri* bigVertData = meshExtData_.getVrtAttri(bigVert->m_id);
					meshExtData_.getMSEdges().push_back(MeshMSEdge(vertData, bigVertData));
					//vertData->msEdges_.push_back(meshExtData_.getLastMSEdge());
					if (bigVertData->type_==VERTEX_REGULAR)
						findMSNeighbors(bigVert->m_id, true);
				}
				if (smallVert!=NULL)
				{
					MeshVrtAttri* smallVertData = meshExtData_.getVrtAttri(smallVert->m_id);
					meshExtData_.getMSEdges().push_back(MeshMSEdge(vertData, smallVertData));
					//vertData->msEdges_.push_back(meshExtData_.getLastMSEdge());
					if (smallVertData->type_==VERTEX_REGULAR)
						findMSNeighbors(smallVert->m_id, false);
				}
			}
		}
		else
		{
			std::vector<HE_vert*> neighbors;
			HE_edge* edge = vert->m_pedge;
			do 
			{
				HE_vert* v = edge->m_pvert;
				neighbors.push_back(v);
				edge = edge->m_ppair->m_pnext;
			} while (edge!=NULL && edge!=vert->m_pedge);
			HE_vert* bigVert, *smallVert;
			getBigSmallVerts(neighbors, bigVert, smallVert);
			HE_vert* nextVert = bAscending ? bigVert : smallVert;
			if (nextVert!=NULL)
			{
				MeshVrtAttri* nextVertData = meshExtData_.getVrtAttri(nextVert->m_id);
				meshExtData_.getMSEdges().push_back(MeshMSEdge(vertData, nextVertData));
				if (nextVertData->type_==VERTEX_REGULAR)
					findMSNeighbors(nextVert->m_id, bAscending);
			}
		}

		// 
	}

	void ZQuadrangulation::findConnectedChains(std::vector<HE_vert*>& inputVertices, std::vector<std::vector<HE_vert*>>& chains)
	{
		// find connected chains
		//std::vector<std::vector<HE_vert*>> highChains;
		int count = 0;
		int i=0;
		int nbSize = inputVertices.size();
		if (nbSize==0) return;
		// start from an unconnected vrt
		while (true)
		{
			if (MeshTools::getEdge(inputVertices[i], inputVertices[(i+1)%nbSize])==NULL)
				break;
			i++; if(i>=nbSize) break;
		}
		chains.push_back(std::vector<HE_vert*>());
		i++; if(i>=nbSize) i-=nbSize;
		while (true)
		{
			std::vector<HE_vert*>& oneChain = chains[chains.size()-1];
			oneChain.push_back(inputVertices[i]);
			count ++;
			i++; if(i>=nbSize) i-=nbSize;
			if (count==nbSize) break;

			if (MeshTools::getEdge(inputVertices[i], inputVertices[(i+1)%nbSize])==NULL)
				chains.push_back(std::vector<HE_vert*>());
		}
	}

	void ZQuadrangulation::getBigSmallVerts(std::vector<HE_vert*>& chain, HE_vert* &bigVert, HE_vert* &smallVert)
	{
		bigVert = smallVert = NULL;
		if (chain.empty()) return;

		float fBig, fSmall;
		fBig = fSmall = vrtScaleField_[chain[0]->m_id];
		bigVert = smallVert = chain[0];
		for (int i=1; i<chain.size(); i++)
		{
			HE_vert* v = chain[i];
			float f = vrtScaleField_[v->m_id];
			if (f>fBig)
			{
				fBig = f;
				bigVert = v;
			}
			if (f<fSmall)
			{
				fSmall = f;
				smallVert = v;
			}
		}
	}
}
#pragma once
#ifndef ZQUADRANGULATION_H_
#define ZQUADRANGULATION_H_

#include "GlobalDefs.h"
#include "../Scene/Mesh3D.h"
#include "../Eigen/Eigen/Dense"
#include "../Eigen/Eigen/Sparse"
#include "../Eigen/Eigen/CholmodSupport"
#include "../Eigen/mkl_addonNew.h"
#include "ZExpMapGenerateWarp.h"

namespace ZMeshSpace
{

	typedef Eigen::SparseMatrix<float> SpMatrix;
	typedef Eigen::CholmodSimplicialLDLT<SpMatrix> SpMatrixSolver;
	typedef Eigen::Triplet<float> SpMatrixTriplet;

	class ZQuadrangulation
	{
		// mesh extra data
		typedef Eigen::Matrix<float, 6, 1> Vec6f;

		// forward class declaration
		class MeshMSEdge;

	public:
		enum MeshVrtType {
			VERTEX_ADDED = -1,
			VERTEX_DEFAULT = 0,
			VERTEX_MAXIMUM = 1,
			VERTEX_MINIMUM = 2,
			VERTEX_SADDLE = 3,
			VERTEX_REGULAR = 4,
		};
		class MeshVrtAttri
		{
		public:
			MeshVrtType type_;
			int id_;
			float area_;
			float area_sqrt_;
			std::vector<MeshMSEdge*> msEdges_;
		};
		class MeshVrtLocalCoord : public MeshVrtAttri
		{
		public:
			Vec3f directionU_;
			Vec3f directionV_;
			std::vector<int> neighborIds_;
			Eigen::MatrixXf neighborTransMat_;
			std::vector<Vec2f> neighborCoords_; // (u, v)
			std::vector<float> neighborDistToCenters_;	// the distance between the neighborhood and the center vertex
			Vec6f operatorMat_;	//(C, Cu, Lv, Cuu, Cuv, Cvv)

			void addNeighborVrt(int id, const Vec3f& newNeighbor)
			{
				neighborIds_.push_back(id);
				neighborCoords_.push_back(Vec2f(newNeighbor.DotProduct(directionU_), newNeighbor.DotProduct(directionV_)));
				neighborDistToCenters_.push_back(newNeighbor.length());
			}

			void addNeighborVrt(int id, const Vec2f& uv)
			{
				neighborIds_.push_back(id);
				neighborCoords_.push_back(uv);
				//std::cout << uv << "\n";
			}

			void computeOperatorMat(const std::vector<float>& scaleField)
			{
				Eigen::MatrixXf constrainedField(neighborIds_.size(), 1);
				for (int i=0; i<neighborIds_.size(); i++)
				{
					constrainedField(i, 0) = scaleField[neighborIds_[i]];
				}
				operatorMat_ = neighborTransMat_*constrainedField;
			}
		};
		class MeshVrtOnEdge : public MeshVrtAttri
		{
		public:
			int vId0_;
			int vId1_;
			float scale_;
		};

		class MeshEdgeExtData
		{
		public:
			int id_;
			float cotAlpha_;
			float laplacianW_;
		};

		class MeshFaceExtData
		{
		public:
			int id_;
			float scaleFactor_;
		};

		class MeshVrt2RingNeighbor
		{
		public:
			std::vector<HE_vert*> neighbors_;

			void addNeighborVrt(HE_vert* v) {
				if (std::find(neighbors_.begin(), neighbors_.end(), v)==neighbors_.end()) {
					neighbors_.push_back(v);
				}
			}
		};

		class MeshMSEdge
		{
		public:
			MeshVrtAttri	*vert1_, *vert2_;
		public:
			MeshMSEdge(MeshVrtAttri* vert1, MeshVrtAttri* vert2) : vert1_(vert1), vert2_(vert2){}

		};

		class MeshExtData
		{
		protected:
			std::vector<MeshVrtLocalCoord> vrtData_;
			std::vector<MeshEdgeExtData> edgeData_;
			std::vector<MeshFaceExtData> faceData_;
			std::vector<MeshVrtOnEdge> vrtExtData_;
			std::vector<MeshVrt2RingNeighbor> vrt2RingNeibors_;
			std::vector<MeshMSEdge> msEdges_;
		public:
			int vrtSize_;

			MeshVrtAttri* getVrtAttri(int id) {
				if (id>=0 && id<vrtSize_) return &vrtData_[id];
				else if (id>=vrtSize_ && id<vrtSize_+vrtExtData_.size()) return &vrtExtData_[id-vrtSize_];
				else return NULL;
			}
			MeshEdgeExtData* getEdgeAttri(int id) {return &edgeData_[id];}
			MeshFaceExtData* getFaceAttri(int id) {return &faceData_[id];}
			MeshFaceExtData* getFaceAttri(HE_face* face) {return face==NULL ? NULL : &faceData_[face->m_id];}
			int getFaceAttriSize() {return faceData_.size();}
			void resizeFaceAttris(int size) {faceData_.resize(size);}
			std::vector<MeshVrt2RingNeighbor>& getVrt2RingNeighbors() {return vrt2RingNeibors_;}
			MeshVrt2RingNeighbor* getVrt2RingNeighbor(int vId) {return &vrt2RingNeibors_[vId];}
			int getNeighborSize(int vId) {return vrt2RingNeibors_[vId].neighbors_.size();}
			HE_vert* getNeighbor(int vId, int idx) {return vrt2RingNeibors_[vId].neighbors_[idx];}
			MeshMSEdge* getMSEdge(int eId) {return &msEdges_[eId];}
			MeshMSEdge* getLastMSEdge() {return msEdges_.empty() ? NULL : &msEdges_[msEdges_.size()-1];}
			std::vector<MeshMSEdge>& getMSEdges() {return msEdges_;}
			const std::vector<MeshMSEdge>& getMSEdges() const {return msEdges_;}

			void clearMSEdges() {msEdges_.clear();}

			void resetNeighbors(int vSize) {
				vrt2RingNeibors_.resize(vSize);
			}

			void resetData(int vSize, int eSize, int fSize) {
				vrtData_.resize(vSize);
				edgeData_.resize(eSize);
				faceData_.resize(fSize);
				vrtSize_ = vSize;
			}

			void addVrtData(MeshVrtType type, float area, int vId0, int vId1, float scale) {
				int id = vrtExtData_.size();
				int realId = vrtData_.size() + id;
				if (type==VERTEX_ADDED)
				{
					vrtExtData_.push_back(MeshVrtOnEdge());
					MeshVrtOnEdge* v = (MeshVrtOnEdge*)getVrtAttri(realId);
					v->id_ = realId;
					v->type_ = type;
					v->vId0_ = vId0;
					v->vId1_ = vId1;
					v->scale_ = scale;
					v->area_ = area;
					v->area_sqrt_ = sqrt(area);
				}
				else if (type==VERTEX_DEFAULT)
				{
					// undo
				}
			}
			void clear();
		};	

		// parameters
		struct Parameters
		{
			float omega1_;
			float omega2_;
			float lambda_;
			float threshold_stopIteration_;
		};
	
	public:
		ZQuadrangulation(Mesh3D* mesh=NULL);
		virtual ~ZQuadrangulation();

		void init(Mesh3D* mesh, const std::vector<Vec3f>& localU, const std::vector<Vec3f>& localV, const std::vector<float>& faceFactors=std::vector<float>());
		void destroy();

		void setLocalUV(const std::vector<Vec3f>& localU, const std::vector<Vec3f>& localV);
		void setLambda(float lambda) {paras_.lambda_ = lambda;}
		void setFaceFactors(const std::vector<float>& faceFactors);

		void go();
		void generateInitField();
		double computeOnce();
		void smoothScaleField();

		void tagVertices();
		void generateMSEdges();

		// getter and setters
		void setVrtScaleFiled(const std::vector<double>& scaleFiled) {
			vrtScaleField_.resize(scaleFiled.size());
			for (int i=0; i<scaleFiled.size(); i++) {
				vrtScaleField_(i) = scaleFiled[i];
			}
		}

		Mesh3D* getMesh() {return pMesh_;}
		//const MeshExtData& getMeshExtData() const {return meshExtData_;}
		int getMSEdgeSize() const {return meshExtData_.getMSEdges().size();}
		MeshMSEdge* getMSEdge(int idx) {return meshExtData_.getMSEdge(idx);}

		Eigen::VectorXf getVrtScaleField() {return vrtScaleField_;}

		// for debug
		void testingEigen();

	private:
		void buildVertex2RingNeighbor();
		void resetVrtAreas();
		void initExtData(const std::vector<Vec3f>& localU, const std::vector<Vec3f>& localV, const std::vector<float>& faceFactors);
		void resetLocalCoordinates(const std::vector<Vec3f>& localU, const std::vector<Vec3f>& localV);
		void resetLaplacianCoord();
		void resetLocalTransformation();
		void initLinearSystem();
		void initLinearSystem2();
		void recomputeQuadraticCoordinate();
		void findMSNeighbors(int vId, bool bAscending);

	private:
		float getLaplacian(int vId, int vId1);
		float getLaplacian(HE_edge* edge);
		float getLaplacianTerm(int vId0, int vId1);
		float getOrientationTerm(int vId0, int vId1);

		float getLaplacianTerm2(int vId0, int vId1);
		float getOrientationTerm2(int vId0, int vId1);
		void findConnectedChains(std::vector<HE_vert*>& inputVertices, std::vector<std::vector<HE_vert*>>& chains);
		void getBigSmallVerts(std::vector<HE_vert*>& chain, HE_vert* &bigVert, HE_vert* &smallVert);

	private:
		Mesh3D* pMesh_;	// only a pointer
		//std::vector<double> vrtScaleField_;
		Eigen::VectorXf vrtScaleField_;
		MeshExtData meshExtData_;
		Parameters paras_;

		ZExpMapGenerateWarp expMapGen_;

		numcNew::SparseSolver solver_;
		numcNew::LeastSquareSparseSolver lsSolver_;

		SpMatrixSolver cholmodSolver_;

	private:
		void saveTransformations(const char* fileName);
	};
}


#endif//ZQUADRANGULATION_H_
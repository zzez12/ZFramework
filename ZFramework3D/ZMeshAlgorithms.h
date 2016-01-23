#pragma once

#include "GlobalDefs.h"
#include "../Scene/Mesh3D.h"
#include "ZPlane.h"
#include "ZMeshSurfacePoint.h"
#include "ZIsoLine.h"
#include "ZMeshProjector.h"
#include "ZEigenData.h"
#include "ZIsoLineContainer.h"
#include "ZQuadrangulation.h"
#include <fstream>

#include "ZMeshFilterManifold.h"

namespace ZMeshSpace
{
	class ZMeshAlgorithms
	{
	private:
		Mesh3D *mesh_;

	public:
		class MeshLaplacianData
		{
		public:
			std::vector<std::pair<int, int>> edges;
			std::vector<float> weights;
		};

		class IsoBand
		{
		public:
			std::vector<int> vrtIdxs;
			double bandMin, bandMax;
		};	

		enum LaplacianType {
			UNIFORM = 0,
			COTANGENT,
		};

		int maxEigens;

		ZMeshAlgorithms(Mesh3D* mesh=NULL);

		void reset();
		void setMesh(Mesh3D* mesh);

		void setLaplacianType(LaplacianType type) {laplacianType_ = type;}
		void computeLaplacian();
		bool saveLaplacian(const char* fileName);
		bool loadEigenVector(const char* fileName);
		bool saveEigenVectors(const char* fileName);
		bool loadEigenVectors(const char* fileName);
		bool computeEigens();

		void setMeshVrtColors();
		void setMeshVrtColors(int eigIdx);

		int getCurrentEigenIdx(){return currentEigenIdx_;}
		void setCurrentEigenIdx(int idx) {currentEigenIdx_ = idx%maxEigens;}

		//void computeIsoLine(HE_vert* vert, double value, ZEigenIsoLine& outLine);
		//void computeIsoLine(ZMeshSurfacePoint& startP, double value, ZEigenIsoLine& outLine);
		void computeIsoLines(int idx, double value);
		void computeProjectionPlanes();

		const ZEigenIsoLine& getIsoLine(int idx);
		int getIsoLineSize();

		void setCurrentIsoLineIdx(int idx) {currentIsoLineIdx_=idx;}
		int getCurrentIsoLineIdx() {return currentIsoLineIdx_;}

		ZEigenData& getEigenData() {return eigenDatas_[getCurrentEigenIdx()];}
		const ZIsoLineContainer& getIsoLineContainer() const {return isoLineContainer_;}
		void setKmeansK(int k) {isoLineContainer_.setKmeansK(k);}

		//ZMeshSpace::ZQuadrangulation& getQuadrangulation() {return meshQuadrangulation;}
		ZMeshSpace::ZQuadrangulation* getQuadrangulation() {return &meshQuadrangulation;}
		//const ZMeshSpace::ZQuadrangulation& getQuadrangulation() const {return meshQuadrangulation;}

		void clearIsoLineData();// {isoLines_.clear();}

		void intersectMeshByPlane(ZPlane* plane);

		bool initQuadrangulation();
		bool runQuadrangulation(float faceScale=20, float lambda=10);

		// manifold smoothing commands
		bool runBilateralFiltering(float sigma_r, float sigma_s);
		bool runManifoldSmooth(float sigma_r, float sigma_s);
		void resetMeshPositions();
		Eigen::MatrixXf computeNewPositions(const Eigen::MatrixXf& newNormals);
		void setRangeLevelWeight(int level, float weight);
		void updateRange();
		void updateVerticePos(const Eigen::MatrixXf& newNormals, const Eigen::MatrixXf& faceCenters);
		void addNoiseToMesh(float sigma);
		void updateVerticesPos();

	public:
		void test0();
		void test1();

	private:
		bool saveEigens(const char* fileName, int eSize, int nCol, double* eigenValues, double* eigenVecs);

	private:
		std::vector<ZEigenData> eigenDatas_;
		LaplacianType laplacianType_;
		int currentEigenIdx_;

		//std::vector<ZEigenIsoLine> isoLines_;
		int currentIsoLineIdx_;	// -1: all; others: one index
		ZIsoLineContainer isoLineContainer_;

		ZMeshProjector meshProjector_;
		ZMeshSpace::ZQuadrangulation meshQuadrangulation;
		Eigen::MatrixXf meshInitPositions_;
		Eigen::MatrixXf meshNewPositions_;
		Eigen::MatrixXf meshNewNormals_;

		ZMeshFilterManifold *pMeshFilterManifold_;

	private:
		float laplacianWeight(HE_edge* he);
		float laplacianWeight(HE_vert* v);

		//HE_edge* findEdge(HE_vert* v0, HE_vert* v1);
		//void findIsoLineInTriangle(HE_edge* e, const std::vector<double>& field, double value, ZMeshSurfacePoint& startP, ZMeshSurfacePoint& findP);
		//Vec3f getSurfacePointPosOnEdge(HE_edge* e, const std::vector<double>& field, double value);

		// return cotA
		static float computeCotangent(float a, float b, float c)
		{
			float ret = 0;
			float cosA = (b*b+c*c-a*a)/2/b/c;
			float sinA = sqrt((double)(1-cosA*cosA));
			ret = cosA / sinA;
			return ret;
		}

		static float computeCotangent(const Vec3f& A, const Vec3f& B, const Vec3f& c)
		{
			Vec3f BA = B-A;
			Vec3f CA = c-A;
			float l = (BA*CA).length();
			float ret = BA.DotProduct(CA)/l;
			return ret;
		}
	};
}
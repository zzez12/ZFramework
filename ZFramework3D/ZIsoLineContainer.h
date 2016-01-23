#pragma once

#include "GlobalDefs.h"
#include "ZIsoLine.h"
#include "ZMeshSurfacePoint.h"
#include "ZEigenData.h"
#include "ANNWarper_N.h"
#include <vector>


namespace ZMeshSpace
{
	typedef ZEigenIsoLine IsoLineData;
	typedef std::vector<IsoLineData> IsoLineDataArray;

	class ZIsoLineContainer
	{
	public:
		struct ZIsoLineParas {
			double gap;	// the minimal distance between the iso-lines
			double rm2Gap;
			double normDGap;
			int k_kmeans;
		};

	public:
		ZIsoLineContainer();
		~ZIsoLineContainer();

		void setMesh(Mesh3D* pMesh);
		void clear();

		void extract(const ZEigenData& eigenData);
		void extract2(const std::vector<ZEigenData>& eigenDatas, int numUsedEigen=10);
		void findNextIsoLine(const std::vector<double>& field, float value);

		void computeIsoLine(const std::vector<double>& field, HE_vert* vert, double value, IsoLineData& outLine);
		void computeIsoLine(const std::vector<double>& field, ZMeshSurfacePoint& startP, double value, IsoLineData& outLine);
		void computeIsoLines(const std::vector<double>& field, double value);
		void computeIsoLines(const std::vector<double>& field, double value, IsoLineDataArray& outLines);

		void computeCuttingPlane(ZPlane *plane);

		void computeProjectionPlanes();

		void clusterPlanes();
		void clusterPlanes2();

		// getter and setters
		const IsoLineDataArray& getIsoLines() const{return isoLines_;}
		IsoLineDataArray& getIsoLines() {return isoLines_;}
		const IsoLineData& getIsoLine(int idx) const {return isoLines_[idx];}
		IsoLineData& getIsoLine(int idx) {return isoLines_[idx];}
		int getIsoLineSize() const {return isoLines_.size();}
		int getKmeansK() const {return para_.k_kmeans;}
		const std::vector<ZPlane>& getKmeansPlanes() const {return kmeansPlanes_;}
		const std::vector<int>& getPlaneList() const {return planeList_;}
		const std::vector<std::vector<int>>& getPlaneLists() const {return planeLists_;}

		void setKmeansK(int k) {para_.k_kmeans=k;} 

		void clearData() {isoLines_.clear();}

	private:
		Mesh3D* mesh_;
		ZIsoLineParas para_;
		IsoLineDataArray isoLines_;
		std::vector<ZPlane> kmeansPlanes_;

		std::vector<int> planeList_;
		std::vector<std::vector<int>> planeLists_;

	private:
		void initParas();

		HE_edge* findEdge(HE_vert* v0, HE_vert* v1);
		void findIsoLineInTriangle(HE_edge* e, const std::vector<double>& field, double value, ZMeshSurfacePoint& startP, ZMeshSurfacePoint& findP);
		Vec3f getSurfacePointPosOnEdge(HE_edge* e, const std::vector<double>& field, double value);

		bool checkAddtionConstraint(const IsoLineData& oneData, const IsoLineDataArray& isos);
		void getDistance(const MATH::vector5f& p, const MATH::vector5f& q, float& pos, float& angle);

		int findStartIsoPlaneId(const std::vector<float>& density, const std::vector<bool>& traversedTag, int clusterId);
		void findIsoPlaneList(int startId, const std::vector<MATH::vector5f>& dataPts, AnnWarper5& annData, 
			float threshold_angle, float threshold_pos, std::vector<bool>& traversedTag, std::vector<int>& planeList);

		void buildAnnData(const std::vector<MATH::vector5f>& dataPts, const std::vector<bool> bTags, AnnWarper5& annData);
	};
}
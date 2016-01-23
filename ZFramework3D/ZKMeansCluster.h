#pragma once
#ifndef ZKMEANSCLUSTER_H_
#define ZKMEANSCLUSTER_H_

#include "GlobalDefs.h"
#include <vector>
#include "vectorN.h"

namespace ZMathAlgorithms
{
	typedef float DataType;
	typedef MATH::vector5f KMeansData;
	typedef KMeansData::DistanceFunc DistanceF;
	//typedef Vec3f KMeansData;
	//typedef std::pair<Vec3f, Vec3f> KMeansData;

	static DataType dataDistance(const KMeansData& p, const KMeansData& q)
	{
		Vec3f np(sin(p[0])*cos(p[1]), sin(p[0])*sin(p[1]), cos(p[0]));
		Vec3f nq(sin(q[0])*cos(q[1]), sin(q[0])*sin(q[1]), cos(q[0]));
		Vec3f vp(p[2], p[3], p[4]);
		Vec3f vq(q[2], q[3], q[4]);
		return acos(np.DotProduct(nq)) + (vp-vq).length();
	}

	class ZKMeansCluster
	{
	public:
		//typedef float (*DistanceFunc)(const KMeansData&, const KMeansData&);
		
		//static float dataDistance(const KMeansData& p, const KMeansData& q);
		//static float dataNewDistance(const KMeansData& p, const KMeansData& q);
		//static float dataNewDistanceWithoutD(const KMeansData& p, const KMeansData& q);


	public:
		ZKMeansCluster();
		~ZKMeansCluster();

		void addDataPoint(const KMeansData& p);
		void setDistanceFunction(DistanceF func);
		std::vector<KMeansData> getMeans(int k);
		std::vector<int> getClusterIds();
		const KMeansData& getDataPoint(int idx) {return dataPts_[idx];}
		DistanceF getDistanceFunction() {return distanceFunc_;}

	private:

		int nk_;
		std::vector<KMeansData> dataPts_;
		std::vector<KMeansData> means_;
		std::vector<int> dataClusterIds_;

		DistanceF distanceFunc_;
	};
}


#endif//ZKMEANSCLUSTER_H_
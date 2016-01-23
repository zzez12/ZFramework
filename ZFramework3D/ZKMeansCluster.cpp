#include "ZKMeansCluster.h"

namespace ZMathAlgorithms
{

	ZKMeansCluster::ZKMeansCluster()
	{
		distanceFunc_ = KMeansData::defaultDistace;
	}

	ZKMeansCluster::~ZKMeansCluster()
	{
		dataPts_.clear();
		means_.clear();
	}

	void ZKMeansCluster::addDataPoint(const KMeansData& p)
	{
		dataPts_.push_back(p);
	}

	void ZKMeansCluster::setDistanceFunction(DistanceF func)
	{
		distanceFunc_ = func;
	}

	std::vector<KMeansData> ZKMeansCluster::getMeans(int k)
	{
		int nDataPts = dataPts_.size();
		int *nearMeanIndex = new int[nDataPts];
		// initialize
		for (int i=0; i<nDataPts; i++)
		{
			nearMeanIndex[i] = -1;
		}

		// seed means
		KMeansData dataLoBound = KMeansData(MATH::CVN_MAX);//KMeansData(FLT_MAX, FLT_MAX, FLT_MAX);
		KMeansData dataHiBound = KMeansData(MATH::CVN_MIN);//KMeansData(-FLT_MAX, -FLT_MAX, -FLT_MAX);

		means_.clear();

		for (int i=0; i<nDataPts; i++)
		{
			dataLoBound = dataLoBound.minimal(dataPts_[i]);
			dataHiBound = dataHiBound.maximal(dataPts_[i]);
		}

		for (int j=0; j<k; j++)
		{
			KMeansData data = KMeansData(MATH::CVN_RAND);
			data = dataLoBound + (dataHiBound - dataLoBound).pointWiseDot(data);
			means_.push_back(data);
		}

		bool bConverged = false;
		while (!bConverged)
		{
			// assignmentt (expectation step
			for (int i=0; i<nDataPts; i++)
			{
				float fMinDist = FLT_MAX;
				int nMinIndex = -1;
				for (int j=0; j<k; j++)
				{
					float fEachDist = distanceFunc_(means_[j], dataPts_[i]);//(means_[j] - dataPts_[i]).length();
					if (fEachDist<fMinDist) 
					{
						fMinDist = fEachDist;
						nMinIndex = j;
					}
				}
				nearMeanIndex[i] = nMinIndex;
			}

			// update (maximization step)
			bConverged = true;
			for (int j=0; j<k; j++)
			{
				KMeansData newMean = KMeansData(MATH::CVN_Default);
				int nPtsInCluster = 0;
				for (int i=0; i<nDataPts; i++)
				{
					if (j==nearMeanIndex[i])
					{
						newMean += dataPts_[i];
						nPtsInCluster++;
					}
				}

				if (nPtsInCluster>0)
				{
					newMean = newMean*(1.f/(float)nPtsInCluster);
					if (distanceFunc_(newMean, means_[j])>1.f)
						bConverged = false;
					means_[j] = newMean;
				}
			}
		}

		delete []nearMeanIndex;
		return means_;
	}

	std::vector<int> ZKMeansCluster::getClusterIds()
	{
		if (means_.empty())
			return dataClusterIds_;

		dataClusterIds_.resize(dataPts_.size(), -1);
		for (int i=0; i<dataPts_.size(); i++)
		{
			// find nearest mean
			float nearestDist = FLT_MAX;
			int nearestId = -1;
			for (int j=0; j<means_.size(); j++)
			{
				float dis = distanceFunc_(means_[j], dataPts_[i]);
				if (dis<nearestDist)
				{
					nearestDist = dis;
					nearestId = j;
				}
			}
			dataClusterIds_[i] = nearestId;
		}
		return dataClusterIds_;
	}

// 	float ZKMeansCluster::dataDistance(const KMeansData& p, const KMeansData& q)
// 	{
// 		return (p-q).length();
// 	}
// 
// 	float ZKMeansCluster::dataNewDistance(const KMeansData& p, const KMeansData& q)
// 	{
// 		Vec3f np(sin(p.x)*cos(p.y), sin(p.x)*sin(p.y), cos(p.x));
// 		Vec3f nq(sin(q.x)*cos(q.y), sin(q.x)*sin(q.y), cos(q.x));
// 		return acos(np.DotProduct(nq)) + abs(p.z - q.z);
// 	}
// 
// 	float ZKMeansCluster::dataNewDistanceWithoutD(const KMeansData& p, const KMeansData& q)
// 	{
// 		Vec3f np(sin(p.x)*cos(p.y), sin(p.x)*sin(p.y), cos(p.x));
// 		Vec3f nq(sin(q.x)*cos(q.y), sin(q.x)*sin(q.y), cos(q.x));
// 		return acos(np.DotProduct(nq));
// 	}

}
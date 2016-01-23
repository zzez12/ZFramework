#include "ZEigenData.h"
#include <algorithm>

namespace ZMeshSpace
{
	void ZEigenData::normalize()
	{
		maxVal = *max_element(eigVec.begin(), eigVec.end());
		minVal = *min_element(eigVec.begin(), eigVec.end());
		normalizeValues.resize(eigVec.size());
		for (int i=0; i<eigVec.size(); i++)
		{
			normalizeValues[i] = (eigVec[i]-minVal)/(maxVal-minVal);
		}
	}
}
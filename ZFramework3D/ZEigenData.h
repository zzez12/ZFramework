#pragma once

#include "GlobalDefs.h"
#include <vector>

namespace ZMeshSpace
{
	class ZEigenData
	{
	public:
		double eigValue;
		std::vector<double> eigVec;
		double maxVal, minVal;
		std::vector<double> normalizeValues;
		//std::vector<IsoBand> isoBands;

		void normalize();
		//void extractIsoBands(Mesh3D* mesh);
	};
}
#pragma once

#include "../Scene/Mesh3D.h"
#include "../Eigen/Eigen/dense"
#include <vector>

using namespace ZMeshSpace;

namespace ZMeshAlgorithms
{

	class WeightedGeodesics
	{
	public:
		typedef float (*DistanceFunc)(float a, float b);
		static float distance(float a, float b);
		static float inverseDistance(float a, float b);
		static float inverseDistanceAscending(float a, float b);
		static float inverseDistanceDescending(float a, float b);

		static float Dijkstra_ShortestPath(Mesh3D* pMesh, const std::vector<float>& vrtWeights,
			int nStart, int nEnd, std::vector<int>& path_list, DistanceFunc distFunc=distance);
		static float Dijkstra_ShortestPath(Mesh3D* pMesh, const Eigen::VectorXd& vrtWeights,
			int nStart, int nEnd, std::vector<int>& path_list, DistanceFunc distFunc=distance);

	};
}
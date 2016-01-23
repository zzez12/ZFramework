#pragma once

#include "ZMeshSurfacePoint.h"
#include "ZPlane.h"

namespace ZMeshSpace
{
	class ZIsoLine
	{
	public:
		double value;
		std::vector<ZMeshSurfacePoint> linePoints;
		ZPlane plane;
		int clusterId;
	};

	class ZEigenIsoLine : public ZIsoLine
	{
	public:
		int eigenIdx;
		bool isGood_;

		ZEigenIsoLine();

		void projectToPlane();
		bool checkGood();

	private:
		bool isPlane();

	private:
		struct Parameters {
			float fDistPlaneThreshold;
		};
		Parameters para_;
	};
}
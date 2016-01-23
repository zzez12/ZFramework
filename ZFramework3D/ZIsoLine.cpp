#include "ZIsoLine.h"
#include <numeric>

namespace ZMeshSpace
{
	ZEigenIsoLine::ZEigenIsoLine()
		: isGood_(false)
	{
		clusterId = 0;
		para_.fDistPlaneThreshold = 1e-4;
	}

	bool ZEigenIsoLine::checkGood()
	{
		if (plane.getProjectedPoints().size()<=3) return false;
		isGood_ = isPlane();
		return isGood_;
	}

	bool ZEigenIsoLine::isPlane()
	{
		const std::vector<Vec3f>& projPs = plane.getProjectedPoints();
		float err=0;
		std::vector<float> dists;
		int size = projPs.size();
		for (int i=0; i<size; i++)
		{
			dists.push_back((projPs[i]-linePoints[i].pos).length());
		}
		float aver = std::accumulate(dists.begin(), dists.end(), 0.f, std::plus<float>());
		aver /= size;
		for (int i=0; i<size; i++)
		{
			float resi = dists[i] - aver;
			dists[i] = resi*resi;
		}
		err = std::accumulate(dists.begin(), dists.end(), 0.f, std::plus<float>());
		err /= size;
		//std::cout << "Plane check: " << err << "\n";
		return err<para_.fDistPlaneThreshold;
	}

	void ZEigenIsoLine::projectToPlane()
	{
		std::vector<Vec3f> points;
		for (int j=0; j<linePoints.size(); j++)
		{
			points.push_back(linePoints[j].pos);
		}
		plane.buildPlane(points);
	}
}
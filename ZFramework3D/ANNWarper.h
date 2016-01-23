#pragma once

#include "GlobalDefs.h"
#include "ANN.h"
#include <vector>

class AnnWarper
{
public:
	AnnWarper();
	~AnnWarper();

public:
	void init(const std::vector<Vec3f>& points);
	void destroy();

	void queryPoint(const Vec3f& pQuery, int nearNum, std::vector<int>& nnIdx, std::vector<float>& nnDists);
	void queryPointIds(const Vec3f& pQuery, std::vector<int>& ids);
	void queryPointPos(const Vec3f& pQuery, std::vector<Vec3f>& poss);

private:
	ANNkd_tree *ann_kdTree_;
	int nPts_;
	int dim_;
	int k_knn_;
	ANNpointArray		dataPts_;		    		// data points
	ANNpoint			queryPt_;			    	// query point

};
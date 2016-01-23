#pragma once
#ifndef ANNWARPER_EIGEN_H_
#define ANNWARPER_EIGEN_H_

#include "ANN.h"
#include "../Eigen/Eigen/Dense"
#include <vector>

class AnnWarper_Eigen
{
public:
	AnnWarper_Eigen() {
		ann_kdTree_ = NULL;
		dataPts_ = NULL;
		k_knn_ = 10;
		dim_ = 3;
		max_query_time_ = 3;
	}
	~AnnWarper_Eigen() {
		destroy();
	}

public:
	void init(const Eigen::MatrixXf& points);
	void init(const Eigen::MatrixXf& points, const std::vector<bool>& bTags);
	void destroy();

	void setFlags(const std::vector<bool>& flags);
	int queryPoint(const Eigen::VectorXf& pQuery, int nearNum, Eigen::VectorXi& nnIdx, Eigen::VectorXf& nnDists);
	int queryPointIds(const Eigen::VectorXf& pQuery, int nearNum, Eigen::VectorXi& nnIdx);
	int queryFRPoint(const Eigen::VectorXf& pQuery, float sqRad, int k, Eigen::VectorXi& nnIdx, Eigen::VectorXf& nnDists);

private:
	ANNkd_tree *ann_kdTree_;
	int nPts_;
	int dim_;
	int k_knn_;
	ANNpointArray	dataPts_;
	ANNpoint		queryPt_;
	std::vector<bool> ptTags_;

	int max_query_time_;
	//Eigen::VectorXi mapIdx_;
};

#endif//ANNWARPER_EIGEN_H_
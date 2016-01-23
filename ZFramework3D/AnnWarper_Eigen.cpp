#include "AnnWarper_Eigen.h"

void AnnWarper_Eigen::destroy()
{
	if (dataPts_)
	{
		annDeallocPts(dataPts_);
		dataPts_ = NULL;
	}
	if (ann_kdTree_)
	{
		delete ann_kdTree_;
		ann_kdTree_ = NULL;
	}
}

void AnnWarper_Eigen::init(const Eigen::MatrixXf& points)
{
	destroy();
	nPts_ = points.rows();
	dim_ = points.cols();
	dataPts_ = annAllocPts(nPts_, dim_);
	//mapIdx_.resize(nPts_, -1);
	ptTags_.resize(nPts_, true);

	// init dataPts
	for (int i=0; i<nPts_; i++)
	{
		for (int j=0; j<dim_; j++)
		{
			dataPts_[i][j] = points(i, j);
		}
		//mapIdx_(i) = i;
	}

	// build kd-tree
	ann_kdTree_ = new ANNkd_tree(dataPts_, nPts_, dim_);
}

void AnnWarper_Eigen::init(const Eigen::MatrixXf& points, const std::vector<bool>& bTags)
{
	init(points);
	ptTags_ = bTags;
}

void AnnWarper_Eigen::setFlags(const std::vector<bool>& flags)
{
	if (ptTags_.size()!=flags.size()) return;
	ptTags_ = flags;
}

int AnnWarper_Eigen::queryPoint(const Eigen::VectorXf& pQuery, int nearNum, Eigen::VectorXi& nnIdxs, Eigen::VectorXf& nnDists)
{
	if (ann_kdTree_==NULL) return -1;
	if (pQuery.size()!=dim_) return -2;
	k_knn_ = nearNum;
	if (k_knn_>nPts_/2) return -3;

	queryPt_ = annAllocPt(dim_);
	for (int i=0; i<dim_; i++)
	{
		queryPt_[i] = pQuery(i);
	}

	int queriedNum = k_knn_;
	bool bFinish = false;
	for (int times = 0; times<max_query_time_; times++)	// query max_query_time_ time to find at most nearNum points
	{
		ANNidxArray			nnIdx;					    // near neighbor indices
		ANNdistArray		dists;		    			// near neighbor distances
		nnIdx = new ANNidx[k_knn_];
		dists = new ANNdist[k_knn_];
		ann_kdTree_->annkSearch(queryPt_, k_knn_, nnIdx, dists);
		int count = 0;
		for (int i=0; i<k_knn_; i++)
		{
			if (ptTags_[nnIdx[i]])
				count++;
		}
		if (count>=nearNum) 
		{
			queriedNum = nearNum;
			nnIdxs.resize(queriedNum);
			nnDists.resize(queriedNum);
			int j=0;
			for (int i=0; i<k_knn_; i++)
			{
				if (ptTags_[nnIdx[i]])
				{
					nnIdxs(j) = nnIdx[i];
					nnDists(j) = dists[j];
					j++;
				}
				if (j>=nearNum) break;
			}
			bFinish = true;
		}
		else if (times==max_query_time_-1)
		{
			queriedNum = count;
			nnIdxs.resize(queriedNum);
			nnDists.resize(queriedNum);
			int j=0;
			for (int i=0; i<k_knn_; i++)
			{
				if (ptTags_[nnIdx[i]])
				{
					nnIdxs(j) = nnIdx[i];
					nnDists(j) = dists[j];
					j++;
				}
			}
		}
		delete [] nnIdx;
		delete [] dists;
		k_knn_ = k_knn_*2;
		if (bFinish) break;
	}
	return queriedNum;
}

int AnnWarper_Eigen::queryPointIds(const Eigen::VectorXf& pQuery, int nearNum, Eigen::VectorXi& nnIdx)
{
	Eigen::VectorXf nnDist;
	return queryPoint(pQuery, nearNum, nnIdx, nnDist);
}

int AnnWarper_Eigen::queryFRPoint(const Eigen::VectorXf& pQuery, float sqRad, int k, Eigen::VectorXi& nnIdxV, Eigen::VectorXf& nnDistsV)
{	
	if (ann_kdTree_==NULL) return -1;
	if (pQuery.size()!=dim_) return -2;

	//!! the query is using L2^2 distance, NOT L2!!
	sqRad = sqRad*sqRad;

	// first query the number of the in-side points
	queryPt_ = annAllocPt(dim_);
	for (int i=0; i<dim_; i++)
	{
		queryPt_[i] = pQuery(i);
	}

	ANNidxArray			nnIdx;					    // near neighbor indices
	ANNdistArray		dists;		    			// near neighbor distances

	int queriedNum = ann_kdTree_->annkFRSearch(queryPt_, sqRad, 0, nnIdx, dists);

	// query the queriedNum points
	nnIdx = new ANNidx[queriedNum];
	dists = new ANNdist[queriedNum];
	ann_kdTree_->annkFRSearch(queryPt_, sqRad, queriedNum, nnIdx, dists);

	int exactQueriedNum = 0;
	for (int i=0; i<queriedNum; i++)
	{
		if (ptTags_[nnIdx[i]])
			exactQueriedNum++;
	}
	nnIdxV.resize(exactQueriedNum);
	nnDistsV.resize(exactQueriedNum);
	int count = 0;
	for (int i=0; i<queriedNum; i++)
	{
		if (ptTags_[nnIdx[i]])
		{
			nnIdxV(count) = nnIdx[i];
			nnDistsV(count) = dists[i];
			count++;
		}
	}

	delete [] nnIdx;
	delete [] dists;
	return exactQueriedNum;
}
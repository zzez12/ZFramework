#include "ANNWarper.h"

AnnWarper::AnnWarper()
{
	ann_kdTree_ = NULL;
	dataPts_ = NULL;
	dim_ = 3;
	k_knn_ = 10;
}

AnnWarper::~AnnWarper()
{
	destroy();
}

void AnnWarper::destroy()
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

void AnnWarper::init(const std::vector<Vec3f>& points)
{
	nPts_ = points.size();
	dataPts_ = annAllocPts(nPts_, dim_);			// allocate data points

	// init dataPts
	for (int i=0; i<nPts_; i++)
	{
		const Vec3f& p = points[i];
		dataPts_[i][0] = p.x;
		dataPts_[i][1] = p.y;
		dataPts_[i][2] = p.z;
	}

	ann_kdTree_ = new ANNkd_tree(					// build search structure
		dataPts_,					// the data points
		nPts_,						// number of points
		dim_);						// dimension of space

}

void AnnWarper::queryPoint(const Vec3f& pQuery, int nearNum, std::vector<int>& nnIdxV, std::vector<float>& nnDistsV)
{
	if (ann_kdTree_==NULL) return;
	k_knn_ = nearNum;

	queryPt_ = annAllocPt(dim_);

	ANNidxArray			nnIdx;					    // near neighbor indices
	ANNdistArray		dists;		    			// near neighbor distances
	nnIdx = new ANNidx[k_knn_];
	dists = new ANNdist[k_knn_];
	queryPt_[0] = pQuery.x;
	queryPt_[1] = pQuery.y;
	queryPt_[2] = pQuery.z;
	ann_kdTree_->annkSearch(queryPt_, k_knn_, nnIdx, dists);
	nnIdxV.clear();
	nnDistsV.clear();
	for (int i=0; i<k_knn_; i++)
	{
		nnIdxV.push_back(nnIdx[i]);
		nnDistsV.push_back(dists[i]);
	}
	annDeallocPt(queryPt_);
	delete [] nnIdx;
	delete [] dists;
}

void AnnWarper::queryPointIds(const Vec3f& pQuery, std::vector<int>& ids)
{
	if (ann_kdTree_==NULL) return;

	queryPt_ = annAllocPt(dim_);

	ANNidxArray			nnIdx;					    // near neighbor indices
	ANNdistArray		dists;		    			// near neighbor distances
	nnIdx = new ANNidx[k_knn_];
	dists = new ANNdist[k_knn_];
	queryPt_[0] = pQuery.x;
	queryPt_[1] = pQuery.y;
	queryPt_[2] = pQuery.z;
	ann_kdTree_->annkSearch(queryPt_, k_knn_, nnIdx, dists);
	ids.clear();
	for (int i=0; i<k_knn_; i++)
	{
		ids.push_back(nnIdx[i]);
	}
	annDeallocPt(queryPt_);
	delete [] nnIdx;
	delete [] dists;
}

void AnnWarper::queryPointPos(const Vec3f& pQuery, std::vector<Vec3f>& poss)
{
	std::vector<int> ids;
	queryPointIds(pQuery, ids);
	poss.resize(ids.size());
	for (int i=0; i<ids.size(); i++)
	{
		int id = ids[i];
		poss[i] = Vec3f(dataPts_[id][0], dataPts_[id][1], dataPts_[id][2]);
	}
}

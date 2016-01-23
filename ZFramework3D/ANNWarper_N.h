#pragma once

#include "GlobalDefs.h"
#include "ANN.h"
#include "vectorN.h"

template <typename T, int N>
class AnnWarper_N
{
	typedef MATH::CVectorN<T, N> DataType;
public:
	AnnWarper_N() {	
		ann_kdTree_ = NULL;
		dataPts_ = NULL;
		dim_ = N;
		k_knn_ = 10;
	}
	~AnnWarper_N(){destroy();}

public:
	void init(const std::vector<DataType>& points);
	void init(const std::vector<DataType>& points, const std::vector<bool> bTags);
	void destroy();

	void queryPoint(const DataType& pQuery, int nearNum, std::vector<int>& nnIdx, std::vector<T>& nnDists);
	void queryPointIds(const DataType& pQuery, std::vector<int>& ids);
	void queryPointPos(const DataType& pQuery, std::vector<DataType>& poss);

private:
	ANNkd_tree *ann_kdTree_;
	int nPts_;
	int dim_;
	int k_knn_;
	ANNpointArray		dataPts_;		    		// data points
	ANNpoint			queryPt_;			    	// query point
	std::vector<int>	mapIdx_;	// a map when using init(points, bTags)
};

// the implements
template <typename T, int N>
void AnnWarper_N<T,N>::destroy()
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
	mapIdx_.clear();
}

template <typename T, int N>
void AnnWarper_N<T,N>::init(const std::vector<DataType>& points)
{
	destroy();
	nPts_ = points.size();
	dataPts_ = annAllocPts(nPts_, dim_);			// allocate data points

	// init dataPts
	for (int i=0; i<nPts_; i++)
	{
		const DataType& p = points[i];
		for (int j=0; j<dim_; j++)
		{
			dataPts_[i][j] = p[j];
		}
		// the map
		mapIdx_.push_back(i);
	}

	ann_kdTree_ = new ANNkd_tree(					// build search structure
		dataPts_,					// the data points
		nPts_,						// number of points
		dim_);						// dimension of space

}

template <typename T, int N>
void AnnWarper_N<T,N>::init(const std::vector<DataType>& points, const std::vector<bool> bTags)
{
	destroy();
	if (bTags.size()!=points.size()) return init(points);

	nPts_ = 0;
	for (int i=0; i<bTags.size(); i++)
	{
		if (bTags[i]==false) 
		{
			mapIdx_.push_back(i);
			nPts_++;
		}
	}
	if (nPts_==0) return;

	//nPts_ = points.size();
	dataPts_ = annAllocPts(nPts_, dim_);			// allocate data points

	// init dataPts
	int nId = 0;
	for (int i=0; i<nPts_; i++)
	{
		if (bTags[i]==true) continue;
		const DataType& p = points[i];
		for (int j=0; j<dim_; j++)
		{
			dataPts_[nId][j] = p[j];
		}
		nId++;
	}

	ann_kdTree_ = new ANNkd_tree(					// build search structure
		dataPts_,					// the data points
		nPts_,						// number of points
		dim_);						// dimension of space

}

template <typename T, int N>
void AnnWarper_N<T,N>::queryPoint(const DataType& pQuery, int nearNum, std::vector<int>& nnIdxV, std::vector<T>& nnDistsV)
{
	if (ann_kdTree_==NULL) return;
	k_knn_ = nearNum;
	if (k_knn_>nPts_/2) return;

	queryPt_ = annAllocPt(dim_);

	ANNidxArray			nnIdx;					    // near neighbor indices
	ANNdistArray		dists;		    			// near neighbor distances
	nnIdx = new ANNidx[k_knn_];
	dists = new ANNdist[k_knn_];
	for (int i=0; i<dim_; i++)
	{
		queryPt_[i] = pQuery[i];
	}
	ann_kdTree_->annkSearch(queryPt_, k_knn_, nnIdx, dists);
	nnIdxV.clear();
	nnDistsV.clear();
	for (int i=0; i<k_knn_; i++)
	{
		//nnIdxV.push_back(nnIdx[i]);
		nnIdxV.push_back(mapIdx_[nnIdx[i]]);
		nnDistsV.push_back(dists[i]);
	}
	annDeallocPt(queryPt_);
	delete [] nnIdx;
	delete [] dists;
}

template <typename T, int N>
void AnnWarper_N<T,N>::queryPointIds(const DataType& pQuery, std::vector<int>& ids)
{
	if (ann_kdTree_==NULL) return;
	if (k_knn_>nPts_/2) return;

	queryPt_ = annAllocPt(dim_);

	ANNidxArray			nnIdx;					    // near neighbor indices
	ANNdistArray		dists;		    			// near neighbor distances
	nnIdx = new ANNidx[k_knn_];
	dists = new ANNdist[k_knn_];
	for (int i=0; i<dim_; i++)
	{
		queryPt_[i] = pQuery[i];
	}
	ann_kdTree_->annkSearch(queryPt_, k_knn_, nnIdx, dists);
	ids.clear();
	for (int i=0; i<k_knn_; i++)
	{
		//ids.push_back(nnIdx[i]);
		ids.push_back(mapIdx_[nnIdx[i]]);
	}
	annDeallocPt(queryPt_);
	delete [] nnIdx;
	delete [] dists;
}

template <typename T, int N>
void AnnWarper_N<T,N>::queryPointPos(const DataType& pQuery, std::vector<DataType>& poss)
{
	std::vector<int> ids;
	queryPointIds(pQuery, ids);
	poss.resize(ids.size());
	for (int i=0; i<ids.size(); i++)
	{
		int id = ids[i];
		//poss[i] = Vec3f(dataPts_[id][0], dataPts_[id][1], dataPts_[id][2]);
		poss[i] = DataType(dataPts_[id]);
	}
}


// typedef
typedef AnnWarper_N<float,3> AnnWarper3;
typedef AnnWarper_N<float,5> AnnWarper5;
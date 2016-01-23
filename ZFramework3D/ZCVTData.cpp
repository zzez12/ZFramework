#include "ZCVTData.h"

namespace ZCVT
{
	bool ZCVTData::init(int nSample) 
	{
		destroy();
		lines_ = Random_Generator::randomNonIntersectLines(nSample);
		nSamples_ = lines_.size();
		lineCells_.resize(nSamples_);
		for (int i=0; i<nSamples_; i++)
			lineCells_[i].id_ = i;
		return nSample==nSamples_;
	}

	void ZCVTData::destroy()
	{
		nSamples_ = 0;
		lines_.clear();
		lineCells_.clear();
		voronoi_.clear();
		mpLine2VoronoiCells_.clear();
	}

	void ZCVTData::setCorresponding(int lineIdx, const std::vector<int>& corresVIdx)
	{
		mpLine2VoronoiCells_[&lines_[lineIdx]] = corresVIdx;
	}

	Vec2f ZCVTData::getSamplePoint(int lineIdx, int lineSampleIdx)
	{
		assert(lineIdx>=0 && lineIdx<lines_.size());
		assert(lineSampleIdx>=0 && lineSampleIdx<=subSize_);

		ZLine& line = lines_[lineIdx];

		return line.p0 + (line.p1-line.p0)*(1.f*lineSampleIdx/subSize_);
	}

	void ZCVTData::getLineEndPoints(int lineIdx, double* x)
	{
		x[0] = lines_[lineIdx].p0.x;
		x[1] = lines_[lineIdx].p0.y;
		x[2] = lines_[lineIdx].p1.x;
		x[3] = lines_[lineIdx].p1.y;
	}

	void ZCVTData::setLineEndPoints(int lineIdx, const Vec2f& p0, const Vec2f& p1)
	{
		lines_[lineIdx].p0 = p0;
		lines_[lineIdx].p1 = p1;
	}

}
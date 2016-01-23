#pragma once
#ifndef ZCVT_DATA_H_H
#define ZCVT_DATA_H_H

#include "GlobalDefs.h"
#include "ZLine.h"
#include "randomGenerator.h"
#include "ZVoronoi.h"
#include <vector>
#include <map>

namespace ZCVT
{
	class ZCVTLineCell
	{
	public:
		int id_;
		std::vector<int> neighborLineIds_;
	};

	class ZCVTData
	{
	private:
		int nSamples_;
		std::vector<ZLine> lines_;

		std::vector<ZCVTLineCell> lineCells_;
		int subSize_;	// each line will be subdivided into (subSize_) size

		ZVoronoi voronoi_;
		std::map<ZLine*, std::vector<int>> mpLine2VoronoiCells_;	// each line correspond to many voronoi cells

	public:
		ZCVTData(int nSamples=0) : nSamples_(nSamples), subSize_(20){}
		~ZCVTData() {destroy();}

		bool init(int nSample);
		void destroy();

		void getLineEndPoints(int lineIdx, double* x);
		void setLineEndPoints(int lineIdx, const Vec2f& p0, const Vec2f& p1);

		int getVoronoiCellIdx(int lineIdx, int lineSampleIdx) {return mpLine2VoronoiCells_[&lines_[lineIdx]][lineSampleIdx];}
		void setCorresponding(int lineIdx, const std::vector<int>& corresVIdx); 

		Vec2f getSamplePoint(int lineIdx, int lineSampleIdx); 
		int getSubSize() {return subSize_;}
		int numSampleSize() {return numLines()*(subSize_+1);}
		int sampleOnLine(int sampleIdx) {return sampleIdx/(subSize_+1);}
		int sampleOnLineInterval(int sampleIdx) {return sampleIdx-sampleOnLine(sampleIdx)*(subSize_+1);}

		int numLines() {return lines_.size();}
		std::vector<ZLine>& getLines() {return lines_;}

		ZVoronoi* getVoronoi() {return &voronoi_;}

		ZLine* getLine(int idx) {return &lines_[idx];}
		ZCVTLineCell* getLineCell(int idx) {return &lineCells_[idx];}
		int getNeighborSize(int lineIdx) {return getLineCell(lineIdx)->neighborLineIds_.size();}
		int getNeighborLine(int lineIdx, int neiIdx) {return getLineCell(lineIdx)->neighborLineIds_[neiIdx];}
	};
}

#endif//ZCVT_DATA_H_H
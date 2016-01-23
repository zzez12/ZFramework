#pragma once

#include "../Scene/Mesh3D.h"
#include "ZMeshAlgorithms.h"
#include "ZCVTAlgorithms.h"
#include <vector>


class ZDataManager
{
private:
	static ZDataManager* dataManager_;	
	ZDataManager(){init();}
	class Cleaner
	{
	public:
		Cleaner();
		~Cleaner();
	};
	void init();

public:
	~ZDataManager();

	static ZDataManager* getDataManager();

	void clear();
	void addMesh(ZMeshSpace::Mesh3D* mesh);
	void clearMeshes();
	void addCVTData(ZCVT::ZCVTData* cvtData);
	void clearCVTData();

	ZMeshSpace::Mesh3D* getMesh();
	ZMeshSpace::ZMeshAlgorithms* getAlgorithmHandler();

	ZCVT::ZCVTData* getCVTData();
	ZCVT::ZCVTAlgorithms* getCVTAlgorithmHandler();

	ZMeshSpace::Mesh3D* getMesh(int idx);
	ZCVT::ZCVTData* getCVTData(int idx);

private:
	std::vector<ZMeshSpace::Mesh3D*> meshes_;
	ZMeshSpace::ZMeshAlgorithms* meshAlgorithms_;

	ZCVT::ZCVTAlgorithms* cvtAlgorithms_;
	std::vector<ZCVT::ZCVTData*> cvtDatas_;
};
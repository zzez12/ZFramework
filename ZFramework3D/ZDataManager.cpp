#include "ZDataManager.h"

ZDataManager* ZDataManager::dataManager_ = NULL;

ZDataManager::Cleaner::Cleaner(){}
ZDataManager::Cleaner::~Cleaner()
{
	if (getDataManager())
	{
		delete getDataManager();
	}
}

ZDataManager::~ZDataManager()
{
	if (!meshes_.empty())
	{
		std::vector<ZMeshSpace::Mesh3D*>::iterator iter = meshes_.begin();
		do
		{
			ZMeshSpace::Mesh3D* mesh = *iter;
			delete mesh;
			iter++;
		} while(iter!=meshes_.end());
		meshes_.clear();
	}


	if (meshAlgorithms_)
	{
		delete meshAlgorithms_;
		meshAlgorithms_ = NULL;
	}
}

void ZDataManager::clear()
{
	clearMeshes();
	clearCVTData();
	SAFE_DELETE(meshAlgorithms_);
	SAFE_DELETE(cvtAlgorithms_);
	init();
}

void ZDataManager::clearMeshes()
{
	if (!meshes_.empty())
	{
		std::vector<ZMeshSpace::Mesh3D*>::iterator iter = meshes_.begin();
		do
		{
			ZMeshSpace::Mesh3D* mesh = *iter;
			delete mesh;
			iter++;
		} while(iter!=meshes_.end());
		meshes_.clear();
	}

}

void ZDataManager::clearCVTData()
{
	if (!cvtDatas_.empty())
	{
		std::vector<ZCVT::ZCVTData*>::iterator iter = cvtDatas_.begin();
		do
		{
			ZCVT::ZCVTData* cvt = *iter;
			delete cvt;
			iter++;
		} while(iter!=cvtDatas_.end());
		cvtDatas_.clear();
	}
}

void ZDataManager::init()
{
	meshAlgorithms_ = new ZMeshSpace::ZMeshAlgorithms;
	cvtAlgorithms_ = new ZCVT::ZCVTAlgorithms;
}

ZDataManager* ZDataManager::getDataManager()
{
	if (dataManager_==NULL)
	{
		dataManager_ = new ZDataManager;
		static Cleaner cleaner;
	}

	return dataManager_;
}

void ZDataManager::addMesh(ZMeshSpace::Mesh3D* mesh)
{
	meshes_.push_back(mesh);
}

void ZDataManager::addCVTData(ZCVT::ZCVTData* cvt)
{
	cvtDatas_.push_back(cvt);
}

ZMeshSpace::Mesh3D* ZDataManager::getMesh(int idx)
{
	return meshes_[idx];
}

ZMeshSpace::Mesh3D* ZDataManager::getMesh()
{
	return meshes_.empty()? NULL : meshes_[0];
}

ZMeshSpace::ZMeshAlgorithms* ZDataManager::getAlgorithmHandler()
{
	return meshAlgorithms_;
}

ZCVT::ZCVTAlgorithms* ZDataManager::getCVTAlgorithmHandler()
{
	return cvtAlgorithms_;
}

ZCVT::ZCVTData* ZDataManager::getCVTData()
{
	return cvtDatas_.empty() ? NULL : cvtDatas_[0];
}

ZCVT::ZCVTData* ZDataManager::getCVTData(int idx)
{
	return cvtDatas_[idx];
}
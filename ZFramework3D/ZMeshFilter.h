#pragma once
#ifndef ZMESHFILTER_H_
#define ZMESHFILTER_H_

#include "GlobalDefs.h"
#include "../Scene/Mesh3D.h"
#include "ZQueryMesh.h"
#include "ZKernelFunc.h"

namespace ZMeshSpace
{
	enum ZMeshFilterType
	{
		MESH_FILTER_TYPE = 0,
		MESH_FILTER_GAUSSIAN = 1,
		MESH_FILTER_GAUSSIAN2,
		MESH_FILTER_WEIGHTED_GAUSSIAN,
		MESH_FILTER_BILATERAL,
		MESH_FILTER_JOINT_BILATERAL,
		MESH_FILTER_MANIFOLD,
	};

	class ZMeshFilterParas
	{
	public:
		float spatial_sigma;
		float range_sigma;
		int tree_height;
	};

	class ZMeshFilterParaNames
	{
	public:
		static const std::string SpatialSigma;// = "SpatialSigma";
		static const std::string RangeSigma;// = "RangeSigma";
	};

	class ZMeshFilter
	{
	public:
		typedef float (*KernelFuncA) (const Eigen::VectorXf& v1, const Eigen::VectorXf& v2, float sigma);
		//typedef float (*KernelFuncA) (const Eigen::Vector3f& v1, const Eigen::Vector3f& v2, float sigma);
		//typedef float (*KernelFuncA) (const Vec3f& v1, const Vec3f& v2, float sigma);

		ZMeshFilter(Mesh3D* mesh, ZMeshFilterType type=MESH_FILTER_TYPE)
			: pMesh_(mesh), type_(type), kernelFuncA_(NULL){queryMeshTool_=new ZQueryMesh(pMesh_);}
		virtual ~ZMeshFilter(){
			SAFE_DELETE(queryMeshTool_);
		}

		virtual bool apply() = 0;
		virtual bool apply(const Eigen::MatrixXf& input, const std::vector<bool>& tags) = 0;
		virtual bool apply(const Eigen::MatrixXf& input, const Eigen::VectorXf& weights, const std::vector<bool>& tags){return false;}
		virtual bool apply(const Eigen::MatrixXf& input, int spatialDim, int rangeDim, 
			const Eigen::VectorXf& weights=Eigen::VectorXf(), const std::vector<bool>& tags=std::vector<bool>()) {return false;}
		virtual bool apply(const Eigen::MatrixXf& input, const Eigen::MatrixXf& joint_input, 
			const Eigen::VectorXf& weights=Eigen::VectorXf(), const std::vector<bool>& tags=std::vector<bool>()){return false;}

		virtual void setPara(const std::string& name, float value){
			if (name.compare(ZMeshFilterParaNames::SpatialSigma)==0)
			{
				filterPara_.spatial_sigma = value;
			}
			else if (name.compare(ZMeshFilterParaNames::RangeSigma)==0)
			{
				filterPara_.range_sigma = value;
			}
		}

		void setParas(ZMeshFilterParas paras) {
			filterPara_.range_sigma = paras.range_sigma;
			filterPara_.spatial_sigma = paras.spatial_sigma;
			filterPara_.tree_height = paras.tree_height;
		}

		void setQueryToolType(QueryMeshType type) {
			queryMeshTool_->setQueryType(type);
		}

		void setRangeDim(int dim) {rangeDim_=dim;}
		void setSpatialDim(int dim) {spatialDim_=dim;}

		int getRangeDim() {return rangeDim_;}
		int getSpatialDim() {return spatialDim_;}
		int getPointSize() {return output_.rows();}

		void setKernelFunc(KernelFuncA func) {kernelFuncA_=func;}
		void setMesh(Mesh3D* mesh) {pMesh_ = mesh;}
		const Eigen::MatrixXf& getResult() {return output_;}

	protected:
		Mesh3D* pMesh_;	// only a pointer
		ZMeshFilterType type_;
		int rangeDim_;
		int spatialDim_;
		ZQueryMesh *queryMeshTool_;

		KernelFuncA kernelFuncA_;
		ZMeshFilterParas filterPara_;
		//std::vector<Eigen::Vector3f> output_;
		Eigen::MatrixXf output_;
	};

}

#endif//ZMESHFILTER_H_
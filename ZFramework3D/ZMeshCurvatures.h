#pragma once
#ifndef ZMESH_CURVATURES_H_
#define ZMESH_CURVATURES_H_

#include "GlobalDefs.h"
#include "../Scene/Mesh3D.h"
#include "../Eigen/Eigen/Dense"

namespace ZMeshSpace
{
	class ZMeshCurvatures
	{
	public:
		ZMeshCurvatures(Mesh3D* pMesh);

		bool apply();

		Eigen::VectorXf getDirectionalCurvatures() {return directionalCurvatures_;}
		const Eigen::VectorXf& getDirectionalCurvatures() const {return directionalCurvatures_;}

		std::vector<Eigen::Vector3f> getPrincipalDirections1() {return principalDirections1_;}
		const std::vector<Eigen::Vector3f> getPrincipalDirections1() const {return principalDirections1_;}
		std::vector<Eigen::Vector3f> getPrincipalDirections2() {return principalDirections2_;}
		const std::vector<Eigen::Vector3f> getPrincipalDirections2() const {return principalDirections2_;}

		// for testing
		void setToMesh();

	private:
		Mesh3D* mesh_;
		std::vector<Eigen::Vector3f> weightedNormals_;	// N_i
		Eigen::VectorXf directionalCurvatures_;		// K_ij
		std::vector<Eigen::Matrix3f> curvatureTensor_;	// M_i
		std::vector<Eigen::Vector3f> principalDirections1_, principalDirections2_;
	};
}

#endif//ZMESH_CURVATURES_H_
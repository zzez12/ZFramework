#pragma once
#ifndef ZFILE_HELPER_H
#define ZFILE_HELPER_H_

#include "../Eigen/Eigen/Dense"
#include <fstream>

namespace ZFileHelper
{
	static void saveEigenMatrixToFile(const char* fileName, const Eigen::MatrixXf& mat)
	{
		std::ofstream fout(fileName);
		fout << mat << "\n";
		fout.close();
	}

	static void saveEigenVectorToFile(const char* fileName, const Eigen::MatrixXf& vec)
	{
		std::ofstream fout(fileName);
		fout << vec << "\n";
		fout.close();
	}
}

#endif//ZFILE_HELPER_H_
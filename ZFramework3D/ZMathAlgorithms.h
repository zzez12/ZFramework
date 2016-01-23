#pragma  once

#include "GlobalDefs.h"
#include "../Math/vector3.h"
#include "../Eigen/svd.h"
#include <vector>
#include <cassert>

using namespace std;

namespace ZMathAlgorithms
{
	static float zMax(float a, float b)
	{
		return a<b ? b : a;
	}

	static void calPCA(const vector<Vec3f>& input, const Vec3f& ori, 
		vector<float>& values, vector<Vec3f>& axes)
	{
		int m = input.size()+1;

		// init the matrix A
		float** A;
		A = new float* [m];
		A[0] = new float [4]; 
		int j = 1;

		for (int i=0; i<input.size(); i++)
		{
			Vec3f p = input.at(i) - ori;
			A[j] = new float [4]; 
			A[j][1] = p.x;
			A[j][2] = p.y;
			A[j][3] = p.z;
			j++;
		}
		assert(j == m);

		// init the matrix V and the vector W
		float **V;
		V=new float*[4];
		V[0]=new float[4];
		V[1]=new float[4];
		V[2]=new float[4];
		V[3]=new float[4];

		float W[4];

		// SVD
		svdcmp(A, m-1, 3, W, V);

		// sort
		int order[3] = {1, 2, 3};

		for (int i=0; i<2; i++)
		{
			if ((W[order[i]]*W[order[i]]) < (W[order[i+1]]*W[order[i+1]]))
			{
				int temp = order[i];
				order[i] = order[i+1]; 
				order[i+1] = temp;
			}
		}

		if ((W[order[0]]*W[order[0]]) < (W[order[1]]*W[order[1]]))
		{
			int temp = order[0];
			order[0] = order[1]; 
			order[1] = temp;
		}

		// save
		values.push_back(W[order[0]]*W[order[0]]); // first 
		values.push_back(W[order[1]]*W[order[1]]); // second
		values.push_back(W[order[2]]*W[order[2]]); // third

		MATH::vector3f axis;     
		axis.x = V[1][order[0]];          // first 
		axis.y = V[2][order[0]];
		axis.z = V[3][order[0]];
		axes.push_back(axis);

		axis.x = V[1][order[1]];          // second
		axis.y = V[2][order[1]];
		axis.z = V[3][order[1]];
		axes.push_back(axis);

		axis.x = V[1][order[2]];          // third
		axis.y = V[2][order[2]];
		axis.z = V[3][order[2]];
		axes.push_back(axis);

		// free memory
		for (int i=0; i<4; i++)     // V
		{
			if (V[i]) delete [] V[i];
		}
		if (V) delete [] V;

		for (int i=0; i<m; i++) // A
		{
			if (A[i]) delete [] A[i];
		}
		if (A) delete [] A;

	}
}
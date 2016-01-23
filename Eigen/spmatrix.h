#ifndef __J_CSRMATRIX_H_
#define __J_CSRMATRIX_H_

#include <vector>
#include <map>
#include <cassert>
#include <algorithm>
using namespace std;

namespace numc
{

template <typename T>
class CSRMatrix
{
public:
	//enum MAT_TYPE{RSS=1, RSP=2, RSI=-2, CSS=3, CHP=4, CHI=-4, CS=6, RU=11, CU=13};
	int m,n;
	bool symmetric;
	vector<T> v;
	vector<int> i, j;
	CSRMatrix<T>():m(0),n(0),symmetric(false) {};


	T getElement(int x, int y) {
		if (x<0||y<0)
			return 0;
		if (symmetric && x>y)
			swap(x, y);
		vector<int>::iterator it = lower_bound(j.begin()+i[x], j.begin()+i[x+1], y);
		if (it==j.begin()+i[x+1] || *it!=y) {
			return 0;
		}
		else {
			return v[it-j.begin()];
		}
	}
	bool ChangeBase(bool oneBase) {
		int ii;
		if (oneBase) {
			assert(this->i[0] == 0);
			if (i[0] != 0)
				return false;

			for (ii=0; ii<i.back(); ii++) {
				j[ii]++;
			}
			for (ii=0; ii<m+1; ii++) {
				i[ii]++;
			}
		}
		else {
			assert(i[0] == 1);
			if (i[0] != 1)
				return false;

			for (ii=0; ii<i.back()-1; ii++) {
				j[ii]--;
			}
			for (ii=0; ii<m+1; ii++) {
				i[ii]--;
			}
		}
		return true;
	}
};

/// Create CSR matrix from vector of rowmap
template <typename T>
void CreateCSRMatrixFromRowMap(CSRMatrix<T> &matC, const vector<map<int, T> >& rowC, const int nCols)
{
	matC.m = rowC.size();
	matC.n = nCols;
	matC.i.clear();
	matC.j.clear();
	matC.v.clear();

	matC.i.reserve(matC.m+1);
	int i,nnz=0;
	for (i=0; i<matC.m; i++) {
		nnz += rowC[i].size();
	}
	matC.j.reserve(nnz);
	matC.v.reserve(nnz);

	// copy rows into matC
	matC.i.push_back(0);
	for (i=0; i<matC.m; i++) {
		matC.i.push_back(matC.i.back());
		for (map<int, T>::const_iterator it=rowC[i].begin(); it!=rowC[i].end(); it++) {
			matC.i.back()++;
			matC.j.push_back(it->first);
			matC.v.push_back(it->second);
		}
	}
}


//////////////////////////////////////////////////////////////////////////
/// Computes the transpose of a matrix.
template <typename T>
void CSRMatrixTranspose(const CSRMatrix<T> &matA, CSRMatrix<T> &matAT)
{
	matAT.m = matA.n;
	matAT.n = matA.m;
	matAT.symmetric = matA.symmetric;

	if (matA.symmetric) {		// symmetric - just copy the matrix
		matAT = matA;
		printf("Symmetric Matrix\n");
		return;
	}

	// non-symmetric matrix -> need to build data structure.
	// we'll go over the columns and build the rows
	vector<map<int, T> > rowC(matA.n);
	for (int i=0; i<matA.m; i++) {
		for (int j=matA.i[i]; j<matA.i[i+1]; j++) {
			rowC[matA.j[j]][i] = matA.v[j];
		}
	}
//	printf("liguo2\n");
	CreateCSRMatrixFromRowMap(matAT, rowC, matA.m);
}

//////////////////////////////////////////////////////////////////////////
// multiplication of sparse matrix
// Assuming nothing about the result (the result is not stored symmetric)
template <typename T>
bool Mul2Matrices(const CSRMatrix<T> &matA, const CSRMatrix<T> &matB, CSRMatrix<T> &matC)
{
	// Compatibility of dimensions
	if (matA.n != matB.m)
		return false;

	// (m x n)*(n x k) = (m x k)
	const int m=matA.m;
	const int n=matA.n;
	const int k=matB.n;

	T aiv, valB;
	int colInd, colB;
	vector<map<int, T> > rowsC(m);
	for (int i=0; i<m; ++i) {					// creating row i of C
		map<int, T> &mapRow2Val = rowsC[i];
		for (int iAi = matA.i[i]; iAi < matA.i[i+1]; ++iAi) {			// travel on ai
			colInd = matA.j[iAi];
			aiv = matA.v[iAi];
			// make aiv*b_{rowInd} and insert into mapRow2Val
			for (int iB=matB.i[colInd]; iB<matB.i[colInd+1]; ++iB) {
				colB=matB.j[iB];
				valB=matB.v[iB];
				// insert valA*aiv into map
				map<int, T>::iterator it = mapRow2Val.find(colB);
				if (it == mapRow2Val.end()) {		// first time
					mapRow2Val[colB] = valB*aiv;
				}
				else {
					it->second = it->second + valB*aiv;
				}
			}
		}
	}// now column i is created

	CreateCSRMatrixFromRowMap(matC, rowsC, k);						// modified by jianwei hu @ 16/09/07
	return true;
}


//////////////////////////////////////////////////////////////////////////
// multiplication of sparse matrix
// The result is symmetric
template <typename T>
bool Mul2MatricesSymmResult(const CSRMatrix<T> &matA, const CSRMatrix<T> &matB, CSRMatrix<T> &matC)
{
	// Compatibility of dimensions
	if (matA.n != matB.m || matA.m != matB.n)
		return false;

	// (m x n)*(n x m) = (m x m)
	const int m=matA.m;
	const int n=matA.n;

	T aiv, valB;
	int colInd, colB;
	vector<map<int, T> > rowsC(m);
	for (int i=0; i<m; ++i) {					// creating row i of C
		map<int, T> &mapRow2Val = rowsC[i];
		for (int iAi = matA.i[i]; iAi < matA.i[i+1]; ++iAi) {			// travel on ai
			colInd = matA.j[iAi];
			aiv = matA.v[iAi];
			// make aiv*b_{colInd} and insert into mapRow2Val
			for (int iB=matB.i[colInd]; iB<matB.i[colInd+1]; ++iB) {
				colB=matB.j[iB];
				if (colB >= i) {
					valB=matB.v[iB];
					// insert valA*aiv into map
					map<int, T>::iterator it = mapRow2Val.find(colB);
					if (it == mapRow2Val.end()) {		// first time
						mapRow2Val[colB] = valB*aiv;
					}
					else {
						it->second = it->second + valB*aiv;
					}
				}
			}
		}
	}// now column i is created

	matC.symmetric = true;
	CreateCSRMatrixFromRowMap(matC, rowsC, m);
	return true;
}

}

#endif		//#def __J_CSRMATRIX_H_
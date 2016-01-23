#ifndef __J_MKL_ADDON_H_
#define __J_MKL_ADDON_H_

#include "spmatrix.h"
#include <mkl.h>
#include <mkl_spblas.h>

//using namespace numc;

namespace numc
{
typedef double mklReal;

bool MklSolveSparseSystem(CSRMatrix<mklReal> &matA, double *const b, double *const x, int nRhs=1);

bool LeastSquareSolver(const CSRMatrix<mklReal> &matA, double *const b, double *const x, int nRhs=1);

class CSparseSolver
{
public:
	CSparseSolver(const CSRMatrix<mklReal> *pMatA=NULL);
	~CSparseSolver();
	bool solve(double *b, double *x, int nRhs=1);
	bool init(const CSRMatrix<mklReal> *pMatA=NULL);
	CSRMatrix<mklReal>* getMatA(){
		return &m_matA;
	}

private:
	CSRMatrix<mklReal> m_matA;

	mklReal *m_pa;
	int *m_pia;
	int *m_pja;

	std::vector<void*> m_ppt;		// Internal solver memory pointer
	std::vector<int> m_piparm;		// Pardiso control parameters
	int m_pmtype;					// Real unsymmetric matrix
	int m_pmaxfct, m_pmnum, m_pphase, m_perror, m_pmsglvl;
	/* Auxiliary variables. */
	double m_pddum;				// Double dummy
	int m_pidum;				// Integer dummy
	int m_pn;
	int m_pRhs;					// nRhs: Number of right hand sides.
};

class CLeastSquareSpareSolver
{
public:
	CLeastSquareSpareSolver(){;}
	~CLeastSquareSpareSolver(){;}
	inline bool solve(double *b, double *x, int nRhs=1)
	{
		/// right side updating
		int m = m_matAT.n;
		int n = m_matAT.m;
		vector<mklReal> ATb(m*nRhs, 0);									//ATb = A' x b
		for (int i=0; i<nRhs; i++) {
			mkl_dcsrgemv("N", &n, &m_matAT.v.front(), &m_matAT.i.front(), &m_matAT.j.front(), b+i*m, &ATb.front()+i*n);
		}

		/// solver starts
		return m_sparseSolver.solve(&ATb.front(), x, nRhs);
	}

	bool init(const CSRMatrix<mklReal> &matA);

private:
	CSRMatrix<mklReal> m_matAT;
	CSparseSolver m_sparseSolver;
};


}

#endif		//#def __J_MKL_ADDON_H_
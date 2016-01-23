#ifndef NUMC_MKL_ADDON_H_
#define NUMC_MKL_ADDON_H_

#include "spmatrixNew.h"
#include <mkl.h>
#include <mkl_spblas.h>

#ifdef P4
#undef P4
#endif


_NUMC_BEGIN

namespace numcNew
{

typedef double mklReal;


inline void MulMatrix2Vector(const CSRMatrix<mklReal> &matA, const double *const b, double *const Ab, int nVect, bool transMatA=false){
	//Ab = A x b
	int n = matA.nRow();

	for(int i=0; i<nVect; i++){
		if( matA.issymm() ){
			// symmetric matrix lower triangle non-zero are saved
			mkl_dcsrsymv("U", &n,
				const_cast<mklReal*>(&matA.mAv[0]), 
				const_cast<int*>(&matA.mAi[0]),
				const_cast<int*>(&matA.mAj[0]),
				const_cast<mklReal*>(b+n*i), Ab+n*i);
			return;
		}
		
		mkl_dcsrgemv(transMatA?"T":"N", &n,
			const_cast<mklReal*>(&matA.mAv[0]), 
			const_cast<int*>(&matA.mAi[0]), 
			const_cast<int*>(&matA.mAj[0]), 
			const_cast<mklReal*>(b), 
			Ab);

	}
	
}


class SparseSolver
{
public:
	SparseSolver(const CSRMatrix<mklReal> *pMatA=NULL):mIsInitiated(false)
	{ clear(); if(pMatA) mMatA = *pMatA; }

	SparseSolver(const RowMat<mklReal>& rMat):mIsInitiated(false)
	{ clear(); mMatA = CSRMatrix<mklReal>(rMat); }

	~SparseSolver() { clear(); }
	void clear();
	bool solve(const double* b, double* x, int nRhs=1, int nItr=0, int nPrecDigit=-1);
	bool init(const CSRMatrix<mklReal> *pMatA=NULL);
	inline bool init(const RowMat<mklReal> &rMatA)
	{ mMatA = CSRMatrix<mklReal>(rMatA);	return init(); }
	inline void setMatrix(const CSRMatrix<mklReal>& mat){
		if( mMatA.empty() ){
			fprintf(stderr, "error!trying to set empty matrix!");
			return;
		}
	
		mMatA = mat;
		pA = &mMatA.mAv[0];
		pIa = &mMatA.mAi[0];
		pJa = &mMatA.mAj[0];
		mMatA.ChangeBase(true);
	}
	inline CSRMatrix<mklReal>& getMatA()	{ return mMatA; }
	inline bool isInitiated() const { return mIsInitiated; }

private:
	CSRMatrix<mklReal> mMatA;
	bool mIsInitiated;

	mklReal *pA;
	int *pIa;
	int *pJa;

	void *pPt[64];		// Internal solver memory pointer
	int pIparm[64];		// Pardiso control parameters
	int pMtype;					// Real unsymmetric matrix
	int pMaxfct, pMnum, pPhase, pError;
	/* Auxiliary variables. */
	double pDdum;				// Double dummy
	int pIdum;				// Integer dummy
	int pN;
	int pRhs;					// nRhs: Number of right hand sides.

public:
	static int pMsglvl;
};

class LeastSquareSparseSolver
{
public:
	LeastSquareSparseSolver()	{}
	~LeastSquareSparseSolver()	{}

	inline bool solve(const double *b, double *x, int nRhs=1)
	{
		/// right side updating
		int m = mMatAT.nCol();
		int n = mMatAT.nRow();
		std::vector<mklReal> ATb(m*nRhs, 0);									//ATb = A' x b
		for (int i=0; i<nRhs; i++) {
			mkl_dcsrgemv("N", &n, &mMatAT.mAv[0], &mMatAT.mAi[0], &mMatAT.mAj[0], (double*)b+i*m, &ATb[0]+i*n);
		}

		/// solver starts
		return mSparseSolver.solve(&ATb[0], x, nRhs);
	}

	bool init(const CSRMatrix<mklReal> &matA);

private:
	CSRMatrix<mklReal> mMatAT;
	SparseSolver mSparseSolver;
};

typedef std::map<int,double> CSRRow;


inline void AddCoeff2CSRRow(CSRRow &row, int idx, double val){
	if( fabs(val)<std::numeric_limits<double>::epsilon() ) return;
	CSRRow::iterator it = row.find(idx);
	if(it == row.end()) row.insert( std::make_pair(idx, val) );
	else it->second += val;
}

inline void Add2RowMat(std::vector<CSRRow>& ma, const std::vector<CSRRow>& mb, float ca, float cb){
	int nrow = ma.size();
	if( nrow != (int)mb.size() ){
		fprintf(stderr, "error matrix!");
		return;
	}
	for(int i=0; i<nrow; i++){
		for(CSRRow::iterator ait=ma[i].begin(); ait!=ma[i].end(); ++ait) ait->second *= ca;
		for(CSRRow::const_iterator bit=mb[i].begin(); bit!=mb[i].end(); ++bit){
			AddCoeff2CSRRow(ma[i], bit->first, bit->second*cb);
		}
	}
}



inline bool LeastSquareSolve(const CSRMatrix<mklReal> &matA, mklReal *const b, mklReal *const x, int nRhs=1)
{
	LeastSquareSparseSolver solver;
	if(solver.init(matA))	return solver.solve(b, x, nRhs);

	return false;
}

inline bool MklSolveSparseSystem(CSRMatrix<mklReal> &matA, mklReal *const b, mklReal *const x, int nRhs=1)
{
	SparseSolver solver;
	if( solver.init(&matA) ) return solver.solve(b, x, nRhs);

	//if(solver.init())	return solver.solve(b, x, nRhs);
	return false;
}
}

_NUMC_END

#endif		//#def NUMC_MKL_ADDON_H
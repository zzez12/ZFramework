#include "mkl_addonNew.h"
#include <algorithm>
#include <ctime>
#include <cstdio>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

_NUMC_BEGIN


template<>
void numcNew::CSRMatrix<double>::MultiVect(double* in, double *out) const
{
	if( !onebase() ) {
		fprintf(stderr, "\nzero-based sparse matrix, change to one-based with changeBase(true) first!");
		return;
	}

	int n = mNCol;
	double *a = const_cast<double*>(&mAv.front());
	int *ai = const_cast<int*>(&mAi.front());
	int *aj = const_cast<int*>(&mAj.front());
	if( issymm() )	mkl_dcsrsymv("U", &n, a, ai, aj, in, out);
	else mkl_dcsrgemv("N", &n, a, ai, aj, in, out);
}

template<>
void numcNew::CSRMatrix<float>::MultiVect(float* in, float *out) const
{
	int n = mNCol;
	float *a = const_cast<float*>(&mAv.front());
	int *ai = const_cast<int*>(&mAi.front());
	int *aj = const_cast<int*>(&mAj.front());
	if( issymm() ){
		fprintf(stderr, "\nnot supported!");
		//mkl_dcsrsymv("U", &n, a, ai, aj, in, out);
	}
	else mkl_scsrgemv("N", &n, a, ai, aj, in, out);
}



//template<>
//void CSRMatrix<float>::MultiVect(float* in, float *out) const
//{
//	if( !onebase() ) {
//		fprintf(stderr, "\nzero-based sparse matrix, change to one-based with changeBase(true) first!");
//		return;
//	}
//
//	int n = mNCol;
//	float *a = const_cast<float*>(&mAv.front());
//	int *ai = const_cast<int*>(&mAi.front());
//	int *aj = const_cast<int*>(&mAj.front());
//	if( issymm() ){
//		fprintf(stderr, "\nnot supported!");
//		//mkl_dcsrsymv("U", &n, a, ai, aj, in, out);
//	}
//	else mkl_scsrgemv("N", &n, a, ai, aj, in, out);
//}



// Multiplies matA by x and stores the result in b. Assumes all memory has been allocated
// and the sizes match; assumes matA is not symmetric!!
// void MulNonSymmMatrixVector(const CSRMatrix<mklReal> &matA, const mklReal *const x, mklReal *const b)
// {
//
// }


void numcNew::SparseSolver::clear()
{
	if( mIsInitiated )	{
		//////////////////////////////////////////////////////////////////////////
		// .. Termination and release of memory
		//////////////////////////////////////////////////////////////////////////
		pPhase = -1; /* Release internal memory. */
		pError = 0;
		PARDISO (pPt, &pMaxfct, &pMnum, &pMtype, &pPhase, &pN, &pDdum, pIa, pJa,
				 &pIdum, &pRhs, pIparm, &pMsglvl, &pDdum, &pDdum, &pError);

		//MKL_FreeBuffers();
	}

	std::fill_n(pPt, 64, (void*)NULL);
	std::fill_n(pIparm, 64, 0);
	pA = NULL;
	pIa = NULL;
	pJa = NULL;
	pRhs = 1;
	mIsInitiated = false;
}

bool numcNew::SparseSolver::init(const CSRMatrix<mklReal> *pMatA)
{
	if(mIsInitiated)	clear();

	if(pMatA)	mMatA = *pMatA;

	if( mMatA.empty() ) return false;

	mMatA.ChangeBase(true);
	pN = mMatA.mAi.size() - 1;
	const int nnz = mMatA.mAi.back() - 1;
	if ( (int)mMatA.mAj.size()!=nnz || (int)mMatA.mAv.size()!=nnz)
		return false;

	pA = &mMatA.mAv.front();
	pIa = &mMatA.mAi.front();
	pJa = &mMatA.mAj.front();

	pMtype = mMatA.mMtype;	// Real unsymmetric matrix
	pIparm[0] = 0;			// No solver default					// revised by jie @ 14/05/2007
// 	pIparm[1] = 2;			// Fill-in reordering from METIS */
// 	pIparm[2] = 1;			// omp_get_max_threads();	/* Numbers of processors, value of OMP_NUM_THREADS */
// 	pIparm[7] = 2;			// Max numbers of iterative refinement steps
// 	pIparm[9] = 13;		// Perturb the pivot elements with 1E-13
// 	pIparm[10] = 1;		// Use nonsymmetric permutation and scaling MPS
// 	pIparm[17] = -1;		// Output: Number of nonzeros in the factor LU
// 	pIparm[18] = -1;		// Output: Mflops for LU factorization
// 	pIparm[19] = 0;		// Output: Numbers of CG Iterations
	pMaxfct = 1;			// Maximum number of numerical factorizations
	pMnum = 1;				// Which factorization to use
	//pMsglvl = 0;			// Print statistical information in file		// changed to static member
	pError = 0;			// Initialize error flag


	//////////////////////////////////////////////////////////////////////////
	// .. Reordering and Symbolic Factorization. This step also allocates
	// all memory that is necessary for the factorization. */
	//////////////////////////////////////////////////////////////////////////
	pPhase = 11;
	PARDISO (pPt, &pMaxfct, &pMnum, &pMtype, &pPhase, &pN, pA, pIa, pJa,
	         &pIdum, &pRhs, pIparm, &pMsglvl, &pDdum, &pDdum, &pError);

	if (pError != 0) {
		printf("\nERROR during symbolic factorization: %d", pError);
		return false;
	}
	//////////////////////////////////////////////////////////////////////////
	// .. Numerical factorization
	//////////////////////////////////////////////////////////////////////////
	pPhase = 22;
	PARDISO (pPt, &pMaxfct, &pMnum, &pMtype, &pPhase, &pN, pA, pIa, pJa,
	         &pIdum, &pRhs, pIparm, &pMsglvl, &pDdum, &pDdum, &pError);
	// !!!!!!!!! don't use random value for pRhs

	if (pError != 0) {
		printf("\nERROR during numerical factorization: %d", pError);
		return false;
	}

	mIsInitiated = true;
	return true;
}


bool numcNew::SparseSolver::solve(const double *b, double *x, int nRhs, int nItr, int nPrecDigit)
{
	//////////////////////////////////////////////////////////////////////////
	// .. Back substitution and iterative refinement
	//////////////////////////////////////////////////////////////////////////
	pPhase = (nItr<0)?23:33;		// nItr < 0 for refactoring when no convergence

 	pIparm[3] = (nPrecDigit<0)?0:( nPrecDigit*10+(mMatA.issymm()?2:1) );		//preconditioned CGS.
 	pIparm[7] = abs(nItr);			// nItr: iteration numbers of substitution

	pError = 0;
	pRhs = nRhs;
	pIparm[26] = 1;
	PARDISO (pPt, &pMaxfct, &pMnum, &pMtype, &pPhase, &pN, pA, pIa, pJa,
		&pIdum, &pRhs, pIparm, &pMsglvl, (double*)b, x, &pError);
	if (pError != 0) {
		printf("\nERROR during solution: %d", pError);
		return false;
	}

	//	std::vector<mklReal> y(pN, 0);
	//	if(mMatA.issymm())	mkl_dcsrsymv("U", &pN, pA, pIa, pJa, x, &y.front());
	//	else mkl_dcsrgemv("N", &pN, pA, pIa, pJa, x, &y.front());
	//for(int i=0; i<pN; i++){ y[i] -= b[i];	}

	//printf("iteration number = %d\n", pIparm[6]);
	return true;
}



//////////////////////////////////////////////////////////////////////////
// Least Square Solver
//////////////////////////////////////////////////////////////////////////
bool numcNew::LeastSquareSparseSolver::init(const CSRMatrix<mklReal> &matA)
{
	CSRMatrix<mklReal> &matATA = mSparseSolver.getMatA();
	CSRMatrix<mklReal> &matAT = mMatAT;

//	printf("\n\n=====\tLeast square solver initiating ...");
	clock_t t1 = clock();

	/// matrix transposing
//	printf("\n\tMatrix transposing ...");
	clock_t t0 = clock();
	CSRMatrixTranspose(matA, matAT);
//	printf(" finished in %d ms", clock()-t0);

	/// matrix multiplication
//	printf("\n\tMatrix multiplication ...");
	t0 = clock();
	Mul2MatricesSymmResult(matAT, matA, matATA);			//ATA = A' x A
//	printf(" finished in %d ms", clock()-t0);
	mSparseSolver.init();

	/// updating index for matrix A': for rightside updating A'B
//	printf("\n\tSparse matrix index updating for MKL ...");
	t0 = clock();
	const int m = matA.nRow();
	const int n = matA.nCol();
	matAT.ChangeBase(true);
//	printf(" finished in %d ms", clock()-t0);

	return true;
}


int numcNew::SparseSolver::pMsglvl = 0;

_NUMC_END
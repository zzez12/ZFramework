#include "mkl_addon.h"
#include <algorithm>
#include <ctime>
#include <cstdio>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

bool numc::LeastSquareSolver(const CSRMatrix<mklReal> &matA, mklReal *const b, mklReal *const x, int nRhs)
{
	CLeastSquareSpareSolver solver;
	if(solver.init(matA)){
		return solver.solve(b, x, nRhs);
	}
	else
		return false;
}

bool numc::MklSolveSparseSystem(CSRMatrix<mklReal> &matA, mklReal *const b, mklReal *const x, int nRhs)
{
	CSparseSolver solver(&matA);

	if(solver.init())
		return solver.solve(b, x, nRhs);
	
	return false;

#if 0
	// nRhs: Number of right hand sides.
	matA.ChangeBase(true);
	vector<mklReal> &a = matA.v;
	vector<int> &ia = matA.i;
	vector<int> &ja = matA.j;

	int n = ia.size() - 1;
	const int nnz = ia.back() - 1;
	if (ja.size()!=nnz || a.size()!=nnz)
		return false;

	int mtype = matA.symmetric?(-2):11;				// Real unsymmetric matrix
	std::vector<void*> pt(64, NULL);		// Internal solver memory pointer
	std::vector<int> iparm(64, 0);			// Pardiso control parameters
	int maxfct, mnum, phase, error, msglvl = 0;
	/* Auxiliary variables. */
	double ddum;			// Double dummy
	int idum;				// Integer dummy
	iparm[0] = 0;			// No solver default					// revised by jie @ 14/05/2007
// 	iparm[1] = 2;			// Fill-in reordering from METIS */
// 	iparm[2] = 1;			// omp_get_max_threads();	/* Numbers of processors, value of OMP_NUM_THREADS */
// 	iparm[7] = 2;			// Max numbers of iterative refinement steps
// 	iparm[9] = 13;			// Perturb the pivot elements with 1E-13
// 	iparm[10] = 1;			// Use nonsymmetric permutation and scaling MPS
// 	iparm[17] = -1;			// Output: Number of nonzeros in the factor LU
// 	iparm[18] = -1;			// Output: Mflops for LU factorization
// 	iparm[19] = 0;			// Output: Numbers of CG Iterations
	maxfct = 1;				// Maximum number of numerical factorizations
	mnum = 1;				// Which factorization to use
	msglvl = 0;				// Print statistical information in file
	error = 0;				// Initialize error flag


	//////////////////////////////////////////////////////////////////////////
	// .. Reordering and Symbolic Factorization. This step also allocates
	// all memory that is necessary for the factorization. */
	//////////////////////////////////////////////////////////////////////////
	phase = 11;
	PARDISO (&pt.front(), &maxfct, &mnum, &mtype, &phase, &n, &a.front(), &ia.front(), &ja.front(),
	         &idum, &nRhs, &iparm.front(), &msglvl, &ddum, &ddum, &error);

	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		return false;
	}
	//////////////////////////////////////////////////////////////////////////
	// .. Numerical factorization
	//////////////////////////////////////////////////////////////////////////
	phase = 22;
	PARDISO (&pt.front(), &maxfct, &mnum, &mtype, &phase, &n, &a.front(), &ia.front(), &ja.front(),
	         &idum, &nRhs, &iparm.front(), &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		return false;
	}

	//////////////////////////////////////////////////////////////////////////
	// .. Back substitution and iterative refinement
	//////////////////////////////////////////////////////////////////////////
	phase = 33;
	PARDISO (&pt.front(), &maxfct, &mnum, &mtype, &phase, &n, &a.front(), &ia.front(), &ja.front(),
	         &idum, &nRhs, &iparm.front(), &msglvl, b, x, &error);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		return false;
	}

	//////////////////////////////////////////////////////////////////////////
	// .. Termination and release of memory
	//////////////////////////////////////////////////////////////////////////
	phase = -1; /* Release internal memory. */
	PARDISO (&pt.front(), &maxfct, &mnum, &mtype, &phase, &n, &ddum, &ia.front(), &ja.front(),
	         &idum, &nRhs, &iparm.front(), &msglvl, &ddum, &ddum, &error);

	//vector<mklReal> y(n, 0);
	//if(matA.symmetric){
	//	mkl_dcsrsymv("U", &n, &a.front(), &ia.front(), &ja.front(), x, &y.front());
	//}
	//mkl_dcsrgemv("N", &n, &a.front(), &ia.front(), &ja.front(), x, &y.front());

	return true;
#endif
}



// Multiplies matA by x and stores the result in b. Assumes all memory has been allocated
// and the sizes match; assumes matA is not symmetric!!
// void MulNonSymmMatrixVector(const CSRMatrix<mklReal> &matA, const mklReal *const x, mklReal *const b)
// {
//
// }


numc::CSparseSolver::CSparseSolver(const CSRMatrix<mklReal> *pMatA)
:m_ppt(64, NULL),m_piparm(64, 0),m_pa(NULL),m_pia(NULL),m_pja(NULL),m_pRhs(1)
{
	if(pMatA)
		m_matA = *pMatA;
}

bool numc::CSparseSolver::init(const CSRMatrix<mklReal> *pMatA)
{
//	printf("liguo1\n");
	if(pMatA){
		m_matA = *pMatA;
	}

	m_matA.ChangeBase(true);
	m_pn = m_matA.i.size() - 1;
	const int nnz = m_matA.i.back() - 1;
	if (m_matA.j.size()!=nnz || m_matA.v.size()!=nnz)
		return false;

	m_pa = &m_matA.v.front();
	m_pia = &m_matA.i.front();
	m_pja = &m_matA.j.front();

	m_pmtype = m_matA.symmetric?(-2):11;		// Real unsymmetric matrix
	m_pmsglvl = 0;
	m_piparm[0] = 0;			// No solver default					// revised by jie @ 14/05/2007
// 	m_piparm[1] = 2;			// Fill-in reordering from METIS */
// 	m_piparm[2] = 1;			// omp_get_max_threads();	/* Numbers of processors, value of OMP_NUM_THREADS */
// 	m_piparm[7] = 2;			// Max numbers of iterative refinement steps
// 	m_piparm[9] = 13;			// Perturb the pivot elements with 1E-13
// 	m_piparm[10] = 1;			// Use nonsymmetric permutation and scaling MPS
// 	m_piparm[17] = -1;			// Output: Number of nonzeros in the factor LU
// 	m_piparm[18] = -1;			// Output: Mflops for LU factorization
// 	m_piparm[19] = 0;			// Output: Numbers of CG Iterations
	m_pmaxfct = 1;				// Maximum number of numerical factorizations
	m_pmnum = 1;				// Which factorization to use
	m_pmsglvl = 0;				// Print statistical information in file
	m_perror = 0;				// Initialize error flag

	//////////////////////////////////////////////////////////////////////////
	// .. Reordering and Symbolic Factorization. This step also allocates
	// all memory that is necessary for the factorization. */
	//////////////////////////////////////////////////////////////////////////
	m_pphase = 11;
	PARDISO (&m_ppt.front(), &m_pmaxfct, &m_pmnum, &m_pmtype, &m_pphase, &m_pn, m_pa, m_pia, m_pja,
	         &m_pidum, &m_pRhs, &m_piparm.front(), &m_pmsglvl, &m_pddum, &m_pddum, &m_perror);

	if (m_perror != 0) {
		printf("\nERROR during symbolic factorization: %d", m_perror);
		return false;
	}
	//////////////////////////////////////////////////////////////////////////
	// .. Numerical factorization
	//////////////////////////////////////////////////////////////////////////
	m_pphase = 22;
	PARDISO (&m_ppt.front(), &m_pmaxfct, &m_pmnum, &m_pmtype, &m_pphase, &m_pn, m_pa, m_pia, m_pja,
	         &m_pidum, &m_pRhs, &m_piparm.front(), &m_pmsglvl, &m_pddum, &m_pddum, &m_perror);
	// !!!!!!!!! don't use random value for m_pRhs

	if (m_perror != 0) {
		printf("\nERROR during numerical factorization: %d", m_perror);
		return false;
	}

	return true;
}

bool numc::CSparseSolver::solve(double *b, double *x, int nRhs)
{

	//////////////////////////////////////////////////////////////////////////
	// .. Back substitution and iterative refinement
	//////////////////////////////////////////////////////////////////////////
	m_pphase = 33;
	m_perror = 0;
	m_pRhs = nRhs;
	PARDISO (&m_ppt.front(), &m_pmaxfct, &m_pmnum, &m_pmtype, &m_pphase, &m_pn, m_pa, m_pia, m_pja,
	         &m_pidum, &m_pRhs, &m_piparm.front(), &m_pmsglvl, b, x, &m_perror);
	if (m_perror != 0) {
		printf("\nERROR during solution: %d", m_perror);
		return false;
	}

	//vector<mklReal> y(n, 0);
	//if(matA.symmetric){
	//	mkl_dcsrsymv("U", &n, a, ia, ja, x, &y.front());
	//}
	//mkl_dcsrgemv("N", &n, a, ia, ja,, x, &y.front());

	return true;
}

numc::CSparseSolver::~CSparseSolver()
{
	//////////////////////////////////////////////////////////////////////////
	// .. Termination and release of memory
	//////////////////////////////////////////////////////////////////////////
	m_pphase = -1; /* Release internal memory. */
	m_perror = 0;
	PARDISO (&m_ppt.front(), &m_pmaxfct, &m_pmnum, &m_pmtype, &m_pphase, &m_pn, &m_pddum, m_pia, m_pja,
	         &m_pidum, &m_pRhs, &m_piparm.front(), &m_pmsglvl, &m_pddum, &m_pddum, &m_perror);
}


bool numc::CLeastSquareSpareSolver::init(const CSRMatrix<mklReal> &matA)
{
	CSRMatrix<mklReal> &matATA = *m_sparseSolver.getMatA();
	CSRMatrix<mklReal> &matAT = m_matAT;

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
	t0 = clock();
//	printf("\nSparse Solver Init...");
	m_sparseSolver.init();
//	printf(" finished in %d ms", clock() - t0);

	/// updating index for matrix A': for rightside updating A'B
//	printf("\n\tSparse matrix index updating for MKL ...");
	t0 = clock();
	const int m = matA.m;
	const int n = matA.n;
	matAT.ChangeBase(true);
//	printf(" finished in %d ms\n", clock()-t0);

	return true;
}

/*
bool numc::CLeastSquareSpareSolver::solve(double *b, double *x, int nRhs)
{
	/// right side updating
	printf("\n\n\tRight side updating ...");
	clock_t t0 = clock();
	int m = m_matAT.n;
	int n = m_matAT.m;
	vector<mklReal> ATb(m*nRhs, 0);									//ATb = A' x b
	for (int i=0; i<nRhs; i++) {
		mkl_dcsrgemv("N", &n, &m_matAT.v.front(), &m_matAT.i.front(), &m_matAT.j.front(), b+i*m, &ATb.front()+i*n);
	}
	printf(" finished in %d ms", clock()-t0);

	/// solver starts
	printf("\n\tSolver starts ...");
	t0 = clock();
	bool solveresult = m_sparseSolver.solve(&ATb.front(), x, nRhs);
	printf(" finished in %d ms", clock()-t0);

	return solveresult;
}
*/
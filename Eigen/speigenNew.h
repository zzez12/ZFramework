#ifndef _EIGEN_H_
#define _EIGEN_H_

#include "mkl_addonNew.h"

#include <ctime>
#include <cstdlib>

_NUMC_BEGIN

namespace numcNew
{

template<typename REAL>
class EigenProblem
{
public:
	EigenProblem(CSRMatrix<REAL> &m):_m(m) {
		if( !m.issymm() ){
			fprintf(stderr, "\nsymmetric matrix required!");
			throw std::exception( "\nsymmetric matrix required!" );
		}
		_m.ChangeBase(true); sprintf(_which, "SM"); }
	~EigenProblem() {};

	void solve(int nev, REAL*, REAL*eigvec=NULL);

	typedef void (* FP1)(int*, char*, int *, char*, int*, REAL*, REAL*, int *, REAL*, int*, int*, int*, REAL*, REAL*, int*, int*);
	typedef void (* FP2)(int*, char*, int*, REAL*, REAL*, int*, REAL*, char*,	int*,  char*, int*, REAL*, REAL*, int*, REAL*, int*, int*, int*, REAL*, REAL*, int*, int*);

	void setParam(int whichEig){
		switch(whichEig){
			case LM: sprintf(_which, "LM");	break;
			case SM: sprintf(_which, "SM"); break;
			case LA: sprintf(_which, "LA"); break;
			case SA: sprintf(_which, "SA");	break;
			case LI: sprintf(_which, "LI"); break;
			case SI: sprintf(_which, "SI"); break;
			default:
				fprintf(stderr, "wrong parameter of which eigen value/vector to calculate!!\n");
				fprintf(stderr, "setting to lm: largest magnitude!!\n");
				sprintf(_which, "LM");
		}
	}
public:
	enum {LM, SM, LA, SA, LI, SI};
	char _which[3];

private:
	static const FP1 saupd;
	static const FP2 seupd;
	//	static void (*const saupd)(int*, char*, int *, char*, int*, REAL*, REAL*, int *, REAL*, int*, int*, int*, REAL*, REAL*, int*, int*);

private:
	CSRMatrix<REAL> &_m;
	static REAL _tol;
	static int _maxiter;
};

#define USING_ARPACK_INTEL 0
#if USING_ARPACK_INTEL
#pragma comment (lib, "ifort/arpack.lib")
#pragma comment (lib, "ifort/ifconsol.lib")
#pragma comment (lib, "ifort/libifcoremt.lib")
#pragma comment (lib, "ifort/libifport.lib")
#pragma comment (lib, "ifort/libmmt.lib")
#pragma comment (lib, "ifort/libirc.lib")
#pragma comment (lib, "ifort/svml_disp.lib")

extern "C" {
	void DSAUPD(int *ido, char *bmat, int *n, char *which,	int *nev, double *tol, double *resid, int *ncv,
		double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info);

	void DSEUPD(int *rvec, char *All, int *select, double *d, double *v, int *ldv, double *sigma, 
		char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v2,
		int *ldv2, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *ierr);

	void SSAUPD(int *ido, char *bmat, int *n, char *which, int *nev, float *tol, float *resid, int *ncv,
		float *v, int *ldv, int *iparam, int *ipntr, float *workd, float *workl, int *lworkl, int *info);

	void SSEUPD(int *rvec, char *All, int *select, float *d, float *v, int *ldv, float *sigma, 
		char *bmat, int *n, char *which, int *nev, float *tol, float *resid, int *ncv, float *v2,
		int *ldv2, int *iparam, int *ipntr, float *workd, float *workl, int *lworkl, int *ierr);
}


template<>
const EigenProblem<double>::FP1 EigenProblem<double>::saupd = &DSAUPD;
template<>
const EigenProblem<double>::FP2 EigenProblem<double>::seupd = &DSEUPD;

template<>
const EigenProblem<float>::FP1 EigenProblem<float>::saupd = &SSAUPD;
template<>
const EigenProblem<float>::FP2 EigenProblem<float>::seupd = &SSEUPD;
#else
#pragma comment (lib, "arpack.lib")
extern "C" {
	void dsaupd_(int *ido, char *bmat, int *n, char *which,	int *nev, double *tol, double *resid, int *ncv,
		double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info);

	void dseupd_(int *rvec, char *All, int *select, double *d, double *v, int *ldv, double *sigma, 
		char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v2,
		int *ldv2, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *ierr);

	void ssaupd_(int *ido, char *bmat, int *n, char *which, int *nev, float *tol, float *resid, int *ncv,
		float *v, int *ldv, int *iparam, int *ipntr, float *workd, float *workl, int *lworkl, int *info);

	void sseupd_(int *rvec, char *All, int *select, float *d, float *v, int *ldv, float *sigma, 
		char *bmat, int *n, char *which, int *nev, float *tol, float *resid, int *ncv, float *v2,
		int *ldv2, int *iparam, int *ipntr, float *workd, float *workl, int *lworkl, int *ierr);
}

template<>
const EigenProblem<double>::FP1 EigenProblem<double>::saupd = &dsaupd_;
template<>
const EigenProblem<double>::FP2 EigenProblem<double>::seupd = &dseupd_;

template<>
const EigenProblem<float>::FP1 EigenProblem<float>::saupd = &ssaupd_;
template<>
const EigenProblem<float>::FP2 EigenProblem<float>::seupd = &sseupd_;

#endif




template<typename REAL>
void EigenProblem<REAL>::solve(int nev, REAL* eigval, REAL* eigvec)
{
	int n = _m.nCol();
	if( n != _m.nRow() ){
		fprintf(stderr, "Wrong dimension of matrix!");
		return;
	}

	int ido = 0;		/* Initialization of the reverse communication parameter. */
	char bmat[2] = "I"; /* Specifies that the right hand side matrix
						should be the identity matrix; this makes
						the problem a standard eigenvalue problem.
						Setting bmat = "G" would have us solve the
						problem Av = lBv (this would involve using
						some other programs from BLAS, however). */

	char *which = _which;
	/*char which[3] = "LM"; /* Ask for the nev eigenvalues of smallest
	magnitude.  The possible options are
	LM: largest magnitude
	SM: smallest magnitude 
	LA: largest real component
	SA: smallest real compoent
	LI: largest imaginary component
	SI: smallest imaginary component */

	REAL tol = _tol;		/* Sets the tolerance; tol<=0 specifies machine precision */

	std::vector<REAL> resid_(n);
	REAL *resid = &resid_.front();

	int ncv = 4*nev;	/* The largest number of basis vectors that will
						be used in the Implicitly Restarted Arnoldi
						Process.  Work per major iteration is
						proportional to N*NCV*NCV. */
	if (ncv>n) ncv = n;

	int ldv = n;

	std::vector<REAL> v_(ldv*ncv);
	REAL *v = &v_.front();

	int iparam[11]; /* An array used to pass information to the routines about their functional modes. */
	iparam[0] = 1;				// Specifies the shift strategy (1->exact)
	iparam[2] = _maxiter>0?_maxiter:3*n;			// Maximum number of iterations
	iparam[6] = 1;			 /* Sets the mode of dsaupd.
							 1 is exact shifting,
							 2 is user-supplied shifts,
							 3 is shift-invert mode,
							 4 is buckling mode,
							 5 is Cayley mode. */

	int ipntr[11]; /* Indicates the locations in the work array workd
				   where the input and output vectors in the
				   callback routine are located. */

	std::vector<REAL> workd_(3*n);
	REAL *workd = &workd_.front();

	std::vector<REAL> workl_( ncv*(ncv*8) );
	REAL *workl = &workl_.front();

	int lworkl = ncv*(ncv+8); /* Length of the workl array */

	int info = 0; /* Passes convergence information out of the iteration routine. */


	/* Here we enter the main loop where the calculations are
	performed.  The communication parameter ido tells us when
	the desired tolerance is reached, and at that point we exit
	and extract the solutions. */

//	fprintf(stdout, "\nsolving eigen problem with dim=%d nev=%d which=%s ...", n, nev, which );
	clock_t t=clock();


	SparseSolver solver;
	bool reverse = false;
	const double eps = .001;
	//if(0){
	if( which[0] == 'S' && which[1] == 'M' ){
		fprintf(stdout, "\nswitch to shift-invert mode for better performance!\n");
		//iparam[0] = 0;
		//iparam[6] = 3;
		which[0] = 'L';
		reverse = true;

		RowMatSym<REAL> rm(_m);
		for(int i = 0; i<n; i++)		rm[i][i] += eps;	/**m.getElementP( i, i )*/

		solver.getMatA() = rm;
		solver.getMatA().mMtype = CSRMatrix<double>::RealSymmIndef;
		solver.init();
	}


	int nit=0;
	do {
		saupd(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv,
			iparam, ipntr, workd, workl, &lworkl, &info);

		nit++;
		if( 0==nit%1000 ) fprintf(stderr, "..");
		if ((ido==1)||(ido==-1)){
			if( !reverse )	_m.MultiVect(workd+ipntr[0]-1, workd+ipntr[1]-1);
			else	solver.solve(workd+ipntr[0]-1, workd+ipntr[1]-1);
		}

	} while ((ido==1)||(ido==-1));

//	fprintf(stdout, " finished in %.3fs after %diters\n", (clock()-t)/1000., nit);

	/* From those results, the eigenvalues and vectors are extracted. */

	if(info<0) {
		fprintf(stderr, "\nError with dsaupd, info = %d", info);
		fprintf(stderr, "\nCheck documentation in dsaupd\n" );
	}else {
		std::vector<int> select_(ncv);
		std::vector<REAL> d_(ncv*2); /* This vector will return the eigenvalues from the second routine, dseupd. */

		int rvec = (NULL!=eigvec); /* Specifies that eigenvectors should not be calculated */
		int *select = &select_.front();
		REAL *d = &d_.front();
		REAL sigma;
		int ierr;

//		fprintf(stdout, "getting eigen values and eigen vectors ... ");
		t=clock();

		seupd(&rvec, "All", select, d, v, &ldv, &sigma, bmat,
			&n, which, &nev, &tol, resid, &ncv, v, &ldv,
			iparam, ipntr, workd, workl, &lworkl, &ierr);

//		fprintf(stdout, "finished in %.3fs\n", (clock()-t)/1000. );

		if (ierr!=0) {
			fprintf(stderr, "\nError with dseupd, info = %d", ierr);
			fprintf(stderr, "\nCheck the documentation of dseupd.\n");
		} else if (info==1) {
			fprintf(stderr, "\nMaximum number of iterations reached.\n");
		} else if (info==3) {
			fprintf(stderr, "\nNo shifts could be applied during implicit\n");
			fprintf(stderr, "\nArnoldi update, try increasing NCV.\n");
		}

		/* Before exiting, we copy the solution information over to
		the arrays of the calling program, then clean up the
		memory used by this routine.  For some reason, when I
		don't find the eigenvectors I need to reverse the order of
		the values. */

		if( !reverse ){
			std::copy(d, d+nev, eigval);
			if( NULL != eigvec )	std::copy(v, v+nev*n, eigvec);
		}
		else{
// 			for(int i=0; i<nev; i++) eigval[nev-i-1] = 1/d[i] - eps;
// 			if( NULL != eigvec ){
// 				for(int i=0; i<nev; i++)	std::copy(v+n*i, v+n*(i+1), eigvec+(nev-i-1)*n);
// 			}
			for (int i=0; i<nev; i++) eigval[i] = 1/d[i] - eps;
			std::copy(v, v+nev*n, eigvec);
		}
	}
}

double EigenProblem<double>::_tol = 0.0;
int EigenProblem<double>::_maxiter = 0;

}

_NUMC_END

#endif

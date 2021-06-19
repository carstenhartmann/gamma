/*  my nonrecursive acor generalized */
// clang -dynamiclib dcoeff.c -o dcoeff.so
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define TAUMAX  2  /*  Compute tau directly only if tau < TAUMAX.  
                       Otherwise compute tau using the pairwise sum series   */
#define WINMULT 5  /*  Compute autocovariances up to lag s = WINMULT*TAU     */
#define MAXLAG  TAUMAX*WINMULT   /*  The autocovariance array is double
                        C[MAXLAG+1] so that C[s] makes sense for s = MAXLAG. */
#define MINFAC   5  /*  Stop and print an error message if the array is
                                      shorter   than MINFAC * MAXLAG.        */

int ondiag(double *C0, double *D0, double X[], int L0, int *reducs){
	double mean = 0.;
	for ( int i = 0; i < L0; i++) mean += X[i];
	mean = mean / L0;
	for ( int i = 0; i <  L0; i++ ) X[i] -= mean;
  int k = 0;
  int L = L0;
	double D;
	while (true) {
		if ( L < MINFAC*MAXLAG ) {
			printf("Acor error 1: The autocorrelation time ");
			printf("is too long relative to the variance.\n");
			return 1; }
		double C[MAXLAG+1]; 
		for ( int s = 0; s <= MAXLAG; s++ )  C[s] = 0.;
		int iMax = L - MAXLAG;
		for ( int i = 0; i < iMax; i++ ) 
			for ( int s = 0; s <= MAXLAG; s++ )
				C[s] += X[i]*X[i+s];
		for ( int s = 0; s <= MAXLAG; s++ ) C[s] = C[s]/iMax;
		D = C[0];
		for ( int s = 1; s <= MAXLAG; s++ ) D += 2*C[s];
		double tau   = D / C[0]; 
		if (k == 0) *C0 = C[0];
		
		// modified to require tau >= 0
		if ( tau >= 0 && tau*WINMULT < MAXLAG ) break;
		k++;
		L = L/2;
	  int j1 = 0;
	  int j2 = 1;
	  for ( int i = 0; i < L; i++ ) {
			X[i] = X[j1] + X[j2];
			j1  += 2;
			j2  += 2; }
	} // end while
	*D0 = pow(0.25, k)*D*L0/(double)L;
  *reducs = k;
	return 0;
}

int offdiag(double *C0, double *D0, double X[], double Y[], int L0,
						int reducs){
	// assumed that X and Y have zero mean and are not overlapping in memory
	double mean = 0.;
	for ( int i = 0; i < L0; i++) mean += X[i];
	mean = mean / L0;
	for ( int i = 0; i <  L0; i++ ) X[i] -= mean;
	mean = 0.;
	for ( int i = 0; i < L0; i++) mean += Y[i];
	mean = mean / L0;
	for ( int i = 0; i <  L0; i++ ) Y[i] -= mean;
  int k = 0;
  int L = L0;
	double D;
	while (true) {
		double C[MAXLAG+1]; 
		for ( int s = 0; s <= MAXLAG; s++ )  C[s] = 0.;
		int iMax = L - MAXLAG;
		for ( int i = 0; i < iMax; i++ ) 
			for ( int s = 0; s <= MAXLAG; s++ )
				C[s] += X[i]*Y[i+s];
		for ( int s = 0; s <= MAXLAG; s++ ) C[s] = C[s]/iMax;
		D = C[0];
		for ( int s = 1; s <= MAXLAG; s++ ) D += 2*C[s];
		if (k == 0) *C0 = C[0];
		
		if ( k == reducs ) break;
		k++;
		L = L/2;
	  int j1 = 0;
	  int j2 = 1;
	  for ( int i = 0; i < L; i++ ) {
			X[i] = X[j1] + X[j2];
			Y[i] = Y[j1] + Y[j2];
			j1  += 2;
			j2  += 2; }
	} // end while
	*D0 = pow(0.25, k)*D*L0/(double)L;
	return 0;
}

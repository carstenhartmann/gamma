/*  my nonrecursive acor */
// clang -c acc.c
// clang -dynamiclib acc.o -o acc.so
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define TAUMAX  2
#define WINMULT 5
#define MAXLAG  TAUMAX*WINMULT
#define MINFAC   5

int acorR(double *sigma, double *tau, double X[], int L0);

int my_acor(double *mean, double *sigma, double *tau, double X[], int L0){
	*mean = 0.;
	for ( int i = 0; i < L0; i++) *mean += X[i];
	*mean = *mean / L0;
	for ( int i = 0; i <  L0; i++ ) X[i] -= *mean;
	return acorR(sigma, tau, X, L0);
}

int acorR(double *sigma, double *tau, double X[], int L0) {
	// assumed that X has zero mean
  double C0 = 0.;
  int iMax = L0 - MAXLAG;
	for ( int i = 0; i < iMax; i++ ) C0 += X[i]*X[i];
	C0 = C0/iMax;
	double D;
  int k = 0;
	int Lk = L0;
	while (true) {
		if ( Lk < MINFAC*MAXLAG ) {
			printf("Acor error 1: The autocorrelation time ");
			printf("is too long relative to the variance.\n");
			return 1; }
		double C[MAXLAG+1]; 
		for ( int s = 0; s <= MAXLAG; s++ )  C[s] = 0.;
		iMax = Lk - MAXLAG;
		for ( int i = 0; i < iMax; i++ ) 
			for ( int s = 0; s <= MAXLAG; s++ )
				C[s] += X[i]*X[i+s];
		for ( int s = 0; s <= MAXLAG; s++ ) C[s] = C[s]/iMax;
		D = C[0];
		for ( int s = 1; s <= MAXLAG; s++ ) D += 2*C[s];
		*tau   = D / C[0]; 
		
		if ( *tau*WINMULT < MAXLAG ) break;
		int Lkh = Lk/2;
	  int j1 = 0;
	  int j2 = 1;
	  for ( int i = 0; i < Lkh; i++ ) {
			X[i] = X[j1] + X[j2];
			j1  += 2;
			j2  += 2; }
		k++;
		Lk = Lkh;
	} // end while
	double D0 = pow(0.25, k)*D*L0/Lk;
	*sigma = sqrt(D0/L0);
  *tau = D0/C0;
	return 0;
}

// propagate.c
// clang -dynamiclib propagate.c -o propagate.so
#include <stdio.h>
#include <stdlib.h>
#include<stdbool.h>
#include <math.h>

int twoGauss(int N, double *q, double *p,
						 double gamma, double d, double dt, int seed){
	double twopi = 8.*atan(1.0);
	srand(seed);
	double expgt = exp(-gamma*dt);
	double sigma = sqrt(1. - expgt*expgt);
  double qn = q[0];
	double Fn = -qn + d*tanh(d*qn);
  double pn = p[0];
  double rn;
  double rnNext = NAN;
  for (int n = 0; n < N; n++){
		q[n] = qn; p[n] = pn;
		pn = pn + 0.5*dt*Fn;
		qn = qn + 0.5*dt*pn;
    rn = rnNext; rnNext = NAN;
		if (isnan(rn)){
			double u1 = ((double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
			double u2 = ((double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
			double sqrt2log = sqrt(-2.*log(u1));
			double twopiu2 = twopi*u2;
			rn = cos(twopiu2)*sqrt2log;
			rnNext = sin(twopiu2)*sqrt2log;}
		pn = expgt*pn + sigma*rn;
		qn = qn + 0.5*dt*pn;
		Fn = -qn + d*tanh(d*qn);
		pn = pn + 0.5*dt*Fn;
	}
	return 0;}

int threeGauss(int N, double *q, double *p,
						 double gamma, double d, double dt, int seed){
  // q = [x1, y1, x2, y2, ..., xN, yN], same for p
  if (N % 2 == 1) printf("Error!\n");
  int Nby2 = N/2;
	double twopi = 8.*atan(1.0);
	srand(seed);
	double expgt = exp(-gamma*dt);
	double sigma = sqrt(1. - expgt*expgt);
  double xn = q[0]; double yn = q[1];
  double pxn = p[0]; double pyn = p[1];
  // U(x, y) = 0.5*(x^2 + y^2) - log(c exp(dx) + 2exp(-dx/2)cosh(sqrt(3)dy/2))
  // where xi = exp(dx/2), eta = exp(sqrt(3)dy/2), and c = exp(-d^2/8)
  // Fx = -x + d(c xi^2 - (1/xi)(eta+1/eta))/(c xi^2 + (2/xi)(eta+1/eta))
  // Fy = -y + sqrt(3)d(1/xi)(eta+1/eta)/c (xi^2 + (2/xi)(eta+1/eta))
	double hd = 0.5*d, r3d = sqrt(3.)*d, hr3d = 0.5*r3d;
  double xi = exp(hd*xn); double eta = exp(hr3d*yn);
  double term = (eta + 1./eta)/xi;
	double denom = xi*xi + 2.*term;
	double Fxn = - xn + d*(xi*xi - term)/denom;
	double Fyn = -yn + r3d*(eta - 1./eta)/xi/denom;
  for (int n = 0; n < Nby2; n++){
	  int nx = 2*n; int ny = 2*n+1;
  	q[nx] = xn; q[ny] = yn;
	  p[nx] = pxn; p[ny] = pyn;
	  pxn = pxn + 0.5*dt*Fxn; pyn = pyn + 0.5*dt*Fyn;
		xn = xn + 0.5*dt*pxn; yn = yn + 0.5*dt*pyn;
		double u1 = ((double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
		double u2 = ((double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
		double sqrt2log = sqrt(-2.*log(u1));
		double twopiu2 = twopi*u2;
		double rxn = cos(twopiu2)*sqrt2log;
		double ryn = sin(twopiu2)*sqrt2log;
		pxn = expgt*pxn + sigma*rxn;
		pyn = expgt*pyn + sigma*ryn;
		xn = xn + 0.5*dt*pxn;	yn = yn + 0.5*dt*pyn;
		xi = exp(hd*xn); eta = exp(hr3d*yn);
		term = (eta + 1./eta)/xi;
		denom = xi*xi + 2.*term;
		Fxn = - xn + d*(xi*xi - term)/denom;
		Fyn = -yn + r3d*(eta - 1./eta)/xi/denom;
		pxn = pxn + 0.5*dt*Fxn; pyn = pyn + 0.5*dt*Fyn;
  }
	return 0;}

int LeMa13(int N, double *q, double *p,
						 double gamma, double dt, int seed){
	double twopi = 8.*atan(1.0);
	srand(seed);
	double expgt = exp(-gamma*dt);
	double sigma = sqrt(1. - expgt*expgt);
  double qn = q[0];
	double Fn = -qn*qn*qn - 5.*cos(1. + 5.*qn);
  double pn = p[0];
  double rn;
  double rnNext = NAN;
  for (int n = 0; n < N; n++){
		q[n] = qn; p[n] = pn;
		pn = pn + 0.5*dt*Fn;
		qn = qn + 0.5*dt*pn;
    rn = rnNext; rnNext = NAN;
		if (isnan(rn)){
			double u1 = ((double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
			double u2 = ((double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
			double sqrt2log = sqrt(-2.*log(u1));
			double twopiu2 = twopi*u2;
			rn = cos(twopiu2)*sqrt2log;
			rnNext = sin(twopiu2)*sqrt2log;}
		pn = expgt*pn + sigma*rn;
		qn = qn + 0.5*dt*pn;
		Fn = -qn*qn*qn - 5.*cos(1. + 5.*qn);
		pn = pn + 0.5*dt*Fn;
	}
	return 0;}

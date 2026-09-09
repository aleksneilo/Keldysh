#include "SFS.h" 
#include <complex>
#include <iostream>
#define pi 3.141592653589793
#define icom (complex<double>(0, 1.))

using namespace std;

double J(double E)
{
	return (tanh(E / 2 / T));
}

/*
void sigma(double hw, complex<double> *Del, complex<double> *Sigma)
{	
	double E;
	int Enum = 20000;

	complex<double> *F, *F1, *G, *G1;
	F = new complex<double>[N];
	F1 = new complex<double>[N];
	G = new complex<double>[N];
	G1 = new complex<double>[N];
	double *sigma1, *sigma2;
	sigma1 = new double[N];
	sigma2 = new double[N];
	for (int i = 0; i < N; i++)
	{
		sigma1[i] = 0;
		sigma2[i] = 0;
	}

	double K1, K2, K3;

	for (int ie = 0; ie < Enum; ie++)
	{
		E = -20 + ie*40. / Enum;

		DOS(E, G, F, Del);
		DOS(E+hw, G1, F1, Del);
		for (int i = 0; i < N; i++)
		{
			K1 = imag(F[i])*imag(F1[i]) + real(G[i])*real(G1[i]);
			K2 = real(F[i])*imag(F1[i]) - imag(G[i])*real(G1[i]);
			K3 = real(F1[i])*imag(F[i]) - real(G[i])*imag(G1[i]);
			sigma1[i] += (J(E + hw) - J(E))*K1 / 2. / hw * 40. / Enum;
			sigma2[i] += (J(E + hw)*K2 + J(E)*K3) / 2. / hw * 40. / Enum;
		}
	}
	for (int i = 0; i < N; i++)
		Sigma[i] = sigma1[i] + icom* sigma2[i];

	delete[] F;
	delete[] F1;
	delete[] G1;
	delete[] G;
	delete[] sigma1;
	delete[] sigma2;
}

complex<double> root_new(complex<double> x, int t)
{
	//if ((real(exp(icom*(arg(x))/2.))>0)&&(imag(exp(icom*(arg(x))/2.))<0))
	if (t == 0)
		return(sqrt(abs(x))*exp(icom*(arg(x)) / 2.));
	else
		return(sqrt(abs(x))*exp(icom*(arg(x) + pi) / 2.));
}

void Z(complex<double> *Z, complex<double> *Sigma, double hw )
{
	double X;
	for (int i = 0; i < N; i++)
	{
		Z[i] = (1. - icom)*root_new(hw / 1.76 / Sigma[i], 0);
		X = -imag(Z[i]);
		if (X<0)
		{
			Z[i] = (1. - icom)*root_new(hw / 1.76 / Sigma[i], 1);
		}
	}

} */

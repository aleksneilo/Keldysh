#include "SFS.h"
#include <complex>
#include <iostream>
#define pi 3.141592653589793
#define icom (complex<double>(0, 1.))

using namespace std;

void DOS(double E, complex<double>* G, complex<double>* G1, complex<double>* F, complex<double>* Del, double* q, double I)
{
	// double  Xi1;
	double  h, dGmax;
	complex<double> t1, t2, w, wm;
	double ro;
    int itermax;
	complex<double> *Fi, *Fi1;
	Fi = new complex<double>[N];
	Fi1 = new complex<double>[N];
	//G1 = new complex<double>[N];
	double delE = 1e-4;
        itermax=250;
        epsG=1e-5;
        delE = 1e-3;

	//h=L_S/(N_S-1);

	if(abs(E)<5)
    {   itermax=1000;
        epsG=5e-5;
        delE = 1e-7;
        for (int i = 0; i < N; i++)
        {
            G[i] = 1.+0.1*icom;
            G1[i] = 1.+0.1*icom;
        }
    }
	 for (int i = 1; i < N - 1; i++)
	if (E > abs(Del[i]))
	{
	G[i] = 1.;
	G1[i] = 1.;
	}
	else
	{
	G[i] = 0.001;
	G1[i] = 0.001;
	}

	w = -icom*E + delE;


	dGmax = 1;
	int iterG = 0;
	while ((dGmax > epsG)&&(iterG<itermax))
	{
		// calc Fi(w), Fi(-w)
		if(iterG>5000)
			iterG++;
		Prog(Fi, G, Del, w,q,I);
		Prog(Fi1, G1, Del, -w,q,I);
		GcalcDOS(G1, &dGmax, Fi1, Fi, Del, -w);
		GcalcDOS(G, &dGmax, Fi, Fi1, Del, w);
		iterG++;
		//cout <<iterG<<"  "<< G[N-1] << " " << dGmax << endl;
	}
	if(iterG>itermax-5) cout<<fixed<<E<<"  "<<iterG<<"  "<<dGmax<<endl;
	for (int i = 0; i < N; i++)
	{
		F[i] = Fi[i] * G[i] / get_wm(i, w);
	}

	delete[] Fi;
	delete[] Fi1;
	//delete[] G1;

}


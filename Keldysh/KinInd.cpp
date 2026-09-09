#include "SFS.h"
#include <complex>
#include <iostream>
#define pi 3.141592653589793
#define icom (complex<double>(0, 1.))

using namespace std;

double KinInd(complex<double>* G, complex<double>* Del, double* q, double I)
{
	// double  Xi1;
	int iterG;
	double w, h, dGmax, KI;
	complex<double> t1, t2, wm, Fx, F1x, Gx, G1x;
	double ro, D, E, Ksi;
	ro = 1; E = 0;
     //Ksi= get_ksi(iN);
	  //D=2*pi*Ksi*Ksi; // Diffusion coefficient for iN point
	complex<double> *Fi, *Fi1, *Ei, *F, *F1, *G1;
	double  *lambda2;
	Fi = new complex<double>[N];
	Fi1 = new complex<double>[N];
	F = new complex<double>[N];
	F1 = new complex<double>[N];
	G1 = new complex<double>[N];
	Ei = new complex<double>[N];
	lambda2 = new double[N];

	//D = 2 * pi * 1 * 1;
	for (int i = 1; i < N - 1; i++)
	{
		lambda2[i] = 0;
		// cout << Del[i] << " ";
	}


	for (int iw = 0; iw<w_obrez; iw++)
	{
		w = pi*T*(2 * iw + 1);
		//  wm=w+ icom*H;
        for (int i=0; i<N; i++)	    G[i]=1.+0.1*icom;
		dGmax = 1; iterG = 0;
		while ((dGmax > epsG)&&(iterG<1000))
		{
			// calc Fi(w), Fi(-w)
			for (int i = 0; i<N; i++)
				G1[i] = -conj(G[i]);
			Prog(Fi, G, Del, w,q,I);
			//Prog(Fi1, G1, Del, -w,q,I);
			 for (int i = 0; i<N; i++)
			Fi1[i] = conj(Fi[i]);
			// recalc G
			Gcalc(G, &dGmax, Fi, Fi1, Del, w,q);
			iterG++;
		}

		// E calculation

		for (int i = 1; i < N - 1; i++)
		{
			wm = get_wm(i, w) + real(w) / abs(w) * 2*pi * get_ksi(i) *get_ksi(i) *G[i]* q[i] * q[i] / 2.;
			F[i] = Fi[i] / sqrt(wm * wm + Fi[i] * conj(Fi1[i]));
			wm = get_wm(i, -w) + real(w) / abs(w) * 2 * pi * get_ksi(i) * get_ksi(i) * G[i] * q[i] * q[i] / 2.;
			F1[i] = Fi1[i] / sqrt(wm * wm + Fi1[i] * conj(Fi[i]));

			ro = Roi[Layer(i)];

			lambda2[i] = lambda2[i] + T / ro * real(F[i] * F[i]);
		}
		//cout << F[1]<<" " << F[1] * F[1] << " " << T / ro << " " << T / ro* real(F[1] * F[1]) << " " << lambda2[1] << endl;
	}
	//for (int i = 1; i < N - 1; i++) lambda2[i] = lambda2[i] + T / ro * real(F[i] * F[i]);

	KI = 0;

	for (int i = 1; i < N - 1; i++)
	{
		h = get_h(i);
		KI = KI + lambda2[i]*h;

		//if ((i - i / 10 * 10) == 0)
			//cout << Del[i] << " " << Fi[i] << " " << F[i] << " " << lambda2[i] << " " << h << " " << KI <<endl;
	}

	for (int i = 1; i < N - 1; i++)
	{	
		q[i] = I /KI;
		
	}
	//for (int i = 1; i < N - 1; i++) cout << i << "  " << T / ro * real(F[i] * F[i]) << "  " << lambda2[i] << "  " <<I<< "  " <<KI<< "  " << q[i] << endl;

	return(1/KI);

	delete[] Fi;
	delete[] Fi1;
	delete[] F;
	delete[] F1;
	delete[] G1;
	delete[] Ei;
	delete[] lambda2;

}//*/


double SuF(complex<double> *G, complex<double> *Del, double* q, double I)
{
	// double  Xi1;
	double w, h, dGmax, KI;
	complex<double> t1, t2, wm, Fx, F1x, Gx, G1x;
	double ro, D, E;
	ro = 1; E = 0;

	complex<double> *Fi, *Fi1, *Ei, *F, *F1, *G1;
	double  *lambda2;
	Fi = new complex<double>[N];
	Fi1 = new complex<double>[N];
	F = new complex<double>[N];
	F1 = new complex<double>[N];
	G1 = new complex<double>[N];
	Ei = new complex<double>[N];
	lambda2 = new double[N];

	D = 2 * pi * 1 * 1;
	for (int i = 1; i < N - 1; i++)
	{
		lambda2[i] = 0;
		// cout << Del[i] << " ";
	}




	for (int iw = 0; iw<w_obrez; iw++)
	{
		w = pi*T*(2 * iw + 1);
		//  wm=w+ icom*H;


		dGmax = 1;
		while (dGmax > epsG)
		{
			// calc Fi(w), Fi(-w)
			for (int i = 0; i<N; i++)
				G1[i] = -conj(G[i]);
			Prog(Fi, G, Del, w,q,I);
			Prog(Fi1, G1, Del, -w, q, I);
			// recalc G
			Gcalc(G, &dGmax, Fi, Fi1, Del, w,q);
		}

		// E calculation

		for (int i = 1; i < N - 1; i++)
		{
			wm = get_wm(i, w);
			F[i] = Fi[i] / sqrt(wm*wm + Fi[i] * conj(Fi1[i]));
			wm = get_wm(i, -w);
			F1[i] = Fi1[i] / sqrt(wm*wm + Fi1[i] * conj(Fi[i]));

			ro = Rbi[Layer(i)];

			lambda2[i] = lambda2[i] + real(F[i] * F[i]);



		}

		//cout << lambda2[1]<<endl;
	}



	return(lambda2[1]);

	delete[] Fi;
	delete[] Fi1;
	delete[] F;
	delete[] F1;
	delete[] G1;
	delete[] Ei;
	delete[] lambda2;
}

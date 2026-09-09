 #include "SFS.h"
 #include <complex>
  #include <iostream>
 #define pi 3.141592653589793
 #define icom (complex<double>(0, 1.))

 using namespace std;


////////// getABC: coefficients of difference method for Usadel equation ///////
//         We use difference scheme. At the nod n equation is: a*F(n-1) - b*F(n) + c* F(n+1) + d = 0
//         Here we set all coefficints a,b,c,d for boundary broblem.
//
////INPUT: iN - number of the nod, for which coefficient are required
//         array G - Normal Green Function from previous step
//         array Del - Pair potential
//         w - Matsubara frequence for certain step
//
// OUTPUT: a, b, c, d - coefficients of equation

void getABC(complex<double> *a, complex<double> *b, complex<double> *c, complex<double> *d, int iN, complex<double> *G, complex<double> *Del, complex<double> w, double* q, double I)
 {
   //   double  Xi1;
      double D, h, Ksi;
      double h_F, h_SL, h_SR, h_SM,alphaTT;
      complex<double> wm, wm1, wm2, G0, dG2_S,  dG2_F;


	  Ksi= get_ksi(iN);
	  D=2.*pi*Ksi*Ksi; // Diffusion coefficient for iN point
	  G0 = w / sqrt(w*w + Del0*Del0);
      wm = get_wm(iN, w) +real(w)/abs(w)*D*q[iN]*G[iN]* q[iN] /2.;
        h= get_h(iN);

	if (MODE == 1)
		 if (iN==0)
		 {  *a=0;
			*b=-1;
			*c=0;
			*d=1.46*exp(icom*0.);
		 }

    /*if (iN == 0)
	 {
		 *a = 0;
		 *b = Gb_SF*Ksi + G0 / G[iN] * h;
		 *c = Gb_SF*Ksi;
		 *d = -G0 / G[iN] * h *wm / w*Del0*exp(icom*Xi1);
	 }*/


	 // Free Boundary
	if (MODE == 0)
	 if (iN == 0)
	 {
		 *a = 0;
		 *b = 1;
		 *c = 1;
		 *d = 0;/*/
		 *a=0;
		 *b=-1;
		 *c=0;
		 *d=Del0;//*/;
	 }
     // F0-F8 - layers
     if ((iN>0)&&(iN<N-1))//((NUM_Tech)*N_Mid-1)))
     {  dG2_S=(G[iN+1]*G[iN+1]-G[iN-1]*G[iN-1])*0.25/G[iN]/G[iN];

        *a=-dG2_S+1.;
        *b=(2./D*(wm)/G[iN]*h*h+2.);
        *c=dG2_S+1.;
        *d=-2./D*(wm)/G[iN]*get_type(iN)*Del[iN]*h*h;
     }
     //for (int ki = 1; ki < (NUM_Tech - 1); ki++)
     //{
         if ((iN == N_S - 1) || (iN == N_S + N_F - 1) || (iN == N_S + N_F + N_S1 - 1) || (iN == N_S + N_F + N_S1 + N_F1 - 1))//||(iN==N_S + N_F + N_S1 + N_F1 + N_N-1))//(ki*N_Mid-1))
         {
             int ki = Layer(iN);
             wm1 = get_wm(iN, w);
             wm2 = get_wm(iN + 1, w);
             *a = -Rbi[ki] / Roi[ki];
             *b = -(Rbi[ki] / Roi[ki] + G[iN + 1] / G[iN] * h);
             //*a=-Gb_FS*Ksi;
             //*b=-(Gb_FS*Ksi + G[iN+1]/G[iN]*h);
             *c = -G[iN + 1] / G[iN] * h * wm1 / wm2;
             *d = 0;
         }
         if ((iN == N_S) || (iN == N_S + N_F) || (iN == N_S + N_F + N_S1) || (iN == N_S + N_F + N_S1 + N_F1))//||(iN==N_S + N_F + N_S1 + N_F1 + N_N))
         {
             int ki = Layer(iN);
             wm1 = get_wm(iN - 1, w);
             wm2 = get_wm(iN, w);
             *a = G[iN - 1] / G[iN] * h * wm2 / wm1;
             *b = Rbi[ki] / Roi[ki] + G[iN - 1] / G[iN] * h;
             *c = Rbi[ki] / Roi[ki];
             //*b=Gb_SF*Ksi + G[iN-1]/G[iN]*h;
             //*c=Gb_SF*Ksi;
             *d = 0;
             //}//*/
         }
     //}
         /*/if (iN == (N_S + N_F + N_S1 + N_F1 + N_S2 - 1))//
         //if(N - N_N - 1)
         {
             *a = -R0AN;
             *b = -(R0AN + G[iN + 1] / G[iN] * h);
             *c = -G[iN + 1] / G[iN] * h;
             *d = 0;
         }

         if (iN == (N_S + N_F + N_S1 + N_F1 + N_S2))//
         //if(iN == N - N_N)
         {
             *a = G[iN - 1] / G[iN] * h;
             *b = R0AN / ro_N + G[iN - 1] / G[iN] * h;
             *c = R0AN / ro_N;
             *d = 0;
         }/*/
     
     /*if (MODE == 1)
		 if (iN==(N-1))
		 {  *a=0;
			*b=-1;
			*c=0;
			*d=Del0*exp(icom*Xi2*pi);
		 }//*/

	  //if (iN==(N-1))
	  //{
		 // *a = -Gb_FS*Ksi;
		 // *b = -(Gb_FS*Ksi + G0 / G[iN] * h);
		 // *c = 0;
		 // *d = G0 / G[iN] * h *wm / w*Del0*exp(icom*Xi2);
	  //}

		 // Free Boundary
	 if (MODE==0)
		 if (iN==(N-1))
		 {  *a=1;
			*b=1;
			*c=0;
			*d=0;
		 }
}




////////// Prog: 3 diagonal sweep method to solve linearized Usadel equation ///////
////INPUT: array G - Normal Green Function
//         array Del - Pair potential
//         w - Matsubara frequence for certain step
// OUTPUT: array Fi - Parametrized Anomalous Green Function

void Prog( complex<double> *Fi, complex<double> *G, complex<double> *Del, complex<double> w, double* q, double I)
 {
        complex<double> a, b, c, d, *p, *qq ;

        p=new complex<double>[N+1];
        qq=new complex<double>[N+1];

        p[0]=0;
        qq[0]=0;

        // Forward sweep

        for (int i = 0; i<N; i++)
        {   // Set difference method coefficients for sweep.
            getABC(&a, &b, &c, &d, i, G, Del, w, q,I);

            p[i+1]=c/(b-a*p[i]);
            qq[i+1]=(a*qq[i]-d)/(b-a*p[i]);

        }

         // Backward sweep

        Fi[N-1]=qq[N];
        for (int i=N-2; i>-1; i--)
        {
            Fi[i]= p[i+1]* Fi[i+1]+qq[i+1];
        }

        delete [] p;
        delete [] qq;

}

////////// Gcalc: calcualtion of Normal Green function ///////
////INPUT: array G - Normal Green Function from previous step
//         arrays Fi1 and Fi - Anomalous Green Functions for positive and negative Matsubara
//         array Del - Pair potential
//         w - Matsubara frequence for certain step
// OUTPUT: array G - Resulting Normal Green Function
//         dGmax - maximal mismatch between old and new G.

void Gcalc( complex<double> *G, double *dGmax, complex<double> *Fi1, complex<double> *Fi2, complex<double> *Del, complex<double> w, double* q)
{   complex<double> wm, Gbuf;
    *dGmax=0;
    //wm=w+ icom*H;
    for (int i=0; i<N; i++)
    {   Gbuf= G[i];
        wm= get_wm(i, w) + real(w) / abs(w) * 2. * pi * get_ksi(i) * get_ksi(i) *G[i]* q[i] * q[i] / 2.;   // function which determines value of modified matsubara from coordinate
        G[i] = 1. * (wm) / sqrt(wm * wm + (Fi1[i] * conj(Fi2[i])));// +*Gbuf;
        if (*dGmax< abs(G[i]-Gbuf))
           *dGmax = abs(G[i]-Gbuf);
        // cout<< i<<" "<<*dGmax<<endl;
    }
    //cout<< w<<" "<<*dGmax<<G[N-100]<<endl;
}


void GcalcDOS(complex<double> *G, double *dGmax, complex<double> *Fi1, complex<double> *Fi2, complex<double> *Del, complex<double> w)
{
	complex<double> wm, Gbuf;
	*dGmax = 0;
	//wm=w+ icom*H;
	for (int i = 0; i<N; i++)
	{
		Gbuf = G[i];
		wm = get_wm(i, w) ;   // function which determines value of modified matsubara from coordinate
		G[i] = (wm) / sqrt(wm*wm + (Fi1[i] * conj(Fi2[i])));
		if (real(G[i]) < 0)
			G[i] = -G[i];
		G[i] = 0.1*G[i] + 0.9* Gbuf;

		if (*dGmax< abs(G[i] - Gbuf))
			*dGmax = abs(G[i] - Gbuf);
		//cout<< G[i]<<" " << Gbuf<<" "<<*dGmax<<" ";
	}

}

#include "SFS.h"
#include <complex>
#define pi 3.141592653589793
#define icom (complex<double>(0, 1.))

// S-I-S-F-S structure

int Layer(int iN)
{   int k;//=(iN/N_Mid);
    if(iN<N_S) k=0;
	/*/if ((iN > N_S - 1) && (iN < N_S + N_F)) k = 1;
	if ((iN > N_S + N_F - 1) && (iN < N_S + N_F + N_S1)) k = 2;
    if((iN>N_S+N_F+N_S1-1)&&(iN<N_S+N_F+N_S1+N_F1)) k=3;
	if ((iN > N_S + N_F + N_S1 + N_F1 - 1) && (iN < N)) k = 4;// _S + N_F + N_S1 + N_F1 + N_S2)) k = 4;
     //if((iN>(N_S + N_F + N_S1 + N_F1 + N_S2-1))&&(iN<N)) k=5;//*/
      /*if((iN>N_S+N_F+2*N_S1+N_F1-1)&&(iN<N_S+2*N_F+2*N_S1+N_F1)) k=5;
      if((iN>N_S+2*N_F+2*N_S1+N_F1-1)&&(iN<N_S+2*N_F+3*N_S1+N_F1)) k=6;
      if((iN>N_S+2*N_F+3*N_S1+N_F1-1)&&(iN<N_S+2*N_F+3*N_S1+2*N_F1)) k=7;//*/
    //if((iN>N_S+N_F+N_F1+N_S1-1)&&(iN<N)) k=4;
    //if((iN>N_S+2*N_F+3*N_S1+2*N_F1-1)&&(iN<N)) k=8;
    return k;
}

complex<double> get_wm(int iN, complex<double> w)
{
	complex<double> wm;
	//if(iN<N_S) il=0;
	//else if(iN<N_S+N_Mid) il=1;//il=(iN/N_Mid);
         //else il=2;
    int il; il=Layer(iN);
	wm=w+icom*(Hi[il]);
	return wm;
}

/*double get_ro(int iN)
{   if ((iN>=(N_SL+N_SM))&&(iN<=(N_SL+N_SM+N_F-1)))
        return ro_F
    else
        return ro_S;
}*/


double get_h(int iN)
{   double h;
	int il;
    h=1./(1.*N_Mid);
	/*il=(iN/N_Mid);
	h=Li[il]/(N_Mid-1);
	if(iN<N_S) h=Li[0]/(N_S-1);
	else h=Li[1]/(N_Mid-1);//*/
    return h;
}

double get_h_out(int iN)
{   double h;
	int il;
    h=0;
	il=(iN/N_Mid);
	h=Li[il]/(N_Mid-1);
	for (int k=0; k<NUM_Tech; k++)
	{
		if (iN==k*N_Mid-1)
			h=0;
	}
    return h;
}

double get_h_E(int iN)
{
	double h;
	int il;
	h = 0;
	il = (iN / N_Mid);
	h = Li[il] / (N_Mid - 1);
	for (int k = 0; k<NUM_Tech; k++)
	{
		if (iN == k*N_Mid - 1)
			h = 0;
	}
	return h;
}

double get_HE(int iN)
{   double HE;
	int il;
    HE=0;
    il=Layer(iN);
	//il=(iN/N_Mid);
	HE=Hi[il];

    return HE;
}

double get_ksi(int iN)
{   double ksi;
	int il;
	//if(iN<N_S) il=0;
	//else if(iN<N_S+N_Mid) il=1;//il=(iN/N_Mid);
         //else il=2;
    il=Layer(iN);
	ksi=Ksii[il];
	//if (ksi>Li[il])	ksi=Li[il];
	return ksi;
}

double get_type(int iN)
{   double s;
	int il;
	//if(iN<N_S) il=0;
	//else if(iN<N_S+N_Mid) il=1;//il=(iN/N_Mid);
         //else il=2;
    il=Layer(iN);
	s=Stype[il];

    return s;
}

double get_tc(int iN)
{
	double s;
	int il;
	//if(iN<N_S) il=0;
	//else if(iN<N_S+N_Mid) il=1;//il=(iN/N_Mid);
         //else il=2;
	il=Layer(iN);
	s = Tci[il];

	return s;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
////      _____    _    _____     _____   _____       ////
////      \  __\  |_|   \  __\   |  __/   \  __\      ////
////       \ \    | |    \ \     | |__     \ \        ////
////      __\ \   | |   __\ \    | __/    __\ \       ////
////      \____\  |_|   \____\   |_|      \____\      ////
////                                                  ////
//////////////////////////////////////////////////////////
////////////////////////by R3ZZ///////////////////////////
//////////////////////////////////////////////////////////

#include <cstring>
#include "stdafx.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#include "SFS.h"
#include <cmath>
#include <iomanip>
#include <thread>
#include <string>
#include <chrono>


#define pi 3.141592653589793
#define icom (complex<double>(0, 1.))

using namespace std;

//////// INITIALIZATION /////////

double Gb_I, Gb_FS, Gb_SF, R0A,  Ksi_F, Ksi_S, Del0, Del0a, Del0b, Xi1, Xi2, T, w, w1, E1, E, ro_F, ro_S,H, alphaT;//
double  L_SL, L_SR, L_S, L_N, L_F, L_S1, L_F1, L_S2, L_s;
int  w_obrez, iter, MODE;
int N_SL, N_SR, N_Mid, N_S, N_N, N_F, N_S1, N_F1, N_S2, N, N_I, N_CPR, NUM_Tech;
double epsG, epsDel, alpha;
double x;
double *Hi, *Ksii, *Roi, *Li, *Stype, *Tci, *Rbi;

int main()
{

///////// VARIABLES //////////////

	// Number of Layers

	NUM_Tech=6;

    //   Number of points in grid for each layer
    N_Mid=1000;// in self-sonst
    N_S = int(5.0 * N_Mid); N_F = int(0.15 * N_Mid); N_S1 = int(0.2 * N_Mid); N_F1 = int(0.25 * N_Mid); N_S2 = int(2.5 * N_Mid); N_I = int(0.01 * N_Mid); N_N = int(0.01 * N_Mid);
    N=N_S+N_F+N_S1+N_F1+N_S2+N_N;//*NUM_Tech;
    L_S=N_S/(1.*N_Mid); L_F=N_F/(1.*N_Mid); L_S1=N_S1/(1.*N_Mid); L_F1=N_F1/(1.*N_Mid); L_S2=N_S2/(1.*N_Mid); L_N=N_N/(1.*N_Mid); L_N=N_N/(1.*N_Mid);
    //cout<<N_S<<" "<<N_F<<" "<<N_S1<<" "<<N_F1<<" "<<N_S2<<" "<<N<<endl;
	N_CPR=1;

    Hi = new double[NUM_Tech];
    Ksii = new double[NUM_Tech]; // Normalized on S
    Roi = new double[NUM_Tech];  // Normalized on S
	Li = new double[NUM_Tech];
	Stype = new double[NUM_Tech];
	Tci = new double[NUM_Tech];
	Rbi = new double[NUM_Tech];
    double *DifDel, DifDelmax, DifDelmax_place, KIDP,KID;
    complex<double> *DelmaxP, *DelmaxAP;
    DifDel = new double[126];
    int dd,findmax,stor, Ns=100000;
	// Initialization of arrays
    DelmaxP = new complex<double>[10*N];
    DelmaxAP = new complex<double>[10*N];

    unsigned int n = std::thread::hardware_concurrency();
    std::cout << n << " concurrent threads are supported.\n";

T=0.5;
    //   Iteration constants
    epsG=1e-9;  // accuracy of normal Green function G iteration loop
    epsDel=5e-6; // accuracy of pair potential \Delta iteration loop
    alpha=0.59; // parameter in Delta loop. Just for increase of convergence
    iter=0; // variable to count number of iterations in Delta loop
    MODE=0;
    //   Material parameters
 Ksi_S=1; // Coherence length of superconductor. Used for normalization of all other lengths and scales
 Ksi_F=2.5;//3; // Coherence length of ferromagnet
 Xi1=0; // Phases of pair potential on the left (Xi1) and right (Xi2) electrodes
 Xi2=-Xi1; // Josephson phase is their difference: Xi = Xi2 - Xi1
 w_obrez=int(30/T); // Number of Matsubara frequencies used in calculation.
 H =100; // Exchange field
 ro_S=1;//520;  // Resistivity of superconductor
 ro_F=1;//440;  // Resistivity of ferromagnet
Gb_FS=0.3; // boundary parameter of Ferromagnet/Superconductor interface. It depends from resistivity of the interface.
R0A = Gb_FS *ro_F* Ksi_F;
Gb_SF= Gb_FS *ro_F* Ksi_F/ ro_S/ Ksi_S; // it also deterimines Superconductor/Ferromagnet interface
Gb_I = 10000;
double DELS, DELPmax,DELAPmax, roo,ksii, dss, Hext, GB_I;
int kk = 0, iter1, itermax, itermax1, culc, dN;
//_____________  Creation of OUTPUT file_______________________//
        string name1("6th(ds)H="),name2("6th(X,ds)H="),name3("2000stDSmax(X,ksi=ro,ds)H="),str1,str2,str3,str4,str5,str6,str7,str8,str9;// cheate one big file
		stringstream s1,s2,s3,s4,s5,s6,s7,s8,s9;
		s1<<H; s2<<T; s3<<L_S1; s4<<L_F; s5<<L_F1; s9<<L_S2; s6<<Gb_FS; s7<<Ksi_F; s8<<ro_F;
		s1>>str1; s2>>str2; s3>>str3; s4>>str4; s5>>str5; s6>>str6; s7>>str7; s8>>str8; s9>>str9;
		name1.append(str1); name1.append("T=");name1.append(str2);name1.append("LS1=");name1.append(str3);name1.append("LF=");name1.append(str4);name1.append("LF1=");name1.append(str5);name1.append("Ls=");name1.append(str9);name1.append("gb=");name1.append(str6);name1.append("KsiF=");name1.append(str7);name1.append("roF=");name1.append(str8);                                            name1.append(".txt");
		name2.append(str1); name2.append("T=");name2.append(str2);name2.append("LS1=");name2.append(str3);name2.append("LF=");name2.append(str4);name2.append("LF1=");name2.append(str5);name2.append("Ls=");name2.append(str9);name2.append("gb=");name2.append(str6);name2.append("KsiF=");name2.append(str7);name2.append("roF=");name2.append(str8);                                            name2.append(".txt");
		name3.append(str1); name3.append("T=");name3.append(str2);name3.append("LS1=");name3.append(str3);name3.append("LF=");name3.append(str4);name3.append("LF1=");name3.append(str5);name3.append("Ls=");name3.append(str9);name3.append("gb=");name3.append(str6);name3.append("KsiF=");name3.append(str7);name3.append("roF=");name3.append(str8);                                            name3.append(".txt");
		const char *file1=name1.c_str();
		const char *file2=name2.c_str();
		const char *file3=name3.c_str();
		ofstream fout1(file1);
		ofstream fout2(file2);
		//ofstream fout3(file3);
double x0,x1,x2; char chh;
		ifstream file ("11e5DELth(X,dF)H=100T=0.5LS1=0LF=0LF1=0Ls=0gb=0.3KsiF=2.5roF=1.txt");
		if (file.is_open()) cout << "good\n";
		else cout << "beed\n\n" << endl;/*
		for(int i=0;i<20;i++)
        { cout<<i;
            for(int j=0;j<7569;j++) {//cout<<j<<endl;
            file>>x1;}
            DifDel[i]=x1; for(int j=0;j<12;j++) file>>chh; file>>x1;
        }
for(int j=0;j<16;j++) fout1<<j<<"  "<< DifDel[j]<<endl;/*/
/*for(T=0.1;T<0.701;T+=0.05)
//if((abs(T-0.5)<1e-5)||(abs(T-0.25)<1e-5)||(abs(T-0.5)<1e-5))//||(abs(T-0.4)<1e-5)||(abs(T-0.5)<1e-5))
//*/
//for(H=5;H<50.46;H+=5)
//if((abs(H-0.2)<1e-5)||(abs(H-0.3)<1e-5))//||(abs(H-50)<1e-5))//||(abs(H-4)<1e-5)||(abs(H-10)<1e-5))//||(abs(H-0.5)<1e-5)||(abs(H-1)<1e-5)||(abs(H-2)<1e-5)||(abs(H-4)<1e-5)||(abs(H-10)<1e-5))
//*////
/*for(N_S1=int(0.2*N_Mid);N_S1<int(2.1*N_Mid);N_S1+=int(0.2*N_Mid))
if((abs(N_S1-int(0.2*N_Mid))<1e-5))//||(abs(N_S1-int(1*N_Mid))<1e-5)||(abs(N_S1-int(2*N_Mid))<1e-5))
{    L_S1=N_S1/(1.*N_Mid);
    DifDelmax_place=1.88;/*/

/*for(Gb_FS=0.3; Gb_FS<4; Gb_FS+=0.05)
if((abs(Gb_FS-0.3)<1e-5))//||(abs(Gb_FS-0.1)<1e-5))//||(abs(Gb_FS-0.2)<1e-5)||(abs(Gb_FS-0.3)<1e-5)||(abs(Gb_FS-0.4)<1e-5)||(abs(Gb_FS-0.5)<1e-5)||(abs(Gb_FS-1)<1e-5))//||(abs(Gb_FS-4)<1e-5)||(abs(Gb_FS-10)<1e-5))
{/*/
    //fout1<<"Gb_FS="<<Gb_FS<<endl;
    //fout2<<"Gb_FS="<<Gb_FS<<endl;
double** StorageDEL = new double *[5];  //ds    DELP    DELAP   ||DELP|-|DELAP||
     for (int i = 0; i < 5; i++) StorageDEL[i] = new double[150];

for(Ksi_F=2.5;Ksi_F<4.01;Ksi_F+=10.1)
//if((abs(Ksi_F-1)<1e-5))//||(abs(Ksi_F-1.5)<1e-5)||(abs(Ksi_F-2)<1e-5)||(abs(Ksi_F-3)<1e-5))//||(abs(Ksi_F-1.5)<1e-5)||(abs(Ksi_F-2)<1e-5))//||(abs(Ksi_F-3.5)<1e-5))//||(abs(Ksi_F-0.75)<1e-5))
{ //fout1<<"xi="<<Ksi_F<<endl;/*
    DifDelmax_place=2.5;
for(ro_F=1;ro_F<3.01;ro_F+=10.1)
//if((abs(ro_F-0.8)<1e-5)||(abs(ro_F-1.4)<1e-5)||(abs(ro_F-1)<1e-5)||(abs(ro_F-1.2)<1e-5))
//if((abs(ro_F-Ksi_F)<1e-5)||((abs(ro_F-1)<1e-5)&&(abs(Ksi_F-10)<1e-5))||((abs(ro_F-10)<1e-5)&&(abs(Ksi_F-1)<1e-5)))
//*/
for(N_F=int(0.15*N_Mid);N_F<int(0.41*N_Mid);N_F+=int(10.02*N_Mid))
//if((abs(L_F-0.5)<1e-5))//||(abs(L_F-0.6)<1e-5)||(abs(L_F-0.7)<1e-5)||(abs(L_F-0.8)<1e-5)||(abs(L_F-0.9)<1e-5)||(abs(L_F-1.1)<1e-5))
{

//for(N_N=int(0.1*N_Mid);N_N<int(2.01*N_Mid);N_N+=int(10.1*N_Mid))
//{
//for(int gg=3; gg<4;gg++)
//for(N_S1=int(0.2*N_Mid);N_S1<int(2.99*N_Mid);N_S1+=int(10.5*N_Mid))
//{   //Gb_I=1;//pow(10,1.0*gg);
//cout<<"Gb_I="<<Gb_I<<endl;

    L_S1=N_S1/(1.*N_Mid);
    N_F1=N_F+int(0.1*N_Mid);
    L_F=N_F/(1.*N_Mid);
    L_F1=N_F1/(1.*N_Mid);//*/
    R0A = Gb_FS *ro_F* Ksi_F;
    kk = 0; iter1 = 0; itermax = 0; itermax1 = 0;  culc = 0; dN = int(0.05 * N_Mid); N_S2 = int((DifDelmax_place + 0.0) * N_Mid); L_S2 = N_S2 / (1. * N_Mid);//+55*4;
DifDelmax=0; stor=0;

for (int i = 0; i < 5; i++) for (int j = 0; j < 150; j++) StorageDEL[i][j]=0.;
while((kk<1))
{
    //if(kk==1) N_S2-=20*4;
    //if(kk==2) N_S2-=10*4;
    //if(k==3) N_S2-=5*4;
findmax=0; DifDelmax=0;
while((findmax==0)&&(L_S2>2.411))//
{
    //N_S2-=dN;
    N=N_S+N_F+N_S1+N_F1+N_S2;//+N_N;//+N_N;//
    //N=N_S+N_F+N_S1;
    L_S2=N_S2/(1.*N_Mid);
    //L_N=N_N/(1.*N_Mid);
    //Ns=0;//N_S+N_F+N_S1+N_F1+N_S2-2;
cout<<N_S<<" "<<N_F<<" "<<N_S1<<" "<<N_F1<<" "<<N_S2<<" "<<N_N<<" "<<N<<endl;
	for (int i = 0; i < NUM_Tech; i++)
	{
		if (i == 0)
		{
			 Hi[i] = 0; Ksii[i] = Ksi_S; Roi[i] = ro_S; Li[i] = L_S; Stype[i] = 1.; Tci[i] =1; Rbi[i] = R0A;
		}
		if ((i == 1))//||(i == 5))
		{
			Hi[i] = H; Ksii[i] = Ksi_F; Roi[i] = ro_F; Li[i] = L_F; Stype[i] = 0.;  Tci[i] = 0.; Rbi[i] = R0A;
		}
		if ((i == 2))//||(i == 4)||(i == 6))
		{
			 Hi[i] = 0; Ksii[i] = Ksi_S; Roi[i] = ro_S; Li[i] = L_S1; Stype[i] = 1.; Tci[i] =1; Rbi[i] = R0A;
		}//*/
		if ((i == 3)||(i == 7))
		{
			Hi[i] = H; Ksii[i] = Ksi_F; Roi[i] = ro_F; Li[i] = L_F1; Stype[i] = 0.;  Tci[i] = 0.; Rbi[i] = R0A;
		}
		if (i == 4)
		{
			 Hi[i] = 0; Ksii[i] = Ksi_S; Roi[i] = ro_S; Li[i] = L_S2; Stype[i] = 1.; Tci[i] =1; Rbi[i] = R0A;
		}
		if (i == 5)
		{
			 Hi[i] = 0; Ksii[i] = Ksi_S; Roi[i] = ro_S; Li[i] = L_N; Stype[i] = 1.; Tci[i] =1; Rbi[i] =  R0A;
		}
	}
	complex<double> *Del, *DelP, *Delbuf;
    double *Is;//, *Is1, *dIdV, ro, dG,x, Ej;
    Del = new complex<double>[N];
    DelP = new complex<double>[N];
    Delbuf = new complex<double>[N];
    Is = new double[N];
    complex<double> *G, *Fi, *Fi_old, *Fi1,*Fi1_old;//, *Del, *dDel,*Delbuf;
    //double *Is, *Is1, *dIdV, ro, dG,x;
    G = new complex<double>[10*N];
    //G1 = new complex<double>[10*N];
	Fi = new complex<double>[10*N];
	Fi_old = new complex<double>[10*N];
	Fi1 = new complex<double>[10*N];
	Fi1_old = new complex<double>[10*N];


    
    Del0 = SelfConsZero();
culc=1;
for (int i=0; i<150; i++)
if(abs(StorageDEL[0][i]-L_S2)<1e-5)
{
    DelP[N-1]=StorageDEL[1][i]; Del[N-1]=StorageDEL[2][i]; iter1=int(StorageDEL[3][i]); iter=int(StorageDEL[4][i]);
    culc=0; // we found the last calculated at the ds lenght
}
if(culc==1)
{
////////////////PARALLEL ORIENTATION///////////

        MODE=0;
    SelfConsParal(G, DelP, 0);     // Self-consistent  calculation of pair potential.
    iter1=iter;
    //KIDP = KinInd(G, DelP);
    //MODE=0;
    //E_IS_calc(Is,G,DelP);
    //fout1<<fixed<<Xi2<<"  "<<L_S2<<"  "<<GB_I<<"  "<<ro_F<<"   "<<Ksi_F<<"   "<<real(DelP[N-N_N-1])<<"  "<<Is[10]<<"  "<<Is[N-N_N-10]<<endl;
    double dGmax;
	for(int iw=0; iw<1; iw++)
    {       dGmax=1;
            w=pi*T*(2.*iw+1.);
	        for (int i=0; i<N; i++) G[i]=1.+0.1*icom;
	        //for (int i=0; i<N; i++)	Del[i]= get_type(i)*Del0;
	        while ((dGmax > epsG))
	        {
				        for (int i=0; i<N; i++) Fi_old[i]= Fi[i];
				        Prog( Fi, G, DelP, w);
                        for (int i=0; i<N; i++)
                        {	G[i]= -conj(G[i]);
					        Fi[i]=0.6*Fi[i]+0.4*Fi_old[i];
				        }
				        for (int i=0; i<N; i++)	Fi1_old[i]= Fi1[i];
            	        for (int i=0; i<N; i++)
				        {
					        G[i]= -conj(G[i]);
					        Fi1[i]=conj(Fi[i]);
				        }
                        Gcalc( G, &dGmax, Fi, Fi1, DelP, w);
	        }//*/
	    //II+=T*Del0*exp(icom*Xi2*pi)/sqrt(w*w+Del0*Del0)*G[N-1]*Fi[N-1]/get_wm(N-1,w);
        for(int n=0; n<N;++n)	fout2<<fixed<<n<<"  "<<real(DelP[n])<<"  "<<imag(DelP[n])<<"  "<<abs(G[n]*Fi[n]/get_wm(n,w))<<"  "<<real(G[n]*Fi[n]/get_wm(n,w))<<"  "<<imag(G[n]*Fi[n]/get_wm(n,w))<<"  "<<arg(G[n]*Fi[n]/get_wm(n,w))<<"  "<<Is[n]<<endl;
        }//for(int n=0; n<N;++n)	fout2<<fixed<<n<<"  "<<real(Del[n])<<endl;/*/
//fout2<<fixed<<Xi2<<"  "<<abs(II)<<"  "<<real(II)<<endl;
//}

//////////////ANTIPARALLEL ORIENTATION///////////////////////////////////
 Hi[1]=-H; //Hi[7]=-H;
    //for (int i=0; i<N; i++)	file>>x0>>Delbuf[i]>>x0>>x0>>x0>>x1>>x0;
    //for (int i=0; i<N; i++)	DelP[i]= get_type(i)*Delbuf[i];
SelfConsParal(G, Del, 0);
//KID = KinInd(G, Del);
	//iterG = 0;
	//Xi2=0.5;
	//E_IS_calc(Is,G,DelP);
	for (int i=0; i<N; i++)	Is[i]=0.;
    for (int iw = 0; iw < 1; iw++)
    {
        dGmax = 1;
        w = pi * T * (2. * iw + 1.);
        for (int i = 0; i < N; i++) G[i] = 1. + 0.1 * icom;
        //for (int i=0; i<N; i++)	Del[i]= get_type(i)*Del0;
        while ((dGmax > epsG))
        {
            for (int i = 0; i < N; i++) Fi_old[i] = Fi[i];
            Prog(Fi, G, Del, w);
            for (int i = 0; i < N; i++)
            {
                G[i] = -conj(G[i]);
                Fi[i] = 0.6 * Fi[i] + 0.4 * Fi_old[i];
            }
            for (int i = 0; i < N; i++)	Fi1_old[i] = Fi1[i];
            for (int i = 0; i < N; i++)
            {
                G[i] = -conj(G[i]);
                Fi1[i] = conj(Fi[i]);
            }
            Gcalc(G, &dGmax, Fi, Fi1, Del, w);
        }//*/
    //II+=T*Del0*exp(icom*Xi2*pi)/sqrt(w*w+Del0*Del0)*G[N-1]*Fi[N-1]/get_wm(N-1,w);
        for (int n = 0; n < N; ++n)	fout2 << fixed << n << "  " << real(Del[n]) << "  " << imag(Del[n]) << "  " << abs(G[n] * Fi[n] / get_wm(n, w)) << "  " << real(G[n] * Fi[n] / get_wm(n, w)) << "  " << imag(G[n] * Fi[n] / get_wm(n, w)) << "  " << arg(G[n] * Fi[n] / get_wm(n, w)) << "  " << Is[n] << endl;
    } 
    //for(int n=0; n<N;++n)	fout2<<fixed<<n<<"  "<<real(Del[n])<<endl;
//E_IS_calc(Is,G,Del);
//for(int i=0; i<N_Mid; i++) Del[i]=Del0+icom*0.;for(int i=N_Mid; i<N; i++) Del[i]=0.+icom*0.;
/*double x1,x2,x3,x4,x5,Delbuf;/*/
//for(int n=0; n<N;++n)	fout2<<fixed<<n<<"  "<<L_F1<<"   "<<Xi2/pi<<"  "<<real(Del[n])<<"   "<<imag(Del[n])<<"   "<<Is[n]<<endl;//Delbuf; Del[n]=Delbuf+icom*0.0;
//}
//*/
//for(int n=0; n<N;++n)	fout2<<fixed<<n<<"  "<<real(Del[n])<<endl;
    //DELS=real(Del[N-1]);
    StorageDEL[0][stor]=L_S2; StorageDEL[1][stor]=real(DelP[N-1]); StorageDEL[2][stor]=real(Del[N-1]); StorageDEL[3][stor]=double(iter1); StorageDEL[4][stor]=double(iter);
    stor++;
}//*/

//fout1<<fixed<<L_F1<<"   "<<Xi2/pi<<"   "<<real(Del[0])<<"   "<<imag(Del[0])<<"   "<<real(Del[N-1])<<"   "<<imag(Del[N-1])<<"   "<<Is[N/2]<<endl;
//fout1<<fixed<<L_S2<<"  "<<L_F<<"  "<<L_S1<<"  "<<ro_F<<"   "<<Ksi_F<<"   "<<real(DelP[N-1])<<"   "<<real(Del[N-1])<<"   "<<abs(abs(real(Del[N-1]))-abs(real(DelP[N-1])))<<"  "<<KIDP<<"  "<<KID<<"  "<<iter1<<"  "<<iter<<"  "<<culc<<endl;

    if(abs(abs(real(Del[N-1]))-abs(real(DelP[N-1])))>DifDelmax)
    {
        DifDelmax=abs(abs(real(Del[N-1]))-abs(real(DelP[N-1])));
        DifDelmax_place=L_S2;
        itermax1=iter1; itermax=iter;
        DELPmax=real(DelP[N-1]);DELAPmax=real(Del[N-1]);
        for(int i=0; i<N; i++) DelmaxP[i]=DelP[i];
        for(int i=0; i<N; i++) DelmaxAP[i]=Del[i];
        //for(int i=0; i<N; i+=100) cout<<DelmaxAP[i]<<"  "Del[i]<<endl;
    }
    else findmax=1;
    findmax=1;  //stop calc, we find max in this k-iteration
    if(kk==0) N_S2-=dN;
    if(kk==1) N_S2-=dN/2;
    if(kk==2) N_S2-=dN/4;
    if(kk==3) N_S2-=dN/8;
    //if(kk==3) N_S2-=1*50;
    //if((kk==2)&&(findmax==1)) { cout<<"k=========="<<kk;
    for(int i=0; i<N; i++) fout2<<fixed<<i<<"  "<<ro_F<<"   "<<Ksi_F<<" "<<real(DelP[i])<<" "<<real(Del[i])<<endl;
    //for(int i=0; i<N; i++) fout3<<fixed<<i<<"  "<<ro_F<<"   "<<Ksi_F<<" "<<real(DelmaxAP[i])<<endl;}//*/
    delete [] G;
    delete [] Del;
    delete [] DelP;
    delete [] Delbuf;
    delete [] Fi;
    delete [] Fi_old ;
    delete [] Fi1;
    delete [] Fi1_old;
}
    if(kk==0) N_S2+=(3*dN);
    if(kk==1) N_S2+=(3*dN)/2;
    if(kk==2) N_S2+=(3*dN)/4;
//if(kk==0) N_S2+=(5*dN)/2;if(kk==1) N_S2+=(5*dN)/8;if(kk==2) N_S2+=(5*dN)/16;
    kk++;//*/
}

//fout2<<fixed<<L_F<<"  "<<L_S1<<"  "<<ro_F<<"   "<<Ksi_F<<" "<<DifDelmax<<"   "<<DifDelmax_place<<"   "<<itermax1<<"   "<<itermax<<endl;
fout2<<fixed<<ro_F<<"   "<<Ksi_F<<" "<<DifDelmax<<"   "<<DifDelmax_place<<"   "<<DELPmax<<"   "<<DELAPmax<<"   "<<itermax1<<"   "<<itermax<<endl;


}}//}//}//}
 delete [] DelmaxP;
 delete [] DelmaxAP;

for (int i = 0; i < 5; i++)
    delete[] StorageDEL[i];
delete[] StorageDEL;

fout1<<"\n"<<"H="<<H<<", T="<<T<<", L_S="<<L_S<<", L_F="<<L_F<<", L_S1="<<L_S1<<", L_F1="<<L_F1<<", L_S2="<<L_S2<<", Gb_FS="<<Gb_FS<<", ksi_F="<<Ksi_F<<", ro_F="<<ro_F<<", R0A="<<R0A<<", Nmid="<<N_Mid<<"\n";
fout2<<"\n"<<"H="<<H<<", T="<<T<<", L_S="<<L_S<<", L_F="<<L_F<<", L_S1="<<L_S1<<", L_F1="<<L_F1<<", L_S2="<<L_S2<<", Gb_FS="<<Gb_FS<<", ksi_F="<<Ksi_F<<", ro_F="<<ro_F<<", R0A="<<R0A<<", Nmid="<<N_Mid<<"\n";



    system("PAUSE");
    return EXIT_SUCCESS;
}

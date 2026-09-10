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
#include "sns/runner.hpp"
#include <cmath>
#include <iomanip>
#include <thread>
#include <string>
#include <chrono>


#define pi 3.141592653589793
#define icom (complex<double>(0, 1.))

using namespace std;

//////// INITIALIZATION /////////

double Gb_I, Gb_FS, Gb_SF, R0A, R0AN, Ksi_F, Ksi_S, Ksi_N, Del0, Del0a, Del0b, Xi1, Xi2, T, w, w1, E1, E, ro_F, ro_N, ro_S, H, alphaT;//
double  L_SL, L_SR, L_S, L_N, L_F, L_S1, L_F1, L_S2, L_s;
int  w_obrez, iter, MODE;
int N_SL, N_SR, N_Mid, N_S, N_N, N_F, N_S1, N_F1, N_S2, N, N_CPR, NUM_Tech;
double epsG, epsDel, alpha;
double x;
double* Hi, * Ksii, * Roi, * Li, * Stype, * Tci, * Rbi;

int run_legacy_calculation();

int main(int argc, char** argv)
{
    // Select false to run the original calculation without command-line arguments.
    const bool use_sns = true;
    if (!use_sns || (argc > 1 && std::string(argv[1]) == "--legacy"))
        return run_legacy_calculation();

    // SNS PARAMETERS: energies in k_B*Tc, lengths in xi_S, phase in pi.
    sns::RunSettings settings;
    settings.physical.Delta = 1.76;
    settings.physical.T = 0.0;
    settings.physical.Ksi_N = 1.0;
    settings.physical.L_N = 1.0;// *std::sqrt(settings.physical.diffusion() / settings.physical.Delta); // L/xi_Delta=1
    settings.physical.ro_N = 1.0;
    settings.physical.area = 1.0;
    settings.physical.Xi = 0.0;

    settings.numerical.NF = 4;
    settings.numerical.Nx = 25;
    settings.numerical.Neps = 16;
    settings.numerical.eta = 0.001 * settings.physical.Delta;
    settings.numerical.mixing = 0.15;
    settings.numerical.tolerance = 1e-8;
    settings.numerical.residual_tolerance = 1e-7;
    settings.numerical.max_iterations = 4000;
    // Write G^R_nm(x) and F^R_nm(x) after EVERY accepted gamma iteration.
    // Use an empty string to disable the output.
    //settings.numerical.iteration_log_path =        "C:/Users/user/source/repos/aleksneilo/keldysh/Ricatti/GF_iterations.txt";

    // spectral, normal, iv, benchmark, convergence, help
    //settings.mode = "spectral";
    //settings.mode = "conductance";

    settings.mode = "convergence";
    settings.voltage = 3.0 * settings.physical.Delta;

    settings.sweep_start = 3.0; // eV/Delta (article voltage)
    settings.sweep_end = 2.7;
    settings.sweep_step = 0.1; // Positive magnitude; direction is automatic.
    settings.conductance_path = "C:/Users/user/source/repos/aleksneilo/keldysh/Ricatti/conductance_L1.txt";
    settings.epsilon = 3.0 * settings.physical.Delta;
    //settings.voltage = 0.0; // Set e.g. 2*Delta for a finite-voltage calculation.
    return sns::run_cli(argc, argv, settings);
}

int run_legacy_calculation()
{

    ///////// VARIABLES //////////////

        // Number of Layers

    NUM_Tech = 1;

    //   Number of points in grid for each layer
    N_Mid = 100;// in self-sonst
    N_N = int(1 * N_Mid); //N_S = int(5 * N_Mid); N_F = int(0.5 * N_Mid); N_S1 = int(0.2 * N_Mid); N_F1 = int(0.6 * N_Mid); N_S2 = int(2.8 * N_Mid); N_N = int(5 * N_Mid);
    N = N_S + N_N;// +N_S1 + N_F1 + N_S2;// +N_N;//*NUM_Tech;
    L_S = N_S / (1. * N_Mid);// L_F = N_F / (1. * N_Mid); L_S1 = N_S1 / (1. * N_Mid); L_F1 = N_F1 / (1. * N_Mid); L_S2 = N_S2 / (1. * N_Mid); L_N = N_N / (1. * N_Mid);
    //cout<<N_S<<" "<<N_F<<" "<<N_S1<<" "<<N_F1<<" "<<N_S2<<" "<<N<<endl;
    N_CPR = 1;

    Hi = new double[NUM_Tech];
    Ksii = new double[NUM_Tech]; // Normalized on S
    Roi = new double[NUM_Tech];  // Normalized on S
    Li = new double[NUM_Tech];
    Stype = new double[NUM_Tech];
    Tci = new double[NUM_Tech];
    Rbi = new double[NUM_Tech];
    double* DifDel, DifDelmax, DifDelmax_place, KIDP, KID;
    complex<double>* DelmaxP, * DelmaxAP, * DelPbuf, * DelAPbuf;
    DifDel = new double[126];
    int dd, findmax, stor, Ns = 1000000;
    // Initialization of arrays
    DelmaxP = new complex<double>[10 * N];
    DelmaxAP = new complex<double>[10 * N];
    DelPbuf = new complex<double>[10 * N];
    DelAPbuf = new complex<double>[10 * N];

    unsigned int n = std::thread::hardware_concurrency();
    std::cout << n << " concurrent threads are supported.\n";

    T = 0.1;
    //   Iteration constants11
    epsG = 1e-9;  // accuracy of normal Green function G iteration loop
    epsDel = 5e-6; // accuracy of pair potential \Delta iteration loop
    alpha = 0.59; // parameter in Delta loop. Just for increase of convergence
    iter = 0; // variable to count number of iterations in Delta loop
    MODE = 0;
    //   Material parameters
    Ksi_S = 1; // Coherence length of superconductor. Used for normalization of all other lengths and scales
    Ksi_F = 1;//3; // Coherence length of ferromagnet
    Ksi_N = 10;
    Xi1 = 0; // Phases of pair potential on the left (Xi1) and right (Xi2) electrodes
    Xi2 = -Xi1; // Josephson phase is their difference: Xi = Xi2 - Xi1
    w_obrez = int(30 / T); // Number of Matsubara frequencies used in calculation.
    H = 10; // Exchange field
    ro_S = 1;//520;  // Resistivity of superconductor
    ro_F = 1;//440;  // Resistivity of ferromagnet
    ro_N = 0.1;
    Gb_FS = 0.3; // boundary parameter of Ferromagnet/Superconductor interface. It depends from resistivity of the interface.
    R0A = Gb_FS * ro_F * Ksi_F;
    Gb_SF = Gb_FS * ro_F * Ksi_F / ro_S / Ksi_S; // it also deterimines Superconductor/Ferromagnet interface

    double DELS = 1, DELPmax, DELAPmax, dqmax, epsq, roo, ksii, dss, Hext, IIP, II, IIPmax, IImax, I;
    int kk, iter1, itermax, itermax1, iterq, culc, dN;
    epsq = 1e-5;

    //_____________  Creation of OUTPUT file_______________________//
    string name1("5stLk(I)H="), name2("5stLk(X,I)H="), name3("2000stDSmax(X,ksi=ro,ds)H="), str1, str2, str3, str4, str5, str6, str7, str8, str9;// cheate one big file
    stringstream s1, s2, s3, s4, s5, s6, s7, s8, s9;
    s1 << H; s2 << T; s3 << L_S1; s4 << L_F; s5 << L_F1; s9 << L_S2; s6 << Gb_FS; s7 << Ksi_F; s8 << ro_F;
    s1 >> str1; s2 >> str2; s3 >> str3; s4 >> str4; s5 >> str5; s6 >> str6; s7 >> str7; s8 >> str8; s9 >> str9;
    name1.append(str1); name1.append("T="); name1.append(str2); name1.append("LS1="); name1.append(str3); name1.append("LF="); name1.append(str4); name1.append("LF1="); name1.append(str5); name1.append("Ls="); name1.append(str9); name1.append("gb="); name1.append(str6); name1.append("KsiF="); name1.append(str7); name1.append("roF="); name1.append(str8);                                            name1.append(".txt");
    name2.append(str1); name2.append("T="); name2.append(str2); name2.append("LS1="); name2.append(str3); name2.append("LF="); name2.append(str4); name2.append("LF1="); name2.append(str5); name2.append("Ls="); name2.append(str9); name2.append("gb="); name2.append(str6); name2.append("KsiF="); name2.append(str7); name2.append("roF="); name2.append(str8);                                            name2.append(".txt");
    name3.append(str1); name3.append("T="); name3.append(str2); name3.append("LS1="); name3.append(str3); name3.append("LF="); name3.append(str4); name3.append("LF1="); name3.append(str5); name3.append("Ls="); name3.append(str9); name3.append("gb="); name3.append(str6); name3.append("KsiF="); name3.append(str7); name3.append("roF="); name3.append(str8);                                            name3.append(".txt");
    const char* file1 = name1.c_str();
    const char* file2 = name2.c_str();
    const char* file3 = name3.c_str();
    ofstream fout1(file1);
    ofstream fout2(file2);
    //ofstream fout3(file3);
    double x0, x1, x2; char chh;
    ifstream file("11e5DELth(X,dF)H=100T=0.5LS1=0LF=0LF1=0Ls=0gb=0.3KsiF=2.5roF=1.txt");
    if (file.is_open()) cout << "good\n";
    else cout << "beed\n\n" << endl;



    //for (double hh = 1.; hh < 1.01; hh += 2.)
    //{
    
                for (int i = 0; i < NUM_Tech; i++)
                {
                    if (i == 1)
                    {
                        Hi[i] = 0; Ksii[i] = Ksi_S; Roi[i] = ro_S; Li[i] = L_S; Stype[i] = 1.; Tci[i] = 1; Rbi[i] = R0A;
                    }
                    if ((i == 0))//||(i == 5))
                    {
                        Hi[i] = 0; Ksii[i] = Ksi_N; Roi[i] = ro_N; Li[i] = L_N; Stype[i] = 0.;  Tci[i] = 0.; Rbi[i] = R0A;
                    }
                    if ((i == 2))//||(i == 4)||(i == 6))
                    {
                        Hi[i] = 0; Ksii[i] = Ksi_S; Roi[i] = ro_S; Li[i] = L_S1; Stype[i] = 1.; Tci[i] = 1; Rbi[i] = R0A;
                    }//*/
                    
                }
                complex<double>* Del, * DelP, * Delbuf, * G, * Fi, * Fi_old, * Fi1, * Fi1_old;
                double* Is, * q, * qP, * qbuf;
                Del = new complex<double>[N];
                DelP = new complex<double>[N];
                Delbuf = new complex<double>[N];
                Is = new double[N];
                q = new double[N];
                qP = new double[N];
                qbuf = new double[N];
                G = new complex<double>[N];
                //G1 = new complex<double>[10*N];
                Fi = new complex<double>[N];
                Fi_old = new complex<double>[N];
                Fi1 = new complex<double>[N];
                Fi1_old = new complex<double>[N];


                Del0 = SelfConsZero();
                
                    
                    findmax = 1;
                    delete[] G;
                    delete[] Del;
                    delete[] DelP;
                    delete[] Delbuf;
                    delete[] Fi;
                    delete[] Fi_old;
                    delete[] Fi1;
                    delete[] Fi1_old;
                    delete[] q;
                    delete[] qP;
                    delete[] qbuf;
        
    delete[] DelmaxP;
    delete[] DelmaxAP;
    delete[] DelPbuf;
    delete[] DelAPbuf;


    fout1 << "\n" << "H=" << H << ", T=" << T << ", L_S=" << L_S << ", L_F=" << L_F << ", L_S1=" << L_S1 << ", L_F1=" << L_F1 << ", L_S2=" << L_S2 << ", Gb_FS=" << Gb_FS << ", ksi_F=" << Ksi_F << ", ro_F=" << ro_F << ", R0A=" << R0A << ", Nmid=" << N_Mid << "\n";
    fout2 << "\n" << "H=" << H << ", T=" << T << ", L_S=" << L_S << ", L_F=" << L_F << ", L_S1=" << L_S1 << ", L_F1=" << L_F1 << ", L_S2=" << L_S2 << ", Gb_FS=" << Gb_FS << ", ksi_F=" << Ksi_F << ", ro_F=" << ro_F << ", R0A=" << R0A << ", Nmid=" << N_Mid << "\n";



    system("PAUSE");
    return EXIT_SUCCESS;
}

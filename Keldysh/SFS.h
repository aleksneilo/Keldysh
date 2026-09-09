#include <complex>

 using namespace std;

 extern int N_SL, N_SR, N_Mid, N_S, N_N, N_F, N_S1, N_F1, N_S2, N, N_I, NUM_Tech;
 extern double epsG, epsDel, alpha;
extern int  w_obrez, iter;
extern int  MODE;  // regime of boundary conditions: 0 - Free Boundary, 1 - Fixed phase

 extern double T, Del0, Del0a, Del0b,Xi1, Xi2;
 extern double L_SL, L_SR, L_S, L_N, L_S1,L_S2, L_F1, Ksi_S, ro_S, Tc;
 extern double L_F, Ksi_F, Ksi_N, ro_F, ro_N, H, alphaT;
 extern double Gb_I, Gb_FS, Gb_SF, R0A, R0AN;
 extern double *Hi, *Ksii, *Roi, *Li, *Stype, *Tci, *Rbi;


 void SelfCons(complex<double> *G, complex<double> *Del, int Initial, double* q, double I);
 double SelfConsZero();
 void getABC(complex<double> *a, complex<double> *b, complex<double> *c, complex<double> *d, int iN, complex<double> *G, complex<double> *Del, complex<double> w, double* q, double I);
 void Prog( complex<double> *Fi, complex<double> *G, complex<double> *Del, complex<double> w, double* q, double I);
 void Gcalc( complex<double> *G, double *dGmax, complex<double> *Fi1, complex<double> *Fi2, complex<double> *Del, complex<double> w, double* q);
 void GcalcDOS(complex<double> *G, double *dGmax, complex<double> *Fi1, complex<double> *Fi2, complex<double> *Del, complex<double> w);
 void getABC_ricatti(complex<double>* a, complex<double>* b, complex<double>* c, complex<double>* d, int iN, complex<double>* gamma, complex<double>* gamma_tilde, complex<double>* Del, complex<double> w, double* q, double I, bool tilde_equation);
 void Prog_ricatti(complex<double>* gamma, complex<double>* gamma_tilde, complex<double>* G, complex<double>* Del, complex<double> w, double* q, double I);
 void Prog_ricatti(complex<double>* gamma, complex<double>* gamma_tilde, complex<double>* G, complex<double>* Del, complex<double> w, double* q, double I, bool use_initial_guess);
 void Gcalc_ricatti(complex<double>* G, double* dGmax, complex<double>* gamma, complex<double>* gamma_tilde);
 void SelfConsParal_ricatti(complex<double>* G, complex<double>* Del, int Initial, double* q, double I);
 void Scalc_ricatti(complex<double>* Del, complex<double>** S1, int w_in, int w_fin, double* q, double I);
 double E_IS_calc(double *Is, complex<double> *G,complex<double> *Del, double* q, double I);
 void DOS(double E, complex<double> *G, complex<double> *G1, complex<double> *Fi, complex<double> *Del, double* q, double I);
  void SelfConsParal(complex<double> *G, complex<double> *Del, int Initial, double* q, double I);
  //void ScalcWarm(complex<double>* Del, complex<double>** S1, vector<vector<complex<double>>>& G_cache, vector<int>& iterG_stat, int w_in, int w_fin, double* q, double I);
  void Scalc(complex<double>* Del, complex<double>** S1, int w_in, int w_fin, double* q, double I);// void Scalc(complex<double>* Del, complex<double>** S1, int& w_in, int& w_fin, double* q, double I);
 int Layer(int iN);
 complex<double> get_wm(int iN, complex<double> w);
 double KinInd( complex<double> *G, complex<double> *Del, double* q, double I);
 double SuF(complex<double> *G, complex<double> *Del, double* q, double I);
 double get_h(int iN);
 double get_HE(int iN);
 double get_h_out(int iN);
 double get_ksi(int iN);
 double get_type(int iN);
 double get_tc(int iN);



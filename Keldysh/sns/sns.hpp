#pragma once
#include "matrix.hpp"
#include <functional>
#include <string>

namespace sns {
// Legacy units: E,Delta,eta in k_B*Tc; T in Tc; voltage means e*V/(k_B*Tc).
// Lengths in xi_S=sqrt(hbar*D_S/(2*pi*k_B*Tc)), rho in rho_S, area in xi_S^2.
// Xi is the phase difference divided by pi, as Xi2-Xi1 in the old main.
struct PhysicalParams {
    double Delta = 1.76, T = 0;
    double Ksi_N = 1, L_N = 1, ro_N = 1, area = 1, Xi = 0;
    double diffusion() const { return 2*std::acos(-1.)*Ksi_N*Ksi_N; }
    double thouless() const { return diffusion()/(L_N*L_N); }
    double conductance() const { return area/(ro_N*L_N); }
    double phase_radians() const { return std::acos(-1.)*Xi; }
};
struct NumericalParams {
    int NF=2, Nx=25, Neps=16, max_iterations=4000;
    double eta=0.00176, mixing=0.15, tolerance=1e-8;
    double residual_tolerance=1e-6, kinetic_tolerance=1e-6;
    int continuation_steps=4;
    int energy_threads=1; // 0: hardware threads minus one; 1: serial; >1: explicit count.
    bool use_anderson=false, anderson_verbose=false;
    int anderson_depth=4, anderson_start=2;
    double anderson_regularization=1e-10, anderson_coefficient_limit=20.0;
    size_t max_block_bytes=1024ULL*1024*1024;
    // Empty disables output. The solver appends; run_cli clears once per run.
    std::string iteration_log_path;
};
using Field = std::vector<Matrix>;
struct PairField { Field gamma, tilde; };
struct Normalizers { Matrix N1,N2; };
struct Green { Matrix R,A,gammaA,tildeA; Normalizers NR,NA; };
struct SpectralSolution { PairField amplitudes; int iterations=0; double change=0,residual=0; };
struct Boundaries { Matrix left,right,tilde_left,tilde_right; };
struct Distribution { Field x,tilde; double residual=0; };
struct KeldyshBlocks { std::vector<Matrix> lower,diagonal,upper,rhs; };
struct SparseEntry { int row,col; Complex value; };
// current = e*rho_S*I_SI/(k_B*Tc*xi_S), not amperes; normal current=conductance()*voltage.
struct CurrentResult { double voltage=0,current=0,conservation_error=0,max_spectral_residual=0,max_kinetic_residual=0; std::vector<double> probe_currents; };
struct ConvergenceResult { std::string parameter; double baseline,refined,relative_change; bool passed; };

void validate(const PhysicalParams&,const NumericalParams&,double voltage);
Matrix make_energy_matrix(double epsilon,double voltage,int NF,const PhysicalParams&);
Complex bulk_bcs_gamma(double energy,double Delta,double eta);
Matrix build_gamma_boundary_left(const Matrix& energy,const PhysicalParams&,double eta);
Matrix build_gamma_boundary_right(const Matrix& energy,double voltage,const PhysicalParams&,double eta);
Boundaries build_gamma_boundaries(const Matrix&,double voltage,const PhysicalParams&,double eta);
PairField initial_gamma_linearized(const Matrix&,const Boundaries&,const PhysicalParams&,const NumericalParams&);
// Callback is evaluated at NEGATIVE physical energies, without clipping or cyclic wrap.
Matrix tilde_gamma(double epsilon,double voltage,int NF,const PhysicalParams&,const std::function<Complex(double,double)>& gamma_at);
Normalizers compute_N1_N2(const Matrix& gamma,const Matrix& tilde);
Green compute_retarded_green_functions(const Matrix& gamma,const Matrix& tilde);
Matrix compute_Q(const Matrix& gamma,const Matrix& tilde);
Matrix compute_gamma_derivative(const Field&,size_t i,double spacing);
Matrix compute_nonlinear_gamma_source(const Matrix& derivative,const Matrix& Q);
std::vector<Complex> solve_tridiagonal_complex(std::vector<Complex> lower,std::vector<Complex> diagonal,std::vector<Complex> upper,std::vector<Complex> rhs);
PairField gamma_picard_step(const PairField&,const Matrix&,const Boundaries&,const PhysicalParams&,const NumericalParams&,double lambda);
SpectralSolution solve_gamma_for_energy(double epsilon,double voltage,const PhysicalParams&,const NumericalParams&,const PairField* initial=nullptr);
void solve_gamma_for_voltage(double voltage,const PhysicalParams&,const NumericalParams&,const std::function<void(double,const SpectralSolution&)>& consume);
Boundaries build_distribution_boundaries(const Matrix&,double voltage,const PhysicalParams&,double eta);
Matrix build_keldysh_green_function(const Matrix& gamma,const Matrix& tilde,const Green&,const Matrix& x,const Matrix& tilde_x);
KeldyshBlocks assemble_keldysh_blocks(const SpectralSolution&,const Matrix&,const PhysicalParams&,const NumericalParams&,const Boundaries&);
std::vector<SparseEntry> assemble_keldysh_sparse(const KeldyshBlocks&);
Distribution solve_distribution_x(const SpectralSolution&,double epsilon,double voltage,const PhysicalParams&,const NumericalParams&);
double compute_spectral_current(const Field& retarded,const Field& advanced,const Field& keldysh,size_t i,double dx);
double integrate_current_over_quasienergy(const std::vector<double>& integrand,double width,const PhysicalParams&);
CurrentResult solve_current_for_voltage(double voltage,const PhysicalParams&,const NumericalParams&,const std::vector<PairField>* initial=nullptr,std::vector<PairField>* solutions=nullptr);
std::vector<CurrentResult> compute_IV_curve(const std::vector<double>& voltages,const PhysicalParams&,const NumericalParams&);
std::vector<ConvergenceResult> check_current_convergence(double voltage,const PhysicalParams&,const NumericalParams&,double tolerance=1e-3);
unsigned energy_worker_count(const NumericalParams&);
std::vector<double> voltage_grid(double start,double end,double step);
std::vector<double> differentiate_current(const std::vector<double>& voltage,const std::vector<double>& current);
} // namespace sns

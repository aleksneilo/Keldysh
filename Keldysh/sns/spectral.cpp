#include "sns.hpp"
#include <limits>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace sns {
void validate(const PhysicalParams& p,const NumericalParams& n,double v) {
    for(double q:{p.Delta,p.T,p.Ksi_N,p.L_N,p.ro_N,p.area,p.Xi,v,n.eta,n.mixing,n.tolerance,n.residual_tolerance,n.kinetic_tolerance})
        if(!std::isfinite(q)) throw std::invalid_argument("parameters must be finite");
    if(p.Delta<0 || p.T<0 || p.Ksi_N<=0 || p.L_N<=0 || p.ro_N<=0 || p.area<=0 || n.eta<=0 || n.NF<0 || n.NF>1000 || n.Nx<5 || n.Neps<2 || n.max_iterations<1 || n.continuation_steps<1 || n.mixing<=0 || n.mixing>1 || n.tolerance<=0 || n.residual_tolerance<=0 || n.kinetic_tolerance<=0)
        throw std::invalid_argument("invalid physical/numerical parameter");
}
Matrix make_energy_matrix(double eps,double v,int nf,const PhysicalParams& p) {
    (void)p; // v already contains the charge and energy-unit conversion.
    if(nf<0 || nf>1000 || !std::isfinite(eps) || !std::isfinite(v)) throw std::invalid_argument("invalid energy grid");
    int k=v==0 ? 1 : 2*nf+1; Matrix e(k,k);
    for(int j=0;j<k;++j) e(j,j)=eps+2*(j-(v==0?0:nf))*v;
    return e;
}
Complex bulk_bcs_gamma(double e,double delta,double eta) {
    if(!std::isfinite(e) || !std::isfinite(delta) || !std::isfinite(eta) || delta<0 || eta<=0) throw std::invalid_argument("invalid BCS parameter");
    if(delta==0) return 0.;
    Complex z(e,eta), root=std::sqrt(Complex(delta*delta)-z*z);
    // Principal sqrt has positive real part; z+i*root is the retarded branch at both signs of E.
    return -delta/(z+Complex(0,1)*root);
}
Matrix build_gamma_boundary_left(const Matrix& e,const PhysicalParams& p,double eta) {
    Matrix b(e.rows,e.rows); for(int j=0;j<e.rows;++j) b(j,j)=bulk_bcs_gamma(e(j,j).real(),p.Delta,eta); return b;
}
Matrix build_gamma_boundary_right(const Matrix& e,double v,const PhysicalParams& p,double eta) {
    if(v<0) throw std::invalid_argument("Cuevas spectral boundary requires V>=0");
    Matrix b(e.rows,e.rows); Complex phase=std::exp(Complex(0,p.phase_radians()));
    if(v==0) { for(int j=0;j<e.rows;++j) b(j,j)=phase*bulk_bcs_gamma(e(j,j).real(),p.Delta,eta); }
    else {
        for(int j=0;j<e.rows;++j) { int m=j+1; if(m<e.rows) b(j,m)=phase*bulk_bcs_gamma(e(j,j).real()+v,p.Delta,eta); }
    }
    return b;
}
Boundaries build_gamma_boundaries(const Matrix& e,double v,const PhysicalParams& p,double eta) {
    // This low-level solver uses the positive-voltage Cuevas gauge; negative I(V) is obtained by electrode symmetry.
    if(v<0) throw std::invalid_argument("spectral solver requires V>=0");
    Boundaries b{build_gamma_boundary_left(e,p,eta),build_gamma_boundary_right(e,v,p,eta),Matrix(e.rows,e.rows),Matrix(e.rows,e.rows)};
    for(int j=0;j<e.rows;++j) {
        b.tilde_left(j,j)=std::conj(bulk_bcs_gamma(-e(j,j).real(),p.Delta,eta));
        int m=v==0 ? j : j-1;
        if(m>=0) b.tilde_right(j,m)=std::exp(Complex(0,-p.phase_radians()))*std::conj(bulk_bcs_gamma(-e(j,j).real()+v,p.Delta,eta));
    } return b;
}
Matrix tilde_gamma(double eps,double v,int nf,const PhysicalParams& p,const std::function<Complex(double,double)>& at) {
    Matrix e=make_energy_matrix(eps,v,nf,p),t(e.rows,e.rows);
    for(int j=0;j<e.rows;++j) for(int i=0;i<e.rows;++i) t(i,j)=std::conj(at(-e(i,i).real(),-e(j,j).real()));
    return t;
}
static Complex sinh_ratio(Complex a,double s) {
    if(std::abs(a)<1e-5) return s*(1.+a*a*(s*s-1.)/6.);
    if(a.real()<0) a=-a;
    if(a.real()<300) return std::sinh(a*s)/std::sinh(a);
    return std::exp(a*(s-1.))*(1.-std::exp(-2.*a*s))/(1.-std::exp(-2.*a));
}
PairField initial_gamma_linearized(const Matrix& e,const Boundaries& b,const PhysicalParams& p,const NumericalParams& n) {
    PairField f{Field(n.Nx,Matrix(e.rows,e.rows)),Field(n.Nx,Matrix(e.rows,e.rows))};
    for(int i=0;i<n.Nx;++i) { double z=double(i)/(n.Nx-1);
        for(int m=0;m<e.rows;++m) for(int j=0;j<e.rows;++j) {
            Complex a=p.L_N*std::sqrt((e(j,j)+e(m,m)+Complex(0,2*n.eta))/Complex(0,p.diffusion()));
            Complex l=sinh_ratio(a,1-z),r=sinh_ratio(a,z);
            f.gamma[i](j,m)=l*b.left(j,m)+r*b.right(j,m);
            f.tilde[i](j,m)=l*b.tilde_left(j,m)+r*b.tilde_right(j,m);
        }
    } return f;
}
Normalizers compute_N1_N2(const Matrix& g,const Matrix& t) {
    Matrix id=Matrix::identity(g.rows);
    return {solve_dense(id+g*t,id),solve_dense(id+t*g,id)};
}
Green compute_retarded_green_functions(const Matrix& g,const Matrix& t) {
    const double pi=std::acos(-1.); Matrix id=Matrix::identity(g.rows);
    Green s; s.NR=compute_N1_N2(g,t); s.gammaA=adjoint(t)*(-1.); s.tildeA=adjoint(g)*(-1.);
    s.NA=compute_N1_N2(s.gammaA,s.tildeA);
    s.R=nambu(s.NR.N1*(id-g*t),2.*(s.NR.N1*g),2.*(s.NR.N2*t),s.NR.N2*(t*g-id))*Complex(0,-pi);
    s.A=nambu((id-s.gammaA*s.tildeA)*s.NA.N1,2.*(s.gammaA*s.NA.N2),2.*(s.tildeA*s.NA.N1),(s.tildeA*s.gammaA-id)*s.NA.N2)*Complex(0,pi);
    return s;
}
Matrix compute_Q(const Matrix& g,const Matrix& t) { return solve_dense(Matrix::identity(g.rows)+t*g,t)*(-2.); }
Matrix compute_gamma_derivative(const Field& f,size_t i,double h) {
    if(f.size()<3 || i>=f.size() || h<=0) throw std::invalid_argument("invalid derivative grid");
    if(i==0) return (-3.*f[0]+4.*f[1]-f[2])*(1./(2*h));
    if(i+1==f.size()) return (3.*f[i]-4.*f[i-1]+f[i-2])*(1./(2*h));
    return (f[i+1]-f[i-1])*(1./(2*h));
}
Matrix compute_nonlinear_gamma_source(const Matrix& d,const Matrix& q) { return d*q*d; }
std::vector<Complex> solve_tridiagonal_complex(std::vector<Complex> a,std::vector<Complex> b,std::vector<Complex> c,std::vector<Complex> d) {
    size_t n=b.size(); if(n==0 || a.size()!=n || c.size()!=n || d.size()!=n) throw std::invalid_argument("tridiagonal sizes");
    double scale=0; for(size_t i=0;i<n;++i) scale=std::max(scale,std::abs(a[i])+std::abs(b[i])+std::abs(c[i]));
    for(size_t i=0;i<n;++i) { if(std::abs(b[i])<=1e-14*scale) throw std::runtime_error("Thomas zero pivot");
        if(i+1<n) { Complex q=a[i+1]/b[i]; b[i+1]-=q*c[i]; d[i+1]-=q*d[i]; }
    }
    d[n-1]/=b[n-1]; for(int i=int(n)-2;i>=0;--i) d[i]=(d[i]-c[i]*d[i+1])/b[i]; return d;
}
static Field frozen_step(const Field& f,const Field& t,const Matrix& e,const Matrix& left,const Matrix& right,const PhysicalParams& p,const NumericalParams& n,double lambda) {
    double h=p.L_N/(n.Nx-1),a=1/(h*h); int q=n.Nx-2,k=e.rows;
    Field source(n.Nx,Matrix(k,k)),out=f;
    for(int i=1;i<n.Nx-1;++i) source[i]=compute_nonlinear_gamma_source(compute_gamma_derivative(f,i,h),compute_Q(f[i],t[i]));
    for(int m=0;m<k;++m) for(int j=0;j<k;++j) {
        Complex c=(e(j,j)+e(m,m)+Complex(0,2*n.eta))/Complex(0,p.diffusion());
        std::vector<Complex> lo(q,a),diag(q,-2*a-c),up(q,a),rhs(q);
        for(int i=0;i<q;++i) rhs[i]=-lambda*source[i+1](j,m);
        rhs.front()-=a*left(j,m);rhs.back()-=a*right(j,m);lo.front()=0.;up.back()=0.;
        auto sol=solve_tridiagonal_complex(lo,diag,up,rhs);
        for(int i=0;i<q;++i) out[i+1](j,m)=sol[i];
    } out.front()=left;out.back()=right; return out;
}
PairField gamma_picard_step(const PairField& f,const Matrix& e,const Boundaries& b,const PhysicalParams& p,const NumericalParams& n,double lambda) {
    return {frozen_step(f.gamma,f.tilde,e,b.left,b.right,p,n,lambda),frozen_step(f.tilde,f.gamma,e,b.tilde_left,b.tilde_right,p,n,lambda)};
}
static double equation_residual(const PairField& f,const Matrix& e,const PhysicalParams& p,const NumericalParams& n,double lambda) {
    double h=p.L_N/(n.Nx-1),r=0.;
    for(int side=0;side<2;++side) { const Field& a=side?f.tilde:f.gamma; const Field& b=side?f.gamma:f.tilde;
        for(int i=1;i<n.Nx-1;++i) {
            Matrix second=(a[i-1]-2.*a[i]+a[i+1])*(1/(h*h));
            Matrix c=compute_nonlinear_gamma_source(compute_gamma_derivative(a,i,h),compute_Q(a[i],b[i]))*lambda;
            Matrix en=(e*a[i]+a[i]*e+Complex(0,2*n.eta)*a[i])* (1./Complex(0,p.diffusion()));
            r=std::max(r,norm(second+c-en)/(1+norm(second)+norm(c)+norm(en)));
            if(!finite(a[i]) || !std::isfinite(r)) return std::numeric_limits<double>::infinity();
        }
    } return r;
}
static void set_boundaries(PairField& f,const Boundaries& b) { f.gamma.front()=b.left;f.gamma.back()=b.right;f.tilde.front()=b.tilde_left;f.tilde.back()=b.tilde_right; }

// Write every Floquet element at every spatial point, after an accepted update.
// G and F denote the electron-electron and electron-hole blocks of G^R.
static void write_GF_iteration(std::ostream& out, const SpectralSolution& s,
    const PhysicalParams& p, double eps, double v, double lambda)
{
    const int nx = static_cast<int>(s.amplitudes.gamma.size());
    const int k = s.amplitudes.gamma.front().rows;
    const int offset = (k - 1) / 2;
    const double pi_value = std::acos(-1.0);
    const Matrix identity = Matrix::identity(k);

    for (int i = 0; i < nx; ++i) {
        const Matrix& gamma = s.amplitudes.gamma[i];
        const Matrix& tilde = s.amplitudes.tilde[i];
        const Matrix product = gamma * tilde;
        const Matrix N1 = solve_dense(identity + product, identity);
        const Matrix G = (N1 * (identity - product)) * Complex(0, -pi_value);
        const Matrix F = (N1 * gamma) * Complex(0, -2 * pi_value);
        const double coordinate = p.L_N * i / (nx - 1);


        //INPUT G(x) F(x) for each iteration
        for (int row = 0; row < k; ++row) {
            for (int col = 0; col < k; ++col) {
                out << s.iterations << '\t' << lambda << '\t'
                    << eps << '\t' << v << '\t' << coordinate << '\t'
                    << row - offset << '\t' << col - offset << '\t'
                    << G(row,col).real() << '\t' << G(row,col).imag() << '\t'
                    << F(row,col).real() << '\t' << F(row,col).imag() << '\t'
                    << s.change << '\t' << s.residual << '\n';
            }
        }
    }
    out << '\n';
    out.flush();
}

SpectralSolution solve_gamma_for_energy(double eps,double v,const PhysicalParams& p,const NumericalParams& n,const PairField* initial) {
    validate(p,n,v); Matrix e=make_energy_matrix(eps,v,n.NF,p); auto b=build_gamma_boundaries(e,v,p,n.eta);
    SpectralSolution s; s.amplitudes=initial?*initial:initial_gamma_linearized(e,b,p,n);
    if(s.amplitudes.gamma.size()!=size_t(n.Nx) || s.amplitudes.tilde.size()!=size_t(n.Nx)) throw std::invalid_argument("initial field spatial size");
    for(const auto* field:{&s.amplitudes.gamma,&s.amplitudes.tilde}) for(const auto& m:*field) if(m.rows!=e.rows || m.cols!=e.rows || !finite(m)) throw std::invalid_argument("initial Floquet shape/value");
    set_boundaries(s.amplitudes,b);
    std::ofstream iteration_output;
    if (!n.iteration_log_path.empty()) {
        iteration_output.exceptions(std::ios::failbit | std::ios::badbit);
        iteration_output.open(n.iteration_log_path, std::ios::app);
        iteration_output << std::setprecision(17)
            << "# BEGIN spectral solve: epsilon=" << eps << " v=" << v
            << " K=" << e.rows << " Nx=" << n.Nx << " Delta=" << p.Delta
            << " D=" << p.diffusion() << " L_N=" << p.L_N
            << " Xi=" << p.Xi << " eta=" << n.eta << '\n'
            << "# iteration lambda epsilon v x n m ReG ImG ReF ImF change residual\n";
    }
    int stages=initial?1:n.continuation_steps;
    for(int stage=1;stage<=stages;++stage) {
        double lambda=initial?1.:double(stage)/stages,mix=n.mixing; bool converged=false;
        double old=equation_residual(s.amplitudes,e,p,n,lambda);
        for(int iter=0;iter<n.max_iterations;++iter) {
            auto target=gamma_picard_step(s.amplitudes,e,b,p,n,lambda); PairField next=s.amplitudes;
            double change=0,r=0; bool accepted=false;
            for(int trial=0;trial<18;++trial) {
                change=0;
                for(int i=1;i<n.Nx-1;++i) for(int side=0;side<2;++side) {
                    const Matrix& a=side?s.amplitudes.tilde[i]:s.amplitudes.gamma[i];
                    const Matrix& t=side?target.tilde[i]:target.gamma[i]; Matrix& out=side?next.tilde[i]:next.gamma[i];
                    out=(1-mix)*a+mix*t; change=std::max(change,norm(out-a)/(1+norm(out)));
                }
                try { r=equation_residual(next,e,p,n,lambda); } catch(const std::runtime_error&) { r=std::numeric_limits<double>::infinity(); }
                if(std::isfinite(r) && (r<=old*1.02 || r<n.residual_tolerance)) { accepted=true;break; }
                mix*=0.5;
            }
            if(!accepted || mix<1e-10) throw std::runtime_error("Picard stalled; increase broadening/grid/continuation or use Newton");
            s.amplitudes=std::move(next);s.change=change;s.residual=r;++s.iterations;old=r;
            
            /*/std::cout
                << "eps = " << eps
                << ", stage = " << stage
                << ", iter = " << iter
                << ", change = " << change
                << ", residual = " << r
                << ", residual_tolerance = " << n.residual_tolerance
                << ", mixing = " << mix
                << std::endl; //*/
            
            mix=std::min(n.mixing,mix*1.1);
            if (iteration_output.is_open())
                write_GF_iteration(iteration_output, s, p, eps, v, lambda);
            if(change<n.tolerance && r<n.residual_tolerance) { converged=true;break; }
        }
        if(!converged) { std::ostringstream msg;msg<<"spectral iteration limit: eps/Delta="<<eps/(p.Delta?p.Delta:1)<<", residual="<<s.residual;throw std::runtime_error(msg.str()); }
    }
    if (iteration_output.is_open())
        iteration_output << "# END converged: iterations=" << s.iterations
            << " residual=" << s.residual << "\n\n";
    return s;
}
void solve_gamma_for_voltage(double v,const PhysicalParams& p,const NumericalParams& n,const std::function<void(double,const SpectralSolution&)>& consume) {
    if(v<=0) throw std::invalid_argument("quasienergy integration requires V>0; use solve_gamma_for_energy at V=0");
    PairField prev; for(int j=0;j<n.Neps;++j) { double eps=2*v*(j+0.5)/n.Neps;
        SpectralSolution s;
        try { s=solve_gamma_for_energy(eps,v,p,n,j?&prev:nullptr); }
        catch(const std::runtime_error&) { s=solve_gamma_for_energy(eps,v,p,n); }
        consume(eps,s);prev=std::move(s.amplitudes);
    }
}
} // namespace sns

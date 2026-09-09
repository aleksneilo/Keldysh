#include "sns.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cstdio>

static void require(bool ok,const char* msg) { if(!ok) throw std::runtime_error(msg); }
static void test_GF_iteration_output(double voltage) {
    using namespace sns;
    PhysicalParams p;
    NumericalParams n;
    n.NF=1; n.Nx=5;
    const double epsilon=.37*p.Delta;
    n.iteration_log_path="sns_GF_test_"+std::to_string(
        std::chrono::high_resolution_clock::now().time_since_epoch().count())+".txt";
    struct Cleanup {
        std::string path;
        ~Cleanup() { std::remove(path.c_str()); }
    } cleanup{n.iteration_log_path};

    auto s=solve_gamma_for_energy(epsilon,voltage,p,n);
    const int k=voltage==0?1:2*n.NF+1,offset=(k-1)/2;
    std::ifstream in(n.iteration_log_path);
    require(bool(in),"iteration log exists");
    size_t rows=0,final_rows=0; bool offdiagonal=false;
    std::string line;
    while(std::getline(in,line)) {
        if(line.empty() || line[0]=='#')continue;
        std::istringstream row(line);
        int iteration,fn,fm;
        double lambda,eps,v,x,gr,gi,fr,fi,change,residual;
        require(bool(row>>iteration>>lambda>>eps>>v>>x>>fn>>fm
            >>gr>>gi>>fr>>fi>>change>>residual),"iteration log columns");
        require(iteration>=1 && iteration<=s.iterations,"accepted iteration numbers");
        require(fn>=-offset && fn<=offset && fm>=-offset && fm<=offset,"all Floquet indices in log");
        require(eps==epsilon && v==voltage,"energy/voltage preserved in log");
        if(iteration==s.iterations) {
            int i=static_cast<int>(std::lround(x*(n.Nx-1)/p.L_N));
            auto green=compute_retarded_green_functions(s.amplitudes.gamma[i],s.amplitudes.tilde[i]);
            require(std::abs(Complex(gr,gi)-green.R(fn+offset,fm+offset))<1e-12,"logged G matches retarded block");
            require(std::abs(Complex(fr,fi)-green.R(fn+offset,k+fm+offset))<1e-12,"logged F matches anomalous block");
            require(std::abs(lambda-1)<1e-14,"converged step was written");
            ++final_rows;
        }
        if(fn!=fm && std::abs(Complex(fr,fi))>1e-8)offdiagonal=true;
        ++rows;
    }
    require(rows==size_t(s.iterations)*n.Nx*k*k,"every accepted iteration includes every matrix element");
    require(final_rows==size_t(n.Nx)*k*k,"final iteration not skipped");
    if(voltage!=0)require(offdiagonal,"finite voltage offdiagonal F output");
    in.close();
    solve_gamma_for_energy(epsilon,voltage,p,n,&s.amplitudes);
    std::ifstream appended(n.iteration_log_path);int begins=0;
    while(std::getline(appended,line))if(line.find("# BEGIN")==0)++begins;
    require(begins==2,"subsequent spectral solve appends instead of overwriting");
}
int main() {
    try {
        test_GF_iteration_output(0);
        test_GF_iteration_output(3.52);
        using namespace sns; PhysicalParams p; p.L_N=std::sqrt(p.diffusion()/p.Delta);NumericalParams n;n.Nx=17;n.Neps=8;n.NF=1;
        double v=p.Delta,pi=std::acos(-1.);
        auto e=make_energy_matrix(.3*p.Delta,v,1,p);require(std::abs(e(0,0).real()/p.Delta+1.7)<1e-14,"sideband energy");
        require(make_energy_matrix(p.Delta,0,10,p).rows==1,"stationary K=1");
        require(bulk_bcs_gamma(0,0,n.eta)==Complex{},"normal BCS");
        for(double en:{-100.,-2.,-.5,0.,.5,2.,100.}) {
            auto g=bulk_bcs_gamma(en*p.Delta,p.Delta,n.eta);require(std::abs(g)<=1.,"retarded branch bounded");
            require(std::abs(std::conj(bulk_bcs_gamma(-en*p.Delta,p.Delta,n.eta))+g)<1e-12,"BCS particle-hole symmetry");
        }
        auto b=build_gamma_boundaries(e,v,p,n.eta);require(b.right(0,1)!=Complex{} && b.right(1,0)==Complex{},"right shift");require(b.tilde_right(1,0)!=Complex{} && b.tilde_right(0,1)==Complex{},"tilde opposite shift");
        auto t=tilde_gamma(.3*p.Delta,v,1,p,[&](double a,double c){return Complex(a/p.Delta,2*c/p.Delta);});
        require(std::abs(t(0,2)-Complex(1.7,4.6))<1e-12,"tilde changes both energies, no transpose");
        Matrix a(2,2);a(0,0)=0.;a(1,0)=2.;a(0,1)=1.;a(1,1)=3.;require(norm(a*solve_dense(a,Matrix::identity(2))-Matrix::identity(2))<1e-12,"pivoted matrix solve");
        auto th=solve_tridiagonal_complex({0.,1.,1.},{2.,2.,2.},{1.,1.,0.},{4.,8.,8.});require(std::abs(th[0]-1.)+std::abs(th[1]-2.)+std::abs(th[2]-3.)<1e-12,"Thomas manufactured solution");
        Matrix g(2,2),tg(2,2);g(0,0)=Complex(.1,.2);g(0,1)=Complex(.2,-.1);g(1,0)=.04;tg(0,0)=Complex(-.2,.1);tg(1,0)=Complex(.1,.1);
        auto gr=compute_retarded_green_functions(g,tg);Matrix id=Matrix::identity(4);
        require(norm(gr.R*gr.R+pi*pi*id)<1e-12,"R normalization noncommuting matrices");require(norm(gr.A*gr.A+pi*pi*id)<1e-12,"A normalization");
        Matrix tau=Matrix::identity(4);tau(2,2)=-1.;tau(3,3)=-1.;require(norm(gr.A-tau*adjoint(gr.R)*tau)<1e-12,"advanced Hermitian symmetry");
        Matrix dx=Matrix::identity(2)*.3,dt=Matrix::identity(2)*(-.2);auto kg=build_keldysh_green_function(g,tg,gr,dx,dt);
        require(norm(gr.R*kg+kg*gr.A)<1e-12,"Keldysh normalization");
        auto initial=initial_gamma_linearized(e,b,p,n);require(norm(initial.gamma.front()-b.left)<1e-14 && norm(initial.gamma.back()-b.right)<1e-14,"initial exact boundaries");
        Field polynomial;for(int i=0;i<5;++i){ Matrix m(1,1);m(0,0)=double(i*i);polynomial.push_back(m); }for(int i=0;i<5;++i)require(std::abs(compute_gamma_derivative(polynomial,i,1)(0,0)-double(2*i))<1e-12,"quadratic derivative including endpoints");
        auto stationary=solve_gamma_for_energy(3*p.Delta,0,p,n);require(stationary.amplitudes.gamma[0].rows==1 && stationary.residual<n.residual_tolerance,"stationary BVP");
        auto equilibrium=solve_distribution_x(stationary,3*p.Delta,0,p,n);
        for(int i=0;i<n.Nx;++i) { auto ss=compute_retarded_green_functions(stationary.amplitudes.gamma[i],stationary.amplitudes.tilde[i]);auto kk=build_keldysh_green_function(stationary.amplitudes.gamma[i],stationary.amplitudes.tilde[i],ss,equilibrium.x[i],equilibrium.tilde[i]); require(norm(kk-(ss.R-ss.A))<.01,"equilibrium fluctuation-dissipation"); }
        auto zero=solve_gamma_for_energy(0,0,p,n);auto zeroGreen=compute_retarded_green_functions(zero.amplitudes.gamma[n.Nx/2],zero.amplitudes.tilde[n.Nx/2]);require(-zeroGreen.R(0,0).imag()/pi<.01,"stationary zero-energy minigap");
        auto positive=solve_gamma_for_energy(.37*p.Delta,v,p,n);auto reflected=solve_gamma_for_energy(-.37*p.Delta,v,p,n);
        for(int i=0;i<n.Nx;++i)for(int col=0;col<3;++col)for(int row=0;row<3;++row)require(std::abs(positive.amplitudes.tilde[i](row,col)-std::conj(reflected.amplitudes.gamma[i](2-row,2-col)))<1e-6,"finite-voltage tilde symmetry across energy sectors");
        int consumed=0;sns::solve_gamma_for_voltage(v,p,n,[&](double eps,const SpectralSolution& solution){require(eps>0 && eps<2*v && solution.residual<n.residual_tolerance,"streamed spectral solve");++consumed;});require(consumed==n.Neps,"quasienergy streaming");
        PhysicalParams phaseTest=p;phaseTest.Xi=.5;
        auto phased=build_gamma_boundaries(e,v,phaseTest,n.eta);
        require(norm(phased.right-Complex(0,1)*b.right)<1e-12,"legacy phase in units of pi");
        require(std::abs(p.diffusion()-2*pi*p.Ksi_N*p.Ksi_N)<1e-12,"legacy diffusion coefficient");
        // Changing the spatial reference rescales Ksi,L together; spectra must coincide.
        PhysicalParams scaled=p;scaled.Ksi_N*=2;scaled.L_N*=2;
        auto scaledSpectrum=solve_gamma_for_energy(.37*p.Delta,v,scaled,n);
        for(int i=0;i<n.Nx;++i)require(norm(scaledSpectrum.amplitudes.gamma[i]-positive.amplitudes.gamma[i])<1e-6,"length-unit covariance of spectral solver");
        NumericalParams regression=n;regression.NF=2;regression.Nx=25;regression.Neps=16;
        auto reference=solve_current_for_voltage(2*p.Delta,p,regression);
        require(std::abs(reference.current/(p.conductance()*p.Delta)-3.03883998768)<1e-6,"normalized current matches previous SI implementation");
        auto sweep=voltage_grid(3,.5,.1);
        require(sweep.size()==26 && sweep.front()==3 && sweep.back()==.5,"descending voltage grid");
        std::vector<double> quadratic,linear;
        for(double u:sweep){quadratic.push_back(u*u+2*u+3);linear.push_back(u);}
        auto dq=differentiate_current(sweep,quadratic),dl=differentiate_current(sweep,linear);
        for(size_t i=0;i<sweep.size();++i) {
            require(std::abs(dq[i]-(2*sweep[i]+2))<1e-11,"descending quadratic derivative including endpoints");
            require(std::abs(dl[i]-1)<1e-12,"normalized normal conductance is one");
        }
        bool badGrid=false;try{voltage_grid(3,.5,.3);}catch(const std::invalid_argument&){badGrid=true;}
        require(badGrid,"reject nonintegral sweep");
        p.Delta=0;n.eta=1e-3;
        auto normal=solve_current_for_voltage(v,p,n);require(std::abs(normal.current/(p.conductance()*v)-1)<1e-10,"normal I=GN V through full solver");require(normal.conservation_error<1e-10,"normal current conservation");
        auto negative=solve_current_for_voltage(-v,p,n);require(std::abs(negative.current/normal.current+1)<1e-12,"odd IV");
        PhysicalParams geometry=p;geometry.ro_N=2.3;geometry.area=3.7;geometry.L_N*=1.6;
        auto scaledNormal=solve_current_for_voltage(v,geometry,n);
        require(std::abs(scaledNormal.current/(geometry.area*v/(geometry.ro_N*geometry.L_N))-1)<1e-10,"normal current with nonunit rho, area and length");
        auto ns=solve_gamma_for_energy(.3*v,v,p,n);auto ne=make_energy_matrix(.3*v,v,n.NF,p);auto nb=build_distribution_boundaries(ne,v,p,n.eta);auto blocks=assemble_keldysh_blocks(ns,ne,p,n,nb);auto sparse=assemble_keldysh_sparse(blocks);require(!sparse.empty(),"sparse block export");
        NumericalParams small=n;small.max_block_bytes=1;bool rejected=false;try{assemble_keldysh_blocks(ns,ne,p,small,nb);}catch(const std::runtime_error&){rejected=true;}require(rejected,"memory guard");
        rejected=false;try{solve_current_for_voltage(0,p,n);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"zero bias is not zero-width IV");
        std::cout<<"PASS: G/F iteration output at zero and finite voltage; BCS, Floquet, tilde, matrix/Thomas solves, R/A/K normalization, boundaries, derivatives, stationary BVP, equilibrium FDT, normal Ohm law, conservation, odd IV, sparse export, validation\n";
        return 0;
    } catch(const std::exception& e) { std::cerr<<"FAIL: "<<e.what()<<'\n';return 1; }
}

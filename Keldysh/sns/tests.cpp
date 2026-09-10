#include "sns.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cstdio>
#include <thread>

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
struct CaptureAndersonLog {
    std::ostringstream text;
    std::streambuf* previous;
    CaptureAndersonLog():previous(std::cerr.rdbuf(text.rdbuf())){}
    ~CaptureAndersonLog(){std::cerr.rdbuf(previous);}
};
static void compare_spectral(const sns::SpectralSolution& a,const sns::SpectralSolution& b) {
    using namespace sns;
    double fields=0,greens=0;
    for(size_t i=0;i<a.amplitudes.gamma.size();++i) {
        fields=std::max(fields,norm(a.amplitudes.gamma[i]-b.amplitudes.gamma[i]));
        fields=std::max(fields,norm(a.amplitudes.tilde[i]-b.amplitudes.tilde[i]));
        auto ga=compute_retarded_green_functions(a.amplitudes.gamma[i],a.amplitudes.tilde[i]);
        auto gb=compute_retarded_green_functions(b.amplitudes.gamma[i],b.amplitudes.tilde[i]);
        greens=std::max(greens,norm(ga.R-gb.R)/(1+norm(ga.R)));
    }
    require(fields<3e-6 && greens<3e-6,"Anderson and Picard same paired amplitudes and GR");
    std::cout<<" differences fields="<<fields<<" relative_GR="<<greens<<'\n';
}
static void test_anderson() {
    using namespace sns;
    PhysicalParams p;p.L_N=std::sqrt(p.diffusion()/p.Delta);
    NumericalParams n;n.NF=2;n.Nx=25;
    NumericalParams aa=n;aa.use_anderson=true;aa.anderson_verbose=true;
    const double v=2*p.Delta;
    std::cout<<"Anderson comparison: epsilon/Delta Picard_iterations AA_iterations Picard_residual AA_residual\n";
    for(double u:{.01,.95,1.,1.05,2.3,3.8}) {
        auto plain=solve_gamma_for_energy(u*p.Delta,v,p,n);
        CaptureAndersonLog log;
        auto fast=solve_gamma_for_energy(u*p.Delta,v,p,aa);
        require(plain.residual<n.residual_tolerance && fast.residual<aa.residual_tolerance,"Anderson physical residual");
        std::cout<<u<<' '<<plain.iterations<<' '<<fast.iterations<<' '<<plain.residual<<' '<<fast.residual;
        compare_spectral(plain,fast);
        for(int stage=1;stage<=4;++stage) {
            std::string marker="stage="+std::to_string(stage)+" iteration=0 ";
            auto pos=log.text.str().find(marker);
            require(pos!=std::string::npos,"continuation stage logged");
            auto line=log.text.str().substr(pos,log.text.str().find('\n',pos)-pos);
            require(line.find("used=0 history=1")!=std::string::npos,"Anderson history reset per stage");
        }
        auto seeded=solve_gamma_for_energy((u+.001)*p.Delta,v,p,n,&plain.amplitudes);
        CaptureAndersonLog seedLog;
        auto seedFast=solve_gamma_for_energy((u+.001)*p.Delta,v,p,aa,&plain.amplitudes);
        compare_spectral(seeded,seedFast);
        require(seedLog.text.str().find("lambda=1 ")!=std::string::npos &&
                seedLog.text.str().find("stage=2 ")==std::string::npos,"initial guess uses one lambda=1 stage");
    }
    aa.anderson_verbose=false;
    auto plainI=solve_current_for_voltage(v,p,n);
    auto fastI=solve_current_for_voltage(v,p,aa);
    double error=std::abs(fastI.current-plainI.current)/(p.conductance()*p.Delta);
    std::cout<<"Anderson normalized current difference="<<error<<'\n';
    require(error<1e-6,"Anderson preserves current with epsilon continuation");
    // Reject every KKT combination deterministically: sum|alpha| >= |sum alpha| = 1.
    aa.anderson_depth=8;aa.mixing=.5;aa.anderson_coefficient_limit=.5;aa.anderson_verbose=true;
    NumericalParams plain=aa;plain.use_anderson=false;plain.anderson_verbose=false;
    auto expected=solve_gamma_for_energy(.95*p.Delta,v,p,plain);
    CaptureAndersonLog log;
    auto rejected=solve_gamma_for_energy(.95*p.Delta,v,p,aa);
    require(log.text.str().find("failure=coefficient-limit")!=std::string::npos,"Anderson coefficient rejection tested");
    require(expected.iterations==rejected.iterations,"rejected Anderson reproduces Picard iterations");
    compare_spectral(expected,rejected);
    aa.anderson_coefficient_limit=20;
    auto aggressive=solve_gamma_for_energy(.95*p.Delta,v,p,aa);
    compare_spectral(expected,aggressive);
    std::cout<<"Aggressive mixing=0.5 iterations Picard="<<expected.iterations<<" AA="<<aggressive.iterations
             <<" residual_fallback_seen="<<(log.text.str().find("failure=residual-fallback")!=std::string::npos)<<'\n';
    aa.mixing=.3;plain=aa;plain.use_anderson=false;plain.anderson_verbose=false;
    auto midPlain=solve_gamma_for_energy(.95*p.Delta,v,p,plain);
    auto midAA=solve_gamma_for_energy(.95*p.Delta,v,p,aa);
    compare_spectral(midPlain,midAA);
    std::cout<<"Separate mixing=0.3 Picard="<<midPlain.iterations<<" AA="<<midAA.iterations<<'\n';
    aa=n;aa.use_anderson=true;
    auto zeroPlain=solve_gamma_for_energy(.37*p.Delta,0,p,n);
    auto zeroAA=solve_gamma_for_energy(.37*p.Delta,0,p,aa);
    compare_spectral(zeroPlain,zeroAA);
    require(zeroAA.residual<n.residual_tolerance,"stationary Anderson residual");
    NumericalParams bad=n;bad.anderson_depth=0;
    bool invalid=false;try{validate(p,bad,v);}catch(const std::invalid_argument&){invalid=true;}
    require(invalid,"invalid Anderson depth rejected");
}
static void test_energy_parallelism() {
    using namespace sns;
    PhysicalParams p;p.L_N=std::sqrt(p.diffusion()/p.Delta);
    NumericalParams serial;serial.NF=2;serial.Nx=25;serial.Neps=32;
    serial.use_anderson=true;serial.residual_tolerance=1e-7;serial.kinetic_tolerance=1e-8;
    NumericalParams parallel=serial;parallel.energy_threads=0;
    unsigned hw=std::thread::hardware_concurrency();
    require(energy_worker_count(parallel)==std::min(unsigned(serial.Neps),hw>1?hw-1:1),"automatic all-but-one threads");
    std::vector<PairField> sequentialFields,parallelFields;
    auto start=std::chrono::steady_clock::now();
    auto a=solve_current_for_voltage(2*p.Delta,p,serial,nullptr,&sequentialFields);
    double serialSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    start=std::chrono::steady_clock::now();
    auto b=solve_current_for_voltage(2*p.Delta,p,parallel,nullptr,&parallelFields);
    double parallelSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    double difference=std::abs(a.current-b.current)/(p.conductance()*p.Delta),fieldError=0;
    require(difference<1e-6,"parallel normalized current matches serial");
    require(parallelFields.size()==size_t(serial.Neps),"all energy solutions saved");
    for(size_t j=0;j<parallelFields.size();++j)
        for(int i=0;i<serial.Nx;++i) {
            fieldError=std::max(fieldError,norm(sequentialFields[j].gamma[i]-parallelFields[j].gamma[i]));
            fieldError=std::max(fieldError,norm(sequentialFields[j].tilde[i]-parallelFields[j].tilde[i]));
        }
    require(fieldError<3e-6,"parallel fields retain energy ordering");
    auto repeated=solve_current_for_voltage(2*p.Delta,p,parallel);
    require(repeated.current==b.current,"parallel deterministic repeated reduction");
    parallel.energy_threads=2;
    auto warmSerial=solve_current_for_voltage(1.9*p.Delta,p,serial,&sequentialFields);
    auto warmParallel=solve_current_for_voltage(1.9*p.Delta,p,parallel,&sequentialFields);
    require(warmSerial.current==warmParallel.current,"same per-energy initial gives exact serial/parallel current");
    std::cout<<"PARALLEL workers="<<energy_worker_count(NumericalParams{parallel})
             <<" auto_workers="<<std::min(unsigned(serial.Neps),hw>1?hw-1:1)
             <<" serial_seconds="<<serialSeconds<<" auto_seconds="<<parallelSeconds
             <<" normalized_current_difference="<<difference<<" field_difference="<<fieldError<<'\n';
    // Worker exceptions return to the caller, with no partial solutions published.
    NumericalParams fail=parallel;fail.max_block_bytes=1;fail.Neps=4;
    std::vector<PairField> saved(1);
    bool caught=false;try{solve_current_for_voltage(2*p.Delta,p,fail,nullptr,&saved);}
    catch(const std::runtime_error&){caught=true;}
    require(caught && saved.size()==1,"worker exception propagated without partial publication");
    p.Delta=0;parallel.NF=1;parallel.Nx=5;parallel.Neps=4;
    parallel.iteration_log_path="sns_parallel_log_"+std::to_string(
        std::chrono::high_resolution_clock::now().time_since_epoch().count())+".txt";
    struct CleanupLog {std::string path;~CleanupLog(){std::remove(path.c_str());}} cleanup{parallel.iteration_log_path};
    auto normal=solve_current_for_voltage(2,p,parallel);
    require(std::abs(normal.current/(p.conductance()*2)-1)<1e-10,"parallel normal Ohm law");
    std::ifstream log(parallel.iteration_log_path);std::string line;
    int open=0,blocks=0;
    while(std::getline(log,line)) {
        if(line.find("# BEGIN")==0){require(open==0,"parallel log blocks do not interleave");++open;++blocks;}
        if(line.find("# END")==0){require(open==1,"parallel log block end");--open;}
    }
    require(open==0 && blocks==parallel.Neps,"all parallel energy logs complete");
    parallel.energy_threads=-1;caught=false;
    try{validate(p,parallel,2);}catch(const std::invalid_argument&){caught=true;}
    require(caught,"negative worker count rejected");
}
int main() {
    try {
        test_energy_parallelism();
        test_anderson();
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

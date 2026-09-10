#include "energy_integration.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
using namespace sns;
static void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
static EnergySample linear_sample(double x,const PairField* seed) {
    EnergySample s;s.integrand={2+x,2+x,2+x};s.spectral_iterations=seed?3:12;
    auto state=std::make_shared<PairField>();state->gamma={Matrix(1,1)};state->gamma[0](0,0)=x;
    s.amplitudes=state;return s;
}
static void tests() {
    PhysicalParams p;p.Delta=.7;
    NumericalParams n;n.adaptive_energy=true;n.energy_base_intervals=32;n.energy_threads=1;
    n.energy_use_anchors=true;n.energy_gap_skip_width=.005;n.eta=.001;
    const auto edges=energy_gap_edges(1,p,n.NF);
    auto evaluate=[&](double x,const PairField* seed,const NumericalParams&) {
        for(double edge:edges)
            require(!(x>edge-n.energy_gap_skip_width && x<edge+n.energy_gap_skip_width),"no spectral solve inside skip zone");
        return linear_sample(x,seed);
    };
    auto one=integrate_energy_adaptive(1,p,n,evaluate,nullptr,true,evaluate);
    require(std::abs(one.values[0]-6)<1e-12,"linear gap-zone integral exact");
    require(one.diagnostics.cold_starts==1 && one.diagnostics.seeded_starts>32,"anchors eliminate interval cold starts");
    require(one.diagnostics.gap_skipped_points>0 && one.diagnostics.gap_interpolated_intervals>0,"gap intervals terminal");
    for(int workers:{4,8,15,31}) {
        n.energy_threads=workers;
        auto other=integrate_energy_adaptive(1,p,n,evaluate,nullptr,true,evaluate);
        require(one.values==other.values,"optimized current exactly thread independent");
        require(one.diagnostics.total_spectral_iterations==other.diagnostics.total_spectral_iterations,"same spectral work for all threads");
        require(one.cache.entries.size()==other.cache.entries.size(),"immutable cache independent of workers");
    }
    p.Delta=0;
    auto normal=integrate_energy_adaptive(1,p,n,[](double x,const PairField* seed,const NumericalParams&){return linear_sample(x,seed);});
    require(normal.diagnostics.gap_interpolated_intervals==0,"normal state skip disabled");
    // Boundary-edge zones cannot have two real interior brackets: retain original quadrature there.
    p.Delta=1;
    auto clipped=integrate_energy_adaptive(1,p,n,[](double x,const PairField* seed,const NumericalParams&){return linear_sample(x,seed);});
    require(clipped.diagnostics.unbracketed_gap_zones>0,"boundary zones not extrapolated");
    require(std::abs(clipped.values[0]-6)<1e-12,"boundary zone integral retained");
    std::cout<<"PASS optimized gap-skip, immutable anchors, cold starts, threads 1/4/8/15/31, normal and boundary zones\n";
}
static void benchmark(int argc,char** argv) {
    std::string mode=argc>2?argv[2]:"optimized";
    int workers=argc>3?std::stoi(argv[3]):31;
    double skip=argc>4?std::stod(argv[4]):5;
    PhysicalParams p;
    NumericalParams n;n.NF=4;n.Nx=49;n.Neps=1024;n.eta=.000176;
    n.use_anderson=true;n.anderson_depth=6;n.mixing=.8;n.tolerance=1e-6;
    n.residual_tolerance=1e-6;n.kinetic_tolerance=1e-6;
    n.adaptive_energy=true;n.energy_threads=workers;
    n.energy_base_intervals=mode=="baseline"?64:32;
    n.energy_integration_tolerance=1e-3;
    n.energy_use_anchors=mode!="baseline";
    n.energy_gap_skip_width=mode=="baseline"?0:skip*n.eta;
    auto r=solve_current_for_voltage(5.28,p,n);
    const auto& d=r.energy;
    double success=double(d.successful_cold_starts+d.successful_seeded_starts);
    std::cout<<std::setprecision(16)
             <<"BENCH mode="<<mode<<" workers="<<workers<<" skip_eta="<<skip
             <<" I_star="<<r.current<<" relative_old_adaptive="<<std::abs(r.current-7.62186764839)/7.62186764839
             <<" conservation="<<r.conservation_error<<" spectral_residual="<<r.max_spectral_residual
             <<" kinetic_residual="<<r.max_kinetic_residual<<" estimated_error="<<d.estimated_error
             <<" gap_indicator="<<d.gap_interpolation_indicator
             <<" seconds="<<d.total_energy_time<<" spectral_seconds="<<d.spectral_time
             <<" kinetic_seconds="<<d.kinetic_time<<" recovery_seconds="<<d.recovery_time
             <<" spectral_solves="<<d.cold_starts+d.seeded_starts
             <<" iterations="<<d.total_spectral_iterations
             <<" average_iterations="<<d.total_spectral_iterations/std::max(1.,success)
             <<" cold="<<d.cold_starts<<" seeded="<<d.seeded_starts<<'\n';
    require(r.conservation_error<1e-3,"benchmark conservation");
    require(r.max_spectral_residual<n.residual_tolerance && r.max_kinetic_residual<n.kinetic_tolerance,"benchmark residuals");
    if(mode!="baseline")require(std::abs(r.current-7.62186764839)/7.62186764839<1e-3,"benchmark old adaptive agreement");
}
int main(int argc,char** argv) {
    try {
        if(argc>1 && std::string(argv[1])=="benchmark")benchmark(argc,argv);
        else tests();
        return 0;
    }catch(const std::exception& ex){std::cerr<<"FAIL optimization: "<<ex.what()<<'\n';return 1;}
}
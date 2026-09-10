#include "energy_integration.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
using namespace sns;
static void check(bool ok,const char* message) {if(!ok)throw std::runtime_error(message);}
static EnergySample sample(double energy,double value) {
    EnergySample s;s.integrand={value,value,value};
    auto amplitudes=std::make_shared<PairField>();
    amplitudes->gamma={Matrix(1,1)};amplitudes->gamma[0](0,0)=energy;
    s.amplitudes=amplitudes;
    return s;
}
static void manufactured() {
    PhysicalParams p;p.Delta=.7;
    NumericalParams n;n.adaptive_energy=true;n.energy_base_intervals=32;
    n.energy_gap_width=.1;n.eta=.001;n.energy_integration_tolerance=1e-7;
    n.energy_max_evaluations=262144;n.energy_threads=1;n.energy_max_refinement=16;
    auto edges=energy_gap_edges(1,p,2);
    check(std::find(edges.begin(),edges.end(),p.Delta)!=edges.end(),"left gap edge included");
    check(std::find(edges.begin(),edges.end(),1-p.Delta)!=edges.end(),"shifted right gap edge included");
    auto sharp=[&](double x,const PairField*,const NumericalParams&){return sample(x,1+1/(1+std::pow((x-p.Delta)/n.eta,2)));};
    auto a=integrate_energy_adaptive(1,p,n,sharp);
    double exact=2+n.eta*(std::atan((2-p.Delta)/n.eta)-std::atan(-p.Delta/n.eta));
    check(std::abs(a.values[0]-exact)<1e-6,"narrow gap peak integral matches analytic integral");
    double weights=0;
    for(auto& node:a.nodes) {
        weights+=node.weight;
        check(node.epsilon>0 && node.epsilon<2 && node.weight>=0,"interior positive weighted quadrature");
        check(node.gap_distance>n.gap_edge_avoidance,"quadrature avoids all gap edges");
    }
    check(std::abs(weights-2)<1e-12,"weights cover whole fundamental interval");
    for(int workers:{4,15,31}) {
        n.energy_threads=workers;
        auto b=integrate_energy_adaptive(1,p,n,sharp);
        check(a.values==b.values && a.nodes.size()==b.nodes.size(),"adaptive exact worker-count determinism");
        for(size_t i=0;i<a.nodes.size();++i)
            check(a.nodes[i].epsilon==b.nodes[i].epsilon && a.nodes[i].weight==b.nodes[i].weight,
                  "adaptive same nodes and weights for all worker counts");
    }
    std::cout<<"MANUFACTURED narrow_peak_error="<<std::abs(a.values[0]-exact)
             <<" points="<<a.nodes.size()<<" workers=1/4/15/31 exact_match=1\n";
    p.Delta=0;n.energy_threads=4;
    auto normal=integrate_energy_adaptive(1,p,n,[](double x,const PairField*,const NumericalParams&){return sample(x,3.);});
    check(energy_gap_edges(1,p,2).empty(),"no gap refinement in normal state");
    check(std::abs(normal.values[0]-6)<1e-12,"constant integral");
    n.energy_max_refinement=0;n.energy_integration_tolerance=1e-16;
    bool rejected=false;
    try{integrate_energy_adaptive(1,p,n,[](double x,const PairField*,const NumericalParams&){return sample(x,x*x);});}
    catch(const std::runtime_error&){rejected=true;}
    check(rejected,"unmet quadrature tolerance cannot silently pass");
}
static void recovery() {
    PhysicalParams p;p.Delta=0;
    NumericalParams n;n.adaptive_energy=true;n.energy_base_intervals=1;
    n.energy_max_evaluations=100000;n.energy_recovery_steps=6;n.energy_recovery_attempts=1024;
    auto reduced=integrate_energy_adaptive(1,p,n,[](double x,const PairField*,const NumericalParams& params){
        if(params.mixing>.1)throw SpectralEnergyFailure("manufactured mixing failure");
        return sample(x,1.);
    });
    check(reduced.diagnostics.recovered_points>0,"mixing recovery exercised");
    EnergyCache initial;initial.voltage=1;
    initial.entries.push_back({.1,*sample(.1,1).amplitudes});
    auto continued=integrate_energy_adaptive(1,p,n,[](double x,const PairField* seed,const NumericalParams&){
        if(!seed || std::abs(seed->gamma[0](0,0).real()-x)>.05)
            throw SpectralEnergyFailure("manufactured continuation failure");
        return sample(x,1.);
    },&initial);
    check(continued.diagnostics.recovered_points>0 && continued.values[0]==2,"intermediate energy recovery");
    p.Delta=.7;n.eta=.01;n.energy_gap_width=.1;n.energy_recovery_steps=2;
    n.energy_interpolation_max_width=.02;
    auto missing=[&](double x,const PairField*,const NumericalParams&){
        if(std::abs(x-p.Delta)<.002)throw SpectralEnergyFailure("manufactured failed gap node");
        return sample(x,1.);
    };
    auto recovered=integrate_energy_adaptive(1,p,n,missing,nullptr,true);
    check(recovered.diagnostics.interpolated_points>0,"integrand interpolation exercised");
    check(std::abs(recovered.values[0]-2)<1e-12,"interpolation preserves constant integral");
    for(auto& entry:recovered.cache.entries)check(std::abs(entry.epsilon-p.Delta)>=.002,"no fake gamma in cache");
    for(int test=0;test<3;++test) {
        NumericalParams denied=n;
        if(test==0)denied.energy_allow_interpolation=false;
        if(test==1)denied.energy_interpolation_max_width=.001;
        bool rejected=false;
        try {
            if(test==2)integrate_energy_adaptive(1,p,denied,[&](double x,const PairField*,const NumericalParams&){
                if(std::abs(x-p.Delta)<.002)throw SpectralEnergyFailure("gap");
                return sample(x,x<p.Delta?-100.:100.);
            });
            else integrate_energy_adaptive(1,p,denied,missing);
        }catch(const std::runtime_error&){rejected=true;}
        check(rejected,"disabled/unbracketed/steep interpolation rejected");
    }
    bool rejected=false;
    try{integrate_energy_adaptive(1,p,n,[](double,const PairField*,const NumericalParams&)->EnergySample{
        throw std::runtime_error("kinetic failure must propagate");
    });}catch(const std::runtime_error& ex){rejected=std::string(ex.what())=="kinetic failure must propagate";}
    check(rejected,"kinetic errors never interpolated");
    std::cout<<"RECOVERY mixing="<<reduced.diagnostics.recovered_points
             <<" continuation="<<continued.diagnostics.recovered_points
             <<" interpolated="<<recovered.diagnostics.interpolated_points
             <<" max_width="<<recovered.diagnostics.max_interpolation_width<<'\n';
}
static void physical() {
    PhysicalParams p;p.Delta=0;
    NumericalParams n;n.adaptive_energy=true;n.energy_base_intervals=32;n.energy_threads=31;
    n.NF=1;n.Nx=5;n.use_anderson=true;n.anderson_verbose=false;
    auto normal=solve_current_for_voltage(3.52,p,n);
    check(std::abs(normal.current/(p.conductance()*3.52)-1)<1e-10,"adaptive full normal Ohm law");
    check(normal.conservation_error<1e-10,"adaptive full normal current conservation");
    p.Delta=1.76;n.Nx=9;n.energy_base_intervals=8;n.energy_threads=1;
    n.eta=.01*p.Delta;n.energy_gap_width=.05*p.Delta;n.energy_integration_tolerance=2e-3;
    n.energy_max_evaluations=65536;n.mixing=.15;
    n.residual_tolerance=1e-7;n.kinetic_tolerance=1e-8;
    EnergyCache cache;
    auto a=solve_current_for_voltage(2*p.Delta,p,n,nullptr,nullptr,nullptr,&cache);
    n.energy_threads=4;
    auto b=solve_current_for_voltage(2*p.Delta,p,n);
    check(a.current==b.current && a.probe_currents==b.probe_currents,"physical adaptive exact thread determinism");
    check(a.conservation_error<1e-3,"superconducting adaptive current conservation");
    auto warm=solve_current_for_voltage(1.9*p.Delta,p,n,nullptr,nullptr,&cache);
    NumericalParams one=n;one.energy_threads=1;
    auto warmOne=solve_current_for_voltage(1.9*p.Delta,p,one,nullptr,nullptr,&cache);
    check(warm.current==warmOne.current,"coordinate-aware voltage continuation deterministic");
    NumericalParams fine=n;fine.energy_integration_tolerance*=.25;fine.energy_base_intervals*=2;
    auto c=solve_current_for_voltage(2*p.Delta,p,fine);
    NumericalParams uniform=n;uniform.adaptive_energy=false;uniform.Neps=1024;
    auto old=solve_current_for_voltage(2*p.Delta,p,uniform);
    double unit=p.conductance()*p.Delta;
    check(std::abs(a.current-c.current)/unit<.01,"physical adaptive tolerance refinement");
    std::cout<<"PHYSICAL I_adaptive="<<a.current/unit<<" I_tighter="<<c.current/unit
             <<" I_uniform1024="<<old.current/unit<<" conservation="<<a.conservation_error
             <<" estimated_Istar_error="<<a.energy.estimated_error<<" cache_size="<<cache.entries.size()<<'\n';
}
int main(int argc,char** argv) {
    try {
        std::cout<<std::setprecision(14);
        if(argc>1 && std::string(argv[1])=="physical")physical();
        else {manufactured();recovery();}
        std::cout<<"PASS adaptive energy integration\n";return 0;
    }catch(const std::exception& ex){std::cerr<<"FAIL adaptive: "<<ex.what()<<'\n';return 1;}
}
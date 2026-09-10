#include "runner.hpp"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

int sns::run_cli(int argc,char** argv,const RunSettings& settings) {
    try {
        sns::PhysicalParams p=settings.physical; sns::NumericalParams n=settings.numerical;
        std::string mode=argc>1?argv[1]:settings.mode;
        
        /*/std::cout
            << "ACTUAL MODE = " << mode
            << ", argc = " << argc
            << std::endl; /*/
        
        double scale=p.Delta;
        if(mode=="help" || mode=="--help") { std::cout<<"SNS reference solver (legacy units: E0=k_B*Tc, length=xi_S; physical defaults are set by the entry point)\n"
            <<"Keldysh spectral [E/E0="<<settings.epsilon<<"] [eV/E0="<<settings.voltage<<"] [NF="<<n.NF<<"] [Nx="<<n.Nx<<"]\n"
            <<"Keldysh iv|normal|convergence [eV/E0="<<settings.voltage<<"] [NF="<<n.NF<<"] [Nx="<<n.Nx<<"] [Neps="<<n.Neps<<"] [eta/E0="<<n.eta<<"] [adaptive|uniform]\n"
            <<"Keldysh conductance [Vstart/Delta=3] [Vend/Delta=0.5] [step/Delta=0.1] [NF] [Nx] [Neps] [output.txt] [adaptive|uniform]\n"
            <<"Keldysh benchmark [NF="<<n.NF<<"] [Nx="<<n.Nx<<"] [Neps="<<n.Neps<<"]\n"
            <<"benchmark runs v/Delta=2,1,0.5 with the configured physical parameters.\n"
            <<"Keldysh --legacy runs the original calculation. No arguments uses main.cpp settings.\n"
            <<"Input energies, voltage and eta use E0, NOT Delta or SI. Current I_star=e*rho_S*I/(E0*xi_S). Finite-bias results require convergence checks.\n";return 0; }
        if (!n.iteration_log_path.empty()) {
            std::ofstream reset(n.iteration_log_path, std::ios::trunc);
            if (!reset) throw std::runtime_error("Cannot create G/F iteration log: " + n.iteration_log_path);
            reset << "# Retarded G and F; Ghat^2=-pi^2. Energies in k_B*Tc, x in xi_S.\n";
            reset.close();
            if (!reset) throw std::runtime_error("Cannot write G/F iteration log: " + n.iteration_log_path);
            std::cerr << "G/F iteration output: " << n.iteration_log_path << '\n';
        }
        std::cout<<std::setprecision(12);
        if(mode=="spectral") {
            double e=argc>2?std::stod(argv[2]):settings.epsilon,v=argc>3?std::stod(argv[3]):settings.voltage;if(argc>4)n.NF=std::stoi(argv[4]);if(argc>5)n.Nx=std::stoi(argv[5]);
            auto s=sns::solve_gamma_for_energy(e,v,p,n);std::cout<<"iterations,residual,change\n"<<s.iterations<<','<<s.residual<<','<<s.change<<'\n';return 0;
        }
        if(mode=="conductance") {
            if(!(scale>0))throw std::invalid_argument("Conductance sweep requires Delta>0");
            double start=argc>2?std::stod(argv[2]):settings.sweep_start;
            double end=argc>3?std::stod(argv[3]):settings.sweep_end;
            double step=argc>4?std::stod(argv[4]):settings.sweep_step;
            if(argc>5)n.NF=std::stoi(argv[5]);
            if(argc>6)n.Nx=std::stoi(argv[6]);
            if(argc>7)n.Neps=std::stoi(argv[7]);
            std::string path=argc>8?argv[8]:settings.conductance_path;
            auto grid=sns::voltage_grid(start,end,step);
            if(argc>9) {
                std::string method=argv[9];
                if(method!="uniform" && method!="adaptive")throw std::invalid_argument("energy method: adaptive|uniform");
                n.adaptive_energy=method=="adaptive";
            }
            sns::validate(p,n,grid.front()*scale);
            std::ofstream out(path,std::ios::trunc);
            if(!out)throw std::runtime_error("Cannot create conductance file: "+path);
            out<<std::setprecision(16)
               <<"# V=eV_SI/Delta_SI; L_N=L/xi_S; I=eI_SI/(G_N*Delta_SI); dI/dV=G/G_N\n"
               <<"# L/xi_Delta="<<p.L_N/std::sqrt(p.diffusion()/scale)
               <<" Delta/E0="<<scale<<" T/Tc="<<p.T<<" Ksi_N="<<p.Ksi_N
               <<" ro_N="<<p.ro_N<<" area="<<p.area<<" Xi="<<p.Xi<<"\n"
               <<"# NF="<<n.NF<<" Nx="<<n.Nx<<" Neps="<<n.Neps<<" eta/Delta="<<n.eta/scale<<"\n"
               <<"# adaptive_energy="<<n.adaptive_energy<<" base_intervals="<<n.energy_base_intervals
               <<" integration_tolerance="<<n.energy_integration_tolerance
               <<" gap_width="<<n.energy_gap_width<<" min_step="<<n.energy_min_step
               <<" gap_avoidance="<<n.gap_edge_avoidance<<"\n"
               <<"# energy_use_anchors="<<n.energy_use_anchors
               <<" gap_skip_width="<<(n.energy_gap_skip_width<0?5*n.eta:n.energy_gap_skip_width)<<"\n"
               <<"# energy_workers="<<sns::energy_worker_count(n)<<"\n"
               <<"# use_anderson="<<n.use_anderson<<" depth="<<n.anderson_depth
               <<" start="<<n.anderson_start<<" regularization="<<n.anderson_regularization
               <<" coefficient_limit="<<n.anderson_coefficient_limit<<" mixing="<<n.mixing<<"\n"
               <<"# Central differences; second-order one-sided endpoints. Check grid convergence.\n"
               <<"# V\tL_N\tI\tdI/dV\n";
            out.flush();
            if(!out)throw std::runtime_error("Cannot write conductance header: "+path);
            std::vector<double> currents;
            std::vector<sns::PairField> previous,next;
            sns::EnergyCache energy_previous,energy_next;
            bool conserved=true;
            size_t written=0;
            std::cerr<<"Conductance output: "<<path<<'\n';
            for(size_t i=0;i<grid.size();++i) {
                auto r=sns::solve_current_for_voltage(grid[i]*scale,p,n,
                    n.adaptive_energy?nullptr:(previous.empty()?nullptr:&previous),
                    n.adaptive_energy?nullptr:&next,
                    energy_previous.entries.empty()?nullptr:&energy_previous,
                    n.adaptive_energy?&energy_next:nullptr);
                energy_previous=std::move(energy_next);
                previous=std::move(next);
                currents.push_back(r.current/(p.conductance()*scale));
                conserved=conserved && r.conservation_error<1e-3;
                std::cerr<<"["<<i+1<<"/"<<grid.size()<<"] V="<<grid[i]<<" I="<<currents.back()
                         <<" conservation_error="<<r.conservation_error<<'\n';
                // Save each row once its neighboring current values are available.
                if(i>=2) {
                    std::vector<double> prefix(grid.begin(),grid.begin()+i+1);
                    auto derivative=sns::differentiate_current(prefix,currents);
                    size_t limit=i+1==grid.size()?i+1:i;
                    for(;written<limit;++written)
                        out<<grid[written]<<'\t'<<p.L_N<<'\t'<<currents[written]<<'\t'<<derivative[written]<<'\n';
                    out.flush();
                    if(!out)throw std::runtime_error("Cannot write conductance data: "+path);
                }
            }
            out<<"# COMPLETE points="<<grid.size()<<" current_conservation_passed="<<conserved<<"\n";
            out.close();
            if(!out)throw std::runtime_error("Cannot finish conductance file: "+path);
            if(!conserved)std::cerr<<"Current conservation failed; refine numerical settings.\n";
            return conserved?0:2;
        }
        if(mode=="benchmark") {
            if(argc>2)n.NF=std::stoi(argv[2]);if(argc>3)n.Nx=std::stoi(argv[3]);if(argc>4)n.Neps=std::stoi(argv[4]);
            auto results=sns::compute_IV_curve({2*scale,scale,.5*scale},p,n);
            bool conserved=true;std::cout<<"eV_over_Delta,eI_over_GNDelta,current_conservation_error\n";for(auto r:results){std::cout<<r.voltage/scale<<','<<r.current/(p.conductance()*scale)<<','<<r.conservation_error<<'\n';conserved=conserved && r.conservation_error<1e-3;}return conserved?0:2;
        }
        if(mode!="iv" && mode!="normal" && mode!="convergence")throw std::invalid_argument("unknown mode; use --help");
        double v=argc>2?std::stod(argv[2]):settings.voltage;

        /*/std::cout
            << "CONVERGENCE VOLTAGE = " << v
            << std::endl;/*/
        
        
        if(argc>3)n.NF=std::stoi(argv[3]);if(argc>4)n.Nx=std::stoi(argv[4]);if(argc>5)n.Neps=std::stoi(argv[5]);if(argc>6)n.eta=std::stod(argv[6]);
        if(argc>7) {
            std::string method=argv[7];
            if(method!="uniform" && method!="adaptive")throw std::invalid_argument("energy method: adaptive|uniform");
            n.adaptive_energy=method=="adaptive";
        }
        if(mode=="normal")p.Delta=0;
        if(mode=="convergence") { auto checks=sns::check_current_convergence(v,p,n);std::cout<<"parameter,baseline_Istar,refined_Istar,relative_change,passed\n";bool passed=true;for(auto c:checks){std::cout<<c.parameter<<','<<c.baseline<<','<<c.refined<<','<<c.relative_change<<','<<c.passed<<'\n';passed=passed&&c.passed;}return passed?0:2; }
        auto r=sns::solve_current_for_voltage(v,p,n);
        std::cout<<"eV_over_E0,eI_over_GNDelta0,I_star,conservation_error,spectral_residual,kinetic_residual\n"<<v<<','<<r.current/(p.conductance()*scale)<<','<<r.current<<','<<r.conservation_error<<','<<r.max_spectral_residual<<','<<r.max_kinetic_residual<<'\n';
        if(r.conservation_error>1e-3) {std::cerr<<"Current conservation is not converged; refine Nx, NF, Neps and eta.\n";return 2;}return 0;
    }catch(const std::exception& e){std::cerr<<"SNS error: "<<e.what()<<'\n';return 1;}
}

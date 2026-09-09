#include "sns.hpp"
#include <iostream>

namespace sns {
std::vector<double> voltage_grid(double start,double end,double step) {
    if(!std::isfinite(start)||!std::isfinite(end)||!std::isfinite(step)||step<=0||start<=0||end<=0)
        throw std::invalid_argument("Sweep requires positive finite voltages and step");
    double intervals=std::abs(end-start)/step, count=std::round(intervals);
    if(count<2 || count>100000 || std::abs(intervals-count)>1e-8)
        throw std::invalid_argument("Sweep must contain at least 3 points and an integral number of steps");
    std::vector<double> grid(static_cast<size_t>(count)+1);
    double h=end>start?step:-step;
    for(size_t i=0;i<grid.size();++i)grid[i]=start+h*i;
    grid.back()=end;
    return grid;
}
std::vector<double> differentiate_current(const std::vector<double>& v,const std::vector<double>& current) {
    if(v.size()<3 || current.size()!=v.size())throw std::invalid_argument("Derivative requires at least 3 paired points");
    const double h=v[1]-v[0];
    if(!std::isfinite(h)||h==0)throw std::invalid_argument("Invalid voltage spacing");
    for(size_t i=0;i<v.size();++i) {
        if(!std::isfinite(v[i])||!std::isfinite(current[i]) ||
           (i && std::abs((v[i]-v[i-1])/h-1)>1e-8))
            throw std::invalid_argument("Derivative requires finite data on a uniform grid");
    }
    std::vector<double> d(v.size());
    d[0]=(-3*current[0]+4*current[1]-current[2])/(2*h);
    for(size_t i=1;i+1<v.size();++i)d[i]=(current[i+1]-current[i-1])/(2*h);
    size_t last=v.size()-1;
    d[last]=(3*current[last]-4*current[last-1]+current[last-2])/(2*h);
    return d;
}
double compute_spectral_current(const Field& r,const Field& a,const Field& k,size_t i,double h) {
    Matrix j=r.at(i)*compute_gamma_derivative(k,i,h)+k.at(i)*compute_gamma_derivative(a,i,h);
    Complex trace=0;int nf=j.rows/2;for(int m=0;m<nf;++m) trace+=j(m,m)-j(nf+m,nf+m);
    return trace.real();
}
double integrate_current_over_quasienergy(const std::vector<double>& values,double width,const PhysicalParams& p) {
    if(values.empty() || width<=0) throw std::invalid_argument("invalid quadrature");
    double sum=0,correction=0;for(double x:values) { if(!std::isfinite(x)) throw std::runtime_error("nonfinite current");double y=x-correction,t=sum+y;correction=(t-sum)-y;sum=t; }
    double pi=std::acos(-1.);
    // Positive terminal current for mu_R-mu_L=eV>0. With x_R=x0(E+eV)
    // the +x matrix-current trace has the opposite sign (see normal-state test).
    return -(p.area/p.ro_N)*sum*width/(values.size()*8*pi*pi);
}
CurrentResult solve_current_for_voltage(double voltage,const PhysicalParams& p,const NumericalParams& n,const std::vector<PairField>* initial,std::vector<PairField>* solutions) {
    validate(p,n,voltage);
    if(voltage==0) throw std::invalid_argument("V=0 is a stationary spectral BVP; dc voltage-state integral has zero-width zone. Use solve_gamma_for_energy.");
    if(voltage<0 && p.Xi!=0) throw std::invalid_argument("negative-voltage symmetry requires zero initial phase");
    double v=std::abs(voltage);std::vector<std::vector<double>> integrands(3);int probes[3]={(n.Nx-1)/4,(n.Nx-1)/2,3*(n.Nx-1)/4};
    CurrentResult result;result.voltage=voltage;PairField previous;
    if(solutions) solutions->clear();
    for(int j=0;j<n.Neps;++j) {
        double eps=2*v*(j+0.5)/n.Neps;const PairField* guess=nullptr;
        if(initial && initial->size()==size_t(n.Neps)) guess=&initial->at(j);else if(j) guess=&previous;
        SpectralSolution s;
        try { s=solve_gamma_for_energy(eps,v,p,n,guess); }
        catch(const std::runtime_error&) { if(!guess) throw; s=solve_gamma_for_energy(eps,v,p,n); }
        auto d=solve_distribution_x(s,eps,v,p,n);Field r,a,k;
        for(int i=0;i<n.Nx;++i) {
            auto g=compute_retarded_green_functions(s.amplitudes.gamma[i],s.amplitudes.tilde[i]);
            r.push_back(g.R);a.push_back(g.A);k.push_back(build_keldysh_green_function(s.amplitudes.gamma[i],s.amplitudes.tilde[i],g,d.x[i],d.tilde[i]));
        }
        for(int probe=0;probe<3;++probe) integrands[probe].push_back(compute_spectral_current(r,a,k,probes[probe],p.L_N/(n.Nx-1)));
        result.max_spectral_residual=std::max(result.max_spectral_residual,s.residual);result.max_kinetic_residual=std::max(result.max_kinetic_residual,d.residual);

       /*/ std::cout
            << "eps = " << eps
            << ", spectral = " << s.residual
            << ", kinetic = " << d.residual
            << ", max kinetic = "
            << result.max_kinetic_residual
            << std::endl; /*/


        if(solutions) solutions->push_back(s.amplitudes);previous=std::move(s.amplitudes);
    }
    for(auto& values:integrands) result.probe_currents.push_back((voltage>0?1.:-1.)*integrate_current_over_quasienergy(values,2*v,p));
    result.current=result.probe_currents[1];double spread=0;for(double x:result.probe_currents)spread=std::max(spread,std::abs(x-result.current));
    result.conservation_error=spread/std::max(std::abs(result.current),1e-12*p.conductance()*v);
    
    /*/std::cout << "\n=== CURRENT CONSERVATION ===\n";
    std::cout<< "V = " << voltage        << std::endl;
    std::cout<< "I(L/4) = " << result.probe_currents[0] << std::endl;
    std::cout<< "I(L/2) = " << result.probe_currents[1] << std::endl;
    std::cout<< "I(3L/4) = " << result.probe_currents[2] << std::endl;
    std::cout<< "conservation_error = " << result.conservation_error << std::endl;/*/

    return result;
}
std::vector<CurrentResult> compute_IV_curve(const std::vector<double>& voltages,const PhysicalParams& p,const NumericalParams& n) {
    std::vector<CurrentResult> out;std::vector<PairField> previous,next;
    for(double v:voltages) { out.push_back(solve_current_for_voltage(v,p,n,previous.empty()?nullptr:&previous,&next));previous=std::move(next); }return out;
}
/*/std::vector<ConvergenceResult> check_current_convergence(double v, const PhysicalParams& p, const NumericalParams& n, double tol) {
    double base=solve_current_for_voltage(v,p,n).current;std::vector<ConvergenceResult> out;
    for(int test=0;test<4;++test) { NumericalParams refined=n;std::string name;
        if(test==0) { refined.NF+=2;name="NF+2"; }
        if(test==1) { refined.Nx=2*n.Nx-1;name="half spatial spacing"; }
        if(test==2) { refined.Neps*=2;name="2*Neps"; }
        if(test==3) { refined.eta*=0.5;name="eta/2"; }
        double val=solve_current_for_voltage(v,p,refined).current;
        double error=std::abs(val-base)/std::max(std::abs(val),1e-12*p.conductance()*std::abs(v));out.push_back({name,base,val,error,error<tol});
    }return out;
}/*/
std::vector<ConvergenceResult>
check_current_convergence(
    double v,
    const PhysicalParams& p,
    const NumericalParams& n,
    double tol)
{
    std::vector<ConvergenceResult> out;


    NumericalParams refined = n;
    //refined.Neps *= 2;
    refined.NF += 2;

    std::cout
        << "\nChecking NF: "
        << n.NF
        << " -> "
        << refined.NF
        << std::endl;
    double base =
        solve_current_for_voltage(v, p, n).current;

    double val =
        solve_current_for_voltage(v, p, refined).current;

    double error =
        std::abs(val - base) /
        std::max(
            std::abs(val),
            1e-12 * p.conductance() * std::abs(v)
        );

    std::cout
        << "I_base = " << base
        << "\nI_refined = " << val
        << "\nrelative change = " << error
        << std::endl;

    /*/std::cout
        << "\nChecking Neps: "
        << n.Neps
        << " -> "
        << refined.Neps
        << std::endl;/*/

    

    out.push_back({
        "2*Neps",
        base,
        val,
        error,
        error < tol
        });

    return out;
}

} // namespace sns

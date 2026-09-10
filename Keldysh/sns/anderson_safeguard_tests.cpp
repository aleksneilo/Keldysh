// Standalone white-box test translation unit. Do not link to sns or compile into Keldysh.
#include "spectral.cpp"
static void check(bool value,const char* message) {
    if(!value)throw std::runtime_error(message);
}
int main() {
    using namespace sns;
    try {
        PhysicalParams p;NumericalParams n;n.Nx=5;n.NF=1;
        auto e=make_energy_matrix(.37*p.Delta,2*p.Delta,n.NF,p);
        auto b=build_gamma_boundaries(e,2*p.Delta,p,n.eta);
        auto x=initial_gamma_linearized(e,b,p,n);
        auto raw=gamma_picard_step(x,e,b,p,n,.25),bad=raw;
        auto vec=flatten_interior(x);
        auto copy=x;
        for(auto& z:vec)z+=Complex(.123,.456);
        set_interior_from_vector(copy,vec);
        check(flatten_interior(copy)==vec,"interior round trip");
        check(norm(copy.gamma.front()-x.gamma.front())==0 &&
              norm(copy.gamma.back()-x.gamma.back())==0 &&
              norm(copy.tilde.front()-x.tilde.front())==0 &&
              norm(copy.tilde.back()-x.tilde.back())==0,"fixed boundaries preserved");
        bool rejected=false;
        try{set_interior_from_vector(copy,{1.});}catch(const std::invalid_argument&){rejected=true;}
        check(rejected,"vector size validation");
        // Orthogonal real/imaginary residuals: alpha=(1/2,1/2).
        AndersonEntry a,c;
        a.image=flatten_interior(raw);c.image=a.image;
        a.residual.resize(a.image.size());c.residual.resize(a.image.size());
        a.residual[0]=1.;c.residual[0]=Complex(0,1);
        for(auto& z:c.image)z+=2.;
        auto combination=raw;
        check(anderson_target({a,c},n,combination)==nullptr,"KKT solve");
        auto actual=flatten_interior(combination);
        for(size_t i=0;i<actual.size();++i)
            check(std::abs(actual[i]-(a.image[i]+1.))<1e-12,"real Hermitian KKT minimizer");
        n.anderson_regularization=0;
        check(std::string(anderson_target({a,a},n,combination))=="singular","singular Gram rejected");
        n.anderson_regularization=1e-10;n.anderson_coefficient_limit=.5;
        check(std::string(anderson_target({a,c},n,combination))=="coefficient-limit","coefficient limit");
        // Force all AA trial steps to fail while the ordinary Picard step remains valid.
        bad.gamma[1](0,0)=Complex(std::numeric_limits<double>::infinity(),0);
        double old=equation_residual(x,e,p,n,.25);
        double mix=n.mixing,change=0,r=0;
        PairField next=x;bool fallback=false;
        check(safeguarded_mixed_step(x,raw,bad,true,e,p,n,.25,old,mix,next,change,r,fallback),
              "bad Anderson falls back successfully");
        check(fallback,"physical-residual fallback exercised");
        double plainMix=n.mixing,plainChange=0,plainResidual=0;
        PairField plainNext=x;bool plainFallback=false;
        check(safeguarded_mixed_step(x,raw,raw,false,e,p,n,.25,old,
              plainMix,plainNext,plainChange,plainResidual,plainFallback),"plain Picard succeeds");
        check(flatten_interior(next)==flatten_interior(plainNext) &&
              mix==plainMix && r==plainResidual && change==plainChange,"fallback restores exact Picard step");
        std::cout<<"PASS: interior packing/boundaries, real-complex KKT, singular/limit rejection, forced residual fallback identical to Picard\n";
        return 0;
    }catch(const std::exception& ex){std::cerr<<"FAIL: "<<ex.what()<<'\n';return 1;}
}
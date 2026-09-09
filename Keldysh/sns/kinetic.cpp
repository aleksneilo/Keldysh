#include "sns.hpp"

namespace sns {
static double equilibrium(double e,double t) { return t==0 ? (e>0?1.:(e<0?-1.:0.)) : std::tanh(e/(2*t)); }
Boundaries build_distribution_boundaries(const Matrix& e,double v,const PhysicalParams& p,double eta) {
    int k=e.rows; Boundaries b{Matrix(k,k),Matrix(k,k),Matrix(k,k),Matrix(k,k)};
    auto x0=[&](double energy) { return equilibrium(energy,p.T)*(1-std::norm(bulk_bcs_gamma(energy,p.Delta,eta))); };
    for(int j=0;j<k;++j) { double en=e(j,j).real();b.left(j,j)=x0(en);b.right(j,j)=x0(en+v);b.tilde_left(j,j)=x0(-en);b.tilde_right(j,j)=x0(-en+v); }
    return b;
}
Matrix build_keldysh_green_function(const Matrix& g,const Matrix& t,const Green& s,const Matrix& x,const Matrix& tx) {
    Matrix a=s.NR.N1*(x+g*tx*s.tildeA)*s.NA.N1;
    Matrix b=s.NR.N1*(x*s.gammaA-g*tx)*s.NA.N2;
    Matrix c=s.NR.N2*(t*x-tx*s.tildeA)*s.NA.N1;
    Matrix d=s.NR.N2*(tx+t*x*s.gammaA)*s.NA.N2;
    return nambu(a,b,c,d)*Complex(0,-2*std::acos(-1.));
}
static Matrix pack(const Matrix& x,const Matrix& t) { Matrix v(int(x.data.size()*2),1);std::copy(x.data.begin(),x.data.end(),v.data.begin());std::copy(t.data.begin(),t.data.end(),v.data.begin()+x.data.size());return v; }
static std::pair<Matrix,Matrix> unpack(const Matrix& v,int k) { Matrix a(k,k),b(k,k);std::copy_n(v.data.begin(),k*k,a.data.begin());std::copy_n(v.data.begin()+k*k,k*k,b.data.begin());return {a,b}; }
struct LocalKinetic {
    Matrix g,t,pg,pt,pa,pta,qr,qt,fa,fta,e;
    Green s;
    double diffusion;
};
// Action of a single neighbour/centre block, with central differences in x/xi_S.
static Matrix kinetic_action(const LocalKinetic& c,const Matrix& x,const Matrix& tx,double second,double first,bool centre) {
    Matrix rx=second*x+first*(c.pg*c.qr*x-x*c.fa*c.pta);
    Matrix rt=second*tx+first*(c.pt*c.qt*tx-tx*c.fta*c.pa);
    if(centre) {
        Matrix k22=(-2.)*(c.s.NR.N2*(tx+c.t*x*c.s.gammaA)*c.s.NA.N2);
        Matrix k11=(-2.)*(c.s.NR.N1*(x+c.g*tx*c.s.tildeA)*c.s.NA.N1);
        rx=rx-c.pg*k22*c.pta-(c.e*x-x*c.e)*(1./Complex(0,c.diffusion));
        rt=rt-c.pt*k11*c.pa-(c.e*tx-tx*c.e)*(1./Complex(0,c.diffusion));
    } return pack(rx,rt);
}
KeldyshBlocks assemble_keldysh_blocks(const SpectralSolution& spectral,const Matrix& e,const PhysicalParams& p,const NumericalParams& n,const Boundaries& bc) {
    int k=e.rows,b=2*k*k,q=n.Nx-2;
    // Stored A,B,C and elimination work: conservative upper bound, checked before allocation.
    long double bytes=static_cast<long double>(q)*b*b*sizeof(Complex)*6;
    if(bytes>n.max_block_bytes) throw std::runtime_error("kinetic dense-block memory budget exceeded; reduce NF or use sparse iterative backend");
    const auto& g=spectral.amplitudes.gamma;const auto& t=spectral.amplitudes.tilde;
    if(g.size()!=size_t(n.Nx) || t.size()!=size_t(n.Nx)) throw std::invalid_argument("spectral grid mismatch");
    std::vector<Green> green;Field ga,ta;
    for(int i=0;i<n.Nx;++i) { green.push_back(compute_retarded_green_functions(g[i],t[i]));ga.push_back(green.back().gammaA);ta.push_back(green.back().tildeA); }
    KeldyshBlocks out; double h=p.L_N/(n.Nx-1),a=1/(h*h);
    for(int i=1;i<n.Nx-1;++i) {
        LocalKinetic c{g[i],t[i],compute_gamma_derivative(g,i,h),compute_gamma_derivative(t,i,h),compute_gamma_derivative(ga,i,h),compute_gamma_derivative(ta,i,h),compute_Q(g[i],t[i]),compute_Q(t[i],g[i]),2.*(green[i].gammaA*green[i].NA.N2),2.*(green[i].tildeA*green[i].NA.N1),e,green[i],p.diffusion()};
        Matrix lo(b,b),di(b,b),up(b,b);
        for(int j=0;j<b;++j) { Matrix x(k,k),tx(k,k); if(j<k*k) x.data[j]=1.;else tx.data[j-k*k]=1.;
            Matrix l=kinetic_action(c,x,tx,a,-1/(2*h),false),d=kinetic_action(c,x,tx,-2*a,0,true),u=kinetic_action(c,x,tx,a,1/(2*h),false);
            for(int row=0;row<b;++row) { lo(row,j)=l(row,0);di(row,j)=d(row,0);up(row,j)=u(row,0); }
        }
        Matrix rhs(b,1);
        if(i==1) rhs=rhs-lo*pack(bc.left,bc.tilde_left);
        if(i==n.Nx-2) rhs=rhs-up*pack(bc.right,bc.tilde_right);
        out.lower.push_back(std::move(lo));out.diagonal.push_back(std::move(di));out.upper.push_back(std::move(up));out.rhs.push_back(std::move(rhs));
    } return out;
}
std::vector<SparseEntry> assemble_keldysh_sparse(const KeldyshBlocks& s) {
    std::vector<SparseEntry> entries;int q=int(s.diagonal.size());if(!q) return entries;int b=s.diagonal[0].rows;
    for(int i=0;i<q;++i) for(int shift=-1;shift<=1;++shift) {
        if(i+shift<0 || i+shift>=q) continue;
        const Matrix& a=shift<0?s.lower[i]:(shift>0?s.upper[i]:s.diagonal[i]);
        for(int c=0;c<b;++c) for(int r=0;r<b;++r) if(a(r,c)!=Complex{}) entries.push_back({i*b+r,(i+shift)*b+c,a(r,c)});
    } return entries;
}
static std::vector<Matrix> solve_blocks(const KeldyshBlocks& s) {
    int q=int(s.diagonal.size()),b=s.diagonal[0].rows; std::vector<Matrix> c(q),d(q);
    for(int i=0;i<q;++i) {
        Matrix pivot=s.diagonal[i],rhs=s.rhs[i];
        if(i) { pivot=pivot-s.lower[i]*c[i-1];rhs=rhs-s.lower[i]*d[i-1]; }
        Matrix combined(b,b+1);set_block(combined,0,0,s.upper[i]);for(int j=0;j<b;++j) combined(j,b)=rhs(j,0);
        Matrix sol=solve_dense(pivot,combined);c[i]=block(sol,0,0,b);d[i]=Matrix(b,1);for(int j=0;j<b;++j)d[i](j,0)=sol(j,b);
    }
    for(int i=q-2;i>=0;--i) d[i]=d[i]-c[i]*d[i+1];return d;
}
Distribution solve_distribution_x(const SpectralSolution& s,double eps,double v,const PhysicalParams& p,const NumericalParams& n) {
    validate(p,n,v);Matrix e=make_energy_matrix(eps,v,n.NF,p);auto bc=build_distribution_boundaries(e,v,p,n.eta);
    auto blocks=assemble_keldysh_blocks(s,e,p,n,bc);auto sol=solve_blocks(blocks);
    Distribution out;out.x.push_back(bc.left);out.tilde.push_back(bc.tilde_left);
    for(size_t i=0;i<sol.size();++i) {
        Matrix action=blocks.diagonal[i]*sol[i];if(i) action=action+blocks.lower[i]*sol[i-1];if(i+1<sol.size()) action=action+blocks.upper[i]*sol[i+1];
        out.residual=std::max(out.residual,norm(action-blocks.rhs[i])/(1+norm(blocks.rhs[i])));
        auto pair=unpack(sol[i],e.rows);out.x.push_back(std::move(pair.first));out.tilde.push_back(std::move(pair.second));
    }
    if(out.residual>n.kinetic_tolerance || !std::isfinite(out.residual)) throw std::runtime_error("kinetic residual exceeds tolerance");
    out.x.push_back(bc.right);out.tilde.push_back(bc.tilde_right);return out;
}
} // namespace sns

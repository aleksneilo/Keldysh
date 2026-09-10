#include "energy_integration.hpp"
#include <map>
#include <thread>
#include <atomic>
#include <limits>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <chrono>

namespace sns {

void validate_energy_parameters(const NumericalParams& n) {
    if(n.energy_base_intervals<1 || n.energy_base_intervals>100000 ||
       n.energy_max_refinement<0 || n.energy_max_refinement>24 ||
       n.energy_recovery_steps<0 || n.energy_recovery_steps>12 ||
       n.energy_max_evaluations<3*n.energy_base_intervals || n.energy_recovery_attempts<6 ||
       (!std::isfinite(n.energy_gap_skip_width) || (n.energy_gap_skip_width<0 && n.energy_gap_skip_width!=-1)) ||
       !std::isfinite(n.energy_refinement_factor) || n.energy_refinement_factor<=1 ||
       !std::isfinite(n.energy_gap_width) || n.energy_gap_width<0 ||
       !std::isfinite(n.energy_min_step) || n.energy_min_step<=0 ||
       !std::isfinite(n.gap_edge_avoidance) || n.gap_edge_avoidance<0 ||
       n.gap_edge_avoidance>n.energy_min_step/4 ||
       !std::isfinite(n.energy_integration_tolerance) || n.energy_integration_tolerance<=0 ||
       !std::isfinite(n.energy_interpolation_max_width) || n.energy_interpolation_max_width<=0 ||
       !std::isfinite(n.energy_interpolation_max_variation) || n.energy_interpolation_max_variation<=0)
        throw std::invalid_argument("invalid adaptive energy parameters (avoidance must be <= min_step/4)");
}
std::vector<double> energy_gap_edges(double v,const PhysicalParams& p,int nf) {
    std::vector<double> edges;
    if(p.Delta==0)return edges;
    // Left electrode E_n, and right electrode E_n +/- v in the existing boundaries.
    for(int j=-nf;j<=nf;++j)for(int shift=-1;shift<=1;++shift)for(int sign:{-1,1}) {
        double x=sign*p.Delta-(2*j+shift)*v;
        if(x>=0 && x<=2*v)edges.push_back(x);
    }
    std::sort(edges.begin(),edges.end());
    edges.erase(std::unique(edges.begin(),edges.end()),edges.end());
    return edges;
}
namespace {
using EnergyClock=std::chrono::steady_clock;
static double seconds_since(EnergyClock::time_point start) {
    return std::chrono::duration<double>(EnergyClock::now()-start).count();
}
static void record_success(EnergyDiagnostics& d,const EnergySample& s,bool seeded,bool recovery,double elapsed) {
    if(seeded){++d.successful_seeded_starts;d.seeded_spectral_iterations+=s.spectral_iterations;}
    else {++d.successful_cold_starts;d.cold_spectral_iterations+=s.spectral_iterations;}
    d.total_spectral_iterations+=s.spectral_iterations;
    d.max_spectral_iterations=std::max(d.max_spectral_iterations,s.spectral_iterations);
    d.spectral_time+=s.spectral_seconds;d.kinetic_time+=s.kinetic_seconds;
    if(recovery)d.recovery_time+=elapsed;
}
static void merge_statistics(EnergyDiagnostics& a,const EnergyDiagnostics& b) {
    a.cold_starts+=b.cold_starts;a.seeded_starts+=b.seeded_starts;
    a.successful_cold_starts+=b.successful_cold_starts;
    a.successful_seeded_starts+=b.successful_seeded_starts;
    a.failed_spectral_solves+=b.failed_spectral_solves;
    a.total_spectral_iterations+=b.total_spectral_iterations;
    a.cold_spectral_iterations+=b.cold_spectral_iterations;
    a.seeded_spectral_iterations+=b.seeded_spectral_iterations;
    a.max_spectral_iterations=std::max(a.max_spectral_iterations,b.max_spectral_iterations);
    a.spectral_time+=b.spectral_time;a.kinetic_time+=b.kinetic_time;a.recovery_time+=b.recovery_time;
    a.gap_skipped_points+=b.gap_skipped_points;
    a.gap_interpolated_intervals+=b.gap_interpolated_intervals;
    a.gap_interpolation_indicator+=b.gap_interpolation_indicator;
}
struct GapZone {double left,right;EnergySample lsample,rsample;bool ready=false;};
static bool in_gap_zone(double x,const std::vector<GapZone>& zones) {
    for(const auto& zone:zones)if(x>zone.left && x<zone.right)return true;
    return false;
}
static const PairField* closest_cache(double x,double v,const EnergyCache* cache) {
    if(!cache || cache->entries.empty())return nullptr;
    double old=x*cache->voltage/v;
    auto right=std::lower_bound(cache->entries.begin(),cache->entries.end(),old,
        [](const EnergyCacheEntry& entry,double value){return entry.epsilon<value;});
    if(right==cache->entries.begin())return &right->amplitudes;
    auto left=std::prev(right);
    if(right==cache->entries.end() || old-left->epsilon<=right->epsilon-old)return &left->amplitudes;
    return &right->amplitudes;
}
template<class Job> void parallel_jobs(size_t count,unsigned workers,const Job& job) {
    std::atomic<size_t> next{0};std::vector<std::exception_ptr> failures(count);
    auto run=[&] {
        for(;;) {
            size_t i=next.fetch_add(1);if(i>=count)return;
            try{job(i);}catch(...){failures[i]=std::current_exception();}
        }
    };
    std::vector<std::thread> threads;
    try{for(unsigned i=1;i<std::min<size_t>(workers,count);++i)threads.emplace_back(run);}
    catch(...){next.store(count);for(auto& thread:threads)thread.join();throw;}
    run();for(auto& thread:threads)thread.join();
    for(auto failure:failures)if(failure)std::rethrow_exception(failure);
}
struct Kahan {
    double sum=0,correction=0;
    void add(double value) { double y=value-correction,t=sum+y;correction=(t-sum)-y;sum=t; }
};
struct Stored {
    EnergySample sample;
    std::string status;
    int level=0;
    double weight=0,uncertainty=0,interpolation_width=0;
};
struct IntervalResult {
    std::map<double,Stored> points;
    size_t leaves=0;
    double error=0;
    EnergyDiagnostics statistics;
};
class IntervalIntegrator {
    double left,right,v,raw_budget,normal_density;
    const PhysicalParams& p;
    const NumericalParams& n;
    const EnergyEvaluator& evaluate;
    const std::vector<double>& gaps;
    const EnergyCache* initial;
    const EnergyCache* anchors;
    const std::vector<GapZone>& zones;
    size_t limit,attempts=0;
    int point_attempts=0;
    struct RecoveryLimit {};
    std::string last_spectral_failure;
    IntervalResult result;

    double distance(double x) const {
        double d=std::numeric_limits<double>::infinity();
        for(double g:gaps)d=std::min(d,std::abs(x-g));
        return d;
    }
    double midpoint(double a,double b) const {
        double x=a+(b-a)/2;
        if(!(x>a && x<b))throw std::runtime_error("energy interval below floating-point resolution");
        if(distance(x)>n.gap_edge_avoidance)return x;
        // Usually gaps are panel boundaries. Handle very close distinct edges explicitly.
        std::vector<double> candidates{a+(b-a)/4,a+3*(b-a)/4};
        for(double g:gaps) {
            candidates.push_back(std::nextafter(g-n.gap_edge_avoidance,a));
            candidates.push_back(std::nextafter(g+n.gap_edge_avoidance,b));
        }
        for(double candidate:candidates)
            if(candidate>a && candidate<b && distance(candidate)>n.gap_edge_avoidance)return candidate;
        throw std::runtime_error("no quadrature node outside gap exclusion; reduce gap_edge_avoidance");
    }
    std::vector<std::pair<double,const PairField*>> neighbors(double x) const {
        std::vector<std::pair<double,const PairField*>> found;
        auto upper=result.points.lower_bound(x);
        for(auto i=upper;i!=result.points.end();++i)
            if(i->second.sample.amplitudes) {found.push_back({i->first,i->second.sample.amplitudes.get()});break;}
        for(auto i=upper;i!=result.points.begin();) {
            --i;
            if(i->second.sample.amplitudes) {found.push_back({i->first,i->second.sample.amplitudes.get()});break;}
        }
        // Prior voltage's energies have explicit coordinates, unlike the old bare PairField vector.
        if(initial && initial->voltage>0 && !initial->entries.empty()) {
            double old_x=x*initial->voltage/v;
            auto i=std::lower_bound(initial->entries.begin(),initial->entries.end(),old_x,
                [](const EnergyCacheEntry& a,double value){return a.epsilon<value;});
            const EnergyCacheEntry* closest=nullptr;
            if(i!=initial->entries.end())closest=&*i;
            if(i!=initial->entries.begin()) {
                auto previous=std::prev(i);
                if(!closest || old_x-previous->epsilon<=closest->epsilon-old_x)closest=&*previous;
            }
            if(closest)found.push_back({closest->epsilon*v/initial->voltage,&closest->amplitudes});
        }
        if(anchors && !anchors->entries.empty()) {
            const PairField* closest=closest_cache(x,v,anchors);
            for(const auto& entry:anchors->entries)
                if(&entry.amplitudes==closest){found.push_back({entry.epsilon,closest});break;}
        }
        // Shared boundary states are immutable throughout parallel integration.
        for(const auto& zone:zones)if(zone.ready) {
            found.push_back({zone.left,zone.lsample.amplitudes.get()});
            found.push_back({zone.right,zone.rsample.amplitudes.get()});
        }
        found.erase(std::remove_if(found.begin(),found.end(),[](const auto& a){return !a.second;}),found.end());
        std::stable_sort(found.begin(),found.end(),[x](const auto& a,const auto& b) {
            double da=std::abs(a.first-x),db=std::abs(b.first-x);
            return da==db?a.first<b.first:da<db;
        });
        std::vector<std::pair<double,const PairField*>> nearest;
        bool have_left=false,have_right=false;
        for(const auto& candidate:found) {
            bool left_side=candidate.first<=x;
            if((left_side && !have_left) || (!left_side && !have_right)) {
                nearest.push_back(candidate);
                if(left_side)have_left=true;else have_right=true;
            }
            if(have_left && have_right)break;
        }
        return nearest;
    }
    bool attempt(double x,const PairField* seed,const NumericalParams& parameters,
                 int level,const char* status) {
        if(in_gap_zone(x,zones))return false;
        if(++point_attempts>n.energy_recovery_attempts)throw RecoveryLimit{};
        if(++attempts>limit)throw std::runtime_error("adaptive energy evaluation budget exhausted in interval");
        auto started=EnergyClock::now();
        if(seed)++result.statistics.seeded_starts;else ++result.statistics.cold_starts;
        try {
            auto sample=evaluate(x,seed,parameters);
            record_success(result.statistics,sample,seed!=nullptr,std::string(status)=="continuation",seconds_since(started));
            for(double value:sample.integrand)
                if(!std::isfinite(value))throw std::runtime_error("nonfinite spectral current");
            if(!std::isfinite(sample.spectral_residual)||!std::isfinite(sample.kinetic_residual))
                throw std::runtime_error("nonfinite energy residual");
            result.points.emplace(x,Stored{std::move(sample),status,level,0,0,0});
            return true;
        }catch(const SpectralEnergyFailure& error) {
            last_spectral_failure=error.what();++result.statistics.failed_spectral_solves;
            result.statistics.spectral_time+=error.spectral_seconds;
            if(std::string(status)=="continuation")result.statistics.recovery_time+=seconds_since(started);
            return false;
        }
    }
    bool solve(double x,int level,bool bridge=false) {
        if(result.points.count(x))return result.points.at(x).sample.amplitudes!=nullptr;
        auto nearby=neighbors(x);
        const PairField* closest=nearby.empty()?nullptr:nearby.front().second;
        const char* direct=bridge?"continuation":(level?"refined":"direct");
        if(attempt(x,closest,n,level,direct))return true;
        for(size_t i=1;i<nearby.size();++i)
            if(attempt(x,nearby[i].second,n,level,"continuation"))return true;
        if(closest && attempt(x,nullptr,n,level,"continuation"))return true;
        NumericalParams reduced=n;
        reduced.mixing=std::min(.15,n.mixing*.5);
        if(attempt(x,closest,reduced,level,"continuation"))return true;
        reduced.continuation_steps=std::max(8,2*n.continuation_steps);
        // With initial != nullptr the existing solver uses a single lambda=1 stage.
        // Therefore the extra homotopy stages must be tried with a fresh initial field.
        return attempt(x,nullptr,reduced,level,"continuation");
    }
    bool bridge(double from,double target,int level,int depth) {
        if(depth==0 || std::abs(target-from)<=n.energy_min_step)return false;
        double middle=midpoint(std::min(from,target),std::max(from,target));
        bool ok=solve(middle,level,true);
        if(!ok)ok=bridge(from,middle,level,depth-1);
        if(!ok)return false;
        if(solve(target,level,true))return true;
        return bridge(middle,target,level,depth-1);
    }
    Stored& point(double x,int level) {
        auto existing=result.points.find(x);
        if(existing!=result.points.end())return existing->second;
        point_attempts=0;
        try {
        if(solve(x,level))return result.points.at(x);
        auto nearby=neighbors(x);
        for(const auto& seed:nearby)
            if(seed.first>=left && seed.first<=right && seed.first!=x &&
               bridge(seed.first,x,level,n.energy_recovery_steps))return result.points.at(x);
        }catch(const RecoveryLimit&) {}
        // Reserve a bounded final attempt to obtain both interpolation neighbors.
        point_attempts=0;
        // Supply real bracketing solutions before considering interpolation. No fake gamma.
        double span=std::min(n.energy_interpolation_max_width,
                            1.8*std::min(x-left,right-x));
        double a=x-span/2,b=x+span/2;
        bool bracket=span>2*n.gap_edge_avoidance && a<x && b>x &&
                     distance(a)>n.gap_edge_avoidance && distance(b)>n.gap_edge_avoidance;
        if(bracket) {
            bool haveLeft=false,haveRight=false;
            try {haveLeft=solve(a,level,true);haveRight=solve(b,level,true);}
            catch(const RecoveryLimit&) {}
            // Retry target once from the new neighbors, within the remaining local budget.
            try {
                if((haveLeft || haveRight) && solve(x,level,true))return result.points.at(x);
            }catch(const RecoveryLimit&) {}
            bracket=haveLeft && haveRight;
        }
        if(n.energy_allow_interpolation && bracket &&
           distance(x)<=std::max(4*n.gap_edge_avoidance,2*n.eta)) {
            const auto& l=result.points.at(a).sample.integrand;
            const auto& r=result.points.at(b).sample.integrand;
            double variation=0,scale=normal_density;
            for(int q=0;q<3;++q) {
                variation=std::max(variation,std::abs(l[q]-r[q]));
                scale=std::max({scale,std::abs(l[q]),std::abs(r[q])});
            }
            if(b-a<=n.energy_interpolation_max_width*(1+1e-12) &&
               variation<=n.energy_interpolation_max_variation*scale) {
                EnergySample sample;
                for(int q=0;q<3;++q)sample.integrand[q]=l[q]+(r[q]-l[q])*(x-a)/(b-a);
                auto inserted=result.points.emplace(x,Stored{sample,"interpolated",level,0,variation,b-a});
                return inserted.first->second;
            }
        }
        std::ostringstream message;
        message<<std::setprecision(17)<<"failed adaptive energy point: epsilon="<<x
               <<" eps/Delta="<<x/(p.Delta?p.Delta:1)<<" gap_distance="<<distance(x)
               <<" level="<<level<<"; recovery exhausted, interpolation unavailable or unsafe; last spectral error: "<<last_spectral_failure;
        throw std::runtime_error(message.str());
    }
    void panel(double a,double b,int level) {
        for(const auto& zone:zones)if(zone.ready && a>=zone.left && b<=zone.right) {
            const double width=zone.right-zone.left;
            double wr=(b-a)*((a-zone.left)+(b-zone.left))/(2*width),wl=(b-a)-wr;
            auto l=result.points.emplace(zone.left,Stored{zone.lsample,"gap-boundary",0,0,0,width}).first;
            auto r=result.points.emplace(zone.right,Stored{zone.rsample,"gap-boundary",0,0,0,width}).first;
            l->second.weight+=wl;r->second.weight+=wr;
            double variation=0;
            for(int q=0;q<3;++q)variation=std::max(variation,
                std::abs(zone.lsample.integrand[q]-zone.rsample.integrand[q]));
            result.statistics.gap_interpolation_indicator+=(b-a)*variation/2;
            result.statistics.gap_skipped_points+=3; // midpoint + two half-midpoints avoided
            ++result.statistics.gap_interpolated_intervals;++result.leaves;
            return;
        }
        const double middle=midpoint(a,b);
        // Monotone l -> m -> r shortens jumps without changing the quadrature formula.
        double l=midpoint(a,middle),r=midpoint(middle,b);
        auto& sl=point(l,level);
        auto& sm=point(middle,level);
        auto& sr=point(r,level);
        double error=0;
    EnergyDiagnostics statistics;
        for(int q=0;q<3;++q) {
            double fine=(middle-a)*sl.sample.integrand[q]+(b-middle)*sr.sample.integrand[q];
            double coarse=(b-a)*sm.sample.integrand[q];
            error=std::max(error,std::abs(fine-coarse)/3);
        }
        // This conservative variation term prevents unreported interpolation uncertainty.
        error+=(middle-a)*sl.uncertainty+(b-middle)*sr.uncertainty+(b-a)*sm.uncertainty;
        double tolerance=raw_budget*(b-a)/(2*v);
        if(error<=tolerance) {
            sl.weight+=middle-a;sr.weight+=b-middle;
            result.error+=error;++result.leaves;return;
        }
        if(level>=n.energy_max_refinement || (b-a)/4<n.energy_min_step) {
            std::ostringstream message;
            message<<"energy quadrature tolerance not reached on ["<<a<<','<<b
                   <<"]: estimate="<<error<<" budget="<<tolerance
                   <<"; increase energy_max_refinement or lower energy_min_step";
            throw std::runtime_error(message.str());
        }
        panel(a,middle,level+1);panel(middle,b,level+1);
    }
public:
    IntervalIntegrator(double a,double b,double voltage,const PhysicalParams& physical,
        const NumericalParams& numerical,const EnergyEvaluator& evaluator,const std::vector<double>& edges,
        const EnergyCache* seeds,size_t budget,const EnergyCache* anchor_cache,const std::vector<GapZone>& gap_zones):
        left(a),right(b),v(voltage),
        raw_budget(numerical.energy_integration_tolerance*8*std::acos(-1.)*std::acos(-1.)*voltage/physical.L_N),
        normal_density(4*std::acos(-1.)*std::acos(-1.)/physical.L_N),
        p(physical),n(numerical),evaluate(evaluator),gaps(edges),initial(seeds),anchors(anchor_cache),zones(gap_zones),limit(budget) {}
    IntervalResult boundary(double x) { if(!point(x,0).sample.amplitudes)throw std::runtime_error("gap boundary requires a real spectral solution");return std::move(result); }
    IntervalResult run() {
        std::vector<double> cuts{left,right};
        auto add=[&](double x){if(x>left && x<right && !in_gap_zone(x,zones))cuts.push_back(x);};
        for(const auto& zone:zones){add(zone.left);add(zone.right);}
        add(v); // zero-temperature distribution step of the shifted electrode
        for(double g:gaps) {
            add(g);
            double width=n.energy_gap_width;
            for(int level=0;level<=n.energy_max_refinement && width>=n.energy_min_step; ++level) {
                add(g-width);add(g+width);
                if(width<=std::max(n.eta,4*n.energy_min_step))break;
                width/=n.energy_refinement_factor;
            }
        }
        std::sort(cuts.begin(),cuts.end());
        const double roundoff=64*std::numeric_limits<double>::epsilon()*std::max(1.,v);
        cuts.erase(std::unique(cuts.begin(),cuts.end(),
            [roundoff](double a,double b){return b-a<=roundoff;}),cuts.end());
        cuts.front()=left;cuts.back()=right;
        for(size_t i=1;i<cuts.size();++i)panel(cuts[i-1],cuts[i],0);
        return std::move(result);
    }
};
} // namespace

EnergyIntegral integrate_energy_adaptive(double v,const PhysicalParams& p,const NumericalParams& n,
    const EnergyEvaluator& evaluate,const EnergyCache* initial,bool keep_solutions,const EnergyEvaluator& anchor_evaluate) {
    const auto total_started=EnergyClock::now();
    validate(p,n,v);validate_energy_parameters(n);
    if(!(v>0))throw std::invalid_argument("adaptive quadrature requires v>0");
    if(initial) {
        double last=-1;
        if(!std::isfinite(initial->voltage)||initial->voltage<=0)throw std::invalid_argument("invalid energy cache voltage");
        for(const auto& entry:initial->entries) {
            if(!std::isfinite(entry.epsilon)||entry.epsilon<=last ||
               entry.epsilon<=0 || entry.epsilon>=2*initial->voltage)
                throw std::invalid_argument("energy cache must be sorted with interior energies");
            last=entry.epsilon;
        }
    }
    int count=n.energy_base_intervals;
    auto gaps=energy_gap_edges(v,p,n.NF);
    std::vector<IntervalResult> intervals(count);

    NumericalParams scheduling=n;scheduling.Neps=count>1?count:2;
    unsigned workers=std::min(unsigned(count),energy_worker_count(scheduling));
    std::cerr<<"Energy integration (adaptive): intervals="<<count<<" workers="<<workers
             <<" tolerance="<<n.energy_integration_tolerance<<" eta="<<n.eta<<'\n';
    const double skip=n.energy_gap_skip_width<0?5*n.eta:n.energy_gap_skip_width;
    std::vector<GapZone> zones;
    size_t unbracketed=0;
    if(skip>0)for(double gap:gaps) {
        if(gap-skip<=0 || gap+skip>=2*v){++unbracketed;continue;}
        if(!zones.empty() && gap-skip<=zones.back().right)zones.back().right=gap+skip;
        else zones.push_back(GapZone{gap-skip,gap+skip,EnergySample{},EnergySample{},false});
    }
    EnergyCache anchors;anchors.voltage=v;
    EnergyDiagnostics anchor_stats;
    if(n.energy_use_anchors) {
        const auto& anchor_solver=anchor_evaluate?anchor_evaluate:evaluate;
        for(int j=0;j<count;++j) {
            double a=2*v*j/count,b=2*v*(j+1)/count;
            double x=a+(b-a)/2;
            if(in_gap_zone(x,zones)) {
                double candidates[2]={a+(b-a)/4,a+3*(b-a)/4};
                bool available=false;
                for(double candidate:candidates)if(!in_gap_zone(candidate,zones)){x=candidate;available=true;break;}
                if(!available)continue;
            }
            double nearest=std::numeric_limits<double>::infinity();
            for(double gap:gaps)nearest=std::min(nearest,std::abs(x-gap));
            if(nearest<=std::max(n.gap_edge_avoidance,skip))continue;
            const PairField* seed=closest_cache(x,v,&anchors);
            if(!seed)seed=closest_cache(x,v,initial);
            if(seed)++anchor_stats.seeded_starts;else ++anchor_stats.cold_starts;
            auto start=EnergyClock::now();
            try {
                auto sample=anchor_solver(x,seed,n);
                record_success(anchor_stats,sample,seed!=nullptr,false,seconds_since(start));
                if(sample.amplitudes)anchors.entries.push_back({x,*sample.amplitudes});
            }catch(const SpectralEnergyFailure& error) {
                ++anchor_stats.failed_spectral_solves;anchor_stats.spectral_time+=error.spectral_seconds;
                // No mutable cross-worker fallback: a missing anchor simply remains absent.
            }
        }
    }
    // Evaluate zone boundaries concurrently from the same immutable anchor snapshot.
    std::vector<IntervalResult> boundary_results(2*zones.size());
    parallel_jobs(boundary_results.size(),workers,[&](size_t i) {
        double x=(i%2)?zones[i/2].right:zones[i/2].left;
        double direction=(i%2)?1.:-1.;
        EnergyCache approach=anchors;
        EnergyDiagnostics approach_stats;
        const auto& anchor_solver=anchor_evaluate?anchor_evaluate:evaluate;
        std::vector<double> positions;
        double distance=n.energy_gap_width;
        for(int level=0;level<n.energy_max_refinement && distance>std::max(skip,n.energy_min_step);++level) {
            double candidate=x+direction*distance;
            if(candidate>0 && candidate<2*v && !in_gap_zone(candidate,zones))positions.push_back(candidate);
            distance/=n.energy_refinement_factor;
        }
        positions.push_back(x);
        for(double candidate:positions) {
            const PairField* seed=closest_cache(candidate,v,&approach);
            if(!seed)seed=closest_cache(candidate,v,initial);
            if(seed)++approach_stats.seeded_starts;else ++approach_stats.cold_starts;
            auto started=EnergyClock::now();
            try {
                auto sample=anchor_solver(candidate,seed,n);
                record_success(approach_stats,sample,seed!=nullptr,false,seconds_since(started));
                if(sample.amplitudes) {
                    auto where=std::lower_bound(approach.entries.begin(),approach.entries.end(),candidate,
                        [](const auto& a,double e){return a.epsilon<e;});
                    if(where==approach.entries.end() || where->epsilon!=candidate)
                        approach.entries.insert(where,{candidate,*sample.amplitudes});
                }
            }catch(const SpectralEnergyFailure& error) {
                ++approach_stats.failed_spectral_solves;approach_stats.spectral_time+=error.spectral_seconds;
                // The established recovery below has access to all successful approach anchors.
                break;
            }
        }
        IntervalIntegrator boundary_task(0,2*v,v,p,n,evaluate,gaps,initial,
            size_t(n.energy_recovery_attempts)*3,&approach,zones);
        boundary_results[i]=boundary_task.boundary(x);
        merge_statistics(boundary_results[i].statistics,approach_stats);
    });
    for(size_t i=0;i<zones.size();++i) {
        zones[i].lsample=boundary_results[2*i].points.at(zones[i].left).sample;
        zones[i].rsample=boundary_results[2*i+1].points.at(zones[i].right).sample;
        zones[i].ready=true;
    }
    parallel_jobs(size_t(count),workers,[&](size_t j) {
        size_t budget=n.energy_max_evaluations/count+(j<size_t(n.energy_max_evaluations%count)?1:0);
        IntervalIntegrator task(2*v*j/count,2*v*(j+1)/count,v,p,n,evaluate,gaps,initial,budget,&anchors,zones);
        intervals[j]=task.run();
    });
    EnergyIntegral out;out.cache.voltage=v;
    out.diagnostics.intervals=count;out.diagnostics.workers=workers;
    out.diagnostics.unbracketed_gap_zones=unbracketed;
    merge_statistics(out.diagnostics,anchor_stats);
    // Boundary evaluations have zero initial weight and are counted once, even if a zone crosses jobs.
    intervals.insert(intervals.end(),std::make_move_iterator(boundary_results.begin()),std::make_move_iterator(boundary_results.end()));
    struct Contribution {double epsilon;std::array<double,3> value;};
    std::vector<Contribution> contributions;
    std::array<Kahan,3> sums;Kahan error;
    for(auto& interval:intervals) {
        merge_statistics(out.diagnostics,interval.statistics);
        out.diagnostics.leaf_intervals+=interval.leaves;error.add(interval.error);
        for(auto& item:interval.points) {
            double x=item.first;auto& sample=item.second;
            double distance=std::numeric_limits<double>::infinity();
            for(double g:gaps)distance=std::min(distance,std::abs(x-g));
            out.nodes.push_back({x,sample.weight,distance,sample.sample.spectral_residual,sample.level,sample.status});
            if(sample.status=="gap-boundary") {
                out.diagnostics.max_interpolation_width=std::max(out.diagnostics.max_interpolation_width,sample.interpolation_width);
            }else if(sample.status=="interpolated") {
                ++out.diagnostics.interpolated_points;
                out.diagnostics.max_interpolation_width=std::max(out.diagnostics.max_interpolation_width,sample.interpolation_width);
            }else if(sample.status=="continuation")++out.diagnostics.recovered_points;
            else ++out.diagnostics.direct_points;
            if(sample.level>0)++out.diagnostics.refined_points;
            if(sample.weight>0)++out.diagnostics.quadrature_points;
            Contribution contribution{x,{}};
            for(int q=0;q<3;++q)contribution.value[q]=sample.weight*sample.sample.integrand[q];
            contributions.push_back(contribution);
            out.spectral_residual=std::max(out.spectral_residual,sample.sample.spectral_residual);
            out.kinetic_residual=std::max(out.kinetic_residual,sample.sample.kinetic_residual);
            if(keep_solutions && sample.weight>0 && sample.sample.amplitudes)
                out.cache.entries.push_back({x,*sample.sample.amplitudes});
        }
    }
    std::stable_sort(contributions.begin(),contributions.end(),[](const auto& a,const auto& b){return a.epsilon<b.epsilon;});
    for(const auto& contribution:contributions)for(int q=0;q<3;++q)sums[q].add(contribution.value[q]);
    std::stable_sort(out.nodes.begin(),out.nodes.end(),[](const auto& a,const auto& b){return a.epsilon<b.epsilon;});
    std::sort(out.cache.entries.begin(),out.cache.entries.end(),[](const auto& a,const auto& b){return a.epsilon<b.epsilon;});
    out.cache.entries.erase(std::unique(out.cache.entries.begin(),out.cache.entries.end(),
        [](const auto& a,const auto& b){return a.epsilon==b.epsilon;}),out.cache.entries.end());
    for(int q=0;q<3;++q)out.values[q]=sums[q].sum;
    out.diagnostics.total_energy_time=seconds_since(total_started);
    out.diagnostics.estimated_error=error.sum; // midpoint estimate outside planned gap zones; wrapper converts units
    if(n.energy_verbose)for(const auto& node:out.nodes)
        std::cerr<<std::setprecision(17)<<"Energy epsilon="<<node.epsilon<<" weight="<<node.weight
                 <<" nearest_gap_distance="<<node.gap_distance<<" refinement_level="<<node.level
                 <<" spectral_residual="<<node.spectral_residual<<" status="<<node.status<<'\n';
    return out;
}

CurrentResult solve_current_adaptive(double voltage,const PhysicalParams& p,const NumericalParams& n,
    const EnergyCache* initial,EnergyCache* solutions) {
    double v=std::abs(voltage);
    auto compute=[&](double eps,const PairField* seed,const NumericalParams& parameters,bool anchor_only) {
        SpectralSolution spectral;
        auto spectral_start=EnergyClock::now();
        try {spectral=solve_gamma_for_energy(eps,v,p,parameters,seed);}
        catch(const std::runtime_error& ex) {
            std::string message=ex.what();
            if(message.find("Picard stalled")==0 || message.find("spectral iteration limit")==0 ||
               message.find("singular/nonfinite linear system")==0 || message.find("nonfinite linear solution")==0)
                throw SpectralEnergyFailure(message,seconds_since(spectral_start));
            throw;
        }
        EnergySample sample;
        sample.spectral_seconds=seconds_since(spectral_start);
        sample.spectral_iterations=static_cast<size_t>(spectral.iterations);
        sample.spectral_residual=spectral.residual;
        if(!anchor_only) {
            auto kinetic_start=EnergyClock::now();
            // Kinetic failures are never eligible for interpolation.
            auto distribution=solve_distribution_x(spectral,eps,v,p,parameters);
            Field r,a,k;
            for(int i=0;i<n.Nx;++i) {
                auto green=compute_retarded_green_functions(spectral.amplitudes.gamma[i],spectral.amplitudes.tilde[i]);
                r.push_back(green.R);a.push_back(green.A);
                k.push_back(build_keldysh_green_function(spectral.amplitudes.gamma[i],spectral.amplitudes.tilde[i],
                    green,distribution.x[i],distribution.tilde[i]));
            }
            int probes[3]={(n.Nx-1)/4,(n.Nx-1)/2,3*(n.Nx-1)/4};
            for(int q=0;q<3;++q)sample.integrand[q]=compute_spectral_current(r,a,k,probes[q],p.L_N/(n.Nx-1));
            sample.kinetic_residual=distribution.residual;
            sample.kinetic_seconds=seconds_since(kinetic_start);
        }
        sample.amplitudes=std::make_shared<PairField>(std::move(spectral.amplitudes));
        return sample;
    };
    EnergyEvaluator evaluate=[&](double x,const PairField* seed,const NumericalParams& params){return compute(x,seed,params,false);};
    EnergyEvaluator anchor=[&](double x,const PairField* seed,const NumericalParams& params){return compute(x,seed,params,true);};
    auto integral=integrate_energy_adaptive(v,p,n,evaluate,initial,solutions!=nullptr,anchor);
    double factor=-(p.area/p.ro_N)/(8*std::acos(-1.)*std::acos(-1.));
    CurrentResult result;result.voltage=voltage;result.energy=integral.diagnostics;
    result.energy.estimated_error*=std::abs(factor);
    result.energy.gap_interpolation_indicator*=std::abs(factor);
    result.max_spectral_residual=integral.spectral_residual;
    result.max_kinetic_residual=integral.kinetic_residual;
    for(double value:integral.values)result.probe_currents.push_back((voltage>0?1.:-1.)*factor*value);
    result.current=result.probe_currents[1];
    double spread=0;
    for(double value:result.probe_currents)spread=std::max(spread,std::abs(value-result.current));
    result.conservation_error=spread/std::max(std::abs(result.current),1e-12*p.conductance()*v);
    if(solutions)*solutions=std::move(integral.cache);
    const auto& d=result.energy;
    std::cerr<<"Energy integration: intervals="<<d.intervals<<" leaf_intervals="<<d.leaf_intervals
             <<" direct_points="<<d.direct_points<<" refined_points="<<d.refined_points
             <<" recovered_points="<<d.recovered_points<<" interpolated_points="<<d.interpolated_points
             <<" quadrature_points="<<d.quadrature_points<<" max_interpolation_width="<<d.max_interpolation_width
             <<" estimated_Istar_error="<<d.estimated_error<<" workers="<<d.workers
             <<" estimated_error_scope="<<(d.gap_interpolated_intervals?"outside_gap_skip_zones":"full_interval")
             <<" gap_skipped_points="<<d.gap_skipped_points<<" gap_interpolated_intervals="<<d.gap_interpolated_intervals
             <<" gap_interpolation_indicator="<<d.gap_interpolation_indicator
             <<" unbracketed_gap_zones="<<d.unbracketed_gap_zones
             <<" cold_starts="<<d.cold_starts<<" seeded_starts="<<d.seeded_starts
             <<" failed_spectral_solves="<<d.failed_spectral_solves
             <<" total_spectral_iterations="<<d.total_spectral_iterations
             <<" average_spectral_iterations="<<double(d.total_spectral_iterations)/std::max<size_t>(1,d.successful_cold_starts+d.successful_seeded_starts)
             <<" max_spectral_iterations="<<d.max_spectral_iterations
             <<" average_cold_iterations="<<double(d.cold_spectral_iterations)/std::max<size_t>(1,d.successful_cold_starts)
             <<" average_seeded_iterations="<<double(d.seeded_spectral_iterations)/std::max<size_t>(1,d.successful_seeded_starts)<<'\n';
    std::cerr<<"Energy timing: total_energy_time="<<d.total_energy_time
             <<" spectral_time="<<d.spectral_time<<" kinetic_time="<<d.kinetic_time
             <<" recovery_time="<<d.recovery_time
             <<" average_energy_point_time="<<(d.spectral_time+d.kinetic_time)/std::max<size_t>(1,d.cold_starts+d.seeded_starts)<<'\n';
    return result;
}
} // namespace sns
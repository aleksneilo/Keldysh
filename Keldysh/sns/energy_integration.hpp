#pragma once
#include "sns.hpp"
#include <array>
#include <memory>


namespace sns {
// Only this error is eligible for spectral recovery. Kinetic/I/O/memory errors propagate.
struct SpectralEnergyFailure : std::runtime_error {
    double spectral_seconds;
    explicit SpectralEnergyFailure(const std::string& message,double seconds=0):
        std::runtime_error(message),spectral_seconds(seconds) {}
};
struct EnergySample {
    std::array<double,3> integrand{};
    std::shared_ptr<const PairField> amplitudes;
    double spectral_residual=0, kinetic_residual=0;
    size_t spectral_iterations=0;
    double spectral_seconds=0,kinetic_seconds=0;
};
using EnergyEvaluator=std::function<EnergySample(double,const PairField*,const NumericalParams&)>;
struct EnergyNode {
    double epsilon=0, weight=0, gap_distance=0, spectral_residual=0;
    int level=0;
    std::string status;
};
struct EnergyIntegral {
    std::array<double,3> values{};
    double spectral_residual=0, kinetic_residual=0;
    EnergyDiagnostics diagnostics;
    std::vector<EnergyNode> nodes;
    EnergyCache cache;
};
void validate_energy_parameters(const NumericalParams&);
std::vector<double> energy_gap_edges(double v,const PhysicalParams&,int NF);
EnergyIntegral integrate_energy_adaptive(double v,const PhysicalParams&,const NumericalParams&,
    const EnergyEvaluator&,const EnergyCache* initial=nullptr,bool keep_solutions=false, const EnergyEvaluator& anchor_evaluate=EnergyEvaluator{});
CurrentResult solve_current_adaptive(double voltage,const PhysicalParams&,const NumericalParams&,
    const EnergyCache* initial=nullptr,EnergyCache* solutions=nullptr);


}
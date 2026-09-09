#pragma once
#include "sns.hpp"

namespace sns {
// Settings for both Keldysh/main.cpp and the optional standalone executable.
// Command-line arguments override the corresponding settings.
struct RunSettings {
    PhysicalParams physical;
    NumericalParams numerical;
    std::string mode = "spectral";
    double epsilon = 5.28; // E/(k_B Tc)
    double voltage = 0;   // eV/(k_B Tc)
    double sweep_start = 3.0, sweep_end = 0.5, sweep_step = 0.1; // eV/Delta
    std::string conductance_path = "conductance_L1.txt";
};
int run_cli(int argc, char** argv, const RunSettings& settings);
}

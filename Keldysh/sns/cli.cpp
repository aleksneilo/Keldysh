#include "runner.hpp"

// Optional standalone entry point. Keldysh.sln uses Keldysh/main.cpp instead.
int main(int argc, char** argv) {
    sns::RunSettings settings;
    settings.physical.L_N = std::sqrt(
        settings.physical.diffusion() / settings.physical.Delta);
    settings.mode = "help";
    settings.epsilon = 3 * settings.physical.Delta;
    settings.voltage = 2 * settings.physical.Delta;
    return sns::run_cli(argc, argv, settings);
}
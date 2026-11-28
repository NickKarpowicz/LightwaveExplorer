#include "Fourier.hpp"
#include "LightwaveExplorerInterfaceClasses.hpp"
#include "vpp/vpp.hpp"
#include <concepts>

template <std::floating_point F>
uint32_t run_fitting(const SimulationParameterSet& launch_parameters){
    return 0u;
}

template <std::floating_point F>
uint32_t solve(const SimulationParameterSet& launch_parameters);

template <std::floating_point F>
uint32_t solve_sequence(const SimulationParameterSet& launch_parameters);

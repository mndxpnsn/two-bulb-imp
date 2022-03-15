//
//  main.cpp
//  two-bulb-imp
//
//  Created by mndx on 12/03/2022.
//  Three-component two-bulb diffusion
//  using implicit time discretization
//

#include <iostream>

#include "lib.hpp"
#include "user_types.h"

int main(int argc, const char * argv[]) {
    
    // Number of grid nodes
    int ng = 10;
    
    // Experimental setup parameters
    e_params_t e_params;
    e_params.V = 5e-4;
    e_params.d = 2e-3;
    e_params.len = 1e-2;
    e_params.A = 0.25 * 3.14 * e_params.d * e_params.d;
    e_params.dz = e_params.len / ng;
    
    // Initial composition bulb 1
    b_data_t bulb_data;
    bulb_data.mol_fracs_bulb1.x1 = 0.501;
    bulb_data.mol_fracs_bulb1.x2 = 0.499;
    bulb_data.mol_fracs_bulb1.x3 = 1 - 0.501 - 0.499;
    
    // Initial composition bulb 2
    bulb_data.mol_fracs_bulb2.x1 = 0.0;
    bulb_data.mol_fracs_bulb2.x2 = 0.501;
    bulb_data.mol_fracs_bulb2.x3 = 1 - 0.0 - 0.501;
    
    // Total concentration
    p_params_t p_params;
    p_params.ct = 1.0;
    
    // Time parameters
    t_params_t t_params;
    t_params.to = 0.0;
    t_params.tf = 5.0;
    t_params.nt = 40;
    t_params.dt = (double) (t_params.tf - t_params.to) / t_params.nt;
    
    // Diffusivities
    p_params.D12 = 8.33e-5 * 3600;
    p_params.D13 = 6.8e-5 * 3600;
    p_params.D23 = 1.68e-5 * 3600;

    // Perform two-bulb diffusion experiment
    compute_bulb_compositions(e_params, p_params, t_params, ng, bulb_data);

    // Print results
    std::cout << "bulb 1 frac 1: " << bulb_data.mol_fracs_bulb1.x1 << std::endl;
    std::cout << "bulb 1 frac 2: " << bulb_data.mol_fracs_bulb1.x2 << std::endl;
    std::cout << "bulb 1 frac 3: " << bulb_data.mol_fracs_bulb1.x3 << std::endl;
    
    std::cout << "bulb 2 frac 1: " << bulb_data.mol_fracs_bulb2.x1 << std::endl;
    std::cout << "bulb 2 frac 2: " << bulb_data.mol_fracs_bulb2.x2 << std::endl;
    std::cout << "bulb 2 frac 3: " << bulb_data.mol_fracs_bulb2.x3 << std::endl;
    
    return 0;
}

//
//  lib.cpp
//  two-bulb-imp
//
//  Created by Derek Harrison on 14/03/2022.
//

#include <stdio.h>
#include <iostream>

#include "lib.hpp"
#include "user_types.h"

double a1(p_params_t & p_params, double x2) {
    double v1 = 1.0 / p_params.D13;
    double v2 = 1.0 / p_params.D12;
    
    return (v2 - v1) * x2 + v1;
}

double a2(p_params_t & p_params, double x1) {
    double v1 = 1.0 / p_params.D13;
    double v2 = 1.0 / p_params.D12;
    
    return x1 * (v1 - v2);
}

double b1(p_params_t & p_params, double x2) {
    double v2 = 1.0 / p_params.D12;
    double v3 = 1.0 / p_params.D23;
    
    return x2 * (v3 - v2);
}

double b2(p_params_t & p_params, double x1) {
    double v2 = 1.0 / p_params.D12;
    double v3 = 1.0 / p_params.D23;
    
    return (v2 - v3) * x1 + v3;
}

double alpha1(p_params_t & p_params, double x1, double x2) {
    double ct = p_params.ct;
    
    double a1_loc = a1(p_params, x2);
    double a2_loc = a2(p_params, x1);
    double b1_loc = b1(p_params, x2);
    double b2_loc = b2(p_params, x1);
    
    return -ct / (a2_loc - a1_loc * b2_loc / b1_loc);
}

double alpha2(p_params_t & p_params, double x1, double x2) {
    double ct = p_params.ct;
    
    double a1_loc = a1(p_params, x2);
    double a2_loc = a2(p_params, x1);
    double b1_loc = b1(p_params, x2);
    double b2_loc = b2(p_params, x1);
    
    return (a1_loc / b1_loc) * ct / (a2_loc - a1_loc * b2_loc / b1_loc);
}

double beta1(p_params_t & p_params, double x1, double x2) {
    double ct = p_params.ct;
    double a1_loc = a1(p_params, x2);
    double a2_loc = a2(p_params, x1);
    double alpha1_loc = alpha1(p_params, x1, x2);
    
    return -ct / a1_loc - a2_loc * alpha1_loc / a1_loc;
}

double beta2(p_params_t & p_params, double x1, double x2) {
    double a1_loc = a1(p_params, x2);
    double a2_loc = a2(p_params, x1);
    double alpha2_loc = alpha2(p_params, x1, x2);
    
    return - a2_loc * alpha2_loc / a1_loc;
}

void compute_bulb_compositions(e_params_t e_params,
                               p_params_t p_params,
                               t_params_t t_params,
                               int ng,
                               b_data_t & bulb_data) {
    
    double V = e_params.V;
    double A = e_params.A;
    double dz = e_params.dz;
    
    // Mole fractions in tube
    node_t * tube_fracs = new node_t[ng];
    node_t * tube_fracs_inter = new node_t[ng];
    node_t * tube_fracs_old = new node_t[ng];
    
    // Bulb data
    b_data_t bulb_data_inter = bulb_data;
    b_data_t bulb_data_old = bulb_data;
    
    // Initialize tube composition
    for(int node = 0; node < ng; ++node) {
        tube_fracs[node].x1 = 1.0 / 3;
        tube_fracs[node].x2 = 1.0 / 3;
        tube_fracs[node].x3 = 1.0 / 3;
        tube_fracs_inter[node].x1 = 1.0 / 3;
        tube_fracs_inter[node].x2 = 1.0 / 3;
        tube_fracs_inter[node].x3 = 1.0 / 3;
        tube_fracs_old[node].x1 = 1.0 / 3;
        tube_fracs_old[node].x2 = 1.0 / 3;
        tube_fracs_old[node].x3 = 1.0 / 3;
    }
    
    // Perform iterations
    int max_out_it = MAX_OUT;
    int max_in_it = MAX_IN;
    
    double t = t_params.to;
    double dt = (t_params.tf - t_params.to) / t_params.nt;
    
    while(t < t_params.tf) {
        
        bulb_data_old = bulb_data;
        
        for(int node = 0; node < ng; ++node) {
            tube_fracs_old[node] = tube_fracs[node];
        }
        
        int out_it = 0;
        while(out_it < max_out_it) {
            
            bulb_data_inter = bulb_data;
            
            for(int node = 0; node < ng; ++node) {
                tube_fracs_inter[node] = tube_fracs[node];
            }
            
            int in_it = 0;
            while(in_it < max_in_it) {
                
                // Bulb 1, component 1
                double x1 = bulb_data_inter.mol_fracs_bulb1.x1;
                double x2 = bulb_data_inter.mol_fracs_bulb1.x2;
                double beta1_loc = beta1(p_params, x1, x2);
                double beta2_loc = beta2(p_params, x1, x2);
                double ap1 = V * p_params.ct / dt - 2 * beta1_loc * A / dz;
                double x11 = tube_fracs_inter[0].x1;
                double x21 = tube_fracs_inter[0].x2;
                double old_term1 = V * p_params.ct * bulb_data_old.mol_fracs_bulb1.x1 / dt;
                
                double x1_b1 = -2 * beta1_loc * x11 * A / dz - 2 * beta2_loc * (x21 - x2) * A / dz + old_term1;
                x1_b1 = x1_b1 / ap1;
                
                bulb_data.mol_fracs_bulb1.x1 = x1_b1;
                
                // Bulb 1, component 2
                double alpha1_loc = alpha1(p_params, x1, x2);
                double alpha2_loc = alpha2(p_params, x1, x2);
                double ap2 = V * p_params.ct / dt - 2 * alpha2_loc * A / dz;
                double old_term2 = V * p_params.ct * bulb_data_old.mol_fracs_bulb1.x2 / dt;

                double x2_b1 = -2 * alpha1_loc / dz * (x11 - x1) * A  - 2 * alpha2_loc * x21 * A / dz + old_term2;
                x2_b1 = x2_b1 / ap2;
                
                bulb_data.mol_fracs_bulb1.x2 = x2_b1;
                bulb_data.mol_fracs_bulb1.x3 = 1.0 - x1_b1 - x2_b1;
                
                // Tube node 0, component 1
                double x1w = bulb_data_inter.mol_fracs_bulb1.x1;
                double x2w = bulb_data_inter.mol_fracs_bulb1.x2;
                double x1e = 0.5 * (tube_fracs_inter[0].x1 + tube_fracs_inter[1].x1);
                double x2e = 0.5 * (tube_fracs_inter[0].x2 + tube_fracs_inter[1].x2);
                
                double x2P = tube_fracs[0].x2;
                double x1E = tube_fracs[1].x1;
                double x2E = tube_fracs[1].x2;
                
                double beta1w = beta1(p_params, x1w, x2w);
                double beta2w = beta2(p_params, x1w, x2w);
                double beta1e = beta1(p_params, x1e, x2e);
                double beta2e = beta2(p_params, x1e, x2e);
                
                double ap1_tube0 = p_params.ct * dz / dt - beta1w / (0.5 * dz) - beta1e / dz;
                
                double old_term_tube1 = p_params.ct * dz * tube_fracs_old[0].x1 / dt;
                
                double x1_tube0 = - beta1w / (0.5 * dz) * x1w + beta2w / (0.5 * dz) * (x2P - x2w) - beta1e / dz * x1E - beta2e / dz * (x2E - x2P) + old_term_tube1;
                
                x1_tube0 = x1_tube0 / ap1_tube0;
                
                tube_fracs[0].x1 = x1_tube0;
                
                // Tube node 0, component 2
                x1E = tube_fracs[1].x1;
                x2E = tube_fracs[1].x2;
                double x1P = tube_fracs_inter[0].x1;
                
                double alpha1w = alpha1(p_params, x1w, x2w);
                double alpha2w = alpha2(p_params, x1w, x2w);
                double alpha1e = alpha1(p_params, x1e, x2e);
                double alpha2e = alpha2(p_params, x1e, x2e);
                
                double ap2_tube0 = p_params.ct * dz / dt - alpha2w / (0.5 * dz) - alpha2e / dz;
                
                double old_term_tube2 = p_params.ct * dz * tube_fracs_old[0].x2 / dt;
                
                double x2_tube0 = alpha1w / (0.5 * dz) * (x1P - x1w) - alpha2w / (0.5 * dz) * x2w - alpha1e / dz * (x1E - x1P) - alpha2e / dz * x2E + old_term_tube2;
                
                x2_tube0 = x2_tube0 / ap2_tube0;
                
                tube_fracs[0].x2 = x2_tube0;
                tube_fracs[0].x3 = 1.0 - x1_tube0 - x2_tube0;
                
                // Tube mid nodes, component 1
                for(int node = 1; node < ng - 1; ++node) {
                    double x1w = 0.5 * (tube_fracs_inter[node - 1].x1 + tube_fracs_inter[node].x1);
                    double x2w = 0.5 * (tube_fracs_inter[node - 1].x2 + tube_fracs_inter[node].x2);
                    double x1e = 0.5 * (tube_fracs_inter[node].x1 + tube_fracs_inter[node + 1].x1);
                    double x2e = 0.5 * (tube_fracs_inter[node].x2 + tube_fracs_inter[node + 1].x2);

                    double x1W = tube_fracs[node - 1].x1;
                    double x2W = tube_fracs[node - 1].x2;
                    double x1E = tube_fracs[node + 1].x1;
                    double x2E = tube_fracs[node + 1].x2;
                    double x2P = tube_fracs[node].x2;

                    double beta1w = beta1(p_params, x1w, x2w);
                    double beta2w = beta2(p_params, x1w, x2w);
                    double beta1e = beta1(p_params, x1e, x2e);
                    double beta2e = beta2(p_params, x1e, x2e);

                    double ap1_tube_node = p_params.ct * dz / dt - beta1w / dz - beta1e / dz;
                    
                    double old_term_tube1 = p_params.ct * dz / dt * tube_fracs_old[node].x1;

                    double x1_tube_node = - beta1w / dz * x1W + beta2w / dz * (x2P - x2W) - beta1e / dz * x1E - beta2e / dz * (x2E - x2P) + old_term_tube1;

                    x1_tube_node = x1_tube_node / ap1_tube_node;

                    tube_fracs[node].x1 = x1_tube_node;
                }

                // Tube mid nodes, component 2
                for(int node = 1; node < ng - 1; ++node) {
                    double x1w = 0.5 * (tube_fracs_inter[node - 1].x1 + tube_fracs_inter[node].x1);
                    double x2w = 0.5 * (tube_fracs_inter[node - 1].x2 + tube_fracs_inter[node].x2);
                    double x1e = 0.5 * (tube_fracs_inter[node].x1 + tube_fracs_inter[node + 1].x1);
                    double x2e = 0.5 * (tube_fracs_inter[node].x2 + tube_fracs_inter[node + 1].x2);

                    double x1W = tube_fracs[node - 1].x1;
                    double x2W = tube_fracs[node - 1].x2;
                    double x1E = tube_fracs[node + 1].x1;
                    double x2E = tube_fracs[node + 1].x2;
                    double x1P = tube_fracs[node].x1;

                    double alpha1w = alpha1(p_params, x1w, x2w);
                    double alpha2w = alpha2(p_params, x1w, x2w);
                    double alpha1e = alpha1(p_params, x1e, x2e);
                    double alpha2e = alpha2(p_params, x1e, x2e);

                    double ap2_tube_node = p_params.ct * dz / dt - alpha2w / dz - alpha2e / dz;
                    
                    double old_term_tube2 = p_params.ct * dz / dt * tube_fracs_old[node].x2;

                    double x2_tube_node = alpha1w / dz * (x1P - x1W) - alpha2w / dz * x2W - alpha1e * (x1E - x1P) / dz - alpha2e / dz * x2E + old_term_tube2;

                    x2_tube_node = x2_tube_node / ap2_tube_node;

                    tube_fracs[node].x2 = x2_tube_node;
                }
                
                // Tube node n, component 1
                x1w = 0.5 * (tube_fracs_inter[ng - 2].x1 + tube_fracs_inter[ng - 1].x1);
                x2w = 0.5 * (tube_fracs_inter[ng - 2].x2 + tube_fracs_inter[ng - 1].x2);
                x1e = bulb_data_inter.mol_fracs_bulb2.x1;
                x2e = bulb_data_inter.mol_fracs_bulb2.x2;
                
                double x1W = tube_fracs[ng - 2].x1;
                double x2W = tube_fracs[ng - 2].x2;
                x1E = x1e;
                x2E = x2e;
                x2P = tube_fracs[ng - 1].x2;
                
                beta1w = beta1(p_params, x1w, x2w);
                beta2w = beta2(p_params, x1w, x2w);
                beta1e = beta1(p_params, x1e, x2e);
                beta2e = beta2(p_params, x1e, x2e);
                
                double ap1_tube_n = p_params.ct * dz / dt - beta1w / dz - beta1e / (0.5 * dz);
                
                double old_term_tube_n_1 = p_params.ct * dz / dt * tube_fracs_old[ng - 1].x1;
                
                double x1_tube_n = - beta1w / dz * x1W + beta2w / dz * (x2P - x2W) - beta1e / (0.5 * dz) * x1E - beta2e / dz * (x2E - x2P) + old_term_tube_n_1;
                
                x1_tube_n = x1_tube_n / ap1_tube_n;
                
                tube_fracs[ng - 1].x1 = x1_tube_n;
                
                // Tube node n, component 2
                x1w = 0.5 * (tube_fracs_inter[ng - 2].x1 + tube_fracs_inter[ng - 1].x1);
                x2w = 0.5 * (tube_fracs_inter[ng - 2].x2 + tube_fracs_inter[ng - 1].x2);
                x1e = bulb_data_inter.mol_fracs_bulb2.x1;
                x2e = bulb_data_inter.mol_fracs_bulb2.x2;
                
                x1W = tube_fracs[ng - 2].x1;
                x2W = tube_fracs[ng - 2].x2;
                x1E = x1e;
                x2E = x2e;
                x1P = tube_fracs[ng - 1].x1;
                x2P = tube_fracs[ng - 1].x2;
                
                alpha1w = alpha1(p_params, x1w, x2w);
                alpha2w = alpha2(p_params, x1w, x2w);
                alpha1e = alpha1(p_params, x1e, x2e);
                alpha2e = alpha2(p_params, x1e, x2e);
                
                double ap2_tube_n = p_params.ct * dz / dt - alpha2w / dz - alpha2e / (0.5 * dz);
                
                double old_term_tube_n_2 = p_params.ct * dz / dt * tube_fracs_old[ng - 1].x2;
                
                double x2_tube_n = alpha1w / dz * (x1P - x1W) - alpha2w / dz * x2W - alpha1e / (0.5 * dz) * (x1E - x1P) - alpha2e / (0.5 * dz) * x2E + old_term_tube_n_2;
                
                x2_tube_n = x2_tube_n / ap2_tube_n;
                
                tube_fracs[ng - 1].x2 = x2_tube_n;
                
                // Set mole fraction component 3
                for(int node = 0; node < ng; ++node) {
                    tube_fracs[node].x3 = 1.0 - tube_fracs[node].x1 - tube_fracs[node].x2;
                }
                
                // Bulb 2, component 1
                x1 = bulb_data_inter.mol_fracs_bulb2.x1;
                x2 = bulb_data_inter.mol_fracs_bulb2.x2;
                beta1_loc = beta1(p_params, x1, x2);
                beta2_loc = beta2(p_params, x1, x2);
                ap1 = V * p_params.ct / dt - 2 * beta1_loc * A / dz;
                double x1n = tube_fracs_inter[ng - 1].x1;
                double x2n = tube_fracs_inter[ng - 1].x2;
                old_term1 = V * p_params.ct * bulb_data_old.mol_fracs_bulb2.x1 / dt;
                
                double x1_b2 = -2 * beta1_loc * x1n * A / dz + 2 * beta2_loc * (x2 - x2n) * A / dz + old_term1;
                x1_b2 = x1_b2 / ap1;
                
                bulb_data.mol_fracs_bulb2.x1 = x1_b2;
                
                // Bulb 2, component 2
                alpha1_loc = alpha1(p_params, x1, x2);
                alpha2_loc = alpha2(p_params, x1, x2);
                ap2 = V * p_params.ct / dt - 2 * alpha2_loc * A / dz;
                old_term2 = V * p_params.ct * bulb_data_old.mol_fracs_bulb2.x2 / dt;

                double x2_b2 = 2 * alpha1_loc / dz * (x1 - x1n) * A  - 2 * alpha2_loc * x2n * A / dz + old_term2;
                x2_b2 = x2_b2 / ap2;
                
                bulb_data.mol_fracs_bulb2.x2 = x2_b2;
                bulb_data.mol_fracs_bulb2.x3 = 1.0 - x1_b2 - x2_b2;
                
                in_it++;
            }
            
            out_it++;
        }
        
        t = t + dt;
    }
    
    delete [] tube_fracs;
    delete [] tube_fracs_inter;
    delete [] tube_fracs_old;
}

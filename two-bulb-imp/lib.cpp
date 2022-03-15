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

void bulb1_c1(c_data_t & comp_data) {
    
    p_params_t p_params = comp_data.p_params;
    e_params_t e_params = comp_data.e_params;
    t_params_t t_params = comp_data.t_params;
    b_data_t bulb_data_old = comp_data.bulb_data_old;
    b_data_t bulb_data_inter = comp_data.bulb_data_inter;
    node_t * tube_fracs_inter = comp_data.tube_fracs_inter;
    
    double V = e_params.V;
    double A = e_params.A;
    double dz = e_params.dz;
    double dt = t_params.dt;
    
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
    
    comp_data.bulb_data.mol_fracs_bulb1.x1 = x1_b1;
}

void bulb1_c2(c_data_t & comp_data) {
    
    p_params_t p_params = comp_data.p_params;
    e_params_t e_params = comp_data.e_params;
    t_params_t t_params = comp_data.t_params;
    b_data_t bulb_data_old = comp_data.bulb_data_old;
    b_data_t bulb_data = comp_data.bulb_data;
    b_data_t bulb_data_inter = comp_data.bulb_data_inter;
    node_t * tube_fracs_inter = comp_data.tube_fracs_inter;
    
    double V = e_params.V;
    double A = e_params.A;
    double dz = e_params.dz;
    double dt = t_params.dt;
    
    // Bulb 1, component 2
    double x1 = bulb_data_inter.mol_fracs_bulb1.x1;
    double x2 = bulb_data_inter.mol_fracs_bulb1.x2;
    double alpha1_loc = alpha1(p_params, x1, x2);
    double alpha2_loc = alpha2(p_params, x1, x2);
    double ap2 = V * p_params.ct / dt - 2 * alpha2_loc * A / dz;
    double x11 = tube_fracs_inter[0].x1;
    double x21 = tube_fracs_inter[0].x2;
    double old_term2 = V * p_params.ct * bulb_data_old.mol_fracs_bulb1.x2 / dt;

    double x2_b1 = -2 * alpha1_loc / dz * (x11 - x1) * A  - 2 * alpha2_loc * x21 * A / dz + old_term2;
    x2_b1 = x2_b1 / ap2;
    
    comp_data.bulb_data.mol_fracs_bulb1.x2 = x2_b1;
    comp_data.bulb_data.mol_fracs_bulb1.x3 = 1.0 - bulb_data.mol_fracs_bulb1.x1 - x2_b1;
}

void bulb2_c1(c_data_t & comp_data) {
    
    int ng = comp_data.ng;
    p_params_t p_params = comp_data.p_params;
    e_params_t e_params = comp_data.e_params;
    t_params_t t_params = comp_data.t_params;
    b_data_t bulb_data_old = comp_data.bulb_data_old;
    b_data_t bulb_data_inter = comp_data.bulb_data_inter;
    node_t * tube_fracs_inter = comp_data.tube_fracs_inter;
    
    double V = e_params.V;
    double A = e_params.A;
    double dz = e_params.dz;
    double dt = t_params.dt;
    
    double x1 = bulb_data_inter.mol_fracs_bulb2.x1;
    double x2 = bulb_data_inter.mol_fracs_bulb2.x2;
    double beta1_loc = beta1(p_params, x1, x2);
    double beta2_loc = beta2(p_params, x1, x2);
    double ap1 = V * p_params.ct / dt - 2 * beta1_loc * A / dz;
    double x1n = tube_fracs_inter[ng - 1].x1;
    double x2n = tube_fracs_inter[ng - 1].x2;
    double old_term1 = V * p_params.ct * bulb_data_old.mol_fracs_bulb2.x1 / dt;
    
    double x1_b2 = -2 * beta1_loc * x1n * A / dz + 2 * beta2_loc * (x2 - x2n) * A / dz + old_term1;
    x1_b2 = x1_b2 / ap1;
    
    comp_data.bulb_data.mol_fracs_bulb2.x1 = x1_b2;
}

void bulb2_c2(c_data_t & comp_data) {
    
    int ng = comp_data.ng;
    p_params_t p_params = comp_data.p_params;
    e_params_t e_params = comp_data.e_params;
    t_params_t t_params = comp_data.t_params;
    b_data_t bulb_data_old = comp_data.bulb_data_old;
    b_data_t bulb_data = comp_data.bulb_data;
    b_data_t bulb_data_inter = comp_data.bulb_data_inter;
    node_t * tube_fracs_inter = comp_data.tube_fracs_inter;
    
    double V = e_params.V;
    double A = e_params.A;
    double dz = e_params.dz;
    double dt = t_params.dt;
    
    double x1 = bulb_data_inter.mol_fracs_bulb2.x1;
    double x2 = bulb_data_inter.mol_fracs_bulb2.x2;
    double alpha1_loc = alpha1(p_params, x1, x2);
    double alpha2_loc = alpha2(p_params, x1, x2);
    double ap2 = V * p_params.ct / dt - 2 * alpha2_loc * A / dz;
    double x1n = tube_fracs_inter[ng - 1].x1;
    double x2n = tube_fracs_inter[ng - 1].x2;
    double old_term2 = V * p_params.ct * bulb_data_old.mol_fracs_bulb2.x2 / dt;

    double x2_b2 = 2 * alpha1_loc / dz * (x1 - x1n) * A  - 2 * alpha2_loc * x2n * A / dz + old_term2;
    x2_b2 = x2_b2 / ap2;
    
    comp_data.bulb_data.mol_fracs_bulb2.x2 = x2_b2;
    comp_data.bulb_data.mol_fracs_bulb2.x3 = 1.0 - bulb_data.mol_fracs_bulb2.x1 - x2_b2;
}

void tube0_c1(c_data_t & comp_data) {
    
    p_params_t p_params = comp_data.p_params;
    e_params_t e_params = comp_data.e_params;
    t_params_t t_params = comp_data.t_params;
    b_data_t bulb_data_inter = comp_data.bulb_data_inter;
    node_t * tube_fracs_inter = comp_data.tube_fracs_inter;
    node_t * tube_fracs = comp_data.tube_fracs;
    node_t * tube_fracs_old = comp_data.tube_fracs_old;
    
    double dz = e_params.dz;
    double dt = t_params.dt;
    
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
    
    comp_data.tube_fracs[0].x1 = x1_tube0;
}

void tube0_c2(c_data_t & comp_data) {
    
    p_params_t p_params = comp_data.p_params;
    e_params_t e_params = comp_data.e_params;
    t_params_t t_params = comp_data.t_params;
    b_data_t bulb_data_inter = comp_data.bulb_data_inter;
    node_t * tube_fracs_inter = comp_data.tube_fracs_inter;
    node_t * tube_fracs = comp_data.tube_fracs;
    node_t * tube_fracs_old = comp_data.tube_fracs_old;
    
    double dz = e_params.dz;
    double dt = t_params.dt;
    
    double x1w = bulb_data_inter.mol_fracs_bulb1.x1;
    double x2w = bulb_data_inter.mol_fracs_bulb1.x2;
    double x1e = 0.5 * (tube_fracs_inter[0].x1 + tube_fracs_inter[1].x1);
    double x2e = 0.5 * (tube_fracs_inter[0].x2 + tube_fracs_inter[1].x2);
    
    double x1E = tube_fracs[1].x1;
    double x2E = tube_fracs[1].x2;
    double x1P = tube_fracs_inter[0].x1;
    
    double alpha1w = alpha1(p_params, x1w, x2w);
    double alpha2w = alpha2(p_params, x1w, x2w);
    double alpha1e = alpha1(p_params, x1e, x2e);
    double alpha2e = alpha2(p_params, x1e, x2e);
    
    double ap2_tube0 = p_params.ct * dz / dt - alpha2w / (0.5 * dz) - alpha2e / dz;
    
    double old_term_tube2 = p_params.ct * dz * tube_fracs_old[0].x2 / dt;
    
    double x2_tube0 = alpha1w / (0.5 * dz) * (x1P - x1w) - alpha2w / (0.5 * dz) * x2w - alpha1e / dz * (x1E - x1P) - alpha2e / dz * x2E + old_term_tube2;
    
    x2_tube0 = x2_tube0 / ap2_tube0;
    
    comp_data.tube_fracs[0].x2 = x2_tube0;
    comp_data.tube_fracs[0].x3 = 1.0 - tube_fracs[0].x1 - x2_tube0;
}

void mid_nodes1(c_data_t & comp_data) {
    
    int ng = comp_data.ng;
    p_params_t p_params = comp_data.p_params;
    e_params_t e_params = comp_data.e_params;
    t_params_t t_params = comp_data.t_params;
    node_t * tube_fracs_inter = comp_data.tube_fracs_inter;
    node_t * tube_fracs = comp_data.tube_fracs;
    node_t * tube_fracs_old = comp_data.tube_fracs_old;
    
    double dz = e_params.dz;
    double dt = t_params.dt;
    
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

        comp_data.tube_fracs[node].x1 = x1_tube_node;
    }
}

void mid_nodes2(c_data_t & comp_data)  {
    
    int ng = comp_data.ng;
    p_params_t p_params = comp_data.p_params;
    e_params_t e_params = comp_data.e_params;
    t_params_t t_params = comp_data.t_params;
    node_t * tube_fracs_inter = comp_data.tube_fracs_inter;
    node_t * tube_fracs = comp_data.tube_fracs;
    node_t * tube_fracs_old = comp_data.tube_fracs_old;
    
    double dz = e_params.dz;
    double dt = t_params.dt;
    
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

        comp_data.tube_fracs[node].x2 = x2_tube_node;
    }
}

void tube_n_c1(c_data_t & comp_data) {
    
    int ng = comp_data.ng;
    p_params_t p_params = comp_data.p_params;
    e_params_t e_params = comp_data.e_params;
    t_params_t t_params = comp_data.t_params;
    node_t * tube_fracs_inter = comp_data.tube_fracs_inter;
    node_t * tube_fracs = comp_data.tube_fracs;
    node_t * tube_fracs_old = comp_data.tube_fracs_old;
    b_data_t bulb_data_inter = comp_data.bulb_data_inter;
    
    double dz = e_params.dz;
    double dt = t_params.dt;
    
    double x1w = 0.5 * (tube_fracs_inter[ng - 2].x1 + tube_fracs_inter[ng - 1].x1);
    double x2w = 0.5 * (tube_fracs_inter[ng - 2].x2 + tube_fracs_inter[ng - 1].x2);
    double x1e = bulb_data_inter.mol_fracs_bulb2.x1;
    double x2e = bulb_data_inter.mol_fracs_bulb2.x2;
    
    double x1W = tube_fracs[ng - 2].x1;
    double x2W = tube_fracs[ng - 2].x2;
    double x1E = x1e;
    double x2E = x2e;
    double x2P = tube_fracs[ng - 1].x2;
    
    double beta1w = beta1(p_params, x1w, x2w);
    double beta2w = beta2(p_params, x1w, x2w);
    double beta1e = beta1(p_params, x1e, x2e);
    double beta2e = beta2(p_params, x1e, x2e);
    
    double ap1_tube_n = p_params.ct * dz / dt - beta1w / dz - beta1e / (0.5 * dz);
    
    double old_term_tube_n_1 = p_params.ct * dz / dt * tube_fracs_old[ng - 1].x1;
    
    double x1_tube_n = - beta1w / dz * x1W + beta2w / dz * (x2P - x2W) - beta1e / (0.5 * dz) * x1E - beta2e / dz * (x2E - x2P) + old_term_tube_n_1;
    
    x1_tube_n = x1_tube_n / ap1_tube_n;
    
    comp_data.tube_fracs[ng - 1].x1 = x1_tube_n;
}

void tube_n_c2(c_data_t & comp_data)  {
    
    int ng = comp_data.ng;
    p_params_t p_params = comp_data.p_params;
    e_params_t e_params = comp_data.e_params;
    t_params_t t_params = comp_data.t_params;
    node_t * tube_fracs_inter = comp_data.tube_fracs_inter;
    node_t * tube_fracs = comp_data.tube_fracs;
    node_t * tube_fracs_old = comp_data.tube_fracs_old;
    b_data_t bulb_data_inter = comp_data.bulb_data_inter;
    
    double dz = e_params.dz;
    double dt = t_params.dt;
    
    double x1w = 0.5 * (tube_fracs_inter[ng - 2].x1 + tube_fracs_inter[ng - 1].x1);
    double x2w = 0.5 * (tube_fracs_inter[ng - 2].x2 + tube_fracs_inter[ng - 1].x2);
    double x1e = bulb_data_inter.mol_fracs_bulb2.x1;
    double x2e = bulb_data_inter.mol_fracs_bulb2.x2;
    
    double x1W = tube_fracs[ng - 2].x1;
    double x2W = tube_fracs[ng - 2].x2;
    double x1E = x1e;
    double x2E = x2e;
    double x1P = tube_fracs[ng - 1].x1;
    
    double alpha1w = alpha1(p_params, x1w, x2w);
    double alpha2w = alpha2(p_params, x1w, x2w);
    double alpha1e = alpha1(p_params, x1e, x2e);
    double alpha2e = alpha2(p_params, x1e, x2e);
    
    double ap2_tube_n = p_params.ct * dz / dt - alpha2w / dz - alpha2e / (0.5 * dz);
    
    double old_term_tube_n_2 = p_params.ct * dz / dt * tube_fracs_old[ng - 1].x2;
    
    double x2_tube_n = alpha1w / dz * (x1P - x1W) - alpha2w / dz * x2W - alpha1e / (0.5 * dz) * (x1E - x1P) - alpha2e / (0.5 * dz) * x2E + old_term_tube_n_2;
    
    x2_tube_n = x2_tube_n / ap2_tube_n;
    
    comp_data.tube_fracs[ng - 1].x2 = x2_tube_n;
}

void set_frac_comp3(c_data_t & comp_data) {
    int ng = comp_data.ng;
    
    for(int node = 0; node < ng; ++node) {
        comp_data.tube_fracs[node].x3 = 1.0 - comp_data.tube_fracs[node].x1 - comp_data.tube_fracs[node].x2;
    }
}

void update_tube_fracs(c_data_t & comp_data) {
    
    // Bulb 1, component 1
    bulb1_c1(comp_data);
    
    // Bulb 1, component 2
    bulb1_c2(comp_data);

    // Tube node 0, component 1
    tube0_c1(comp_data);
    
    // Tube node 0, component 2
    tube0_c2(comp_data);
    
    // Tube mid nodes, component 1
    mid_nodes1(comp_data);
    
    // Tube mid nodes, component 2
    mid_nodes2(comp_data);
    
    // Tube node n, component 1
    tube_n_c1(comp_data);
    
    // Tube node n, component 2
    tube_n_c2(comp_data);
    
    // Set mole fraction component 3
    set_frac_comp3(comp_data);
    
    // Bulb 2, component 1
    bulb2_c1(comp_data);
    
    // Bulb 2, component 2
    bulb2_c2(comp_data);
}

void compute_bulb_compositions(e_params_t e_params,
                               p_params_t p_params,
                               t_params_t t_params,
                               int ng,
                               b_data_t & bulb_data) {
    
    // Organize data
    c_data_t comp_data;
    
    // Allocate data for tube composition
    comp_data.tube_fracs = new node_t[ng];
    comp_data.tube_fracs_inter = new node_t[ng];
    comp_data.tube_fracs_old = new node_t[ng];
    
    // Initialize tube composition
    for(int node = 0; node < ng; ++node) {
        comp_data.tube_fracs[node].x1 = 1.0 / 3;
        comp_data.tube_fracs[node].x2 = 1.0 / 3;
        comp_data.tube_fracs[node].x3 = 1.0 / 3;
        comp_data.tube_fracs_inter[node].x1 = 1.0 / 3;
        comp_data.tube_fracs_inter[node].x2 = 1.0 / 3;
        comp_data.tube_fracs_inter[node].x3 = 1.0 / 3;
        comp_data.tube_fracs_old[node].x1 = 1.0 / 3;
        comp_data.tube_fracs_old[node].x2 = 1.0 / 3;
        comp_data.tube_fracs_old[node].x3 = 1.0 / 3;
    }
    
    // Perform iterations
    int max_out_it = MAX_OUT;
    int max_in_it = MAX_IN;
    
    double t = t_params.to;
    double dt = (t_params.tf - t_params.to) / t_params.nt;
    
    // Organize data continued
    comp_data.ng = ng;
    comp_data.p_params = p_params;
    comp_data.e_params = e_params;
    comp_data.t_params = t_params;
    comp_data.bulb_data = bulb_data;
    comp_data.bulb_data_old = bulb_data;
    comp_data.bulb_data_inter = bulb_data;
    
    // Compute composition
    while(t < t_params.tf) {
        
        comp_data.bulb_data_old = comp_data.bulb_data;
        
        for(int node = 0; node < ng; ++node)
            comp_data.tube_fracs_old[node] = comp_data.tube_fracs[node];
        
        // Outer Gauss-Seidel iterations
        int out_it = 0;
        while(out_it < max_out_it) {
            
            comp_data.bulb_data_inter = comp_data.bulb_data;
            
            for(int node = 0; node < ng; ++node)
                comp_data.tube_fracs_inter[node] = comp_data.tube_fracs[node];
            
            // Inner Gauss-Seidel iterations
            int in_it = 0;
            while(in_it < max_in_it) {

                update_tube_fracs(comp_data);
                
                in_it++;
            }
            
            out_it++;
        }
        
        t = t + dt;
    }
    
    // Set bulb data
    bulb_data = comp_data.bulb_data;
    
    delete [] comp_data.tube_fracs;
    delete [] comp_data.tube_fracs_inter;
    delete [] comp_data.tube_fracs_old;
}

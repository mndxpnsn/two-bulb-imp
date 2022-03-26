//
//  lib.cpp
//  two-bulb-imp
//
//  Created by dwh on 14/03/2022.
//

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>

#include "lib.hpp"
#include "user_types.h"

double ** anb_coeff1;
double * ap_coeff1;
double ** anb_coeff2;
double * ap_coeff2;

bool is_stable = true;
bool l_one1 = false;
bool l_one2 = false;

void check_stability(int ng) {

    int n = ng;

    for(int node = 0; node < n + 2; ++node) {
        int abs_sum1 = 0.0;
        int abs_sum2 = 0.0;

        for(int j = 0; j < 3; ++j) {
            abs_sum1 = abs_sum1 + fabs(anb_coeff1[node][j]);
            abs_sum2 = abs_sum2 + fabs(anb_coeff2[node][j]);
        }

        double ratio1 = abs_sum1 / fabs(ap_coeff1[node]);
        bool is_stable_loc1 = ratio1 <= 1;

        double ratio2 = abs_sum2 / fabs(ap_coeff2[node]);
        bool is_stable_loc2 = ratio2 <= 1;

        if(!is_stable_loc1) { is_stable = false; }
        if(!is_stable_loc2) { is_stable = false; }

        if(ratio1 < 1) { l_one1 = true; }
        if(ratio2 < 1) { l_one2 = true; }
    }

    is_stable = is_stable && l_one1 && l_one2;
}

void print_stability() {
    std::cout << "is_stable: " << is_stable << std::endl;
}

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

    double V = comp_data.e_params.V;
    double A = comp_data.e_params.A;
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    double x1 = comp_data.bulb_data.mol_fracs_bulb1.x1;
    double x2 = comp_data.bulb_data.mol_fracs_bulb1.x2;
    
    double beta1_loc = beta1(comp_data.p_params, x1, x2);
    double beta2_loc = beta2(comp_data.p_params, x1, x2);
    
    double ap1 = V * comp_data.p_params.ct / dt - 2 * beta1_loc * A / dz;
    
    double x11 = comp_data.tube_fracs[0].x1;
    double x21 = comp_data.tube_fracs[0].x2;
    
    double old_term1 = V * comp_data.p_params.ct * comp_data.bulb_data_old.mol_fracs_bulb1.x1 / dt;
    
    double x1_b1 = -2 * beta1_loc * x11 * A / dz - 2 * beta2_loc * (x21 - x2) * A / dz + old_term1;
    
    x1_b1 = x1_b1 / ap1;
    
    // Store coefficients for analysis
    ap_coeff1[0] = ap1;
    anb_coeff1[0][0] = -2 * beta1_loc * 1.0 * A / dz;
    anb_coeff1[0][1] = 2 * beta2_loc * 1.0 * A / dz;
    anb_coeff1[0][2] = 2 * beta2_loc * 1.0 * A / dz;

    comp_data.bulb_data.mol_fracs_bulb1.x1 = x1_b1;
}

void bulb1_c2(c_data_t & comp_data) {

    double V = comp_data.e_params.V;
    double A = comp_data.e_params.A;
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    double x1 = comp_data.bulb_data.mol_fracs_bulb1.x1;
    double x2 = comp_data.bulb_data.mol_fracs_bulb1.x2;
    
    double alpha1_loc = alpha1(comp_data.p_params, x1, x2);
    double alpha2_loc = alpha2(comp_data.p_params, x1, x2);
    
    double ap2 = V * comp_data.p_params.ct / dt - 2 * alpha2_loc * A / dz;
    
    double x11 = comp_data.tube_fracs[0].x1;
    double x21 = comp_data.tube_fracs[0].x2;
    
    double old_term2 = V * comp_data.p_params.ct * comp_data.bulb_data_old.mol_fracs_bulb1.x2 / dt;

    double x2_b1 = -2 * alpha1_loc / dz * (x11 - x1) * A  - 2 * alpha2_loc * x21 * A / dz + old_term2;
    
    x2_b1 = x2_b1 / ap2;

    // Store coefficients for analysis
    ap_coeff2[0] = ap2;
    anb_coeff2[0][0] = -2 * alpha1_loc / dz * 1.0 * A;
    anb_coeff2[0][1] = -2 * alpha1_loc / dz * 1.0 * A;
    anb_coeff2[0][2] = 2 * alpha2_loc * 1.0 * A / dz;

    comp_data.bulb_data.mol_fracs_bulb1.x2 = x2_b1;
    comp_data.bulb_data.mol_fracs_bulb1.x3 = 1.0 - comp_data.bulb_data.mol_fracs_bulb1.x1 - x2_b1;
}

void bulb2_c1(c_data_t & comp_data) {
    
    int ng = comp_data.ng;
    double V = comp_data.e_params.V;
    double A = comp_data.e_params.A;
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    double x1 = comp_data.bulb_data.mol_fracs_bulb2.x1;
    double x2 = comp_data.bulb_data.mol_fracs_bulb2.x2;
    
    double beta1_loc = beta1(comp_data.p_params, x1, x2);
    double beta2_loc = beta2(comp_data.p_params, x1, x2);
    
    double ap1 = V * comp_data.p_params.ct / dt - 2 * beta1_loc * A / dz;
    
    double x1n = comp_data.tube_fracs[ng - 1].x1;
    double x2n = comp_data.tube_fracs[ng - 1].x2;
    
    double old_term1 = V * comp_data.p_params.ct * comp_data.bulb_data_old.mol_fracs_bulb2.x1 / dt;
    
    double x1_b2 = -2 * beta1_loc * x1n * A / dz + 2 * beta2_loc * (x2 - x2n) * A / dz + old_term1;
    
    x1_b2 = x1_b2 / ap1;
    
    // Store coefficients for analysis
    ap_coeff1[ng + 1] = ap1;
    anb_coeff1[ng + 1][0] = -2 * beta1_loc * 1.0 * A / dz;
    anb_coeff1[ng + 1][1] = 2 * beta2_loc * 1.0 * A / dz;
    anb_coeff1[ng + 1][2] = 2 * beta2_loc * 1.0 * A / dz;

    comp_data.bulb_data.mol_fracs_bulb2.x1 = x1_b2;
}

void bulb2_c2(c_data_t & comp_data) {
    
    int ng = comp_data.ng;
    double V = comp_data.e_params.V;
    double A = comp_data.e_params.A;
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    double x1 = comp_data.bulb_data.mol_fracs_bulb2.x1;
    double x2 = comp_data.bulb_data.mol_fracs_bulb2.x2;
    
    double alpha1_loc = alpha1(comp_data.p_params, x1, x2);
    double alpha2_loc = alpha2(comp_data.p_params, x1, x2);
    
    double ap2 = V * comp_data.p_params.ct / dt - 2 * alpha2_loc * A / dz;
    
    double x1n = comp_data.tube_fracs[ng - 1].x1;
    double x2n = comp_data.tube_fracs[ng - 1].x2;
    
    double old_term2 = V * comp_data.p_params.ct * comp_data.bulb_data_old.mol_fracs_bulb2.x2 / dt;

    double x2_b2 = 2 * alpha1_loc / dz * (x1 - x1n) * A  - 2 * alpha2_loc * x2n * A / dz + old_term2;
    
    x2_b2 = x2_b2 / ap2;

    // Store coefficients for analysis
    ap_coeff2[ng + 1] = ap2;
    anb_coeff2[ng + 1][0] = -2 * alpha1_loc / dz * 1.0 * A;
    anb_coeff2[ng + 1][1] = -2 * alpha1_loc / dz * 1.0 * A;
    anb_coeff2[ng + 1][2] = 2 * alpha2_loc * 1.0 * A / dz;

    comp_data.bulb_data.mol_fracs_bulb2.x2 = x2_b2;
    comp_data.bulb_data.mol_fracs_bulb2.x3 = 1.0 - comp_data.bulb_data.mol_fracs_bulb2.x1 - x2_b2;
}

void tube0_c1(c_data_t & comp_data) {

    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    double x1w = comp_data.bulb_data_inter.mol_fracs_bulb1.x1;
    double x2w = comp_data.bulb_data_inter.mol_fracs_bulb1.x2;
    double x1e = 0.5 * (comp_data.tube_fracs_inter[0].x1 + comp_data.tube_fracs_inter[1].x1);
    double x2e = 0.5 * (comp_data.tube_fracs_inter[0].x2 + comp_data.tube_fracs_inter[1].x2);
    
    double x2P = comp_data.tube_fracs[0].x2;
    double x1E = comp_data.tube_fracs[1].x1;
    double x2E = comp_data.tube_fracs[1].x2;
    
    double beta1w = beta1(comp_data.p_params, x1w, x2w);
    double beta2w = beta2(comp_data.p_params, x1w, x2w);
    double beta1e = beta1(comp_data.p_params, x1e, x2e);
    double beta2e = beta2(comp_data.p_params, x1e, x2e);
    
    double ap1_tube0 = comp_data.p_params.ct * dz / dt - beta1w / (0.5 * dz) - beta1e / dz;
    
    double old_term_tube1 = comp_data.p_params.ct * dz * comp_data.tube_fracs_old[0].x1 / dt;
    
    double x1_tube0 = - beta1w / (0.5 * dz) * x1w + beta2w / (0.5 * dz) * (x2P - x2w) - beta1e / dz * x1E - beta2e / dz * (x2E - x2P) + old_term_tube1;
    
    x1_tube0 = x1_tube0 / ap1_tube0;

    // Store coefficients for analysis
    ap_coeff1[1] = ap1_tube0;
    anb_coeff1[1][0] = beta1e / dz;
    anb_coeff1[1][1] = beta1w / (0.5 * dz);
    anb_coeff1[1][2] = 0.0;

    comp_data.tube_fracs[0].x1 = x1_tube0;
}

void tube0_c2(c_data_t & comp_data) {
    
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    double x1w = comp_data.bulb_data_inter.mol_fracs_bulb1.x1;
    double x2w = comp_data.bulb_data_inter.mol_fracs_bulb1.x2;
    double x1e = 0.5 * (comp_data.tube_fracs_inter[0].x1 + comp_data.tube_fracs_inter[1].x1);
    double x2e = 0.5 * (comp_data.tube_fracs_inter[0].x2 + comp_data.tube_fracs_inter[1].x2);
    
    double x1E = comp_data.tube_fracs[1].x1;
    double x2E = comp_data.tube_fracs[1].x2;
    double x1P = comp_data.tube_fracs_inter[0].x1;
    
    double alpha1w = alpha1(comp_data.p_params, x1w, x2w);
    double alpha2w = alpha2(comp_data.p_params, x1w, x2w);
    double alpha1e = alpha1(comp_data.p_params, x1e, x2e);
    double alpha2e = alpha2(comp_data.p_params, x1e, x2e);
    
    double ap2_tube0 = comp_data.p_params.ct * dz / dt - alpha2w / (0.5 * dz) - alpha2e / dz;
    
    double old_term_tube2 = comp_data.p_params.ct * dz * comp_data.tube_fracs_old[0].x2 / dt;
    
    double x2_tube0 = alpha1w / (0.5 * dz) * (x1P - x1w) - alpha2w / (0.5 * dz) * x2w - alpha1e / dz * (x1E - x1P) - alpha2e / dz * x2E + old_term_tube2;
    
    x2_tube0 = x2_tube0 / ap2_tube0;

    // Store coefficients for analysis
    ap_coeff2[1] = ap2_tube0;
    anb_coeff2[1][0] = alpha2e / dz;
    anb_coeff2[1][1] = alpha2w / (0.5 * dz);
    anb_coeff2[1][2] = 0.0;

    comp_data.tube_fracs[0].x2 = x2_tube0;
    comp_data.tube_fracs[0].x3 = 1.0 - comp_data.tube_fracs[0].x1 - x2_tube0;
}

void mid_nodes1(c_data_t & comp_data) {
    
    int ng = comp_data.ng;
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    for(int node = 1; node < ng - 1; ++node) {
        double x1w = 0.5 * (comp_data.tube_fracs_inter[node - 1].x1 + comp_data.tube_fracs_inter[node].x1);
        double x2w = 0.5 * (comp_data.tube_fracs_inter[node - 1].x2 + comp_data.tube_fracs_inter[node].x2);
        double x1e = 0.5 * (comp_data.tube_fracs_inter[node].x1 + comp_data.tube_fracs_inter[node + 1].x1);
        double x2e = 0.5 * (comp_data.tube_fracs_inter[node].x2 + comp_data.tube_fracs_inter[node + 1].x2);

        double x1W = comp_data.tube_fracs[node - 1].x1;
        double x2W = comp_data.tube_fracs[node - 1].x2;
        double x1E = comp_data.tube_fracs[node + 1].x1;
        double x2E = comp_data.tube_fracs[node + 1].x2;
        double x2P = comp_data.tube_fracs[node].x2;

        double beta1w = beta1(comp_data.p_params, x1w, x2w);
        double beta2w = beta2(comp_data.p_params, x1w, x2w);
        double beta1e = beta1(comp_data.p_params, x1e, x2e);
        double beta2e = beta2(comp_data.p_params, x1e, x2e);

        double ap1_tube_node = comp_data.p_params.ct * dz / dt - beta1w / dz - beta1e / dz;
        
        double old_term_tube1 = comp_data.p_params.ct * dz / dt * comp_data.tube_fracs_old[node].x1;

        double x1_tube_node = - beta1w / dz * x1W + beta2w / dz * (x2P - x2W) - beta1e / dz * x1E - beta2e / dz * (x2E - x2P) + old_term_tube1;

        x1_tube_node = x1_tube_node / ap1_tube_node;

        // Store coefficients for analysis
        ap_coeff1[node + 1] = ap1_tube_node;
        anb_coeff1[node + 1][0] = beta1w / dz;
        anb_coeff1[node + 1][1] = beta1e / dz;
        anb_coeff1[node + 1][2] = 0.0;

        comp_data.tube_fracs[node].x1 = x1_tube_node;
    }
}

void mid_nodes2(c_data_t & comp_data)  {
    
    int ng = comp_data.ng;
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    for(int node = 1; node < ng - 1; ++node) {
        double x1w = 0.5 * (comp_data.tube_fracs_inter[node - 1].x1 + comp_data.tube_fracs_inter[node].x1);
        double x2w = 0.5 * (comp_data.tube_fracs_inter[node - 1].x2 + comp_data.tube_fracs_inter[node].x2);
        double x1e = 0.5 * (comp_data.tube_fracs_inter[node].x1 + comp_data.tube_fracs_inter[node + 1].x1);
        double x2e = 0.5 * (comp_data.tube_fracs_inter[node].x2 + comp_data.tube_fracs_inter[node + 1].x2);

        double x1W = comp_data.tube_fracs[node - 1].x1;
        double x2W = comp_data.tube_fracs[node - 1].x2;
        double x1E = comp_data.tube_fracs[node + 1].x1;
        double x2E = comp_data.tube_fracs[node + 1].x2;
        double x1P = comp_data.tube_fracs[node].x1;

        double alpha1w = alpha1(comp_data.p_params, x1w, x2w);
        double alpha2w = alpha2(comp_data.p_params, x1w, x2w);
        double alpha1e = alpha1(comp_data.p_params, x1e, x2e);
        double alpha2e = alpha2(comp_data.p_params, x1e, x2e);

        double ap2_tube_node = comp_data.p_params.ct * dz / dt - alpha2w / dz - alpha2e / dz;
        
        double old_term_tube2 = comp_data.p_params.ct * dz / dt * comp_data.tube_fracs_old[node].x2;

        double x2_tube_node = alpha1w / dz * (x1P - x1W) - alpha2w / dz * x2W - alpha1e * (x1E - x1P) / dz - alpha2e / dz * x2E + old_term_tube2;

        x2_tube_node = x2_tube_node / ap2_tube_node;

        // Store coefficients for analysis
        ap_coeff2[node + 1] = ap2_tube_node;
        anb_coeff2[node + 1][0] = alpha2w / dz;
        anb_coeff2[node + 1][1] = alpha2e / dz;
        anb_coeff2[node + 1][2] = 0.0;

        comp_data.tube_fracs[node].x2 = x2_tube_node;
    }
}

void tube_n_c1(c_data_t & comp_data) {
    
    int ng = comp_data.ng;
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    double x1w = 0.5 * (comp_data.tube_fracs_inter[ng - 2].x1 + comp_data.tube_fracs_inter[ng - 1].x1);
    double x2w = 0.5 * (comp_data.tube_fracs_inter[ng - 2].x2 + comp_data.tube_fracs_inter[ng - 1].x2);
    double x1e = comp_data.bulb_data_inter.mol_fracs_bulb2.x1;
    double x2e = comp_data.bulb_data_inter.mol_fracs_bulb2.x2;
    
    double x1W = comp_data.tube_fracs[ng - 2].x1;
    double x2W = comp_data.tube_fracs[ng - 2].x2;
    double x1E = x1e;
    double x2E = x2e;
    double x2P = comp_data.tube_fracs[ng - 1].x2;
    
    double beta1w = beta1(comp_data.p_params, x1w, x2w);
    double beta2w = beta2(comp_data.p_params, x1w, x2w);
    double beta1e = beta1(comp_data.p_params, x1e, x2e);
    double beta2e = beta2(comp_data.p_params, x1e, x2e);
    
    double ap1_tube_n = comp_data.p_params.ct * dz / dt - beta1w / dz - beta1e / (0.5 * dz);
    
    double old_term_tube_n_1 = comp_data.p_params.ct * dz / dt * comp_data.tube_fracs_old[ng - 1].x1;
    
    double x1_tube_n = - beta1w / dz * x1W + beta2w / dz * (x2P - x2W) - beta1e / (0.5 * dz) * x1E - beta2e / dz * (x2E - x2P) + old_term_tube_n_1;
    
    x1_tube_n = x1_tube_n / ap1_tube_n;

    // Store coefficients for analysis
    ap_coeff1[ng] = ap1_tube_n;
    anb_coeff1[ng][0] = beta1w / dz;
    anb_coeff1[ng][1] = beta1e / (0.5 * dz);
    anb_coeff1[ng][2] = 0.0;

    comp_data.tube_fracs[ng - 1].x1 = x1_tube_n;
}

void tube_n_c2(c_data_t & comp_data)  {
    
    int ng = comp_data.ng;
    double dz = comp_data.e_params.dz;
    double dt = comp_data.t_params.dt;
    
    double x1w = 0.5 * (comp_data.tube_fracs_inter[ng - 2].x1 + comp_data.tube_fracs_inter[ng - 1].x1);
    double x2w = 0.5 * (comp_data.tube_fracs_inter[ng - 2].x2 + comp_data.tube_fracs_inter[ng - 1].x2);
    double x1e = comp_data.bulb_data_inter.mol_fracs_bulb2.x1;
    double x2e = comp_data.bulb_data_inter.mol_fracs_bulb2.x2;
    
    double x1W = comp_data.tube_fracs[ng - 2].x1;
    double x2W = comp_data.tube_fracs[ng - 2].x2;
    double x1E = x1e;
    double x2E = x2e;
    double x1P = comp_data.tube_fracs[ng - 1].x1;
    
    double alpha1w = alpha1(comp_data.p_params, x1w, x2w);
    double alpha2w = alpha2(comp_data.p_params, x1w, x2w);
    double alpha1e = alpha1(comp_data.p_params, x1e, x2e);
    double alpha2e = alpha2(comp_data.p_params, x1e, x2e);
    
    double ap2_tube_n = comp_data.p_params.ct * dz / dt - alpha2w / dz - alpha2e / (0.5 * dz);
    
    double old_term_tube_n_2 = comp_data.p_params.ct * dz / dt * comp_data.tube_fracs_old[ng - 1].x2;
    
    double x2_tube_n = alpha1w / dz * (x1P - x1W) - alpha2w / dz * x2W - alpha1e / (0.5 * dz) * (x1E - x1P) - alpha2e / (0.5 * dz) * x2E + old_term_tube_n_2;
    
    x2_tube_n = x2_tube_n / ap2_tube_n;

    // Store coefficients for analysis
    ap_coeff2[ng] = ap2_tube_n;
    anb_coeff2[ng][0] = alpha2w / dz;
    anb_coeff2[ng][1] = alpha2e / (0.5 * dz);
    anb_coeff2[ng][2] = 0.0;

    comp_data.tube_fracs[ng - 1].x2 = x2_tube_n;
}

void set_frac_comp3(c_data_t & comp_data) {
    int ng = comp_data.ng;
    
    for(int node = 0; node < ng; ++node) {
        comp_data.tube_fracs[node].x3 = 1.0 - comp_data.tube_fracs[node].x1 - comp_data.tube_fracs[node].x2;
    }
}

void update_composition_estimate(c_data_t & comp_data) {
    
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

void init_stab_data(int ng) {

    anb_coeff1 = new double * [ng + 2];
    ap_coeff1 = new double[ng + 2];
    anb_coeff2 = new double * [ng + 2];
    ap_coeff2 = new double[ng + 2];

    for(int i = 0; i < ng + 2; ++i) {
        anb_coeff1[i] = new double[3];
        anb_coeff2[i] = new double[3];
        
        for(int j = 0; j < 3; ++j) {
            anb_coeff1[i][j] = SMALL_NUM;
            anb_coeff2[i][j] = SMALL_NUM;
        }
        
        ap_coeff1[i] = SMALL_NUM;
        ap_coeff2[i] = SMALL_NUM;
    }
}

void delete_stab_data(int ng) {

    for(int i = 0; i < ng + 2; ++i) {
        delete [] anb_coeff1[i];
        delete [] anb_coeff2[i];
    }

    delete [] ap_coeff1;
    delete [] ap_coeff2;
    delete [] anb_coeff1;
    delete [] anb_coeff2;
}

void compute_bulb_compositions(e_params_t e_params,
                               p_params_t p_params,
                               t_params_t t_params,
                               int ng,
                               b_data_t & bulb_data) {
    
    // Organize data
    c_data_t comp_data;
    comp_data.ng = ng;

    comp_data.p_params = p_params;
    comp_data.e_params = e_params;
    comp_data.t_params = t_params;

    comp_data.bulb_data = bulb_data;
    comp_data.bulb_data_old = bulb_data;
    comp_data.bulb_data_inter = bulb_data;
    
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
    
    // Initialize stability check data
    init_stab_data(ng);

    // Compute composition at t = tf
    int max_out_it = MAX_OUT;
    int max_in_it = MAX_IN;
    
    double t = t_params.to;
    double dt = (t_params.tf - t_params.to) / t_params.nt;

    // Loop to time t = tf
    while(t < t_params.tf) {
        
        // Update bulb composition of previous time step
        comp_data.bulb_data_old = comp_data.bulb_data;
        
        // Update tube composition of previous time step
        for(int node = 0; node < ng; ++node)
            comp_data.tube_fracs_old[node] = comp_data.tube_fracs[node];
        
        // Outer Gauss-Seidel iterations
        int out_it = 0;
        while(out_it < max_out_it) {

            // Update intermediate bulb composition
            comp_data.bulb_data_inter = comp_data.bulb_data;
            
            // Update intermediate tube composition
            for(int node = 0; node < ng; ++node)
                comp_data.tube_fracs_inter[node] = comp_data.tube_fracs[node];
            
            // Inner Gauss-Seidel iterations
            int in_it = 0;
            while(in_it < max_in_it) {

                // Update estimates bulb and tube composition
                update_composition_estimate(comp_data);
                
                in_it++;
            }
            
            // Check stability
            check_stability(ng);
            
            out_it++;
        }
        
        t = t + dt;
    }
    
    // Set bulb data
    bulb_data = comp_data.bulb_data;
    
    // Deallocate data
    delete_stab_data(ng);

    delete [] comp_data.tube_fracs;
    delete [] comp_data.tube_fracs_inter;
    delete [] comp_data.tube_fracs_old;
}

//
//  main.cpp
//  two-bulb-imp
//
//  Created by mndx on 12/03/2022.
//  Three-component two-bulb diffusion
//  using implicit time discretization
//

#include <iostream>

typedef struct node_data {
    double x1;
    double x2;
    double x3;
} node_t;

typedef struct physical_params {
    double ct;
    double D12;
    double D13;
    double D23;
} p_params_t;

typedef struct bulb_data {
    node_t mol_fracs_bulb1;
    node_t mol_fracs_bulb2;
} b_data_t;

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

int main(int argc, const char * argv[]) {
    
    // Number of grid nodes
    int ng = 10;
    
    // Experimental setup parameters
    double V = 5e-4;
    double d = 2e-3;
    double len = 1e-2;
    double A = 0.25 * 3.14 * d * d;
    double dz = len / ng;
    
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
    double to = 0.0;
    double tf = 5.0;
    int nt = 20;
    double dt = (double) (tf - to) / nt;
    
    // Diffusivities
    p_params.D12 = 8.33e-5 * 3600;
    p_params.D13 = 6.8e-5 * 3600;
    p_params.D23 = 1.68e-5 * 3600;
    
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
    int max_out_it = 300;
    int max_in_it = 300;
    
    double t = to;
    
    while(t < tf) {
        
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
                
                double x2P = tube_fracs_inter[0].x2;
                double x1E = tube_fracs[1].x1;
                double x2E = tube_fracs_inter[1].x2;
                
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
                    double x2W = tube_fracs_inter[node - 1].x2;
                    double x1E = tube_fracs[node + 1].x1;
                    double x2E = tube_fracs_inter[node + 1].x2;
                    double x2P = tube_fracs_inter[node].x2;

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

                    double x1W = tube_fracs_inter[node - 1].x1;
                    double x2W = tube_fracs[node - 1].x2;
                    double x1E = tube_fracs_inter[node + 1].x1;
                    double x2E = tube_fracs[node + 1].x2;
                    double x1P = tube_fracs_inter[node].x1;

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
                double x2W = tube_fracs_inter[ng - 2].x2;
                x1E = x1e;
                x2E = x2e;
                x2P = tube_fracs_inter[ng - 1].x2;
                
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
                
                x1W = tube_fracs_inter[ng - 2].x1;
                x2W = tube_fracs[ng - 2].x2;
                x1E = x1e;
                x2E = x2e;
                x1P = tube_fracs_inter[ng - 1].x1;
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
    
    // Check Bulb 1 computation
    double x1_loc = bulb_data.mol_fracs_bulb1.x1;
    double x2_loc = bulb_data.mol_fracs_bulb1.x2;
    double a1_loc = alpha1(p_params, x1_loc, x2_loc);
    double a2_loc = alpha2(p_params, x1_loc, x2_loc);
    double b1_loc = beta1(p_params, x1_loc, x2_loc);
    double b2_loc = beta2(p_params, x1_loc, x2_loc);
    
    double x1P = tube_fracs[0].x1;
    double x1W = x1_loc;
    
    double x2P = tube_fracs[0].x2;
    double x2W = x2_loc;
    
    double J1_0 = b1_loc * (x1P - x1W) / (0.5 * dz) + b2_loc * (x2P - x2W) / (0.5 * dz);
    double J2_0 = a1_loc * (x1P - x1W) / (0.5 * dz) + a2_loc * (x2P - x2W) / (0.5 * dz);
    
    double var1c1 = - J1_0 * A;
    double var2c1 = V * p_params.ct / dt * (bulb_data.mol_fracs_bulb1.x1 - bulb_data_old.mol_fracs_bulb1.x1);
    
    double var1c2 = - J2_0 * A;
    double var2c2 = V * p_params.ct / dt * (bulb_data.mol_fracs_bulb1.x2 - bulb_data_old.mol_fracs_bulb1.x2);
    
    double ratio1 = var1c1 / var2c1;
    double ratio2 = var1c2 / var2c2;
    
    std::cout << "ratio 1 bulb 1: " << ratio1 << std::endl;
    std::cout << "ratio 2 bulb 1: " << ratio2 << std::endl;
    
    // Check first node
    
    // Check mid tube nodes
    for(int node = 1; node < ng - 1; ++node) {
        double x1_loc_z = 0.5 * (tube_fracs[node - 1].x1 + tube_fracs[node].x1);
        double x2_loc_z = 0.5 * (tube_fracs[node - 1].x2 + tube_fracs[node].x2);
        double b1_loc_z = beta1(p_params, x1_loc_z, x2_loc_z);
        double b2_loc_z = beta2(p_params, x1_loc_z, x2_loc_z);
     
        double x1_loc_z_dz = 0.5 * (tube_fracs[node + 1].x1 + tube_fracs[node].x1);
        double x2_loc_z_dz = 0.5 * (tube_fracs[node + 1].x2 + tube_fracs[node].x2);
        double b1_loc_z_dz = beta1(p_params, x1_loc_z_dz, x2_loc_z_dz);
        double b2_loc_z_dz = beta2(p_params, x1_loc_z_dz, x2_loc_z_dz);
        
        double x1P = tube_fracs[node].x1;
        double x1W = tube_fracs[node - 1].x1;
        double x1E = tube_fracs[node + 1].x1;
        
        double x2P = tube_fracs[node].x2;
        double x2W = tube_fracs[node - 1].x2;
        double x2E = tube_fracs[node + 1].x2;
        
        double J1_z = b1_loc_z * (x1P - x1W) / dz + b2_loc_z * (x2P - x2W) / dz;
        double J1_z_dz = b1_loc_z_dz * (x1E - x1P) / dz + b2_loc_z_dz * (x2E - x2P) / dz;
        
        double var1 = p_params.ct * (tube_fracs[node].x1 - tube_fracs_old[node].x1) / dt;
        double var2 = (J1_z - J1_z_dz) / dz;
        
        double res = var1 / var2;
        
        std::cout << "ratio: " << res << std::endl;
        
    }
    
    // Last node
    
    // Check Bulb 2 computation
    x1_loc = bulb_data.mol_fracs_bulb2.x1;
    x2_loc = bulb_data.mol_fracs_bulb2.x2;
    a1_loc = alpha1(p_params, x1_loc, x2_loc);
    a2_loc = alpha2(p_params, x1_loc, x2_loc);
    b1_loc = beta1(p_params, x1_loc, x2_loc);
    b2_loc = beta2(p_params, x1_loc, x2_loc);
    
    x1P = tube_fracs[ng - 1].x1;
    x2P = tube_fracs[ng - 1].x2;
    
    double x1e = x1_loc;
    double x2e = x2_loc;
    
    double J1_n = b1_loc * (x1e - x1P) / (0.5 * dz) + b2_loc * (x2e - x2P) / (0.5 * dz);
    double J2_n = a1_loc * (x1e - x1P) / (0.5 * dz) + a2_loc * (x2e - x2P) / (0.5 * dz);
    
    var1c1 = J1_n * A;
    var2c1 = V * p_params.ct / dt * (bulb_data.mol_fracs_bulb2.x1 - bulb_data_old.mol_fracs_bulb2.x1);
    
    var1c2 = J2_n * A;
    var2c2 = V * p_params.ct / dt * (bulb_data.mol_fracs_bulb2.x2 - bulb_data_old.mol_fracs_bulb2.x2);
    
    ratio1 = var1c1 / var2c1;
    ratio2 = var1c2 / var2c2;
    
    std::cout << "ratio 1 bulb 2: " << ratio1 << std::endl;
    std::cout << "ratio 2 bulb 2: " << ratio2 << std::endl;
    
    // Print stuff
    for(int node = 0; node < ng; ++node) {
        std::cout << tube_fracs[node].x1 << " ";
        std::cout << tube_fracs[node].x2 << " ";
        std::cout << tube_fracs[node].x3;
        std::cout << std::endl;
    }
    
    std::cout << "bulb 1 frac 1: " << bulb_data.mol_fracs_bulb1.x1 << std::endl;
    std::cout << "bulb 1 frac 2: " << bulb_data.mol_fracs_bulb1.x2 << std::endl;
    std::cout << "bulb 1 frac 3: " << bulb_data.mol_fracs_bulb1.x3 << std::endl;
    
    std::cout << "bulb 2 frac 1: " << bulb_data.mol_fracs_bulb2.x1 << std::endl;
    std::cout << "bulb 2 frac 2: " << bulb_data.mol_fracs_bulb2.x2 << std::endl;
    std::cout << "bulb 2 frac 3: " << bulb_data.mol_fracs_bulb2.x3 << std::endl;
    
    // Testing
//    double res1 = 2 * beta1(p_params, 0.1, 0.7) / dz;
//    double res2 = p_params.ct * dz / dt;
//    double res3 = 2 * beta2(p_params, 0.1, 0.7) / dz;
//    double res4 = -beta2(p_params, 0.7, 0.1) / dz;
//
//    double res5 = 2 * alpha1(p_params, 0.1, 0.7) / dz;
//    double res6 = p_params.ct * dz / dt;
//    double res7 = 2 * alpha2(p_params, 0.1, 0.7) / dz;
//    double res8 = -alpha2(p_params, 0.7, 0.1) / dz;
//
//    std::cout << "2 * beta1 / dz: " << res1 << std::endl;
//    std::cout << "ct * dz / dt: " << res2 << std::endl;
//    std::cout << "2 * beta2 / dz: " << res3 << std::endl;
//    std::cout << "-beta2 / dz: " << res4 << std::endl;
//
//    std::cout << "2 * alpha1 / dz: " << res5 << std::endl;
//    std::cout << "2 * alpha2 / dz: " << res7 << std::endl;
//    std::cout << "-alpha2 / dz: " << res8 << std::endl;
    
    std::cout << "done" << std::endl;
    
    return 0;
}

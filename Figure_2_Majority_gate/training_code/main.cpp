//
//  main.cpp
//  Fast_PAOA
//
//  Created by Abdelrahman Said Abdelrahman on February 4, 2025.
//




#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "PAOA_Circuit.h"
#include "load_mat.h"
#include "omp.h"
#include <math.h>
#include "nlopt.h"
#include "nlopt.hpp"
#include <random>
#include <fstream>
#include <tuple>
#include "matio.h"



// Struct to hold optimization parameters
struct OptimizationParams {
    const Eigen::MatrixXd* J_cost; // pointer to cost matrix for the instance
    int num_layers;
    int num_nodes;
    int num_experiments;
    std::vector<double>* cost_history; // Pointer to store cost history

    // OptimizationParams(const Eigen::MatrixXd* jc, int nl, int nexp, int nsl)
    //     : J_cost(jc), num_layers(nl), num_experiments(nexp), num_sweeps_layer(nsl) {}
};


// Objective function for NLopt
double objective_function(const std::vector<double> &beta, std::vector<double> &grad, void *data) {
    auto* opt_params = static_cast<OptimizationParams*>(data);
    const Eigen::MatrixXd &J_cost = *(opt_params->J_cost);
    int num_layers = opt_params->num_layers;
    int num_nodes = opt_params->num_nodes;
    int num_experiments = opt_params->num_experiments;

    // Run experiments with the given beta
    double cost = run_experiments(beta, J_cost, num_layers, num_experiments, num_nodes);
    // Store the cost in the history vector
    opt_params->cost_history->push_back(cost);
    // If printing in a parallel context, consider protecting with #pragma omp critical
    std::cout << "Current Negative Log-likelihood cost: " << cost << std::endl;

    return cost;
}

void save_cost_history(const std::vector<double>& cost_history, const std::string& filename) {
    mat_t* matfp = Mat_CreateVer(filename.c_str(), NULL, MAT_FT_MAT5);
    if (!matfp) {
        std::cerr << "Error creating MAT file " << filename << std::endl;
        return;
    }

    size_t num_elements = cost_history.size();
    double* data = const_cast<double*>(cost_history.data());

    matvar_t* matvar = Mat_VarCreate("cost_history", MAT_C_DOUBLE, MAT_T_DOUBLE, 1, &num_elements, data, 0);
    if (!matvar) {
        std::cerr << "Error creating variable for MAT file" << std::endl;
        Mat_Close(matfp);
        return;
    }

    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
    Mat_VarFree(matvar);
    Mat_Close(matfp);
}



int main() {
    const int num_experiments = int(1e7); // number of experiments
    std::vector<int> num_layers_vec = {2}; // number of layers
    int num_nodes = 4; // Number of nodes in the graph
    int num_instances = 1;
    for (int instance = 1; instance <= num_instances; ++instance) {
        // Load J that represents the graph weights.
        std::string filename =  "J_weight.mat"; //std::to_string(instance) +
        std::string varname = "J";
        Eigen::MatrixXd J_cost = loadMatFile(filename, varname);

        // Check if loading was successful
        if (J_cost.size() == 0) {
            std::cerr << "Error: Loaded matrix J_cost is empty for file " << filename << std::endl;
            continue;
        }

        // prepare text files for saving optimized beta values
        std::string beta_filename = "optimized_beta_instance_" + std::to_string(instance) + ".txt";
        std::ofstream beta_file(beta_filename);
        if (!beta_file.is_open()) {
            std::cerr << "Failed to open file for writing optimized beta values: " << beta_filename << std::endl;
            continue;
        }

        // Initialize cost history vector
        std::vector<double> cost_history;

        for (int nl : num_layers_vec) {
            int num_layers = nl;

            // Initialize beta
            std::vector<double> beta(num_nodes*num_layers, 1.0);

            // Set up optimizer
            nlopt::opt opt(nlopt::LN_COBYLA, beta.size());
            OptimizationParams opt_params = {&J_cost, num_layers,num_nodes, num_experiments, &cost_history};
            opt.set_min_objective(objective_function, &opt_params);
            opt.set_xtol_abs(1e-6);
            opt.set_maxeval(7000); // Maximum number of evaluations
            opt.set_initial_step(0.1);

            double minf;
            try {
                nlopt::result result = opt.optimize(beta, minf);

                // Output the results
                std::cout << "Optimization result: " << result << std::endl;
                std::cout << "Instance " << instance << ", num_layers " << num_layers << std::endl;
                std::cout << "Minimum cost: " << minf << std::endl;
                // Save cost history to a .mat file
                save_cost_history(cost_history, "cost_history.mat");
                
            

                

                // Save the optimized beta values
                for (double b : beta) {
                    beta_file << b << " ";
                }
                beta_file << "\n";
            }
            catch (std::exception& e) {
                std::cerr << "Optimization failed for instance " << instance << ", num_layers " << num_layers << ": " << e.what() << std::endl;
            }
        }

        // Close the text file and beta file
        beta_file.close();
       
    }

    return 0;
}
    



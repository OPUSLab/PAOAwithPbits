//
//  main.cpp
//  Fast_PAOA
//
//  Created by Abdelrahman Said Abdelrahman on February 4, 2025.
//
// Building the code : 
// ******************* with parallelization *******************




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
    const Eigen::VectorXi* rank;    // <<< added
    int num_layers;
    int num_experiments;
    int num_sweeps_layer;

    OptimizationParams(const Eigen::MatrixXd* jc, const Eigen::VectorXi* ri, int nl, int nexp, int nsl)
        : J_cost(jc), 
          rank(ri), 
          num_layers(nl), 
          num_experiments(nexp), 
          num_sweeps_layer(nsl) {}
};


// Objective function for NLopt
double objective_function(const std::vector<double> &beta, std::vector<double> &grad, void *data) {
    auto* opt_params = static_cast<OptimizationParams*>(data);
    const Eigen::MatrixXd &J_cost = *(opt_params->J_cost);
    const Eigen::VectorXi &rank = *(opt_params->rank);   // <<< added
    int num_layers = opt_params->num_layers;
    int num_experiments = opt_params->num_experiments;
    int num_sweeps_layer = opt_params->num_sweeps_layer;

    // Run experiments with the given beta
    double average_cost = run_experiments(beta, J_cost,rank, num_layers, num_experiments, num_sweeps_layer);

    // If printing in a parallel context, consider protecting with #pragma omp critical
    std::cout << "Current average cost: " << average_cost << std::endl;

    return average_cost;
}



int main() {
    const int num_experiments = 1e5; // 1e6
    std::vector<int> num_layers_vec = {17};
    int num_nodes = 50; // Number of nodes in the graph
    int num_instances = 50;
    int num_sweeps_layer = 1;
    // int num_runs = 100;
    for (int instance = 1; instance <= num_instances; ++instance) {
        // Load J_cost from the file named "1.mat", "2.mat", etc.
        std::string filename = std::to_string(instance) + ".mat";
        std::string varname = "J";
        Eigen::MatrixXd J_cost = loadMatFile(filename, varname);
        // Check if loading was successful
        if (J_cost.size() == 0) {
            std::cerr << "Error: Loaded matrix J_cost is empty for file " << filename << std::endl;
            continue;
        }

        // load the sorted indices
        std::string filename2 = "indices_instance_" + std::to_string(instance) + ".mat";
        std::string varname2 = "sorted_indices";
        Eigen::VectorXd sorted_indices = loadMatFile(filename2, varname2);
        // Check if loading was successful
        if (sorted_indices.size() == 0) {
            std::cerr << "Error: Loaded vector sorted_indices is empty for file " << filename << std::endl;
            continue;
        }

        Eigen::VectorXi rank(J_cost.rows());
        for (int r = 0; r < sorted_indices.size(); ++r) {
            rank[ static_cast<int>(sorted_indices[r]) ] = r;   // rank[i] = position in list
        }

        // prepare text files for saving optimized beta values
        std::string beta_filename = "optimized_beta_" + std::to_string(instance) + ".txt";
        std::ofstream beta_file(beta_filename);
        if (!beta_file.is_open()) {
            std::cerr << "Failed to open file for writing optimized beta values: " << beta_filename << std::endl;
            continue;
        }

        // Prepare text file for saving minf values
        std::string text_filename = "avg_energy_per_spin_instance_" + std::to_string(instance) + ".txt";
        std::ofstream text_file(text_filename);
        if (!text_file.is_open()) {
            std::cerr << "Error creating text file " << text_filename << std::endl;
            continue;
        }

        for (int nl : num_layers_vec) {
            int num_layers = nl;

            // Initialize beta
            // std::vector<double> beta(2*num_layers, 0.1);
            double start = std::log(0.005);
            double end = std::log(0.12);
            std::vector<double> beta(2*num_layers);

            for (int i = 0; i < num_layers; ++i) {
                double t = static_cast<double>(i) / (num_layers - 1); // normalized [0,1]
                double value = std::exp(start + t * (end - start));
                beta[i] = value;
                beta[num_layers + i] = value;     // second schedule (concatenated)
            }

            // Set up optimizer
            nlopt::opt opt(nlopt::LN_COBYLA, beta.size());
            OptimizationParams opt_params(&J_cost, &rank, num_layers, num_experiments, num_sweeps_layer);
            opt.set_min_objective(objective_function, &opt_params);
            opt.set_xtol_abs(1e-5);
            opt.set_maxeval(5000); // Maximum number of evaluations

            double minf;
            try {
                nlopt::result result = opt.optimize(beta, minf);

                std::cout << " instance " << instance << ", num_layers " << num_layers << std::endl;
                std::cout << "Optimization result: " << result << std::endl;
                std::cout << "Minimum cost: " << minf << std::endl;

                // Save minf value to the text file
                text_file << minf/num_nodes << std::endl;

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
        text_file.close();
    }

    return 0;
}
    



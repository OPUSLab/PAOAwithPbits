//
//  PAOA_circuit.cpp
//  Fast_PAOA
//
//  Created by Abdelrahman Said Abdelrahman on February 4, 2025.
//


#include "omp.h"
#include "PAOA_Circuit.h"
#include <random>
#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <thread>




// Function to simulate the PAOA Circuit
double PAOA_Circuit(const std::vector<double>& beta, const Eigen::MatrixXd& J_cost, int num_layers,int num_sweeps_layer, std::mt19937& gen) {
    
    
    
    // Initialize the state vector 'm' with all ones
    Eigen::VectorXd m = Eigen::VectorXd::Constant(J_cost.rows(), 1);

    // Uniform distribution between -1 and 1
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    // beta is a flat std::vector<double> of length N * p

    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > BetaMat(beta.data(), 2, num_layers);
    // Now BetaMat is an 2 x p matrix where
    //  - row 0 = beta[0..p-1]
    //  - row 1 = beta[p..2p-1]
    


    int num_sweeps = num_layers * num_sweeps_layer;
    int idx = -1;
    double beta_value= 0;
    // Galuber Criterion
    for (int i = 0; i < num_sweeps; i++) {
        if (i % num_sweeps_layer == 0) {          
            idx++;
        }
        for (int j = 0; j < J_cost.rows(); j++) {
            if (j<J_cost.rows()/2){
                beta_value = BetaMat(0,idx);
            }
            else {
                beta_value = BetaMat(1,idx);    
            }
            double I = beta_value * (J_cost.row(j).dot(m));
            m[j] = std::tanh(I) > dis(gen) ? 1 : -1; 
        }
    }

    // Compute the cost using J_cost
    double cost = -0.5 * m.dot(J_cost * m)/sqrt(J_cost.rows());
    return cost;
}



double run_experiments(const std::vector<double>& beta, const Eigen::MatrixXd& J_cost, int num_layers, int num_experiments, int num_sweeps_layer) {
    double global_sum = 0.0;

    // OpenMP parallelization with reduction to avoid race conditions
    #pragma omp parallel 
    {
        // Thread-local random number generator
        std::random_device rd;
        int seed = rd() + omp_get_thread_num();
        std::mt19937 gen(seed);
        std::uniform_real_distribution<> dis(-1.0, 1.0);

        #pragma omp for reduction(+:global_sum)
        for (int i = 0; i < num_experiments; i++) {
            double cost = PAOA_Circuit(beta, J_cost, num_layers,num_sweeps_layer, gen);
            global_sum += cost;
        }
    }
    
    double average_cost = global_sum /num_experiments;
    return average_cost;
}

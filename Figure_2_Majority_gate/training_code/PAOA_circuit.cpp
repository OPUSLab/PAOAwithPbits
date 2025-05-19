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

// Target states: Majority gate truth table
const std::vector<int> target_states = {0, 2, 4, 7, 8, 11, 13, 15};

// Function to simulate the PAOA Circuit
int PAOA_Circuit(const std::vector<double>& beta, const Eigen::MatrixXd& J_cost, int num_layers, int num_nodes, std::mt19937& gen) {
    
    // Initialize the state vector 'm' with random spins (-1 or +1)
    Eigen::VectorXd m = ((Eigen::VectorXd::Random(num_nodes).array() > 0).cast<double>() * 2.0 - 1.0).matrix();
    
    
    // Uniform distribution between -1 and 1
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    // beta is a flat std::vector<double> of length N * p

    // Map the data from beta into an Eigen matrix (N x p) in row-major order
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > BetaMat(beta.data(), J_cost.rows(), num_layers);
    // Now BetaMat is an N x p matrix where
    //  - row 0 = beta[0..p-1]
    //  - row 1 = beta[p..2p-1]
    //  - ...
    //  - row N-1 = beta[(N-1)*p..N*p - 1]


    // Galuber Criterion
    for (int i = 0; i < num_layers; i++) {
        for (int j = 0; j < J_cost.rows(); j++) {
            double I = BetaMat(j,i) * (J_cost.row(j).dot(m));
            m[j] = std::tanh(I) > dis(gen) ? 1 : -1; 
        }
    }

     // Convert bipolar vector 'm' to binary and then to decimal
     int decimal_state = 0;
     for (int i = 0; i < num_nodes; ++i) {
         int bit = (m[i] == 1.0) ? 1 : 0;
         decimal_state |= (bit << (num_nodes - i - 1));
     }
 
     return decimal_state;
}




double run_experiments(const std::vector<double>& beta, const Eigen::MatrixXd& J_cost, int num_layers, int num_experiments, int num_nodes) {
    std::vector<int> state_counts(1 << num_nodes, 0); // Initialize counts for all possible states

    // OpenMP parallelization
    #pragma omp parallel
    {

        // Thread-local random number generator
        std::random_device rd;
        std::mt19937 gen(rd());

        std::vector<int> local_state_counts(state_counts.size(), 0);

        #pragma omp for
        for (int i = 0; i < num_experiments; i++) {
            int state = PAOA_Circuit(beta, J_cost, num_layers, num_nodes, gen);
            local_state_counts[state]++;
        }

        // Combine local counts into global counts
        #pragma omp critical
        {
            for (size_t i = 0; i < state_counts.size(); ++i) {
                state_counts[i] += local_state_counts[i];
            }
        }
    }

    // Compute probabilities of target states
    std::vector<double> prob_v;
    for (int s : target_states) {
        double prob = static_cast<double>(state_counts[s]) / num_experiments;
        prob_v.push_back(prob);
    }

    // Compute negative log-likelihood cost
    double cost = 0.0;
    for (double p : prob_v) {
        if (p > 0) {
            cost -= std::log(p);
        } else {
            cost -= std::log(1e-10); // Small value to avoid log(0)
        }
    }

    return cost;
    
}


#ifndef PAOA_CIRCUIT_H
#define PAOA_CIRCUIT_H

#include <vector>
#include <random>
#include "Eigen/Sparse"
#include "Eigen/Dense"

// Eigen::MatrixXd
int PAOA_Circuit(const std::vector<double>& beta, const Eigen::MatrixXd& J_cost, int num_layers, int num_nodes, std::mt19937& gen);
double run_experiments(const std::vector<double>& beta, const Eigen::MatrixXd& J_cost, int num_layers, int num_experiments, int num_nodes);

#endif // PAOA_CIRCUIT_H

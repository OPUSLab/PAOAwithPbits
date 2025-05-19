#ifndef LOAD_MAT_H
#define LOAD_MAT_H

#include "Eigen/Dense"
#include <string>

Eigen::MatrixXd loadMatFile(const std::string& filename, const std::string& varname);

#endif // LOAD_MAT_H

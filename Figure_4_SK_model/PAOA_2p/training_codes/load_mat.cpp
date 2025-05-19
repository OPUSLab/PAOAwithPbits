#include "load_mat.h"
#include <matio.h>
#include <iostream>

Eigen::MatrixXd loadMatFile(const std::string& filename, const std::string& varname) {
    mat_t *matfp;
    matvar_t *matvar;

    // Open MAT file
    matfp = Mat_Open(filename.c_str(), MAT_ACC_RDONLY);
    if (matfp == NULL) {
        std::cerr << "Error opening MAT file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    // Read the variable
    matvar = Mat_VarRead(matfp, varname.c_str());
    if (matvar == NULL) {
        std::cerr << "Error reading variable: " << varname << std::endl;
        Mat_Close(matfp);
        exit(EXIT_FAILURE);
    }

    // Check if the variable is of type double
    if (matvar->data_type != MAT_T_DOUBLE) {
        std::cerr << "Variable is not of type double: " << varname << std::endl;
        Mat_VarFree(matvar);
        Mat_Close(matfp);
        exit(EXIT_FAILURE);
    }

    // Get dimensions of the variable
    size_t rows = matvar->dims[0];
    size_t cols = matvar->dims[1];

    // Create an Eigen matrix to hold the data
    Eigen::MatrixXd matData(rows, cols);

    // Copy the data to the Eigen matrix
    double *data = static_cast<double*>(matvar->data);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            matData(i, j) = data[i + j * rows];
        }
    }

    // Free the MAT variable and close the MAT file
    Mat_VarFree(matvar);
    Mat_Close(matfp);

    return matData;
}

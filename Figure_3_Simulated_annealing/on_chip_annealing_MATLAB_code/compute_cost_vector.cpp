// compute_cost_vector.cpp

#include "mex.h"
#include <Eigen/Dense>
#include <omp.h>

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    // Input validation
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("compute_cost_vector:nrhs",
                          "Four inputs required: s (state matrix), J (coupling matrix), num_pbits, num_replicas.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("compute_cost_vector:nlhs",
                          "One output required: cost vector.");
    }
    // Set number of threads to maximum available
    omp_set_num_threads(omp_get_max_threads());
    // Get input matrices and parameters
    double *s_ptr = mxGetPr(prhs[0]);
    mwSize num_instances = mxGetM(prhs[0]); // Number of instances
    mwSize s_cols = mxGetN(prhs[0]);        // Should be num_pbits * num_replicas

    double *J_ptr = mxGetPr(prhs[1]);

    double num_pbits_dbl = mxGetScalar(prhs[2]);
    int num_pbits = static_cast<int>(num_pbits_dbl);

    double num_replicas_dbl = mxGetScalar(prhs[3]);
    int num_replicas = static_cast<int>(num_replicas_dbl);

    // Validate sizes
    if (s_cols != num_pbits * num_replicas) {
        mexErrMsgIdAndTxt("compute_cost_vector:sSize",
                          "Number of columns in s must be num_pbits * num_replicas.");
    }

    mwSize J_rows = mxGetM(prhs[1]);
    mwSize J_cols = mxGetN(prhs[1]);
    if (J_rows != num_pbits || J_cols != num_pbits) {
        mexErrMsgIdAndTxt("compute_cost_vector:JSize",
                          "J must be a square matrix of size num_pbits x num_pbits.");
    }

    // Map data to Eigen matrices
    Map<MatrixXd> s_mat(s_ptr, num_instances, num_pbits * num_replicas);
    Map<MatrixXd> J(J_ptr, num_pbits, num_pbits);

    // Prepare output vector
    mwSize total_costs = num_instances * num_replicas;
    plhs[0] = mxCreateDoubleMatrix(total_costs, 1, mxREAL);
    double *cost_ptr = mxGetPr(plhs[0]);

    // Compute the costs
    #pragma omp parallel for schedule(static)
    for (mwSize i = 0; i < num_instances; ++i) {
        RowVectorXd s_row = s_mat.row(i); // Size: num_pbits * num_replicas
        for (int r = 0; r < num_replicas; ++r) {
            // Extract and convert m for this replica
            RowVectorXd m = s_row.segment(r * num_pbits, num_pbits);
            m = 2.0 * m.array() - 1.0;

            // Compute energy: E = -0.5 * m * J * m^T
            double energy = -0.5 * m * J * m.transpose();

            // Store energy
            mwSize cost_idx = i * num_replicas + r;
            cost_ptr[cost_idx] = energy;
        }
    }
}

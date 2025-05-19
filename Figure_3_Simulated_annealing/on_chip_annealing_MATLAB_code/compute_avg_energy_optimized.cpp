// compute_avg_energy_optimized.cpp

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
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("compute_avg_energy:nrhs",
                          "Two inputs required: M (states) and J (coupling matrix).");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("compute_avg_energy:nlhs",
                          "One output required: average energy.");
    }
    omp_set_num_threads(omp_get_max_threads());

    // Get input matrices
    double *M_ptr = mxGetPr(prhs[0]);
    mwSize N = mxGetM(prhs[0]); // Number of states
    mwSize P = mxGetN(prhs[0]); // Number of p-bits

    double *J_ptr = mxGetPr(prhs[1]);
    if (mxGetM(prhs[1]) != P || mxGetN(prhs[1]) != P) {
        mexErrMsgIdAndTxt("compute_avg_energy:JSize",
                          "J must be a square matrix of size P x P.");
    }

    // Map MATLAB data to Eigen matrices
    Map<MatrixXd> M_matlab(M_ptr, N, P); // MATLAB uses column-major order
    Map<MatrixXd> J(J_ptr, P, P);

    // Convert M from binary (0,1) to bipolar (-1,1)
    // We can do this in-place to save memory
    Eigen::MatrixXd M = 2.0 * M_matlab.array() - 1.0;

    // Compute the energies in parallel
    double sumE = 0.0;

    // Use OpenMP for parallelization
    #pragma omp parallel
    {
        // Each thread has its own partial sum
        double partial_sumE = 0.0;

        #pragma omp for nowait schedule(static)
        for (mwSize i = 0; i < N; ++i) {
            // Access the i-th state vector
            RowVectorXd m = M.row(i);

            // Compute the energy: E_i = m * J * m^T
            // Since J is symmetric, we can use efficient matrix multiplication
            double energy = -0.5*m * (J * m.transpose());

            partial_sumE += energy;
        }

        // Safely accumulate the partial sums
        #pragma omp atomic
        sumE += partial_sumE;
    }

    // Calculate average energy
    double avgE = sumE / static_cast<double>(N);

    // Set output
    plhs[0] = mxCreateDoubleScalar(avgE);
}
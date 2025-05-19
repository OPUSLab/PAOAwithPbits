#include "mex.h"
#include <stdint.h>

// The MEX function that will replace int2bit/de2bi
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Ensure correct number of input/output arguments
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:bit_conversion_mex:nrhs", "Two inputs required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("MyToolbox:bit_conversion_mex:nlhs", "Only one output allowed.");
    }

    // Get the input data and its size
    uint32_t *s_part = (uint32_t *)mxGetData(prhs[0]);
    mwSize num_elements = mxGetNumberOfElements(prhs[0]);

    // Get msbfirst flag
    bool msbfirst = mxGetScalar(prhs[1]);

    // Prepare the output array
    mwSize num_bits = 32;
    mwSize out_size[2] = {num_elements * num_bits, 1};
    plhs[0] = mxCreateLogicalMatrix(out_size[0], 1);
    mxLogical *out_bits = mxGetLogicals(plhs[0]);

    // Iterate over each uint32 value and convert it to bits
    for (mwSize i = 0; i < num_elements; i++) {
        uint32_t value = s_part[i];
        if (msbfirst) {
            // MSB first
            for (int j = 0; j < 32; j++) {
                out_bits[i * 32 + j] = (value >> (31 - j)) & 1;
            }
        } else {
            // LSB first
            for (int j = 0; j < 32; j++) {
                out_bits[i * 32 + j] = (value >> j) & 1;
            }
        }
    }
}

extern "C" {

    __global__ void setBC(double *var, const unsigned int *idx, const unsigned int size, const double val)
    {
        const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

        if ( i >= size) return;
        var[idx[i]] = val;
    }

    __global__ void copyBC(double *var, const unsigned int *idx1, const unsigned int *idx2, const unsigned int size)
    {
        const unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

        if ( i >= size) return;
        var[idx1[i]] = var[idx2[i]];
    }



}
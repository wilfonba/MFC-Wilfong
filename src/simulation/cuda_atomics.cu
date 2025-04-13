// cuda_atomics.cu
#include <cuda.h>
#include <cuda_fp16.h>

extern "C" {
    __device__ void atomicAdd_half(__half* address, float val) {
        atomicAdd(address, __float2half(val));
        return;
    }
}

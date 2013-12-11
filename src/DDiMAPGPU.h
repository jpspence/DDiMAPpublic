#ifndef DDIMAPGPU_H

#define DDIMAPGPU_H

#include "DDiMAP-lib.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>

__global__ void convert_kernel(Read *bam_data);
int run();

#endif

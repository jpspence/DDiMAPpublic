//============================================================================
// Name        : DDiMAP.cpp
// Author      : Androwis Abumoussa
// Version     :
// Copyright   : All Rights Reserved
// Description : DDiMAP in C
//============================================================================
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>
#include "DDiMAP-lib.h"
#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <getopt.h>
#include <time.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

using namespace BamTools;
using namespace std;

// Default file.
// char *file = "data/Burack_128F/128F_Gen1_Frag_WithBcl2Sanger_sorted.bam";
string file  = "data/128test_Gen1_example_sorted.bam";

/******************************************************************************
 *                                                                        GPU
 ******************************************************************************/

__device__ long long stringToUINT64GPU( char *s) 
{

        long long a = 1;
        long long c = 2;
        long long g = 3;
        long long t = 4;
        long long dash = 7;

        long long temp = 0;
        for( int i = 0; i < 17; i++){
                temp+= (s[i] == 'A') ? a         << (3*i) : 0;
                temp+= (s[i] == 'C') ? c         << (3*i) : 0;
                temp+= (s[i] == 'G') ? g         << (3*i) : 0;
                temp+= (s[i] == 'T') ? t         << (3*i) : 0;
                temp+= (s[i] == '-') ? dash << (3*i) : 0;
        }

        return temp;
}

__global__ void convert_kernel(Read *bam_data)
{
        int idx = blockIdx.x * blockDim.x + threadIdx.x;

        // Read in a single read from global memory
        Read ba = bam_data[idx];

        char left[17];
        char right[17];

        for(int i =0; i<17; i++){
                left[i]  = ba.sequence[i];
                right[i] = ba.sequence[i+17];
        }

        ba.left_sequence_half  = stringToUINT64GPU(left);
        ba.right_sequence_half = stringToUINT64GPU(right);

        // Save the read to global memory.
        bam_data[idx] = ba;

}

//void convertReads(){
//        
//        BamReader *br = new BamReader();
//        br->Open(file);
//
//        BamAlignment ba;
//        int counter = 0;
//        ba
//        while(br->GetNextAlignment(ba)){
//                read( ba, 34);
//                counter++;
//        }
//        
//        cudaStream_t stream1;
//        cudaError_t result;
//        
//        result = cudaStreamCreate(&stream1);
//        result = cudaMemcpyAsync(d_a, a, N, cudaMemcpyHostToDevice, stream1);
//        result = increment<<<1,N,0,stream1>>>(d_a);
//        result = cudaStreamDestroy(stream1);
//        cudaDeviceSynchronize();
//}


/******************************************************************************
 *                                                                        CPU
 ******************************************************************************/

long n = 1024 * 1024;

int correct_output( Read *gpu)
{
        int length = 34;
        
        BamReader *br = new BamReader();
        br->Open(file);
        BamAlignment ba;
        ba.Position = -1;

        for (int i = 0; i < n; i++){
                while(ba.Position < 0)
                  br->GetNextAlignment(ba);
                int offset    = (ba.IsReverseStrand()) ? ba.AlignedBases.length() - length : 0 ;
                string word   = ba.AlignedBases.substr(offset, length);
                Read bam = buildRead(word, 34);

                if (  bam.left_sequence_half  != gpu[i].left_sequence_half || 
                         bam.right_sequence_half  != gpu[i].right_sequence_half)
                {
                        cout << "Error : "<< i << endl;
                        cout << "GPU left = " << gpu[i].left_sequence_half << " | GPU Right = " << gpu[i].right_sequence_half << endl;
                        cout << "CPU left = " << bam.left_sequence_half << " | CPU Right = " << bam.right_sequence_half << endl;
                        return 0;
                }
        }

        return 1;
}

int main (int argc, char **argv) {


        // ------------------------------------------------------------------------
        // Parameters
        // ------------------------------------------------------------------------
        int c;

        static struct option long_options[] = {
                        {"file",         0, 0, 'f'},
                        {NULL,                 0, NULL, 0}
        };

        int option_index = 0;
        while ((c = getopt_long(argc, argv, "f:", long_options, &option_index)) != -1) {

                switch (c) {
                case 'f':
                        printf ("Parsing file :  %s \n",optarg);
                        file = optarg;
                        break;
                default:
                        printf ("?? getopt returned character code 0%o ??\n", c);
                }
        }
        if (optind < argc) {
                printf ("non-option ARGV-elements: ");
                while (optind < argc)
                        printf ("%s ", argv[optind++]);
                printf ("\n");
        }

        // ------------------------------------------------------------------------
        // Setup
        // ------------------------------------------------------------------------

        int devID;
        cudaDeviceProp deviceProps;

        printf("[%s] - Starting...\n", argv[0]);

        // This will pick the best possible CUDA capable device
        devID = findCudaDevice(argc, (const char **)argv);

        // get device name
        checkCudaErrors(cudaGetDeviceProperties(&deviceProps, devID));
        printf("CUDA device [%s]\n", deviceProps.name);

        long alignmentBytes = n * sizeof(Read);
        long aBytes = n * sizeof(BamAlignment);

        // allocate host memory
        Read *a = 0;
        checkCudaErrors(cudaMallocHost((void **)&a, alignmentBytes));
        memset(a, 0, alignmentBytes);

        BamAlignment *alignments = 0;
        checkCudaErrors(cudaMallocHost((void **)&alignments, aBytes));
        memset(alignments, 0, aBytes);

        // allocate device memory
        Read *d_alignments=0;
        checkCudaErrors(cudaMalloc((void **)&d_alignments, alignmentBytes));
        checkCudaErrors(cudaMemset(d_alignments, 0, alignmentBytes));

        // set kernel launch configuration
        dim3 threads = dim3(512, 1);
        dim3 blocks  = dim3(n / threads.x, 1);

        // create cuda event handles
        cudaEvent_t start, stop;
        checkCudaErrors(cudaEventCreate(&start));
        checkCudaErrors(cudaEventCreate(&stop));

        StopWatchInterface *timer = NULL;
        sdkCreateTimer(&timer);
        sdkResetTimer(&timer);
        checkCudaErrors(cudaDeviceSynchronize());

        float gpu_time = 0.0f;

        // ------------------------------------------------------------------------
        // ASYNC DDiMAP Kernel Execution
        // ------------------------------------------------------------------------

        // Read the bamfile
        BamReader *br = new BamReader();
        br->Open(file);
        BamAlignment ba;
        ba.Position = -1;
        int counter = 0;
	int length = 34;
        while(counter < n ){
                while(ba.Position < 0)
                        br->GetNextAlignment(ba);
                int offset    = (ba.IsReverseStrand()) ? ba.AlignedBases.length() - length : 0 ;
                string word   = ba.AlignedBases.substr(offset, length);
                a[counter] = convert(word,34);
                counter++;
        }
        br->Close();
        
        // asynchronously issue work to the GPU (all to stream 0)
        sdkStartTimer(&timer);
        cudaEventRecord(start, 0);

        cudaMemcpyAsync(d_alignments, a, alignmentBytes, cudaMemcpyHostToDevice, 0);
        convert_kernel<<<blocks, threads, 0, 0>>>(d_alignments);
        cudaMemcpyAsync(a, d_alignments, alignmentBytes, cudaMemcpyDeviceToHost, 0);
        cudaEventRecord(stop, 0);
        sdkStopTimer(&timer);


        // have CPU do some work while waiting for stage 1 to finish
        unsigned long int counter2=0;
        while (cudaEventQuery(stop) == cudaErrorNotReady)
        {
                counter2++;
        }
        checkCudaErrors(cudaEventElapsedTime(&gpu_time, start, stop));

        // ------------------------------------------------------------------------
        // Check Correctness. 
        // ------------------------------------------------------------------------

        // print the cpu and gpu times
        printf("time spent executing by the GPU: %.2f (ms) \n", gpu_time);
        printf("time spent by CPU in CUDA calls: %.2f (ms) \n", sdkGetTimerValue(&timer));
        printf("CPU executed %lu iterations while waiting for GPU to finish\n", counter2);

        // check the output for correctness
        bool bFinalResults = (bool) correct_output(a);


        // ------------------------------------------------------------------------
        // End. 
        // ------------------------------------------------------------------------
        // release resources
        checkCudaErrors(cudaEventDestroy(start));
        checkCudaErrors(cudaEventDestroy(stop));
        checkCudaErrors(cudaFreeHost(a));
        checkCudaErrors(cudaFree(d_alignments));

        cudaDeviceReset();

        exit(bFinalResults ? EXIT_SUCCESS : EXIT_FAILURE);
}

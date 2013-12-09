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

using namespace std;

// Default file.
// char *file = "data/Burack_128F/128F_Gen1_Frag_WithBcl2Sanger_sorted.bam";
string file  = "data/128test_Gen1_example_sorted.bam";

/******************************************************************************
 *									GPU
 ******************************************************************************/

//__device__ long long stringToUINT64( string s) 
//{
//
//	long long a = 1;
//	long long c = 2;
//	long long g = 3;
//	long long t = 4;
//	long long dash = 7;
//	
//	long long temp = 0;
//	for( int i = 0; i < s.length(); i++){
//		temp+= (s[i] == 'A') ? a 	<< (3*i) : 0;
//		temp+= (s[i] == 'C') ? c 	<< (3*i) : 0;
//		temp+= (s[i] == 'G') ? g 	<< (3*i) : 0;
//		temp+= (s[i] == 'T') ? t 	<< (3*i) : 0;
//		temp+= (s[i] == '-') ? dash << (3*i) : 0;
//	}
//	
//	return temp;
//}

__global__ void increment_kernel(int *g_data, int inc_value)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	g_data[idx] = g_data[idx] + inc_value;
}
//
//__global__ void read(BamAlignment ba)
//{
//	
//	unsigned int tidx = threadIdx.x + blockDim.x*blockIdx.x;
//
//	// Read in a single read from global memory
//
//	if(ba.Position > 0)
//	{		
//		int name 	  = ba.RefID;
//		int offset    = (ba.IsReverseStrand()) ? ba.AlignedBases.length() - length : 0 ;
//		int position  = ba.Position + offset;
//		string word   = ba.AlignedBases.substr(offset, length);
//
//		Read r;
//		r.count = 1;
//		r.verification_flags = 0;
//		r.left_sequence_half  = stringToUINT64(word.substr(0, length/2));
//		r.right_sequence_half = stringToUINT64(word.substr(lenghth/2, length/2));
//		
//	}
//	
//	// Save the read to global memory.
//}

//void convertReads(){
//	
//	BamReader *br = new BamReader();
//	br->Open(file);
//
//	BamAlignment ba;
//	int counter = 0;
//	ba
//	while(br->GetNextAlignment(ba)){
//		read( ba, 34);
//		counter++;
//	}
//	
//	cudaStream_t stream1;
//	cudaError_t result;
//	
//	result = cudaStreamCreate(&stream1);
//	result = cudaMemcpyAsync(d_a, a, N, cudaMemcpyHostToDevice, stream1);
//	result = increment<<<1,N,0,stream1>>>(d_a);
//	result = cudaStreamDestroy(stream1);
//	cudaDeviceSynchronize();
//}



int correct_output(int *data, const int n, const int x)
{
	for (int i = 0; i < n; i++)
		if (data[i] != x)
		{
			printf("Error! data[%d] = %d, ref = %d\n", i, data[i], x);
			return 0;
		}

	return 1;
}



/******************************************************************************
 *									CPU
 ******************************************************************************/

int main (int argc, char **argv) {


	// ------------------------------------------------------------------------
	// Parameters
	// ------------------------------------------------------------------------
	int c;

	static struct option long_options[] = {
			{"file", 	0, 0, 'f'},
			{NULL, 		0, NULL, 0}
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
	// DDiMAP
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// End. 
	// ------------------------------------------------------------------------
	
	int devID;
	cudaDeviceProp deviceProps;

	printf("[%s] - Starting...\n", argv[0]);

	// This will pick the best possible CUDA capable device
	devID = findCudaDevice(argc, (const char **)argv);

	// get device name
	checkCudaErrors(cudaGetDeviceProperties(&deviceProps, devID));
	printf("CUDA device [%s]\n", deviceProps.name);

	int n = 16 * 1024 * 1024;
	int nbytes = n * sizeof(int);
	int value = 26;

	// allocate host memory
	int *a = 0;
	checkCudaErrors(cudaMallocHost((void **)&a, nbytes));
	memset(a, 0, nbytes);

	// allocate device memory
	int *d_a=0;
	checkCudaErrors(cudaMalloc((void **)&d_a, nbytes));
	checkCudaErrors(cudaMemset(d_a, 255, nbytes));

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

	// asynchronously issue work to the GPU (all to stream 0)
	sdkStartTimer(&timer);
	cudaEventRecord(start, 0);
	cudaMemcpyAsync(d_a, a, nbytes, cudaMemcpyHostToDevice, 0);
	increment_kernel<<<blocks, threads, 0, 0>>>(d_a, value);
	cudaMemcpyAsync(a, d_a, nbytes, cudaMemcpyDeviceToHost, 0);
	cudaEventRecord(stop, 0);
	sdkStopTimer(&timer);

	// have CPU do some work while waiting for stage 1 to finish
	unsigned long int counter=0;

	while (cudaEventQuery(stop) == cudaErrorNotReady)
	{
		counter++;
	}
	checkCudaErrors(cudaEventElapsedTime(&gpu_time, start, stop));

	// print the cpu and gpu times
	printf("time spent executing by the GPU: %.2f\n", gpu_time);
	printf("time spent by CPU in CUDA calls: %.2f\n", sdkGetTimerValue(&timer));
	printf("CPU executed %lu iterations while waiting for GPU to finish\n", counter);

	// check the output for correctness
	bool bFinalResults = (bool)correct_output(a, n, value);

	// release resources
	checkCudaErrors(cudaEventDestroy(start));
	checkCudaErrors(cudaEventDestroy(stop));
	checkCudaErrors(cudaFreeHost(a));
	checkCudaErrors(cudaFree(d_a));

	cudaDeviceReset();

	exit(bFinalResults ? EXIT_SUCCESS : EXIT_FAILURE);
}

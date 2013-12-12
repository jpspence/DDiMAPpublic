//============================================================================
// Name        : DDiMAP.cpp
// Author      : Androwis Abumoussa
// Version     :
// Copyright   : All Rights Reserved
// Description : DDiMAP in C
//============================================================================
#include "DDiMAPGPU.h"
#include <getopt.h>
#include <time.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>

// Default file.
// char *file = "data/Burack_128F/128F_Gen1_Frag_WithBcl2Sanger_sorted.bam";
string file  = "data/128test_Gen1_example_sorted.bam";

/******************************************************************************
 *									GPU
 ******************************************************************************/

__device__ long long stringToUINT64GPU( char *s) 
{

	long long temp = 0;

	for( int i = 0; i < 17; i++){
		temp+= (s[i] == 'A') ? 1 << (3*i) : 0;
		temp+= (s[i] == 'C') ? 2 << (3*i) : 0;
		temp+= (s[i] == 'G') ? 3 << (3*i) : 0;
		temp+= (s[i] == 'T') ? 4 << (3*i) : 0;
		temp+= (s[i] == '-') ? 7 << (3*i) : 0;
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


/******************************************************************************
 *									CPU
 ******************************************************************************/

int check_output( Read *gpu)
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
		Read bam = buildRead(word);

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
	// Read the BAM file
	// We're going to do a simple map/reduce on this data to prep data for GPU.
	// ------------------------------------------------------------------------
	readFile(file, convert);
	
	/****************************************************************************
	 * GPU Setup 
	 ***************************************************************************/

	// ------------------------------------------------------------------------
	// GPU Initialization
	// ------------------------------------------------------------------------
	int devID;
	cudaDeviceProp deviceProps;

	// This will pick the best possible CUDA capable device
	devID = findCudaDevice(argc, (const char **)argv);

	// get device name
	checkCudaErrors(cudaGetDeviceProperties(&deviceProps, devID));
	printf("I'm using the CUDA device [%s]\n", deviceProps.name);

	// --- Choose which GPU to run on, change this on a multi-GPU system.
	cudaError_t cudaStatus;
	cudaStatus = cudaSetDevice(devID);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	// set kernel launch configuration
	dim3 threads = dim3(512, 1);
	dim3 blocks  = dim3(n / threads.x, 1);

	// ------------------------------------------------------------------------
	// Memory Allocation
	// ------------------------------------------------------------------------

	const long alignmentBytes = n * sizeof(Read);
	const long aBytes = n * sizeof(BamAlignment);

	// allocate CPU memory
	Read *a = 0;
	checkCudaErrors(cudaMallocHost((void **)&a, alignmentBytes));
	memset(a, 0, alignmentBytes);

	BamAlignment *alignments = 0;
	checkCudaErrors(cudaMallocHost((void **)&alignments, aBytes));
	memset(alignments, 0, aBytes);

	// allocate GPU memory
	Read *d_alignments=0;
	checkCudaErrors(cudaMalloc((void **)&d_alignments, alignmentBytes));
	checkCudaErrors(cudaMemset(d_alignments, 0, alignmentBytes));

	// ------------------------------------------------------------------------
	// Setup Streams & Timers
	// ------------------------------------------------------------------------
	float ms;  				// total elapsed time in milliseconds
	float gpu_time = 0.0f;	// total time on GPU

	// --- create cuda timers
	cudaEvent_t start, stop;
	cudaEvent_t time1, time2, time3, time4;

	checkCudaErrors(cudaEventCreate(&start));
	checkCudaErrors(cudaEventCreate(&stop));
	checkCudaErrors(cudaEventCreate(&time1));
	checkCudaErrors(cudaEventCreate(&time2));
	checkCudaErrors(cudaEventCreate(&time3));
	checkCudaErrors(cudaEventCreate(&time4));

	StopWatchInterface *timer = NULL;
	sdkCreateTimer(&timer);
	sdkResetTimer(&timer);

	// --- Create & configure CUDA streams
	const int blockSize = 128, nStreams = 4;
	const int streamSize = n / nStreams;
	const int streamBytes = streamSize * sizeof(Read);
	
	cudaStream_t stream[nStreams];
	for (int i = 0; i < nStreams; ++i)
		checkCudaErrors( cudaStreamCreate(&stream[i]) );

	checkCudaErrors(cudaDeviceSynchronize());


	// ------------------------------------------------------------------------
	// Convert BAM Alignments to binary representations
	// ------------------------------------------------------------------------
	cudaEventRecord(time1, 0);

	cudaEventRecord(time2, 0);

	sdkStartTimer(&timer);
	cudaEventRecord(start, 0);

	cudaMemcpyAsync(d_alignments, a, alignmentBytes, cudaMemcpyHostToDevice);

	cudaEventRecord(time3, 0);

	convert_kernel<<<streamSize/blockSize, blockSize, 0, 0>>>(d_alignments);

	cudaEventRecord(time4, 0);

	cudaMemcpyAsync(a, d_alignments, alignmentBytes, cudaMemcpyDeviceToHost);

	cudaEventRecord(stop, 0);
	sdkStopTimer(&timer);

	// have CPU do some work while waiting for stage 1 to finish
	unsigned long int counter2=0;
	while ( cudaEventQuery(stop) == cudaErrorNotReady)
		counter2++;

	checkCudaErrors(cudaEventElapsedTime(&gpu_time, start, stop));

	// print the cpu and gpu times
	printf("time spent executing by the GPU: %.2f (ms) \n", gpu_time);
	printf("time spent by CPU in CUDA calls: %.2f (ms) \n", sdkGetTimerValue(&timer));
	printf("CPU executed %lu iterations while waiting for GPU to finish\n", counter2);

	// ------------------------------------------------------------------------
	// Check Correctness. 
	// ------------------------------------------------------------------------
	bool bFinalResults = (bool) check_output(a);

	// ------------------------------------------------------------------------
	// Verify the reads
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// Print out the new calls
	// ------------------------------------------------------------------------


	// ------------------------------------------------------------------------
	// End. 
	// ------------------------------------------------------------------------
	// release resources
	for (int i = 0; i < nStreams; ++i)
		checkCudaErrors( cudaStreamDestroy(stream[i]) );

	checkCudaErrors(cudaEventDestroy(start));
	checkCudaErrors(cudaEventDestroy(stop));
	checkCudaErrors(cudaFreeHost(a));
	checkCudaErrors(cudaFree(d_alignments));

	cudaDeviceReset();

	exit(bFinalResults ? EXIT_SUCCESS : EXIT_FAILURE);
}

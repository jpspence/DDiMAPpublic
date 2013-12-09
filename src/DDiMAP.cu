//============================================================================
// Name        : DDiMAP.cpp
// Author      : Androwis Abumoussa
// Version     :
// Copyright   : All Rights Reserved
// Description : DDiMAP in C
//============================================================================

#include "DDiMAP-lib.h"
#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <getopt.h>
#include <sys/_types/_clock_t.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#if defined (_WIN32)
#include <Windows.h>
#endif

using namespace BamTools;
using namespace std;

// Default file.
// char *file = "data/Burack_128F/128F_Gen1_Frag_WithBcl2Sanger_sorted.bam";
string file  = "data/128test_Gen1_example_sorted.bam";

/******************************************************************************
*									GPU
******************************************************************************/
__global__ void read()
{
	// Read in a single read from global memory
	
	if(ba.Position > 0)
	{		
		int name 	  =  ba.RefID;
		int offset    = (ba.IsReverseStrand()) ? ba.AlignedBases.length() - length : 0 ;
		int position  = ba.Position + offset;
		string word   = ba.AlignedBases.substr(offset, length);

		// Increment counter for the observed sequence
		int count = reads[name][position][word];
		reads[name][position][word] = ( count ) ? count + 1 : 1;

	}
}

__global__ void verify()
{
	
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
#if defined (_WIN32)
	double start = GetTickCount();
#else
	clock_t t;
	t = clock();
#endif

	BamReader *br = new BamReader();
	br->Open(file);
	
	BamAlignment ba;
	while(br->GetNextAlignment(ba)){
		read( ba, 34);	
	}

	kernel<<<1,1>>>();
	iterate( print );

#if defined (_WIN32)
	double stop = GetTickCount();
	printf("It took me  %6g s ", (stop - start)/1000 );
#else
	t = clock() - t;
	printf ("It took me %lu ticks (%f seconds).\n",t, ((float)t)/CLOCKS_PER_SEC);
#endif


	// ------------------------------------------------------------------------
	// End. 
	// ------------------------------------------------------------------------
	return EXIT_SUCCESS;

}



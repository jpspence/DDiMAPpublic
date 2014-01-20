//============================================================================
// Name        : DDiMAP.cpp
// Author      : Androwis Abumoussa
// Version     :
// Copyright   : All Rights Reserved
// Description : DDiMAP in C
//============================================================================

#include <getopt.h>
#include <time.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "DDiMAP-lib.h"


// Default file.
// char *file = "data/Burack_128F/128F_Gen1_Frag_WithBcl2Sanger_sorted.bam";
string file  = "data/128test_Gen1_example_sorted.bam";
char  *fasta = "data/128test_Gen1_example.fa";

int main (int argc, char **argv)
{
	// ------------------------------------------------------------------------
	// Parameter Parsing
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

	clock_t t;

	t = clock();
	int total = 0, unique;
	unique = readFile(file, fasta, 34, buildRead);
	t = clock() - t;
	printf ("It took me %lu ticks (%f seconds) to read %d | %d reads from BAM file.\n",
			t, ((float)t)/CLOCKS_PER_SEC, unique, total);


	t = clock();
	int verified = iterate(verify);
	t = clock() - t;

	//	int printed = iterate(print);
	printf ("It took me %lu ticks (%f seconds) to verify %d | %d.\n",
			t, ((float)t)/CLOCKS_PER_SEC , verified, unique);

	// ------------------------------------------------------------------------
	// End. 
	// ------------------------------------------------------------------------
	return EXIT_SUCCESS;

}

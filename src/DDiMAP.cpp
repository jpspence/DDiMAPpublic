//============================================================================
// Name        : DDiMAP.cpp
// Author      : Androwis Abumoussa
// Version     :
// Copyright   : All Rights Reserved
// Description : DDiMAP in C
//============================================================================

#include <stdio.h>
#include <ctype.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include <unordered_map>
#include <time.h>

#include "../include/bamtools/src/api/BamReader.h"
#include "../include/bamtools/src/api/BamWriter.h"

using namespace BamTools;
using namespace std;


// Default file.
char *file = "data/128test_Gen1_example_sorted.bam";



void read(){

	BamReader *br = new BamReader();
	BamAlignment ba;
	unordered_map<int, unordered_map<string, int> > reads;

	//	br->Open("data/Burack_128F/128F_Gen1_Frag_WithBcl2Sanger_sorted.bam");
	br->Open(file);

	long actual_count = 0;
	while(br->GetNextAlignment(ba))
		if(ba.Position > 0)
		{

			string word;
			int position = ba.Position;
			int offset	 = 0;

			if(ba.IsReverseStrand())
				offset = ba.AlignedBases.length()-34;

			position += offset;
			word = ba.AlignedBases.substr(offset, 34);


			// Increment counter for the observed sequence
			int count = reads[position][word];
			reads[position][word] = ( count ) ? count + 1 : 1;
			actual_count++;
		}


	long count 	 = 0;
	long uniques = 0;

	// Let's print out the map
	for (unordered_map<int, unordered_map<string, int> >::iterator It = reads.begin(); It != reads.end(); ++It)
	{
		unordered_map<string, int>::iterator ROA = (*It).second.begin();
		for (; ROA != (*It).second.end(); ++ROA)
		{
			// Print [ROA][COUNT][Seq]
//			if((*ROA).second > 40)
//				cout << " [" << (*It).first << "][" << (*ROA).second << "] "  << (*ROA).first <<endl;

			count += (*ROA).second;
			uniques++;

		}
	}

//	cout << " We read : " << count << endl;
//	cout << " with " << uniques << " | " << actual_count << endl;
//	puts("Boom.");

}


int main (int argc, char **argv) {

	// ------------------------------------------------------------------------
	// Paramaters
	// ------------------------------------------------------------------------
	int c;
	int digit_optind = 0;
	int aopt = 0, bopt = 0;
	char *copt = 0, *dopt = 0;
	static struct option long_options[] = {
			{"file", 	0, 0, 'f'},
			{NULL, 		0, NULL, 0}
	};

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "f:", long_options, &option_index)) != -1) {

		int this_option_optind = optind ? optind : 1;

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
	  int f;
	  t = clock();
	  read();
	  t = clock() - t;
	  printf ("It took me %d ticks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);


	return EXIT_SUCCESS;

}

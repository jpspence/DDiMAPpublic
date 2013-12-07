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
#include <map>
#include <time.h>

#include "../include/bamtools/src/api/BamReader.h"
#include "../include/bamtools/src/api/BamWriter.h"

using namespace BamTools;
using namespace std;


// Default file.
// char *file = "data/Burack_128F/128F_Gen1_Frag_WithBcl2Sanger_sorted.bam";
string file  = "data/128test_Gen1_example_sorted.bam";

void func ( map<int , map<int, map<string, int> > > reads,  void (*f)(int, int, string, int) );
void func ( map<int , map<int, map<string, int> > > reads,  void (*f)(int, int, string, int) )
{

	long count 	 = 0;
	long uniques = 0;

	map<int , map<int, map<string, int> > >::iterator genes = reads.begin();
	for(; genes != reads.end(); ++genes)
	{
		map<int, map<string, int> >::iterator ROAs = (*genes).second.begin();
		for (; ROAs != (*genes).second.end(); ++ROAs)
		{
			map<string, int>::iterator sequences = (*ROAs).second.begin();
			for (; sequences != (*ROAs).second.end(); ++sequences)
			{
				(*f)((*genes).first, (*ROAs).first, (*sequences).first, (*sequences).second);
				count += (*sequences).second;
				uniques++;

			}
		}
	}

	cout << " We read : " << count << endl;
	cout << " with " << uniques  << endl;
	puts("Boom.");

}


void print (int gene, int roa, string seq, int count){

	// Print [gene][ROA][COUNT][Seq]
	if(count > 2 && gene < 10)
		cout << gene << " [" << roa << "][" << count << "] "  << seq << endl;

}


void read(){

	BamReader *br = new BamReader();
	BamAlignment ba;

	// [ gene-name [ roa [ seq count ] ] ]
	map<int , map<int, map<string, int> > > reads;
	br->Open(file);

	while(br->GetNextAlignment(ba))
		if(ba.Position > 0)
		{

			string word;
			int name 	 =  ba.RefID;
			int position = ba.Position;
			int offset	 = 0;

			if(ba.IsReverseStrand())
				offset = ba.AlignedBases.length()-34;

			position += offset;
			word = ba.AlignedBases.substr(offset, 34);


			// Increment counter for the observed sequence
			int count = reads[name][position][word];
			reads[name][position][word] = ( count ) ? count + 1 : 1;
		}

	func(reads, print);
}


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
	clock_t t;
	t = clock();
	read();
	t = clock() - t;
	printf ("It took me %d ticks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);


	// ------------------------------------------------------------------------
	// End. 
	// ------------------------------------------------------------------------
	return EXIT_SUCCESS;

}

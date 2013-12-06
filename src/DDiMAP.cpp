//============================================================================
// Name        : DDiMAP.cpp
// Author      : Androwis Abumoussa
// Version     :
// Copyright   : All Rights Reserved
// Description : Hello World in C, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>
#include "../include/bamtools/src/api/BamReader.h"
#include "../include/bamtools/src/api/BamWriter.h"

using namespace BamTools;
using namespace std;

int main(void) {

	BamReader *br = new BamReader();
	BamAlignment ba;
	unordered_map<int, unordered_map<string, int> > reads;

	//	br->Open("data/Burack_128F/128F_Gen1_Frag_WithBcl2Sanger_sorted.bam");
	br->Open("data/128test_Gen1_example_sorted.bam");


	long actual_count = 0;
	while(br->GetNextAlignment(ba))
		if(ba.Position > 0)
		{
			// Increment counter for the observed sequence
			int count = reads[ba.Position][ba.AlignedBases];
			reads[ba.Position][ba.AlignedBases] = ( count ) ? count + 1 : 1;
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
			if((*ROA).second > 1)
				// Print both the key and the value of that key
				// First = Key (Int), Second = Value of that key (String)
				cout << " [" << (*It).first << "][" << (*ROA).second << "] "  << (*ROA).first <<endl;

			count += (*ROA).second;
			uniques++;

		}
	}

	cout << " We read : " << count << endl;
	cout << " with " << uniques << " | " << actual_count << endl;
	puts("Boom.");
	return EXIT_SUCCESS;

}

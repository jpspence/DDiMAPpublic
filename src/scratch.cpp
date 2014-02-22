
//============================================================================
// Name        : DDiMAP-lib.c
// Author      : Androwis Abumoussa
// Version     :
// Copyright   : All Rights Reserved
// Description : DDiMAP library used by other programs.
//============================================================================

// Common & Custom Junctions are causing a lot of trouble for me.
// Still need to implement frags

// TODO : Compare against NCBI only

#include "DDiMAP-lib.h"

#include <api/BamAux.h>
#include <api/SamHeader.h>
#include <api/SamSequence.h>
#include <api/SamSequenceDictionary.h>
#include <kseq.h>
#include <zlib.h>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <regex>
#include <iomanip>

#define THRESHOLD   2
#define PPM         0.00075
#define PPM         0.01
#define ROA_LENGTH  34

// ----------------------------------------------------------------------------
// Iterator Functions
// ----------------------------------------------------------------------------


//	SNV calls | Only in histograms.

//  Both non-zero --> call 1.
//  One  non-zero --> % coverage
//  ow. check  non-verified histogram > 10%
//	----------------------------------------


//	An SNV is called for any base at a position when:
//
//		A letter that does not match the reference sequence
//
//		- AND -
//
//		a) A letter appears in the verified histograms for both start positions.
//
//		b) A letter appears in the verified histogram for one of the start positions at a frequency in excess of a first set minimum threshold.
//		--> the frequency threshold we typically use is 0.3%.  The computed
//		frequency for each covering collection is the ratio of the verified
//		histogram count to the coverage at that location, namely the sum of
//		all letter counts in that covering collection histogram.
//
//		c) A letter appears in either all-words histogram at a frequency in excess of a second set minimum threshold.
//		--> the frequency threshold we typically use is 10%.  The computed
//		 frequency for each covering collection is the ratio of the all-words
//		 histogram count to the coverage at that location, namely the sum of
//		 all letter counts in that covering collection histogram.
//

//       	-----------------------------------------------------------
// track 0  [0               16|17               33|34              49]
//          -------------------------------------------------------------------
// track 1	         [8               24|25      X        41|42              57]
// 		             ----------------------------------------------------------------------
// track 0 	                   [17               33|34               50|51              66]
// 		                       ---------------------------------------------------------------------
// track 1	    		                [25               41|42               58|59              74]
// 		                                ------------------------------------------------------------

int callSNVs()
{

	int count = iterate(reduceSNVs);

	cout << "StartPos \t EndPos \t #reads \t #refreads \t #NRwords \t Avg#Diffs \t RMSDiffs \t #topNRdiffs \t #topNRreads \t #mutantPos" << endl;

	map<string , map<int, map<uint64_t, map <int , int > > > >::iterator genes = SNVs.begin();
	for(; genes != SNVs.end(); ++genes)
	{
		map<int, map<uint64_t, map <int , int > > >::iterator positions = (*genes).second.begin();
		for (; positions != (*genes).second.end(); ++positions)
		{
			map<uint64_t, map <int , int > >::iterator sequences = (*positions).second.begin();
			for (; sequences != (*positions).second.end(); ++sequences)
			{
				map <int , int >::iterator refIDs = (*sequences).second.begin();
				for (; refIDs != (*sequences).second.end(); ++refIDs)
				{
					if((*refIDs).second > 0 ){
						// #refreads	#NRwords	Avg#Diffs	RMSDiffs	#topNRdiffs	#topNRreads	#mutantPos
						cout << std::setw(20) << genes_names[(*refIDs).first] << " ["; // Gene Name
						cout << std::setw(4)  << (*positions).first << " - "; // StartPos
						cout << std::setw(4)  << (*positions).first + 17	<< "] "; // EndPos
						cout << std::setw(3)  << (*refIDs).second << " "; // Count
						cout << UINT64ToStringCompare((*sequences).first, references[genes_names[(*refIDs).first]][(*positions).first] ) << '\t';
						cout << CalculateGC((*sequences).first) << endl;
					}
				}

			}
		}
	}

	cout << "I called " << count << " SNVs.";

	return 0;
}


// -------------------------

int reduceSNVs(string gene, int position, string seq, Read& read)
{

	count = 0;

	map<string , map<int, map<string, Read> > >::iterator genes = reads.begin();
	for(; genes != reads.end(); ++genes)
	{
		map<int, map<string, Read> >::iterator positions = (*genes).second.begin();
		for (; positions != (*genes).second.end(); ++positions)
		{

			map<string, Read>::iterator  left_sequences = (*positions).second.begin();
			++positions;
			map<string, Read>::iterator       sequences = (*positions).second.begin();
			++positions;
			map<string, Read>::iterator right_sequences = (*positions).second.begin();

			for (; sequences != (*positions).second.end(); ++sequences)
			{
				// Check for SNVs
				Read read = (*sequences).second;
				int  pos  = (*positions).first-17;

				if( not read.matches_reference() and read.is_right_left_verified() )
				{


					// TODO make sure this hasn't been counted
					if( SNVs[gene][pos+17][read.right_sequence_half][read.RefID] )
						SNVs[gene][pos+17][read.right_sequence_half][read.RefID]+= read.count;

					// Begin verifying it against the reads
					else
					{
						for (the left_sequences())
						{
							if(not read.matches_ref_on_left() && not hasDash(read.left_sequence_half) )
							{

								// check the left-front
								if ( not (read.track_verification_flags & 0b1) == 0b1){

								}
								// check the left-rear
								if ( not (read.track_verification_flags & 0b10) == 0b10){

								}
							}

							// check the right-front
							if (not read.matches_ref_on_right() and
									not hasDash(read.right_sequence_half and
											not (read.track_verification_flags & 0b100) == 0b100 )
							{

							}
						}

						for (the right_sequences()){

							// check the left-rear
							if (not read.matches_ref_on_left() and
									not hasDash(read.left_sequence_half) and
									not (read.track_verification_flags & 0b10) == 0b10)
							{

							}
							if (not read.matches_ref_on_right() and not hasDash(read.right_sequence_half)
							{
								// check the right-front
								if ( not (read.track_verification_flags & 0b100) == 0b100){

								}
								// check the right-rear
								if ( not (read.track_verification_flags & 0b1000) == 0b1000){

								}
							}
						}


					}
				}
			}
		}
		--positions;
		--positions;
	}
}

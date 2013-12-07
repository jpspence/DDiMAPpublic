//============================================================================
// Name        : DDiMAP-lib.c
// Author      : Androwis Abumoussa
// Version     :
// Copyright   : All Rights Reserved
// Description : DDiMAP library used by other programs.
//============================================================================

#include "DDiMAP-lib.h"
// [ gene-name [ roa [ seq count ] ] ]
map<int , map<int, map<string, int> > > reads;

void iterate ( void (*f)(int, int, string, int) )
{
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
			}
		}
	}
}


void print (int gene, int roa, string seq, int count)
{
	// Print [gene][ROA][COUNT][Seq]
	if( gene > 1 && count > 3 )
		cout << gene << " [" << roa << "][" << count << "] "  << seq << endl;
}


void read( BamAlignment ba, int length )
{
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


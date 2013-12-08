//============================================================================
// Name        : DDiMAP-lib.c
// Author      : Androwis Abumoussa
// Version     :
// Copyright   : All Rights Reserved
// Description : DDiMAP library used by other programs.
//============================================================================

#include "DDiMAP-lib.h"

// [ gene-name [ roa [ seq count ] ] ]
map<int , map<int, map<string, Read> > > reads;

void iterate ( void (*f)(int, int, string, Read) )
{
	map<int , map<int, map<string, Read> > >::iterator genes = reads.begin();
	for(; genes != reads.end(); ++genes)
	{
		map<int, map<string, Read> >::iterator ROAs = (*genes).second.begin();
		for (; ROAs != (*genes).second.end(); ++ROAs)
		{
			map<string, Read>::iterator sequences = (*ROAs).second.begin();
			for (; sequences != (*ROAs).second.end(); ++sequences)
			{
				(*f)((*genes).first, (*ROAs).first, (*sequences).first, (*sequences).second);
			}
		}
	}
}


void print (int gene, int roa, string seq, Read read)
{
	// Print [gene][ROA][COUNT][Seq]
	if( gene > 1 && read.count > 10 )
		cout << gene << " [" << roa << "][" << read.count << "] "  << seq << endl;
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

		if(reads[name][position][word].count)
			reads[name][position][word].count+=1;
		else
			reads[name][position][word] = buildRead(ba, word);

	}
}

string readToString(Read r){
	string s = " test ";
	return s;
}


// Cipher is given by :
uint64_t a    = 0b00000001;
uint64_t c    = 0b00000010;
uint64_t g    = 0b00000011;
uint64_t t 	  = 0b00000100;
uint64_t dash = 0b00000111;

uint64_t charToUINT64(char ch)
{
	switch(ch){
	case 'A':
	case 'a':
		return a;
	case 'T':
	case 't':
		return t;
	case 'C':
	case 'c':
		return c;
	case 'G':
	case 'g':
		return g;
	case '-':
		return dash;
	default :
		return 0;
	}
}

uint64_t stringToUINT64(string s)
{
	uint64_t temp = 0;
	 for ( int i = 0 ; i < s.length()  ;  i++)
	 {
	         temp += charToUINT64(s[i]) << (3 * i);
	  }

     return temp;

}

Read buildRead(BamAlignment ba, string word)
{
	Read r;
	r.count = 1;

    int half = word.length()/2;
	string  left_half = word.substr(0, half);
	string right_half = word.substr(half , half);

	r.left_sequence_half = stringToUINT64( left_half );
	r.right_sequence_half = stringToUINT64( right_half );

	return r;
}

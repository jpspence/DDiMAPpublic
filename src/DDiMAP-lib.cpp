//============================================================================
// Name        : DDiMAP-lib.c
// Author      : Androwis Abumoussa
// Version     :
// Copyright   : All Rights Reserved
// Description : DDiMAP library used by other programs.
//============================================================================

#include "DDiMAP-lib.h"

#include <api/BamAux.h>
#include <api/SamHeader.h>
//#include <api/SamSequence.h>
#include <api/SamSequenceDictionary.h>
#include <cstdint>
#include <cstring>

// [ gene-name [ roa [ seq count ] ] ]
map<int , map<int, map<string, Read> > > reads;
map<int, string > genes;

// ----------------------------------------------------------------------------
// Convenience methods.
// ----------------------------------------------------------------------------

string createWordString(BamAlignment &ba, int length, int &position)
{
	int offset    = (ba.IsReverseStrand()) ? ba.AlignedBases.length() - length : 0 ;
	position = ba.Position + offset;
	string word   = ba.AlignedBases.substr(offset, length);
	return word;

}

const char *createWordArray(BamAlignment &ba, int length, int &position)
{
	string word = createWordString(ba, length, position);
	return word.c_str();

}

// ----------------------------------------------------------------------------
// BAM --> Reads
// ----------------------------------------------------------------------------

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
		temp += charToUINT64(s[i]) << (3 * i);
	return temp;

}

// Build Read converts sequence to INT
Read buildRead( string &word, int length)
{
	Read r;
	r.count = 1;
	r.verification_flags = 0;
	r.left_sequence_half  = stringToUINT64( word.substr(0, length/2));
	r.right_sequence_half = stringToUINT64( word.substr(length/2 , length/2) );
	return r;
}

// Convert does NOT convert sequence to INT
Read convert(string &word, int length)
{

	Read r;
	r.count = 1;
	r.verification_flags = 0;
	memcpy(r.sequence, word.c_str(), length*sizeof(char));
	return r;

}

// ----------------------------------------------------------------------------
// Reading Files
// ----------------------------------------------------------------------------


int reduce( BamAlignment &ba, int length, Read (*f)(string &, int) )
{

//	cout << " Tag Data : " << ba.TagData << endl;
	if(ba.Position > 0)
	{
		int position;
		string word   = createWordString(ba, length, position);
//		string fullname = genes[ba.RefID];
		int name   = ba.RefID;// fullname.substr(0,fullname.find_first_of("_"));

		// Increment counter for the observed sequence
		if(reads[name][position][word].count)
			reads[name][position][word].count+=1;
		else {
			Read r;
			r = f(word, length);

			//TODO: Find out if this matches the reference sequence
			if(ba.CigarData[0].Length == 50)
				r.verification_flags = r.verification_flags | 0b00000100 ;
			reads[name][position][word] = r;
			return 1;
		}
	}

	return 0;
}

// Returns the number of unique reads in the file.
int readFile(string file, int length, Read (*f)(string &, int))
{
	BamReader *bamreader = new BamReader();
	bamreader->Open(file);

	// Read the header file
	int i =0;
	int size = bamreader->GetHeader().Sequences.Size();
	SamSequenceIterator seqs = bamreader->GetHeader().Sequences.Begin() ;
	for( int j=0; j< size; j++){

		genes[i] = (*seqs).Name;
		cout << i << "th gene = "<< (*seqs).Name << " | " << size << endl;
		i++;seqs++;
	}

	// Begin the alignment search
	BamAlignment ba;
	int counter = 0;
	while(bamreader->GetNextAlignment(ba))
		counter += reduce(ba, length, f);
	bamreader->Close();
	return counter;
}


// ----------------------------------------------------------------------------
// Iterators
// ----------------------------------------------------------------------------

void iterateAndSet( Read reads_array[])
{
	long count = 0;
	map<int , map<int, map<string, Read> > >::iterator genes = reads.begin();
	for(; genes != reads.end(); ++genes)
	{
		map<int, map<string, Read> >::iterator ROAs = (*genes).second.begin();
		for (; ROAs != (*genes).second.end(); ++ROAs)
		{
			map<string, Read>::iterator sequences = (*ROAs).second.begin();
			for (; sequences != (*ROAs).second.end(); ++sequences)
			{
				reads_array[count] = (*sequences).second;
				count++;
			}
		}
	}

}

int iterate ( int (*f)(int, int, string, Read) )
{
	long count = 0;
	map<int , map<int, map<string, Read> > >::iterator genes = reads.begin();
	for(; genes != reads.end(); ++genes)
	{
		map<int, map<string, Read> >::iterator ROAs = (*genes).second.begin();
		for (; ROAs != (*genes).second.end(); ++ROAs)
		{
			map<string, Read>::iterator sequences = (*ROAs).second.begin();
			for (; sequences != (*ROAs).second.end(); ++sequences)
			{
				count += (*f)((*genes).first, (*ROAs).first, (*sequences).first, (*sequences).second);
			}
		}
	}

	return count;
}


int count (int gene, int roa, string seq, Read read)
{ return 1; }


int print (int gene, int roa, string seq, Read read)
{
	// Print [gene][ROA][COUNT][Seq]
	if( read.count > 100 )
		cout << gene << " [" << roa << "][" << read.count << "] "  << seq << endl;

	return 1;
}

int verify ( int gene, int roa, string seq, Read read)
{
	map<string, Read> roaVerifier;

	// Verify the left
	if((roaVerifier = reads[gene][roa - seq.length()/2 ]).size() > 0){
		map<string, Read>::iterator sequences = roaVerifier.begin();
		for (; sequences != roaVerifier.end(); ++sequences)
		{
			if(read.left_sequence_half == (*sequences).second.right_sequence_half){
//				cout << gene << " : " << roa << " :                  " << seq <<endl;
//				cout << gene << " : " << roa - seq.length()/2  << " : " << (*sequences).first << endl;
//				cout << endl;
				read.verification_flags = read.verification_flags | 0b00000001;
				break;
			}
		}
	}

	// Verify the right
	if((roaVerifier = reads[gene][roa + seq.length()/2 ]).size() > 0){
		map<string, Read>::iterator sequences = roaVerifier.begin();
		for (; sequences != roaVerifier.end(); ++sequences)
		{
			if(read.right_sequence_half == (*sequences).second.left_sequence_half ){
//				cout << gene << " : " << roa << " : " << seq <<endl;
//				cout << gene << " : " << roa - seq.length()/2  << " :                  " << (*sequences).first << endl;
//				cout << endl;
				read.verification_flags  = read.verification_flags | 0b00000010;
				break;
			}
		}
	}

	if( (read.verification_flags & 0b00000011) == 0b00000011)
		return 1;

	return 0;

}

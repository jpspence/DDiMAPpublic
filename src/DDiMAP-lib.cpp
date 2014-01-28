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
#include <api/SamSequence.h>
#include <api/SamSequenceDictionary.h>
#include <kseq.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <string>

#define THRESHOLD 2

// Declare the type of file handler for fasta reader
KSEQ_INIT(gzFile, gzread)

// [ gene-name [ roa [ seq count ] ] ]
map<string , map<int, map<string, Read> > > reads;
map<int, string > genes;
map<int, string > genes_names;

// Histogram Data Structures
map<string , map<int, map<int, int> > > verified_histogram;
map<string , map<int, map<int, int> > > threshold_histogram;

// Map <Reference Name -> < position -> sequence>
map<string, map<int, uint64_t > > references;

// TODO: Fix the binary flags for unchanged refs in both constructor & verifications
// TODO: Reduce the printed or verified mutations
// TODO: Print out into Fasta File

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

char UINT64ToChar(uint64_t ch)
{
	if( ch == a)
		return 'A';
	if( ch == t)
		return 'T';
	if (ch == c)
		return 'C';
	if (ch == g)
		return 'G';
	if (ch == dash)
		return '-';
	return '.';
}

string UINT64ToString(uint64_t s)
{

	std::stringstream temp;

	while(s!=0){
		temp << UINT64ToChar( s & 0b00000111 );
		s = s >> 3;
	}
	return temp.str();

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


int counters =0;
int reduce( BamAlignment &ba, int length, Read (*f)(string &, int) )
{

	//	cout << " Tag Data : " << ba.TagData << endl;
	if(ba.Position > 0)
	{
		int position;
		string word   = createWordString(ba, length, position);
		if(position > 2000){
			cout << "wtf position : " << position << endl;
			return 0;
		}
		string name   = genes[ba.RefID];

		// Increment counter for the observed sequence
		if(reads[name][position][word].count)
			reads[name][position][word].count+=1;
		else {
			Read r;
			r = f(word, length);

			// Find out if this matches the reference sequence
			if(ba.CigarData[0].Length == 50){

				bool changed = false;
				r.verification_flags = r.verification_flags | 0b00000100 ;

				// Check to see that the sequence doesn't match a reference
				if(references[genes_names[ba.RefID]][position] == r.left_sequence_half)
					r.verification_flags = r.verification_flags | 0b00001000 ;
				//				else if (references[genes_names[ba.RefID]][position] != 0 ){
				//					cout << "GENE : " << genes_names[ba.RefID] << " @ Position : "<<position << " L "<< endl;
				//					cout << "read : " << UINT64ToString(r.left_sequence_half) 							<< " : " << r.left_sequence_half << endl ;
				//					cout << "refs : " << UINT64ToString(references[genes_names[ba.RefID]][position])	<< " : " << references[genes_names[ba.RefID]][position]  << endl;
				//					cout << counters << "WRONG " <<endl;
				//					changed = true;
				//					cout << endl;
				//				}
				if(references[genes_names[ba.RefID]][position+length/2] == r.right_sequence_half)
					r.verification_flags = r.verification_flags | 0b00010000 ;
				//				else if (references[genes_names[ba.RefID]][position+length/2] != 0 ){
				//					cout << "GENE : " << genes_names[ba.RefID] << " @ Position : "<<position+length/2 << " R "<< endl;
				//					cout << "read : " << UINT64ToString(r.right_sequence_half) 									<< " : " << r.right_sequence_half << endl ;
				//					cout << "refs : " << UINT64ToString(references[genes_names[ba.RefID]][position+length/2])	<< " : " << references[genes_names[ba.RefID]][position+length/2]  << endl;
				//					cout << counters << "WRONG " <<endl;
				//					changed = true;
				//					cout << endl;
				//				}

				if(changed) counters++;
			}

			reads[name][position][word] = r;
			return 1;
		}
	}

	return 0;
}

// Returns the number of unique reads in the file.
int readFile(string file, char *fasta, int length, Read (*f)(string &, int))
{

	// Read the fasta file
	gzFile fp;
	kseq_t *seq;
	int n = 0, slen = 0, qlen = 0;
	FILE *fast = fopen(fasta,"r");
	fp = gzdopen(fileno(fast), "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0){
		++n, slen += seq->seq.l, qlen += seq->qual.l;

		string s = seq->seq.s;
		s.erase( std::remove_if( s.begin(), s.end(), ::isspace ), s.end() );

		map<int, uint64_t> reference;
		for(int j= 0; j< s.length()-length; j++){
			reference[j] = stringToUINT64(s.substr(j, length/2));
		}
		references[seq->name.s] = reference;

	}
	printf("I read %d sequences \t of size %d \t Quality scores %d\n", n, slen, qlen);
	kseq_destroy(seq);
	gzclose(fp);


	// Read the bamfile
	BamReader *bamreader = new BamReader();
	bamreader->Open(file);

	// --- Read the header file and assign the gene ID to the names
	int i =0;
	int size = bamreader->GetConstSamHeader().Sequences.Size();
	SamSequenceIterator seqs = bamreader->GetHeader().Sequences.Begin() ;
	for( int j=0; j< size; j++){
		genes[i] = (*seqs).Name.substr(0,(*seqs).Name.find_first_of("_"));
		genes_names[i] = (*seqs).Name;
		i++;seqs++;
	}

	// --- Begin the alignment search
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
	map<string , map<int, map<string, Read> > >::iterator genes = reads.begin();
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

int iterate ( int (*f)(string, int, string, Read) )
{
	long count = 0;
	map<string , map<int, map<string, Read> > >::iterator genes = reads.begin();
	for(; genes != reads.end(); ++genes)
	{
		map<int, map<string, Read> >::iterator positions = (*genes).second.begin();
		for (; positions != (*genes).second.end(); ++positions)
		{
			map<string, Read>::iterator sequences = (*positions).second.begin();
			for (; sequences != (*positions).second.end(); ++sequences)
			{
				count += (*f)((*genes).first, (*positions).first, (*sequences).first, (*sequences).second);
			}
		}
	}

	return count;
}

// ----------------------------------------------------------------------------
// Iterator Functions
// ----------------------------------------------------------------------------


int count (string gene, int position, string seq, Read read)
{
	if(position > 2000){ return 0;}
	return 1;
}


void printFasta()
{

}
int print (string gene, int position, string seq, Read read)
{
	if(position > 2000){ return 0;}

	// Print [gene][ROA][COUNT][Seq]
	if( read.count > 10000 )
		cout << gene << " [" << position << "][" << read.count << "] "  << seq << endl;

	return 1;
}

int verify ( string gene, int position, string seq, Read read)
{
	if(position > 2000){ return 0;}

	map<string, Read> roaVerifier;

	// Verify the left
	if((roaVerifier = reads[gene][position - seq.length()/2 ]).size() > 0){
		map<string, Read>::iterator sequences = roaVerifier.begin();
		for (; sequences != roaVerifier.end(); ++sequences)
		{
			if(read.left_sequence_half == (*sequences).second.right_sequence_half &&
					(*sequences).second.count > THRESHOLD
			){
				//				cout << gene << " : " << roa << " :                  " << seq <<endl;
				//				cout << gene << " : " << roa - seq.length()/2  << " : " << (*sequences).first << endl;
				//				cout << endl;
				read.verification_flags = read.verification_flags | 0b00000001;
				break;
			}
		}
	}

	// Verify the frequency threshold for that ROA
	if((roaVerifier = reads[gene][position]).size() > 0){
		map<string, Read>::iterator sequences = roaVerifier.begin();
		double count = 0;
		for (; sequences != roaVerifier.end(); ++sequences)
			count += (*sequences).second.count;
		if( (double) read.count / count  > 0.00075)
			read.verification_flags = read.verification_flags | 0b00100000;

	}

	// Verify the right
	if((roaVerifier = reads[gene][position + seq.length()/2 ]).size() > 0){
		map<string, Read>::iterator sequences = roaVerifier.begin();
		for (; sequences != roaVerifier.end(); ++sequences)
		{
			if(read.right_sequence_half == (*sequences).second.left_sequence_half &&
					(*sequences).second.count > THRESHOLD
			){
				//				cout << gene << " : " << roa << " : " << seq <<endl;
				//				cout << gene << " : " << roa - seq.length()/2  << " :                  " << (*sequences).first << endl;
				//				cout << endl;
				read.verification_flags  = read.verification_flags | 0b00000010;
				break;
			}
		}
	}


	//	 Print the verified words that are not reference.
	//	if( (read.verification_flags & 0b00000111) == 0b00000111 && // Left and Right Verified
	//			(read.verification_flags & 0b00011000) != 0b00011000 // Is not an exact match to reference
	//	)
	//		cout << gene << " : " << roa << " : " << seq << " is verified and unique." << endl;

	return ( read.is_right_left_verified() 	&&
			 read.is_above_threshold() 		&&
			 not read.matches_ref_on_left_and_right()
	);
}


/*

	SNV Calling and Coverage Profiling
	----------------------------------

	The final analysis phase produces
	1) SNV calls
	2) Associated SNV Frequencies
	3) Coverage profiles.

	We start with the complete collection of all ROA dictionaries.
	- We accept only the words that pass the bidirectional frequency test
		(at least 2x in each direction &&
		at least 750 ppm in each direction relative to the total coverage of all words in the ROA).

	- We choose two start positions, each of which is used to define a covering dictionary collection.
		For instance, a start position of 9,  ROA size of 34, overlap of 17
		uses the collection of ROAs starting at 9, 26, 43, … as corresponding
		covering collections.  These are the ones used for verification.

	Within each covering collection, we calculate two histograms for each location
	in the underlying reference sequence
	- the number of occurrences of each letter (A, C, G, T, -).
	- Once for words that pass the frequency threshold tests
 	- separately for the verified words.

	TODO: NOTE:
	1) An SNV is called for any base at a position for which one of these rules is followed:

	A letter that does not match the reference sequence
	- AND -

	a) A letter appears in the verified histograms for both start positions.

	b) A letter appears in the verified histogram for one of the start positions at a frequency in excess of a first set minimum threshold.
	--> the frequency threshold we typically use is 0.3%.  The computed
	frequency for each covering collection is the ratio of the verified
	histogram count to the coverage at that location, namely the sum of
	all letter counts in that covering collection histogram.

	c) A letter appears in either all-words histogram at a frequency in excess of a second set minimum threshold.
	--> the frequency threshold we typically use is 10%.  The computed
	 frequency for each covering collection is the ratio of the all-words
	 histogram count to the coverage at that location, namely the sum of
	 all letter counts in that covering collection histogram.


	2) Once an SNV is called, a frequency is computed.
	This combines the results from the two start position results.
	This can be done by adding the histograms for the two start positions together and computing frequencies using the combined counts.  For type 1 calls, this is easy.  For type 2 calls, I have been only including the SNV counts from the covering collection in which it is called in the numerator but have been dividing by the total coverage of the two start position verified histograms.  This is probably not the best way, but the paper was written using data generated this way.  I'd like to put in a flag to compute it this way or by combining counts from both verified histograms in the numerator as well.  For type 3 calls, I have been including SNV counts from both all-words covering collections in the numerator and dividing by the total coverage of the two start position all-words histograms.

	The final SNV calling output should be a csv file for each reference sequence containing the following columns of data (and suggested column headers):

	Gene (reference sequence base name like Bcl2, KS, etc)
	SpecimenID (in a run, there can be multiple specimens that are called separately and we could combine them if we processed saved dictionary collections from them all at once)
	Location (1-based location within the reference sequence)
	RefBase (the base in the reference sequence at that location)
	CalledBase (the called SNV base or D for a deletion.  A hyphen “-” is typically used but spreadsheet programs can turn this into a number 0 as it treats entries with a minus and no letters as a number)
	CallReason (rule number 1, 2, or 3)
	MaxFreq (the larger of the two verified frequencies for rules 1 or 2, the larger of the two all-words frequencies for rule 3 – these are the numbers used for the calling rule)
	Count (the numerator of the frequency calculation)
	Coverage (the denominator of the frequency calculation)
	CalcFreq (the ratio of Count/Coverage)
	RawCoverage (the coverage from the all words histograms, averaged at each position between the two start positions – as an int, chopped not rounded)

	A table of control parameters used to do the analysis should be put as a header in the csv file, with a description in the first column and a value in a second column

	For example:
	ControlParameter,ParameterValue
	ROAsize,34
	ROAoverlap,17
	AcceptanceCountThreshold,2
	AcceptanceFrequencyThreshold,0.000750
	StartPosition#1,1
	StartPosition#2,9
	VerifiedThreshold,0.003
	AllWordsThreshold, 0.1

	This would be followed by a blank line, and then a line with the set of
	column headers for the data and then the data.  Note that if there are
	multiple SNVs at a given location, they each get their own line in the file

 */
int buildHistograms(string gene, int position, string seq, Read read)
{
	if(position > 2000){
		cout << "can't do ROA " << position << endl;
		return 0;
	}

	if( not read.matches_ref_on_left_and_right() &&
		(read.is_right_left_verified() || read.is_above_threshold()))
	{
		// Create a histogram for the left half
		uint64_t s = read.left_sequence_half;
		int i = 0;
		while(s!=0){
			if( read.is_right_left_verified() ) // If the read is verified on both sides
				verified_histogram[gene][position + i][ s & 0b00000111] += read.count;
			if( read.is_above_threshold()) // Threshold
				threshold_histogram[gene][position+i][ s & 0b00000111]  += read.count;
			s = s >> 3; i++;
		}

		// Then the right
		s = read.right_sequence_half;
		while(s!=0){
			if( read.is_right_left_verified() ) // If the read is verified on both sides
				verified_histogram[gene][position + i][ s & 0b00000111] += read.count;
			if( read.is_above_threshold() ) // Threshold
				threshold_histogram[gene][position+i][ s & 0b00000111]  += read.count;
			s = s >> 3; i++;
		}

	}
	return 0;
}

void printHistograms()
{
	map<string , map<int, map<int, int> > >::iterator genes = verified_histogram.begin();
	for(; genes != verified_histogram.end(); ++genes)
	{
		string filename = "/Users/androwis/Desktop/"+(*genes).first+"AT.txt";
		string filenameC = "/Users/androwis/Desktop/"+(*genes).first+"CG.txt";
		ofstream f, fc;

		f.open(filename.c_str());
		fc.open(filenameC.c_str());

		map<int, map<int, int> >::iterator positions = (*genes).second.begin();
		for (; positions != (*genes).second.end(); ++positions)
		{
			map<int, int> counts = (*positions).second;
			double total = counts[1]+counts[2]+counts[3]+counts[4];
			f <<  (*positions).first << "\t" << ((double) (counts[1]+ counts[4]) / total) << "\n";
			fc <<  (*positions).first << "\t" << ((double) (counts[2]+ counts[3]) / total) << "\n";

		}

		f.close();fc.close();
	}
}

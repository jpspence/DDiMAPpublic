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
#include <zlib.h>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <string>

#define THRESHOLD   2
#define PPM         0.00075
#define ROA_LENGTH  34

// Declare the type of file handler for fasta reader
KSEQ_INIT(gzFile, gzread)

// [ gene-name [ roa [ seq count ] ] ]
map<string , map<int, map<string, Read> > > reads;
map<int, string > genes;
map<int, string > genes_names;

// Histogram Data Structures
map<string , map<int, map<int, int> > > verified_histogram_0;
map<string , map<int, map<int, int> > > verified_histogram_1;

map<string , map<int, map<int, int> > > threshold_histogram_0;
map<string , map<int, map<int, int> > > threshold_histogram_1;

// Map <Reference Name -> < position -> sequence>
map<string, map<int, uint64_t > > references;


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
	case 'a': return a;

	case 'T':
	case 't': return t;

	case 'C':
	case 'c':	return c;

	case 'G':
	case 'g':	return g;

	case '-': return dash;

	default : return 0;
	}
}

char UINT64ToChar(uint64_t ch)
{
	if( ch == a) return 'A';
	if( ch == t) return 'T';
	if (ch == c) return 'C';
	if (ch == g) return 'G';
	if (ch == dash) return '-';
	return '\0';
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
// Reading Files : Convenience Functions
// ----------------------------------------------------------------------------

// We assign each read to only 1 Track
string createWordString(BamAlignment &ba, int length, int &position)
{
	int to_roa    = ( length/2 - ba.Position ) % (length/2);
	if(to_roa < 0) to_roa+=length/2;
	int offset    = (ba.IsReverseStrand()) ? ba.AlignedBases.length() - ( (ba.Position + ba.AlignedBases.length())%(length/2)) - length  : to_roa ;

	// Place the read in one of two tracks.
	if(ba.IsReverseStrand() && ba.AlignedBases.length() - offset > 8 )
		offset+=8;
	else if((to_roa % length/2) > 8)
		offset-=8;

	position = ba.Position + offset;

	//cout << (ba.IsReverseStrand() ? "reverse " : "forward ") << "[" << ba.Position << " - " << (ba.Position+ba.AlignedBases.length()) << "] @" << position << endl;
	string word   = ba.AlignedBases.substr(offset, length);
	return word;
}

const char *createWordArray(BamAlignment &ba, int length, int &position)
{
	string word = createWordString(ba, length, position);
	return word.c_str();
}

// ----------------------------------------------------------------------------
// Reading Files
// ----------------------------------------------------------------------------

int reduce( BamAlignment &ba, int length, Read (*f)(string &, int) )
{
	if(ba.Position > 0)
	{
		int position;
		string word   = createWordString(ba, length, position);
		string name   = genes[ba.RefID];

		// Increment counter for the observed sequence
		if(reads[name][position][word].count)
			reads[name][position][word].count+=1;

		// Create a new read for the position on this track
		else {
			Read r;
			r = f(word, length);

			if(ba.CigarData[0].Length == 50)
				r.set_no_indels();

			if(references[genes_names[ba.RefID]][position] == r.left_sequence_half)
				r.set_matches_ref_on_left();

			if(references[genes_names[ba.RefID]][position+length/2] == r.right_sequence_half)
				r.set_matches_ref_on_right();

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

int iterate ( int (*f)(string, int, string, Read&) )
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


int count (string gene, int position, string seq, Read& read)
{
	return 1;
}


void printFasta()
{

}
int print (string gene, int position, string seq, Read& read)
{

	if( read.count > 10000 )
		cout << gene << " [" << position << "][" << read.count << "] "  << seq << endl;

	return 1;
}

int verify ( string gene, int position, string seq, Read& read)
{


	map<string, Read> roaVerifier;

	// Verify the left
	if((roaVerifier = reads[gene][position - seq.length()/2 ]).size() > 0)
	{
		map<string, Read>::iterator sequences = roaVerifier.begin();
		for (; sequences != roaVerifier.end(); ++sequences)
			if(read.left_sequence_half == (*sequences).second.right_sequence_half &&
					(*sequences).second.count > THRESHOLD )
			{
				read.set_left_verified();
				break;
			}
	}

	// Verify the frequency
	double count = 0;
	if((roaVerifier = reads[gene][position]).size() > 0)
	{
		map<string, Read>::iterator sequences = roaVerifier.begin();
		for (; sequences != roaVerifier.end(); ++sequences)
			count+=(*sequences).second.count;
	}
	if( (double)read.count / count > PPM)
		read.set_above_ppm_threshold();


	// Verify the right
	if((roaVerifier = reads[gene][position + seq.length()/2 ]).size() > 0)
	{
		map<string, Read>::iterator sequences = roaVerifier.begin();
		for (; sequences != roaVerifier.end(); ++sequences)
			if(read.right_sequence_half == (*sequences).second.left_sequence_half &&
					(*sequences).second.count > THRESHOLD )
			{
				read.set_right_verified();
				break;
			}
	}

	if( read.is_right_left_verified() &&
			not read.matches_reference() )
		return 1;

	return 0;

}

//	SNV calls
//	----------------------------------
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

int callSNVs( string gene, int position, string seq, Read& read)
{
	unsigned int verified = 0;
	if( not read.matches_reference() and read.is_right_left_verified() ){

		if(not read.matches_ref_on_right()){

			// Check at offset of - 8.
			map<string, Read>::iterator sequence = reads[gene][position+8].begin();
			for(; sequence != reads[gene][position+8].end(); ++sequence)
			{

				if(read.right_half_matches_track_left_offset((*sequence).second.left_sequence_half))
					verified = verified | 0b100;

				if(read.right_half_matches_track_right_offset((*sequence).second.right_sequence_half))
					verified = verified | 0b1000;

				if((verified & 0b1100) == 0b1100)
					break;
			}

			// Check at offset of - 8 - 1/2 ROA.
			if( not ((verified & 0b100) == 0b100)  && position > 7)
			{

				map<string, Read>::iterator sequence = reads[gene][position - 8].begin();
				for(; sequence != reads[gene][position-8].end(); ++sequence)

					if(read.right_half_matches_track_left_offset2((*sequence).second.left_sequence_half))
					{
						verified = verified | 0b100;
						break;
					}
			}

			//TODO: I think there's a one offset error w/ either backward or forward reads.
			if( not ((verified & 0b1000) == 0b1000))
			{
				map<string, Read>::iterator sequence = reads[gene][position + 25].begin();
				for(; sequence != reads[gene][position + 25].end(); ++sequence)
					if(read.right_half_matches_track_right_offset_2((*sequence).second.left_sequence_half))
					{
						verified = verified | 0b1000;
						break;
					}
			}



		}

		if(not read.matches_ref_on_left() && position > 7){

			// Check at offset of - 8.
			map<string, Read>::iterator sequence = reads[gene][position-8].begin();
			for(; sequence != reads[gene][position-8].end(); ++sequence)
			{

				if(read.matches_track_left_offset((*sequence).second.left_sequence_half))
					verified = verified | 0b1;

				if(read.matches_track_right_offset((*sequence).second.right_sequence_half))
					verified = verified | 0b10;

				if((verified & 0b11) == 0b11)
					break;

			}

			// Check at offset of - 8 - 1/2 ROA.
			if( not ((verified & 0b1) == 1)  && position > 24)
			{
				map<string, Read>::iterator sequence = reads[gene][position- 25].begin();
				for(; sequence != reads[gene][position-25].end(); ++sequence)

					if(read.matches_track_left_offset((*sequence).second.right_sequence_half))
					{
						verified = verified | 0b1;
						break;
					}
			}

			// Check at offset of - 8 + 1/2 ROA.
			if( not ((verified & 0b10) == 0b10))
			{

				map<string, Read>::iterator sequence = reads[gene][position + 8].begin();
				for(; sequence != reads[gene][position + 8].end(); ++sequence)
					if(read.matches_track_right_offset_2((*sequence).second.left_sequence_half))
					{
						verified = verified | 0b10;
						break;
					}
			}
		}

		if((verified & 0b11) == 0b11){ }
		if((verified & 0b1100) == 0b1100) {}

	}
	return 0;
}



int buildHistograms(string gene, int position, string seq, Read& read)
{

	map<string, map <int, map<int,int> > > * verified;
	map<string, map <int, map<int,int> > > * frequency;

	if( position % ROA_LENGTH == 0 || position % ROA_LENGTH == ROA_LENGTH/2){
		verified  =  &verified_histogram_0;
		frequency = &threshold_histogram_0;
	} else {
		verified  =  &verified_histogram_1;
		frequency = &threshold_histogram_1;
	}

	if( ( read.is_right_left_verified() || read.is_above_ppm() ) &&
			not read.matches_reference() )
	{

		// Create a histogram for the left half
		uint64_t s = read.left_sequence_half;
		int i = 0;
		while(s!=0)
		{
			if( read.is_right_left_verified() )
				(*verified)[gene][position + i][ s & 0b00000111] += read.count;
			if( read.is_above_ppm() )
				(*frequency)[gene][position + i][ s & 0b00000111] += read.count;
			s = s>>3;i++;
		}

		// Then the right
		s = read.right_sequence_half;
		while(s!=0)
		{
			if( read.is_right_left_verified() )
				(*verified)[gene][position + i][ s & 0b00000111] += read.count;
			if( read.is_above_ppm() )
				(*frequency)[gene][position + i][ s & 0b00000111] += read.count;
			s=s>>3;i++;
		}

	}
	return 0;
}

void printHistograms()
{
	map<string , map<int, map<int, int> > >::iterator genes = verified_histogram_0.begin();
	for(; genes != verified_histogram_0.end(); ++genes)
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
			double total = counts[1]+counts[2]+counts[3]+counts[4] +
					verified_histogram_1[(*genes).first][(*positions).first][1] +
					verified_histogram_1[(*genes).first][(*positions).first][2] +
					verified_histogram_1[(*genes).first][(*positions).first][3] +
					verified_histogram_1[(*genes).first][(*positions).first][4];
			f <<  (*positions).first << "\t" << ((double) (counts[1]+ counts[4]+
					verified_histogram_1[(*genes).first][(*positions).first][1] +
					verified_histogram_1[(*genes).first][(*positions).first][4] ) / total) << "\n";
			fc <<  (*positions).first << "\t" << ((double) (counts[2]+ counts[3] +
					verified_histogram_1[(*genes).first][(*positions).first][2] +
					verified_histogram_1[(*genes).first][(*positions).first][3] ) / total) << "\n";

		}

		f.close();fc.close();
	}
}

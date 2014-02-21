
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
#include <ostream>

#define TEST        0
#define THRESHOLD   2
#define PPM         0.01
#define ROA_LENGTH  34

int max_refid = 0;
int tracks [2] = {0, 8};


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

// SNV Calling Data Structures
//   GENE		Position  SEQ         BA.refID Count
map<string , map<int, map<uint64_t, map <int , int > > > > SNVs;

// Map <Reference Name -> < position -> sequence>
map<string, map<int, uint64_t > > references;



// ----------------------------------------------------------------------------
// BAM --> Reads
// ----------------------------------------------------------------------------

// Convert does NOT convert sequence to INT
Read convert(string &word, int length)
{
	Read r;
	r.forward_count = 2;
	r.reverse_count = 2;
	r.verification_flags = 0;
	memcpy(r.sequence, word.c_str(), length*sizeof(char));
	return r;
}

// Build Read converts sequence to INT
Read buildRead( string &word, int length)
{
	Read r;
	r = convert(word, length);
	r.left_sequence_half  = stringToUINT64( word.substr(0, length/2));
	r.right_sequence_half = stringToUINT64( word.substr(length/2 , length/2) );
	return r;
}

// ----------------------------------------------------------------------------
// Reading Files : Convenience Functions
// ----------------------------------------------------------------------------

string createWordString(BamAlignment &ba, int length, int &position, int track)
{
	int to_roa    = ( length/2 - ba.Position ) % (length/2);

	if(to_roa < 0) to_roa+=length/2;
	int offset    = (ba.IsReverseStrand()) ? ba.AlignedBases.length() - ( (ba.Position + ba.AlignedBases.length())%(length/2)) - length  : to_roa ;


	if(ba.IsReverseStrand() )
		offset+=track;
	else
		offset-=(length/2-track);

	if(offset < 0) offset+=length/2;
	if( offset + length > ba.AlignedBases.length()) offset-=length/2;


	position = ba.Position + offset;
	string word = ba.AlignedBases.substr(offset, length);

	// ensuring this is correct.
	if(TEST){
		// Check that there isn't a more approriate ROA
		if(ba.IsReverseStrand()){
			for(int i = ba.Position+ba.AlignedBases.length(); i > ba.Position; i--)
				if(i%(length/2) == track and i > position+length){
					cout << " Problem with the reverse. Read [" << ba.Position << " - " << (ba.Position + ba.AlignedBases.length() ) << "] " << endl;
					cout << " ROA should end at " << i << " not " << position + length << endl;
					cout << ((ba.IsReverseStrand()) ? " <<<  "  :" >>> ") << " track : " << track << " Starting : " << ba.Position << endl;
					cout << ba.AlignedBases << endl;
					for(int i =0; i < offset; i++)
						cout << " ";
					cout << word << endl;
				}
		}
		else
		{
			for(int i = ba.Position; i < ba.Position + ba.AlignedBases.length(); i++)
				if(i%(length/2) == track and i < position){
					cout << " Problem with the forward. Read [" << ba.Position << " - " << (ba.Position + ba.AlignedBases.length()) << "] " << endl;
					cout << " ROA should start at " << i << " not! " << position << endl;
					cout << ((ba.IsReverseStrand()) ? " <<<  "  :" >>> ") << " track : " << track << " Starting : " << ba.Position << endl;
					cout << ba.AlignedBases << endl;
					for(int i =0; i < offset; i++)
						cout << " ";
					cout << word << endl;
				}
		}
	}
	return word;
}

const char *createWordArray(BamAlignment &ba, int length, int &position, int track)
{
	string word = createWordString(ba, length, position, track);
	return word.c_str();
}

// ----------------------------------------------------------------------------
// DDiMAP
// ----------------------------------------------------------------------------

int cntr = 0;
int reduce( BamAlignment &ba, int length, Read (*f)(string &, int) )
{
	int uniques = 0;

	if( references[genes[ba.RefID]].size() > 0)
	{
		for(int track : tracks){

			if(TEST){
				cntr++;
				cout << "I'm Reduce()ing " << cntr << endl;
			}

			int position;
			string word   = createWordString(ba, length, position, track);

			if(word.size() > 0){

				string name   = genes[ba.RefID];

				// Increment counter for the observed sequence
				if( reads[name][position][word].total_count() )
					if(ba.IsReverseStrand())
						reads[name][position][word].reverse_count+=1;
					else
						reads[name][position][word].forward_count+=1;

				// Create a new read for the position on this track
				else {
					Read r;
					r = f(word, length);
					if(ba.IsReverseStrand())
						++r.reverse_count;
					else
						++r.forward_count;

					r.RefID = ba.RefID;

					if(ba.CigarData[0].Length == 50)
						r.set_no_indels();

					// Check NCBI
					if(references[name][position] == r.left_sequence_half)
						r.set_matches_ref_on_left();

					if(references[name][position+length/2] == r.right_sequence_half)
						r.set_matches_ref_on_right();

					reads[name][position][word] = r;

					uniques++;
				}
			}
		}
	}
	return uniques;
}

// ----------------------------------------------------------------------------
// Reading Files
// ----------------------------------------------------------------------------

// ref_id --> offset
map < int, int> 	frag_offset;

// Returns the number of unique reads in the file.
int readFile(string file, char *fasta, int length, Read (*f)(string &, int))
{

	// Read in the NCBI sequences from Fasta File, assign appropriate offsets
	gzFile fp;
	kseq_t *seq;
	int n = 0;
	FILE *fast = fopen(fasta,"r");
	fp = gzdopen(fileno(fast), "r");

	regex e ("[^a-zA-Z0-9\\-]+");
	regex frag ("[fF][rR][aA][gG][0-9]+_([0-9]*)");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0){

		string seq_name = seq->name.s;
		string clean = std::regex_replace (seq_name,e,"_");

		string s = seq->seq.s;
		s.erase( std::remove_if( s.begin(), s.end(), ::isspace ), s.end() );

		std::smatch m;
		std::regex_search(clean, m, frag);

		if(m.size()>1)
			frag_offset[n] = atoi(m[1].str().c_str());
		else
			frag_offset[n] = 0;

		if(clean.find("NCBI") != -1){
			map<int, uint64_t> reference;
			for(int j= 0; j< s.length()-length; j++){
				reference[j] = stringToUINT64(s.substr(j, length/2));
			}
			references[clean] = reference;
		}

		genes[n] = clean;
		++n;
	}

	printf("I read %d sequences from the fasta file \n",n);
	kseq_destroy(seq);
	gzclose(fp);


	// Read the bamfile
	BamReader *bamreader = new BamReader();
	bamreader->Open(file);

	// --- Begin the alignment search
	BamAlignment ba;
	int counter = 0, total = 0;
	while(bamreader->GetNextAlignment(ba)){
		(&ba)->Position += frag_offset[ba.RefID];
		counter += reduce(ba, length, f);
		if(TEST && references[genes[ba.RefID]].size())
			total++;
	}

	if(TEST)
		cout << "I read a total of " << total << " reads that matched NCBI"<< endl;
	bamreader->Close();

	return counter;
}

// ----------------------------------------------------------------------------
// Iterators
// ----------------------------------------------------------------------------

int iterate ( int (*f)(string, int, string, Read&) )
{
	long count = 0;

	for(auto genes = reads.begin(); genes != reads.end(); ++genes)
		for (auto positions = (*genes).second.begin(); positions != (*genes).second.end(); ++positions)
			for (auto sequences = (*positions).second.begin(); sequences != (*positions).second.end(); ++sequences)
				count += (*f)((*genes).first, (*positions).first, (*sequences).first, (*sequences).second);

	return count;
}

// ----------------------------------------------------------------------------
// Iterator Functions
// ----------------------------------------------------------------------------

int count (string gene, int position, string seq, Read& read)
{ return 1; }

map<string, int> frag_counts;

ofstream fasta_file;

int count_verified (string gene, int position, string seq, Read& read)
{
	int frags = 0;
	if(read.is_right_left_verified() and
	   position % 17 == 0 and
	   not read.matches_reference() and
	   not hasDash(read.left_sequence_half) and
	   not hasDash(read.right_sequence_half)){

		map<string, Read> left_track = reads[gene][position  - ROA_LENGTH/2];
		map<string, Read> right_track = reads[gene][position + ROA_LENGTH/2];

		for (auto left = left_track.begin(); left != left_track.end(); ++left)
			if( (*left).second.is_above_ppm() &&
				not hasDash((*left).second.left_sequence_half) &&
					(*left).second.right_sequence_half == read.left_sequence_half)
				for(auto right = right_track.begin(); right != right_track.end(); ++right)
					if((*right).second.is_above_ppm() &&
					   (*right).second.left_sequence_half == read.right_sequence_half&&
					   not hasDash((*right).second.right_sequence_half))
					{
						if(frag_counts[gene])
							++frag_counts[gene];
						else
							frag_counts[gene] = 1;
						frags++;
						fasta_file << ">" << gene << "_Frag_" << (position-17) << endl;;
						fasta_file << UINT64ToString((*left).second.left_sequence_half);
						fasta_file << UINT64ToString((*left).second.right_sequence_half);
						fasta_file << UINT64ToString((*right).second.left_sequence_half);
						fasta_file << UINT64ToString((*right).second.right_sequence_half);
						fasta_file << endl;

					}

	}
	return frags;
}

int printFasta()
{
	fasta_file.open ("/Users/androwis/Desktop/fasta.fa");
	int frags = iterate(count_verified);
	fasta_file.close();

	return frags;
}


void frequency_filter(string gene, int position)
{
	map<string, Read> roaVerifier;
	if((roaVerifier = reads[gene][position]).size() > 0)
	{

		// Count the number of reads at each location.
		int forward_count = 0;
		int reverse_count = 0;

		for (auto sequences = roaVerifier.begin(); sequences != roaVerifier.end(); ++sequences)
		{
			forward_count+=(*sequences).second.forward_count;
			reverse_count+=(*sequences).second.reverse_count;
		}

		// Apply frequency filter
		int forward_threshold = ((2 > PPM * forward_count) ? 2 : PPM * forward_count);
		int reverse_threshold = ((2 > PPM * reverse_count) ? 2 : PPM * reverse_count);

		for (auto sequences = roaVerifier.begin(); sequences != roaVerifier.end(); ++sequences)
			if( (*sequences).second.forward_count >= forward_threshold &&
				(*sequences).second.reverse_count >= reverse_threshold )
				reads[gene][position][(*sequences).first].set_above_ppm_threshold();

	}
}

void verify( map<string, Read> *left_track, map<string, Read> *right_track)
{
	if((*left_track).size() > 0 && (*right_track).size() > 0)
		for (auto left = (*left_track).begin(); left != (*left_track).end(); ++left)
			if((*left).second.is_above_ppm())
				for(auto right = (*right_track).begin(); right != (*right_track).end(); ++right)
					if((*right).second.is_above_ppm())
						if((*left).second.right_sequence_half == (*right).second.left_sequence_half)
						{
							((*left_track)[(*left).first]).set_right_verified();
							((*right_track)[(*right).first]).set_left_verified();
						}
}

void callSNVs()
{
}

void sequential()
{

	// Frequency Filter.
	for(auto genes = reads.begin(); genes != reads.end(); ++genes)
		for (auto positions = (*genes).second.begin() ; positions != (*genes).second.end(); ++positions)
			frequency_filter((*genes).first, (*positions).first);


	// Verify
	for(auto genes2 = reads.begin(); genes2 != reads.end(); ++genes2)
	{
		auto positions = (*genes2).second.begin();
		positions++; positions++;
		for (; positions != --(--((*genes2).second.end())); ++positions)
		{

			int size = sizeof(tracks)/sizeof(tracks[0]);

			for(int i = 0; i < size; i ++)
				--positions;

			map<string, Read> * left_track = &reads[(*genes2).first][(*positions).first];

			for(int i = 0; i < 2; i ++)
				++positions;

			map<string, Read> * track = &reads[(*genes2).first][(*positions).first];

			verify(left_track, track);

			for(int i = 0; i < size; i ++)
				++positions;

			map<string, Read> * right_track = &reads[(*genes2).first][(*positions).first];

			for(int i = 0; i < size; i ++)
				--positions;

			verify(track, right_track);

		}
	}

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
				(*verified)[gene][position + i][  s & 0b00000111] += read.total_count();
			if( read.is_above_ppm() )
				(*frequency)[gene][position + i][ s & 0b00000111] += read.total_count();
			s = s>>3;i++;
		}

		// Then the right
		s = read.right_sequence_half;
		while(s!=0)
		{
			if( read.is_right_left_verified() )
				(*verified)[gene][position + i][ s & 0b00000111] += read.total_count();
			if( read.is_above_ppm() )
				(*frequency)[gene][position + i][ s & 0b00000111] += read.total_count();
			s=s>>3;i++;
		}

	}
	return 0;
}

void printHistograms()
{
	for(auto frag = frag_counts.begin(); frag !=frag_counts.end(); ++frag)
		cout << "I have " << (*frag).second << " frags for "<< (*frag).first << endl;

	for(auto genes = verified_histogram_0.begin(); genes != verified_histogram_0.end(); ++genes)
	{
		string filename = "/Users/androwis/Desktop/"+(*genes).first+"AT.txt";
		string filenameC = "/Users/androwis/Desktop/"+(*genes).first+"CG.txt";
		ofstream f, fc;

		f.open(filename.c_str());
		fc.open(filenameC.c_str());

		for (auto positions = (*genes).second.begin(); positions != (*genes).second.end(); ++positions)
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

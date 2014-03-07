
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

#define TEST 0

int ROA_LENGTH = 34;

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

map<string , map<int, map<int, int> > > ppm_histogram_0;
map<string , map<int, map<int, int> > > ppm_histogram_1;

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
	r.forward_count = 0;
	r.reverse_count = 0;
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
	string read;
	int i = 0;
	for(auto element = ba.CigarData.begin(); element < ba.CigarData.end(); element++)
	{
		//		if(genes[ba.RefID] == "CD83" && ba.Position > 688 && ba.Position < 730)
		//			cout << (*element).Length << (*element).Type << " being processed on CD83 @ "<< ba.Position << endl;

		if((*element).Type != 'I'){
			read.append( ba.AlignedBases.substr(i, (*element).Length));
		}
		//			else if(genes[ba.RefID] == "CD83" && ba.Position > 688 && ba.Position < 730)
		//				cout << "I'm getting rid of an insertion on " << ba.AlignedBases << " @ "<< ba.Position << endl;
		i+=(*element).Length;

	}

	//	if(genes[ba.RefID] == "CD83" && ba.Position > 688 && ba.Position < 730){
	//		cout << ba.AlignedBases << " Original"<< endl;
	//		cout << read << " processed" << endl;
	//	}

	int to_roa    = ( length/2 - ba.Position ) % (length/2);
	if(to_roa < 0) to_roa+=length/2;

	int offset    = (ba.IsReverseStrand()) ? read.length() - ( (ba.Position + read.length())%(length/2)) - length  : to_roa ;

	if(ba.IsReverseStrand() )
		offset+=track;
	else
		offset-=(length/2-track);

	if(offset < 0) offset+=length/2;
	if( offset + length > read.length()) offset-=length/2;

	position = ba.Position + offset;

	if(offset < 0 || offset + length > read.length()){
		//		cout << endl;
		//		cout << track << " dropping a read ["<< position << " , " << position+read.length() << "]" << endl;
		//		cout << " bc offset " << offset<< " & " << read.length() << " on is Reverse " <<ba.IsReverseStrand() << endl;
		return "";
	}

	string word = read.substr(offset, length);

	//	// ensuring this is correct.
	//	if(TEST){
	//		// Check that there isn't a more approriate ROA
	//		if(ba.IsReverseStrand()){
	//			for(int i = ba.Position+ba.AlignedBases.length(); i > ba.Position; i--)
	//				if(i%(length/2) == track and i > position+length){
	//					cout << " Problem with the reverse. Read [" << ba.Position << " - " << (ba.Position + ba.AlignedBases.length() ) << "] " << endl;
	//					cout << " ROA should end at " << i << " not " << position + length << endl;
	//					cout << ((ba.IsReverseStrand()) ? " <<<  "  :" >>> ") << " track : " << track << " Starting : " << ba.Position << endl;
	//					cout << ba.AlignedBases << endl;
	//					for(int i =0; i < offset; i++)
	//						cout << " ";
	//					cout << word << endl;
	//				}
	//		}
	//		else
	//		{
	//			for(int i = ba.Position; i < ba.Position + ba.AlignedBases.length(); i++)
	//				if(i%(length/2) == track and i < position){
	//					cout << " Problem with the forward. Read [" << ba.Position << " - " << (ba.Position + ba.AlignedBases.length()) << "] " << endl;
	//					cout << " ROA should start at " << i << " not! " << position << endl;
	//					cout << ((ba.IsReverseStrand()) ? " <<<  "  :" >>> ") << " track : " << track << " Starting : " << ba.Position << endl;
	//					cout << ba.AlignedBases << endl;
	//					for(int i =0; i < offset; i++)
	//						cout << " ";
	//					cout << word << endl;
	//				}
	//		}
	//	}
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

	bool hasDeletion  = false;
	bool hasInsertion = false;

	for(auto element = ba.CigarData.begin(); element < ba.CigarData.end(); element++){
		if((*element).Type == 'D') hasDeletion = true;
		else if((*element).Type == 'I')hasInsertion = true;

	}

	if( references[genes[ba.RefID]].size() > 0 && not (hasInsertion && hasDeletion) )
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
			frag_offset[n] = (atoi(m[1].str().c_str()) - 1);
		else
			frag_offset[n] = 0;

		if(clean.find("NCBI") != -1){
			clean = clean.substr(0, clean.find_first_of("_"));

			map<int, uint64_t> reference;
			for(int j= 0; j< s.length()-length; j++){
				reference[j] = stringToUINT64(s.substr(j, length/2));

				// Add reference to ROAs
				if(j % 17 == 0  || j % 17 == 8){
					Read r;
					string seq = s.substr(j, length);
					r = buildRead( seq , length);
					r.set_matches_ref_on_right();
					r.set_matches_ref_on_left();
					r.RefID = n;
					reads[clean][j][r.sequence] = r;
				}

			}
			references[clean] = reference;
		}
		else{
			clean = clean.substr(0, clean.find_first_of("_"));
		}
		genes[n] = clean;
		++n;
	}

	//	printf("I read %d sequences from the fasta file \n",n);
	kseq_destroy(seq);
	gzclose(fp);


	// Read the bamfile
	BamReader *bamreader = new BamReader();
	bamreader->Open(file);

	// --- Begin the alignment search
	BamAlignment ba;
	int counter = 0, total = 0;
	while(bamreader->GetNextAlignment(ba)){
		(&ba)->Position += (frag_offset[ba.RefID]);
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
{
	if( false && position > 990 && read.is_right_left_verified()){
		cout << gene << "  " <<  position <<" ";
		cout << UINT64ToStringCompare(read.left_sequence_half, references[gene][position]) ;
		cout << UINT64ToStringCompare(read.right_sequence_half, references[gene][position+17]) ;
		cout <<  " " << read.forward_count << " | "<< read.reverse_count << " PPM : "<< read.is_above_ppm() << " ver-R: " <<read.is_right_left_verified() << " " << endl;
	}
	return 1;
}

map<string, int> frag_counts;

ofstream fasta_file;

//OK, so here's the 10% frag making rule one more time:
//
//If over 10%  && not right/left verified at 1%,
//then extend in any unverified direction using ref sequence.
//
//If a core word is over 1% and fully verified,
// then
//    if there are any unverified verifiers in either direction,
//       replace the unverified half of the verifier with ref sequence and treat the synthetic sequence as if it is verified
//    extend using all combinations of verified verifiers in each direction

// At this point, I know that
// 1) not read.matches_reference() and
// 2) read.is_above_non_verified() OR (read.is_right_left_verified() and read.is_above_frag())
// 3) not hasDash(read.left_sequence_half)
// 4) not hasDash(read.right_sequence_half)

map<string, Read> returnMatches(string gene, int position, int offset, Read& read)
{

	map<string, Read> matches;
	map<string, Read> roa = reads[gene][position  + offset];
	bool add_ref = false;

	if( offset < 0 and read.is_left_verified() )
	{
		for (auto seq = roa.begin(); seq != roa.end(); ++seq)
			if( (*seq).second.right_sequence_half == read.left_sequence_half
					and not hasDash((*seq).second.left_sequence_half))
				if((*seq).second.is_above_frag() || (*seq).second.matches_reference())
					if((*seq).second.is_right_left_verified())
						matches[(*seq).first] = (*seq).second;
					else if (read.is_above_non_verified())
						add_ref = true;


	} else if(offset > 0 and read.is_right_verified() )
		for (auto seq = roa.begin(); seq != roa.end(); ++seq)
			if((*seq).second.left_sequence_half == read.right_sequence_half
					and not hasDash((*seq).second.right_sequence_half))
				if((*seq).second.is_above_frag() || (*seq).second.matches_reference())
					if( (*seq).second.is_right_left_verified())
						matches[(*seq).first] = (*seq).second;
					else if (read.is_above_non_verified())
						add_ref = true;

	// Append reference if no matches.
	if((matches.size() == 0 and read.is_above_non_verified()) || add_ref)
	{
		Read ref;
		ref.left_sequence_half = read.right_sequence_half;
		ref.right_sequence_half = read.left_sequence_half;
		uint64_t reference = references[gene][position + ( (offset > 0) ? 2 : 1) * offset];
		if(reference){
			if(offset < 0) // add the left reference
				ref.left_sequence_half = reference;
			else 		   // add the right reference
				ref.right_sequence_half = reference;
			matches[UINT64ToString(ref.left_sequence_half)+UINT64ToString(ref.right_sequence_half)] = ref;
		}
	}

	return matches;
}

int generateFrags (string gene, int position, string seq, Read& read)
{
	int frags = 0;
	if(     not read.matches_reference() and
			(  read.is_above_non_verified() ||
					(read.is_right_left_verified() and read.is_above_frag())
			) and
			not hasDash(read.left_sequence_half) and
			not hasDash(read.right_sequence_half))
	{

		map<string, Read> left_track = returnMatches(gene, position, (-1*ROA_LENGTH/2), read);
		map<string, Read> right_track = returnMatches(gene, position,    ROA_LENGTH/2,  read);

		if(left_track.size() && right_track.size())
			for (auto left = left_track.begin(); left != left_track.end(); ++left)
				for(auto right = right_track.begin(); right != right_track.end(); ++right)
				{
					if(frag_counts[gene])
						++frag_counts[gene];
					else
						frag_counts[gene] = 1;
					frags++;

					fasta_file << ">" << gene << "_Frag_" << (position-17) << "@" <<read.total_count()<< endl;;
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
	int frags = iterate(generateFrags);
	fasta_file.close();

	return frags;
}

void frequency_filter(string gene, int position, int threshold, double ppm, double frag, double non_verified)
{
	map<string, Read> roaVerifier;
	if((roaVerifier = reads[gene][position]).size() > 0)
	{

		// Count the number of reads at each location.
		int total = 0;
		for (auto sequences = roaVerifier.begin(); sequences != roaVerifier.end(); ++sequences)
			total+=(*sequences).second.total_count();

		double value = 0;
		double  ppm_threshold = (( threshold > (value = ppm * (double)total)) ? threshold : value);
		double frag_threshold = (( threshold > (value = frag * (double)total)) ? threshold : value);

		for (auto sequences = roaVerifier.begin(); sequences != roaVerifier.end(); ++sequences){
			Read read = (*sequences).second;
			if( ( read.forward_count >= ppm_threshold && read.reverse_count >= ppm_threshold ) || read.matches_reference())
				reads[gene][position][(*sequences).first].set_above_ppm_threshold();
			else
				total -= read.total_count();
		}
		double   nv_threshold = (( threshold > (value = non_verified * (double)total)) ? threshold : value);

		// apply frequency filters.
		for (auto sequences = roaVerifier.begin(); sequences != roaVerifier.end(); ++sequences)
		{
			Read read = (*sequences).second;
			if( (read.forward_count >= frag_threshold && read.reverse_count >= frag_threshold) || read.matches_reference() )
				if(read.forward_count + read.reverse_count >= nv_threshold )
					reads[gene][position][(*sequences).first].set_above_nv_threshold();
				else
					reads[gene][position][(*sequences).first].set_above_frag_threshold();
		}

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



void sequential(int threshold, double ppm, double frag, double non_verified)
{

	// Frequency Filter.
	for(auto genes = reads.begin(); genes != reads.end(); ++genes)
		for (auto positions = (*genes).second.begin() ; positions != (*genes).second.end(); ++positions)
			frequency_filter((*genes).first, (*positions).first, threshold, ppm, frag, non_verified);


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
	map<string, map <int, map<int,int> > > * ppm;

	if( position % ROA_LENGTH/2 == 0 ){
		verified  =  &verified_histogram_0;
		ppm = &ppm_histogram_0;
	} else {
		verified  =  &verified_histogram_1;
		ppm = &ppm_histogram_1;
	}

	if( read.is_right_left_verified() || read.is_above_ppm() )
	{

		// Create a histogram for the left half
		uint64_t s = read.left_sequence_half;
		int i = 0;
		while(s!=0)
		{
			if( read.is_right_left_verified() )
				(*verified)[gene][position + i][  s & 0b00000111] += read.total_count();
			if( read.is_above_ppm() )
				(*ppm)[gene][position + i][ s & 0b00000111] += read.total_count();
			s = s>>3;i++;
		}

		// Then the right
		s = read.right_sequence_half;
		while(s!=0)
		{
			if( read.is_right_left_verified() )
				(*verified)[gene][position + i][ s & 0b00000111] += read.total_count();
			if( read.is_above_ppm() )
				(*ppm)[gene][position + i][ s & 0b00000111] += read.total_count();
			s=s>>3;i++;
		}

	}
	return 0;
}

int printSNV(int reason, string gene, int pos, int i, uint64_t ref, double freq )
{
	cout << reason <<" Reference " << gene << " @ " << pos << " is : " << UINT64ToString(ref)<< endl;
	uint64_t convert (i);
	cout << UINT64ToString(ref & 0b111) << " -> " << UINT64ToString(convert) << " : ";
	cout << freq << "\n";
	return 1;
}

void callSNVs(double snv_verified_threshold, double snv_total_threshold)
{
	int snvs = 0;

	for(auto genes = verified_histogram_0.begin(); genes != verified_histogram_0.end(); ++genes)

		for (auto positions = (*genes).second.begin(); positions != (*genes).second.end(); ++positions)
		{

			map<int, int> verified_counts = (*positions).second;
			map<int, int> verified_counts2 = verified_histogram_1[(*genes).first][(*positions).first];

			double verified_total  	= 0;
			double verified_total2 	= 0;
			double ppm_total 		= 0;
			double ppm_total2		= 0;

			for(int i = 1; i < 5; i++)
			{
				verified_total 	+= verified_counts[i];
				verified_total2 += verified_counts2[i];
				ppm_total  		+= ppm_histogram_0[(*genes).first][(*positions).first][i];
				ppm_total2  	+= ppm_histogram_1[(*genes).first][(*positions).first][i];
			}

			double verified = verified_total + verified_total2;
			double ppm = ppm_total + ppm_total2;

			double freq;
			for (int i = 1; i < 5; i++)
				if((references[(*genes).first][(*positions).first] & 0b111) != i)
				{
					uint64_t ref = references[(*genes).first][(*positions).first];

					// --- Call type #1
					// If the reads are in both verified histograms.
					if( (freq = ((double) (verified_counts[i]+ verified_counts2[i]) / verified)) and verified_counts[i] and verified_counts2[i])
						snvs += printSNV(1, (*genes).first,(*positions).first, i, ref, freq );

					// --- Call type #2
					// If the reads are only in one histogram
					else if( (freq = ((double) verified_counts[i]) / verified_total) > snv_verified_threshold )
						snvs += printSNV(2, (*genes).first,(*positions).first, i, ref, freq );

					else if( (freq = ((double) verified_counts2[i]) / verified_total2) > snv_verified_threshold )
						snvs += printSNV(2, (*genes).first,(*positions).first, i, ref, freq );

					// # Call type 3
					// If the reads exceed a 3rd threshold in either track
					else if( (freq = ((double) ppm_histogram_0[(*genes).first][(*positions).first][i]) / ppm_total) > snv_total_threshold)
						snvs += printSNV(3, (*genes).first,(*positions).first, i, ref, freq );

					else if( (freq = ((double) ppm_histogram_1[(*genes).first][(*positions).first][i]) / ppm_total2) > snv_total_threshold)
						snvs += printSNV(3, (*genes).first,(*positions).first, i, ref, freq );

				}

		}
	cout << " I read " << snvs << " SNVs";
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



void check(int total, double ppm, double frag,  int position,int i, string name, string seq, Read read, int above)
{
	//	 Check the left half.
	if( position == i &&
			read.left_sequence_half == stringToUINT64(seq.substr(0,17)) &&
			read.right_sequence_half == stringToUINT64(seq.substr(17,17))){
		cout << read.sequence;
		for (int j = 0; j<34; j++)
			cout << " ";
		cout << " (f " <<read.forward_count <<"-" <<read.matches_ref_on_left()<< " + " <<read.reverse_count << "-"<<read.matches_ref_on_right()<<" / " << ((above > 0) ?  total : ((above == 0) ? frag :  ppm))  << ")" << (( above > 9 )? " > 10% ": (( above > 0) ? " > 1%" : (( above == 0) ? " > ppm" :  " :  NOT PPM"))) << endl ;
	}

	//	 Check the middle
	if( position == i+17
			&& read.left_sequence_half == stringToUINT64(seq.substr(17,17))
			&& read.right_sequence_half == stringToUINT64(seq.substr(34,17))
	){
		for (int j = 0; j<17; j++)
			cout << " ";
		cout << read.sequence;
		for (int j = 0; j<17; j++)
			cout << " ";
		cout << " (f " <<read.forward_count <<"-" <<read.matches_ref_on_left()<< " + " <<read.reverse_count << "-"<<read.matches_ref_on_right()<<" / " << total  << ")" << (( above > 9 )? " > 10% ": (( above > 0) ? " > 1%" : (( above == 0) ? " > ppm" : "  NOT PPM"))) << endl ;
	}

	// Check the right half
	if( position == i+34
			&& read.left_sequence_half == stringToUINT64(seq.substr(34,17))
			&& read.right_sequence_half == stringToUINT64(seq.substr(51,17))
	){
		for (int j = 0; j<34; j++)
			cout << " ";
		cout << read.sequence;
		cout << " (f " <<read.forward_count <<"-" <<read.matches_ref_on_left()<< " + " <<read.reverse_count << "-"<<read.matches_ref_on_right()<<" / " << ((above > 0) ?  total : ((above == 0) ? frag :  ppm))  << ")" << (( above > 9 )? " > 10% ": (( above > 0) ? " > 1%" : (( above == 0) ? " > ppm" :  " :  NOT PPM"))) << endl ;
	}

}

void check_frequency_filter(string gene, int position, string name, string sequence)
{
	int VERIFY_THRESHOLD  =  2;
	double PPM = 0.00075;
	double FRAG_THRESHOLD = .01;
	double NON_VERIFIED_THRESHOLD = .1;

	map<string, Read> roaVerifier;
	for(int i = 0; i < 3; i++)
	if((roaVerifier = reads[gene][position + i * 17]).size() > 0)
	{
		// Count the number of reads at each location.
		int total = 0;
		for (auto sequences = roaVerifier.begin(); sequences != roaVerifier.end(); ++sequences)
			total+=(*sequences).second.total_count();

		double value = 0;
		double ppm_threshold = (( VERIFY_THRESHOLD > (value = PPM * (double)total)) ? VERIFY_THRESHOLD : value);
		double frag_threshold = (( VERIFY_THRESHOLD > (value = FRAG_THRESHOLD * (double)total)) ? VERIFY_THRESHOLD : value);

		for (auto sequences = roaVerifier.begin(); sequences != roaVerifier.end(); ++sequences){
			Read read = (*sequences).second;
			if( ( read.forward_count >= ppm_threshold && read.reverse_count >= ppm_threshold ) || read.matches_reference()){ }
			else
				total -= read.total_count();
		}

		double   nv_threshold = (( VERIFY_THRESHOLD > (value = NON_VERIFIED_THRESHOLD * (double)total)) ? VERIFY_THRESHOLD : value);

		// apply frequency filters.
		for (auto sequences = roaVerifier.begin(); sequences != roaVerifier.end(); ++sequences)
		{
			Read read = (*sequences).second;
			if( (read.forward_count >= frag_threshold && read.reverse_count >= frag_threshold) || read.matches_reference() )
				if(read.forward_count + read.reverse_count >= nv_threshold )
					check(total,ppm_threshold, frag_threshold,  position+i*17,  position , name , sequence, read, 10);
				else
					check(total,ppm_threshold, frag_threshold, position+i*17,  position , name , sequence, read, 1);

			else if( ( read.forward_count >= ppm_threshold && read.reverse_count >= ppm_threshold ) || read.matches_reference()){
				check(total, ppm_threshold, frag_threshold,  position+i*17,  position , name , sequence, read, 0);
			} else {
				check(total,ppm_threshold, frag_threshold,  position+i*17,  position , name , sequence, read, -1);
			}
		}

	}
}



void test()
{

	map <string, string> johns, johns_alph;
	map <string, string> mine, mine_alph;

	char  *fasta = "/Users/androwis/Downloads/TestCaseFixedHigh.fa";
	char  *fasta2 = "/Users/androwis/Desktop/fasta.fa";


	// Read in the NCBI sequences from Fasta File, assign appropriate offsets
	gzFile fp, fp2;
	kseq_t *seq, *seq2;
	int n = 0;
	FILE *fast = fopen(fasta,"r");
	fp = gzdopen(fileno(fast), "r");

	regex e ("[^a-zA-Z0-9\\-]+");
	regex frag ("[fF][rR][aA][gG][0-9]+_([0-9]*)");

	// Read Johns
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0){

		string seq_name = seq->name.s;
		string clean = std::regex_replace (seq_name,e,"_");

		string s = seq->seq.s;
		s.erase( std::remove_if( s.begin(), s.end(), ::isspace ), s.end() );

		std::smatch m;
		std::regex_search(clean, m, frag);

		if(m.size()>1)
			frag_offset[n] = (atoi(m[1].str().c_str()) - 1);
		else
			frag_offset[n] = 0;

		if(clean.find("Frag") && clean.find("NCBI") == -1 && clean.find("Junction") == -1)
			johns[s] = clean;
	}

	FILE *fast2 = fopen(fasta2,"r");
	fp2 = gzdopen(fileno(fast2), "r");
	seq2 = kseq_init(fp2);
	while (kseq_read(seq2) >= 0){

		string seq_name = seq2->name.s;
		string clean = std::regex_replace (seq_name,e,"_");

		string s = seq2->seq.s;
		s.erase( std::remove_if( s.begin(), s.end(), ::isspace ), s.end() );

		std::smatch m;
		std::regex_search(clean, m, frag);

		if(m.size()>1)
			frag_offset[n] = (atoi(m[1].str().c_str()) - 1);
		else
			frag_offset[n] = 0;

		if(clean.find("Frag") && clean.find("NCBI") == -1 && clean.find("Junction") == -1)
			mine[s] = clean;
	}

	int a =0,b=0;

	for(auto my = mine.begin(); my != mine.end(); ++my)
		if( johns[(*my).first].length() == 0)
			mine_alph[(*my).second] = (*my).first;

	for(auto my = mine_alph.begin(); my!=mine_alph.end(); ++my){
		int pos = (*my).first.find_first_of("_");
		string gene = (*my).first.substr(0,pos);
		int position = stoi((*my).first.substr(pos+6, (*my).first.find_last_of("_") - pos - 5));

		// Check
		cout <<  endl << (*my).second<< " "<< (*my).first << " <-- me " << gene << " : "<< (a++)<< endl;
		check_frequency_filter(gene, position,(*my).first , (*my).second);
	}

	cout << "-------- Johns --------- "<< endl;

	for(auto his = johns.begin(); his != johns.end(); ++his)
		if( mine[(*his).first].length() == 0 )
			johns_alph[(*his).second] = (*his).first;

	for(auto his = johns_alph.begin(); his!=johns_alph.end(); ++his){
		int pos = (*his).first.find_first_of("_");
			string gene = (*his).first.substr(0,pos);
			int position = stoi((*his).first.substr(pos+6, (*his).first.find_last_of("_") - pos - 5));

			// Check
			cout << endl<< (*his).second << " " << (*his).first << " <-- John " << (b++)<< endl;
			check_frequency_filter(gene, position - 1, (*his).first , (*his).second);
	}

	cout << " I had " << a <<" / " << mine.size() << " unique | John had "<< b <<" / " << johns.size()<< " unique"<< endl;



}

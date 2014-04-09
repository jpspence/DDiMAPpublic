
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
#include "DDiMAP-test.h"
#include <cctype>
#define TEST 0

int ROA_LENGTH;
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
		if((*element).Type != 'H' && (*element).Type != 'S'){
			if((*element).Type != 'I')
				read.append( ba.AlignedBases.substr(i, (*element).Length));
			i+=(*element).Length;
		}
	}


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

	if(offset < 0 || offset + length > read.length())
		return "";


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
int indel = 0;
int noref = 0;

int reduce( BamAlignment &ba, int length, bool dropID, Read (*f)(string &, int) )
{
	int uniques = 0;

	bool hasDeletion  = false;
	bool hasInsertion = false;
	for(auto element = ba.CigarData.begin(); element < ba.CigarData.end(); element++)
	{
		if((*element).Type == 'D') hasDeletion = true;
		else if((*element).Type == 'I')hasInsertion = true;
	}

	if( references[genes[ba.RefID]].size() > 0 && not (hasInsertion && hasDeletion && dropID) )
	{
		if(TEST){
			cntr++;
			cout << "I'm Reduce()ing : " << " | " << cntr  << endl;
		}
		for(int track : tracks){

			int position;
			string word   = createWordString(ba, length, position, track);

			if(word.size() > 0){

				string name   = genes[ba.RefID];

				// Increment counter for the observed sequence
				if( reads[name][position][word].total_count() )
					reads[name][position][word].set_indels(hasDeletion, hasInsertion, ba.IsReverseStrand());

				// Create a new read for the position on this track
				else {
					Read r;
					r = f(word, length);

					// TODO : handle this r.RefID = ba.RefID;

					r.set_indels(hasDeletion,hasInsertion, ba.IsReverseStrand());

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
int readFile(string file, string fasta, int roa_length, bool dropID, Read (*f)(string &, int))
{
	// Set the global variable
	ROA_LENGTH = roa_length;

	// Read in the NCBI sequences from Fasta File, assign appropriate offsets
	gzFile fp;
	kseq_t *seq;
	int n = 0;
	FILE *fast = fopen(fasta.c_str(),"r");
	fp = gzdopen(fileno(fast), "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0){

		string seq_name = seq->name.s;
		string s = seq->seq.s;

		for (size_t i = 0; i < s.length(); ++i)
			if (s[i]!='a' && s[i]!='A' &&
					s[i]!='c' && s[i]!='C' &&
					s[i]!='t' && s[i]!='T' &&
					s[i]!='g' && s[i]!='G')
				s.erase(i, 1);

		if(seq_name.find("Frag")!=-1)
		{
			//			cout << "Sequence name : "<< seq_name << endl;
			int loc = seq_name.find("Frag");
			string frag = seq_name.substr(loc, seq_name.length()-loc);
			loc = frag.find_first_of("_")+1;
			string locations = frag.substr(loc, frag.length()-loc);
			//			cout << "Number : " << locations.substr(0,locations.find_first_of("_")) << endl;
			frag_offset[n] = atoi(locations.substr(0,locations.find_first_of("_")).c_str()) - 1;
		}
		else
			frag_offset[n] = 0;

		int loc = (seq_name.find_first_of("_") == -1) ? seq_name.length() : seq_name.find_first_of("_");

		// If this includes NCBI
		if(seq_name.find("Frag") == -1 && seq_name.find("Junction") == -1){
			seq_name = seq_name.substr(0, loc);
			map<int, uint64_t> reference;
			for(int j= 0; j< s.length()-ROA_LENGTH; j++){

				reference[j] = stringToUINT64(s.substr(j, ROA_LENGTH/2));
				// Add reference to ROAs
				if(j % (ROA_LENGTH/2) == 0  || j % (ROA_LENGTH/2) == 8){
					Read r;
					string seq = s.substr(j, ROA_LENGTH);
					r = buildRead( seq , ROA_LENGTH);
					r.set_matches_ref_on_right();
					r.set_matches_ref_on_left();
					r.RefID = n;
					reads[seq_name][j][r.sequence] = r;
				}

			}
			references[seq_name] = reference;
		}
		else{
			seq_name = seq_name.substr(0, loc);
		}

		genes[n] = seq_name;
		++n;
	}

	kseq_destroy(seq);
	gzclose(fp);
	if(TEST)
		cout << "I read a total of " << n << " reads from the fasta file"<< endl;



	// Read the bamfile
	BamReader *bamreader = new BamReader();
	bamreader->Open(file);

	// --- Begin the alignment search
	BamAlignment ba;
	int counter = 0, total = 0;
	while(bamreader->GetNextAlignment(ba)){
		(&ba)->Position += (frag_offset[ba.RefID]);
		counter += reduce(ba, ROA_LENGTH, dropID, f);
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
	if( TEST && position > 990 && read.is_right_left_verified()){
		cout << gene << "  " <<  position <<" ";
		cout << UINT64ToStringCompare(read.left_sequence_half, references[gene][position]) ;
		cout << UINT64ToStringCompare(read.right_sequence_half, references[gene][position+(ROA_LENGTH/2)]) ;
		cout <<  " " << read.forward_count << " | "<< read.reverse_count << " PPM : "<< read.is_above_ppm() << " ver-R: " <<read.is_right_left_verified() << " " << endl;
	}
	return 1;
}

map<string, int> frag_counts;

ofstream fasta_file;

map<string, Read> returnMatches(string gene, int position, int offset, Read& read)
														{

	map<string, Read> matches;
	map<string, Read> roa = reads[gene][position  + offset];
	bool add_ref = false;

	if( (offset < 0 and read.is_left_verified())
			|| (offset > 0 and read.is_right_verified()) )
		for (auto seq = roa.begin(); seq != roa.end(); ++seq)
			if( (offset < 0
					and (*seq).second.right_sequence_half == read.left_sequence_half
					and not hasDash((*seq).second.left_sequence_half))
					|| (offset > 0
							and (*seq).second.left_sequence_half == read.right_sequence_half
							and not hasDash((*seq).second.right_sequence_half) )
			)
				if((*seq).second.is_above_frag() || (*seq).second.matches_reference())
					if((*seq).second.is_right_left_verified_at_frag())
						matches[(*seq).first] = (*seq).second;


	// Append reference if no matches.
	if(matches.size() == 0)
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
					(read.is_right_left_verified_at_frag() and read.is_above_frag())
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

					fasta_file << ">" << gene << "_Frag_" << (position-(ROA_LENGTH/2)+1) << "_" << frag_counts[gene] << endl;;
					fasta_file << UINT64ToString((*left).second.left_sequence_half);
					fasta_file << UINT64ToString((*left).second.right_sequence_half);
					fasta_file << UINT64ToString((*right).second.left_sequence_half);
					fasta_file << UINT64ToString((*right).second.right_sequence_half);
					fasta_file << endl;

				}
	}
	return frags;
}

int printFasta(string output)
{
	fasta_file.open (output + "/fasta.fa");
	int frags = iterate(generateFrags);
	fasta_file.close();
	return frags;
}


void frequency_filter(string gene, int position, int threshold, double ppm, double frag, double non_verified, bool testing, string name, string sequence, int test_position)
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
					if(testing)
						check(total, ppm_threshold, frag_threshold, gene, test_position,  position , name , sequence, read, 10);
					else
						reads[gene][position][(*sequences).first].set_above_nv_threshold();
				else
					if(testing)
						check(total, ppm_threshold, frag_threshold, gene, test_position,  position , name , sequence, read, 1);
					else
						reads[gene][position][(*sequences).first].set_above_frag_threshold();

			else if( testing && (( read.forward_count >= ppm_threshold && read.reverse_count >= ppm_threshold ) || read.matches_reference()))
				check(total, ppm_threshold, frag_threshold, gene,  test_position,  position , name , sequence, read, 0);
			else if(testing)
				check(total, ppm_threshold, frag_threshold, gene,  test_position,  position , name , sequence, read, -1);
		}

	}
}

void verify( map<string, Read> *left_track, map<string, Read> *right_track, bool testing, bool is_right)
{
	if((*left_track).size() > 0 && (*right_track).size() > 0)
		for (auto left = (*left_track).begin(); left != (*left_track).end(); ++left)
			if((*left).second.is_above_ppm())
				for(auto right = (*right_track).begin(); right != (*right_track).end(); ++right)
					if((*right).second.is_above_ppm())
						if((*left).second.right_sequence_half == (*right).second.left_sequence_half)
						{
							if(testing){
								Read read = (is_right) ? (*right).second : (*left).second;
								if(is_right)
									for(int i = 0; i <68; i++)
										cout << " ";
								cout << read.sequence << " (f " <<read.forward_count <<"-" <<read.matches_ref_on_left()<< " + " <<read.reverse_count << "-"<<read.matches_ref_on_right()<< ")" << (( read.is_above_non_verified() )? " > 10% ": (( read.is_above_frag()) ? " > 1% " : (( read.is_above_ppm()) ? " > ppm " :  " :  NOT PPM "))) << " vflags:" << read.is_left_verified() << ":"<<read.is_right_verified()<< endl;
							} else {

								if((*left).second.is_above_frag())
									((*right_track)[(*right).first]).set_left_verified_at_frag();
								else
									((*right_track)[(*right).first]).set_left_verified();

								if ((*right).second.is_above_frag())
									((*left_track)[(*left).first]).set_right_verified_at_frag();
								else
									((*left_track)[(*left).first]).set_right_verified();
							}
						}
}

void check_verify ( Read r, bool is_right, string gene, int position)
{
	map<string, Read> left_track, right_track;
	if(is_right)
	{
		left_track[r.sequence] = r;
		right_track = reads[gene][position];
	} else {
		left_track = reads[gene][position];
		right_track[r.sequence] = r;
	}
	verify(&left_track, &right_track, true, is_right);
}


void sequential(int threshold, double ppm, double frag, double non_verified)
{

	// Frequency Filter.
	for(auto genes = reads.begin(); genes != reads.end(); ++genes)
		for (auto positions = (*genes).second.begin() ; positions != (*genes).second.end(); ++positions)
			frequency_filter((*genes).first, (*positions).first, threshold, ppm, frag, non_verified, 0, "", "", 0);


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

			verify(left_track, track, false, false);

			for(int i = 0; i < size; i ++)
				++positions;

			map<string, Read> * right_track = &reads[(*genes2).first][(*positions).first];

			for(int i = 0; i < size; i ++)
				--positions;

			verify(track, right_track, false, false);

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

ofstream snv_file;
ofstream coverage_file;
ofstream dictionary_file;

int callSNV(int reason, string gene, int pos, int i, uint64_t ref_left,uint64_t ref,uint64_t ref_right, double freq, double cov )
{
	uint64_t convert (i);

	snv_file << setw(10) << gene <<" , ";
	snv_file << reason <<" ,";
	snv_file << setw(5) << (pos+1) <<" , ";
	snv_file << UINT64ToString(ref & 0b111) << " , ";
	snv_file << UINT64ToString(convert) << " , ";
	snv_file << setw(10) << freq << " , ";
	snv_file << UINT64ToStringCompare(ref_left, 0);
	snv_file << UINT64ToString(convert);
	snv_file << UINT64ToStringCompare(ref_right, 0) << " , ";
	snv_file << cov;
	snv_file << endl;
	return 1;
}

void callSNVs(double snv_verified_threshold, double snv_total_threshold, string output )
{
	int snvs = 0;

	snv_file.open (output+"snv.csv");
	snv_file << "Gene, CallReason, Loc, RefBase, CallBase, Freq, LocalSeq, Coverage"<< endl;

	coverage_file.open (output+"coverage.csv");
	coverage_file << "Gene, Loc, Coverage "<< endl;

	for(auto genes = verified_histogram_0.begin(); genes != verified_histogram_0.end(); ++genes)

		for (auto positions = (*genes).second.begin(); positions != (*genes).second.end(); ++positions)
		{

			map<int, int> verified_counts = (*positions).second;
			map<int, int> verified_counts2 = verified_histogram_1[(*genes).first][(*positions).first];

			double verified_total  	= 0;
			double verified_total2 	= 0;
			double ppm_total 		= 0;
			double ppm_total2		= 0;

			for(int i = 1; i < 6; i++)
			{
				verified_total 	+= verified_counts[i];
				verified_total2 += verified_counts2[i];
				ppm_total  		+= ppm_histogram_0[(*genes).first][(*positions).first][i];
				ppm_total2  	+= ppm_histogram_1[(*genes).first][(*positions).first][i];
			}

			double verified = verified_total + verified_total2;
			double ppm = ppm_total + ppm_total2;

			coverage_file << (*genes).first <<","<< (*positions).first << ","<< (ppm / 2) << endl;

			double freq;
			for (int i = 1; i < 6; i++)
				if((references[(*genes).first][(*positions).first] & 0b111) != i)
				{
					uint64_t ref       = references[(*genes).first][(*positions).first];
					uint64_t ref_left  = references[(*genes).first][(*positions).first-ROA_LENGTH/2];
					uint64_t ref_right = references[(*genes).first][(*positions).first+1];

					if(ref){
						// --- Call type #1
						// If the reads are in both verified histograms.
						if( (freq = ((double) (verified_counts[i]+ verified_counts2[i]) / verified)) and verified_counts[i] > 0 and verified_counts2[i] > 0)
							snvs += callSNV(1, (*genes).first,(*positions).first, i, ref_left, ref,ref_right, freq, (verified/2) );

						// --- Call type #2
						// If the reads are only in one histogram
						else if( (freq = ((double) verified_counts[i]) / verified_total) > snv_verified_threshold )
							snvs += callSNV(2, (*genes).first,(*positions).first, i, ref_left, ref,ref_right, freq, verified_total );

						else if( (freq = ((double) verified_counts2[i]) / verified_total2) > snv_verified_threshold )
							snvs += callSNV(2, (*genes).first,(*positions).first, i, ref_left, ref,ref_right, freq, verified_total2 );

						// # Call type 3
						// If the reads exceed a 3rd threshold in either track
						else if( (freq = ((double) (ppm_histogram_0[(*genes).first][(*positions).first][i]) + ppm_histogram_1[(*genes).first][(*positions).first][i]) / ppm ) > snv_total_threshold)
							snvs += callSNV(3, (*genes).first,(*positions).first, i, ref_left, ref,ref_right, freq, (ppm/2) );

					}

				}

		}
	snv_file.close();
	coverage_file.close();
	cout << " I read " << snvs << " SNVs \n";
}

int printDictionaries (string gene, int position, string seq, Read& read)
{
	if(read.is_above_ppm())
	{

		int roa_coverage = 0; // PPM
		for (auto sequences = reads[gene][position].begin(); sequences != reads[gene][position].end(); ++sequences)
			if((*sequences).second.is_above_ppm())
				roa_coverage+=(*sequences).second.total_count();

		int n_diffs = countDifferences(read.left_sequence_half, references[gene][position]) + countDifferences(read.right_sequence_half, references[gene][position+ROA_LENGTH/2]);

		dictionary_file << gene << ","<< position <<", "<< UINT64ToStringCompare(read.left_sequence_half, references[gene][position]);
		dictionary_file << UINT64ToStringCompare(read.right_sequence_half, references[gene][position+ROA_LENGTH/2]) << ", " << roa_coverage <<", ";
		dictionary_file << n_diffs << ", "<< read.is_right_left_verified_at_frag();
		dictionary_file <<"," << read.is_right_left_verified_at_frag() << ", ";
		dictionary_file << read.is_left_verified() << ", " << read.is_right_verified();
		dictionary_file << ", " << read.total_count() << ", " << read.forward_count;
		dictionary_file << ", " << read.reverse_count << ", " << read.cigar_counts[0];
		dictionary_file << ", "<< read.cigar_counts[2]<< ", " << read.cigar_counts[1];
		dictionary_file << ", "<< read.cigar_counts[3]<< ", ";
		dictionary_file << ", Frag1,...,Frag" << endl;
		return 1;
	}
	return 0;
}

void printDicitonaries(string output)
{
	int words = 0;
	dictionary_file.open(output+"dictionary.csv");
	dictionary_file << "Gene, ROAstart, Sequence, ROAcover, Ndiffs, LVerPct, RVerPct, LVerPPM, RVerPPM, Total, Fwd, Rev, NoIndel, DelOnly, InsOnly, InsAndDel, Ref, Frag1,...,Frag" << endl;
	words += iterate(printDictionaries);
	dictionary_file.close();

}

void printHistograms(string output)
{
	for(auto frag = frag_counts.begin(); frag !=frag_counts.end(); ++frag)
		cout << "I have " << (*frag).second << " frags for "<< (*frag).first << endl;

	for(auto genes = verified_histogram_0.begin(); genes != verified_histogram_0.end(); ++genes)
	{
		string filename = output+(*genes).first+"AT.txt";
		string filenameC = output+(*genes).first+"CG.txt";
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

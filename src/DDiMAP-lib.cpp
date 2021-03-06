
//============================================================================
// Name        : DDiMAP-lib.c
// Author      : Androwis Abumoussa
// Version     :
// Copyright   : All Rights Reserved
// Description : DDiMAP library used by other programs.
//============================================================================

#include "DDiMAP-lib.h"
#include <bitset>
#include <math.h>
#include <cctype>
#define TEST 0

int FASTA_ENTRIES = 0;
int ROA_LENGTH;
int DICTIONARY_LEVELS = 0;
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
// Map < Gene --> < Position --> < SEQ --> < BA.refID -> Count> > >
map<string , map<int, map<uint64_t, map <int , int > > > > SNVs;

// Map <Reference Name -> < position -> sequence>
map<string, map<int, std::bitset<Read::half_length> > > references;

std::bitset<Read::half_length> mask (string("111"));


// ----------------------------------------------------------------------------
// Templated code.  I know this is a sin.


// Cipher is given by :
uint64_t a    = 0b001;
uint64_t c    = 0b010;
uint64_t g    = 0b011;
uint64_t t 	  = 0b100;
uint64_t dash = 0b101;

template <int length>
bitset<length> generateMask(uint64_t mask)
{
	return bitset<length> (mask);
}

template <int length>
bitset<length> charToBitset(char ch)
{
	switch(ch){
	case 'A':
	case 'a': return generateMask<length>(a);
	case 'T':
	case 't': return generateMask<length>(t);
	case 'C':
	case 'c': return generateMask<length>(c);
	case 'G':
	case 'g': return generateMask<length>(g);
	case '-': return generateMask<length>(dash);
	default : return generateMask<length>(0);
	}
}

template <int length>
char BitsetToChar(std::bitset<length> ch, bool upper_case)
{
	if(upper_case){
		if( ch == a) return 'A';
		if( ch == t) return 'T';
		if (ch == c) return 'C';
		if (ch == g) return 'G';
	} else {
		if( ch == a) return 'a';
		if( ch == t) return 't';
		if (ch == c) return 'c';
		if (ch == g) return 'g';
	}
	if (ch == dash) return '-';
	return '\0';
}

template <int length>
int countDifferences(std::bitset<length> s, std::bitset<length> t)
{
	int diffs = 0;
	while(s.count()!=0)
	{
		diffs += ((s & generateMask<length>(0b111)) != (t & generateMask<length>(0b111)));
		s = s >> 3; t = t >> 3;
	}
	return diffs;

}

template <int length>
string BitsetToStringCompare(std::bitset<length> s, std::bitset<length> ref)
{
	std::stringstream temp;
	while(s.count()!=0){
		temp << BitsetToChar<length>( (s & generateMask<length>(0b111)) , ( (s & generateMask<length>(0b111)) == (ref & generateMask<length>(0b111))) );
		s = s >> 3; ref = ref >> 3;
	}
	return temp.str();
}

template <int length>
string BitsetToString(std::bitset<length> s)
{ return BitsetToStringCompare<length>(s,s); }


template <int length>
std::bitset<length> stringToBitset(string s)
{
	std::bitset<length> temp = 0;
	for ( int i = 0 ; i < s.length();  i++)
		temp |= (charToBitset<length>(s[i]) << (3 * i) );
	return temp;
}

template <int length>
bool hasDash(std::bitset<length> seq)
{
	while(seq.count()!=0){
		if((seq & generateMask<length>(0b111)) == dash)
			return true;
		seq = seq >> 3;
	}
	return false;
}

template <int length>
float CalculateGC(std::bitset<length> seq)
{
	float gc = 0.0, total = 0.0;
	while(seq.count!=0){
		total += 1;
		if((seq & generateMask<length>(0b111)) == g || (seq & generateMask<length>(0b111)) ==c )
			gc+=1;
		seq = seq >> 3;
	}
	return gc / total * 100.0;
}

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
	r.left_sequence_half  = stringToBitset<Read::half_length>( word.substr(0, length/2));
	r.right_sequence_half = stringToBitset<Read::half_length>( word.substr(length/2 , length/2) );
	return r;
}

// ----------------------------------------------------------------------------
// Reading Files : Convenience Functions
// ----------------------------------------------------------------------------
struct listOfWords
{ vector<string> words;
  vector<int> offset;
};

listOfWords createWordString(BamAlignment &ba, int length, int &position, int track)
{

	listOfWords words;
	string read;

	// Process the CIGAR String
	int i = 0;
	for(auto element = ba.CigarData.begin(); element < ba.CigarData.end(); element++)
	{
		if((*element).Type != 'H' && (*element).Type != 'S'){
			if((*element).Type != 'I'){
				read.append( ba.AlignedBases.substr(i, (*element).Length));
				}
			i+=(*element).Length;
		}
	}


	int to_roa    = ( length/2 - ba.Position) % (length/2);
	if(to_roa < 0) to_roa+=length/2;

	int offset    = (ba.IsReverseStrand()) ? ( read.length() - (ba.Position + read.length())%(length/2))%length  : to_roa ;  // changed offset computation for reverse strand

	if(ba.IsReverseStrand() )
		{
		offset+=track;
                }
	else
		{
		offset-=(length/2-track);
	}

	if(offset < 0) offset+=length/2;
	if( offset + length > read.length()) offset-=length/2;

	position = ba.Position + offset;

	//if(ba.IsReverseStrand() )
	//	{
	//	 cout << "reverseStrand position= "<< ba.Position << " track= " << track << " readLength= " << read.length() << " offset= " << offset << " loop limit = " << read.length() - offset - length << endl;
        //        }
	

	if(offset < 0 || offset + length > read.length())
		return words;


	for(int j = 0; j < read.length() - offset - length; j+=length)  //  changed loop limit to use read.length() rather than i (JPS 7/16/2014)
	{
	string word = read.substr(offset+j, length);
	// check to see if the word has any N characters in it or is too short - if not, add it
	if (word.find("N") ==  word.npos && word.length() == length)   //  added check on words shorter than length (JPS 7/16/2014)
	{
		words.offset.push_back(j);
		words.words.push_back(word);
	}
	//	// ensuring this is correct.
	//	if(TEST){
	//		// Check that there isn't a more appropriate ROA
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
	}
	return words;
}

// This needs to be fixed.
const char *createWordArray(BamAlignment &ba, int length, int &position, int track)
{
	listOfWords words = createWordString(ba, length, position, track);
	return words.words[0].c_str();
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
			listOfWords words   = createWordString(ba, length, position, track);

			for(int i = 0; i < words.words.size(); i++)
				if(words.words[i].size() > 0)
			{

				string name   = genes[ba.RefID];

				// Increment counter for the observed sequence
				if( reads[name][position+words.offset[i]][words.words[i]].total_count() )
					reads[name][position+words.offset[i]][words.words[i]].set_data(hasDeletion, hasInsertion, ba.IsReverseStrand(), ba.RefID);

				// Create a new read for the position on this track
				else {
					Read r;
					r = f(words.words[i], length);
					r.set_data(hasDeletion,hasInsertion, ba.IsReverseStrand(), ba.RefID);

					// Check NCBI
					if(references[name][position+words.offset[i]] == r.left_sequence_half)
						r.set_matches_ref_on_left();

					if(references[name][position+words.offset[i]+length/2] == r.right_sequence_half)
						r.set_matches_ref_on_right();

					reads[name][position+words.offset[i]][words.words[i]] = r;

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
	tracks[1] = floor(roa_length/4);
	cout << "Setting tracks[1]= " << tracks[1] << endl;

	// Read in the NCBI sequences from Fasta File, assign appropriate offsets
	gzFile fp;
	kseq_t *seq;
	FILE *fast = fopen(fasta.c_str(),"r");
	fp = gzdopen(fileno(fast), "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0){

		string seq_name = seq->name.s;
		string s = seq->seq.s;

		genes_names[FASTA_ENTRIES] = seq_name;

		for (size_t i = 0; i < s.length(); ++i)
			if (s[i]!='a' && s[i]!='A' &&
					s[i]!='c' && s[i]!='C' &&
					s[i]!='t' && s[i]!='T' &&
					s[i]!='g' && s[i]!='G')
				s.erase(i, 1);

		if(seq_name.find("Frag")!=-1)
		{
			int loc = seq_name.find("Frag");
			string frag = seq_name.substr(loc, seq_name.length()-loc);
			loc = frag.find_first_of("_")+1;
			string locations = frag.substr(loc, frag.length()-loc);
			frag_offset[FASTA_ENTRIES] = atoi(locations.substr(0,locations.find_first_of("_")).c_str()) - 1;
		}
		else
			frag_offset[FASTA_ENTRIES] = 0;

		int loc = (seq_name.find_first_of("_") == -1) ? seq_name.length() : seq_name.find_first_of("_");

		// If this is not a fragment or junction sequence add to references
		if(seq_name.find("Frag") == -1 && seq_name.find("Junction") == -1){
			seq_name = seq_name.substr(0, loc);
			map<int, std::bitset<Read::half_length> > reference;
			for(int j= 0; j< s.length()-ROA_LENGTH/2; j++){

				reference[j] = stringToBitset<Read::half_length>(s.substr(j, ROA_LENGTH/2));
				// Add reference to ROAs
				if((j % (ROA_LENGTH/2) == tracks[0]  || j % (ROA_LENGTH/2) == tracks[1] ) && (j + ROA_LENGTH) < s.length() ){
					Read r;
					string seq = s.substr(j, ROA_LENGTH);
					r = buildRead( seq , ROA_LENGTH);
					r.set_matches_ref_on_right();
					r.set_matches_ref_on_left();
					r.RefID[FASTA_ENTRIES]++;
					reads[seq_name][j][r.sequence] = r;
				}

			}
			references[seq_name] = reference;
		}
		else{
			seq_name = seq_name.substr(0, loc);
		}

		genes[FASTA_ENTRIES] = seq_name;
		++FASTA_ENTRIES;
	}

	kseq_destroy(seq);
	gzclose(fp);
	if(TEST)
		cout << "I read a total of " << FASTA_ENTRIES << " reads from the fasta file"<< endl;



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
	if( TEST && read.is_right_left_verified()){
		cout << gene << "  " <<  position <<" ";
		cout << BitsetToStringCompare<Read::half_length>(read.left_sequence_half, references[gene][position]) ;
		cout << BitsetToStringCompare<Read::half_length>(read.right_sequence_half, references[gene][position+(ROA_LENGTH/2)]) ;
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

	if( (offset < 0 and read.is_left_verified()) || (offset > 0 and read.is_right_verified()) )
		for (auto seq = roa.begin(); seq != roa.end(); ++seq)
			if( (offset < 0
					and (*seq).second.right_sequence_half == read.left_sequence_half
					and not hasDash<Read::half_length>((*seq).second.left_sequence_half))
					|| (offset > 0
							and (*seq).second.left_sequence_half == read.right_sequence_half
							and not hasDash<Read::half_length>((*seq).second.right_sequence_half) )
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
		bitset<Read::half_length> reference = references[gene][position + ( (offset > 0) ? 2 : 1) * offset];
		if(reference.count()){
			if(offset < 0) // add the left reference
				ref.left_sequence_half = reference;
			else 		   // add the right reference
				ref.right_sequence_half = reference;
			matches[BitsetToString<Read::half_length>(ref.left_sequence_half)+BitsetToString<Read::half_length>(ref.right_sequence_half)] = ref;
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
			not hasDash<Read::half_length>(read.left_sequence_half) and
			not hasDash<Read::half_length>(read.right_sequence_half))
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
					fasta_file << BitsetToString<Read::half_length>((*left).second.left_sequence_half);
					fasta_file << BitsetToString<Read::half_length>((*left).second.right_sequence_half);
					fasta_file << BitsetToString<Read::half_length>((*right).second.left_sequence_half);
					fasta_file << BitsetToString<Read::half_length>((*right).second.right_sequence_half);
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

void check(int total, double ppm, double frag, string gene,  int position,int i, string name, string seq, Read read, int above)
{
	//	 Check the left half.
	if( position == i &&
			read.left_sequence_half == stringToBitset<Read::half_length>(seq.substr(0,ROA_LENGTH/2)) &&
			read.right_sequence_half == stringToBitset<Read::half_length>(seq.substr(ROA_LENGTH/2,ROA_LENGTH/2))){
		check_verify(read, false, gene, (position-ROA_LENGTH/2));
		for (int j = 0; j<ROA_LENGTH/2; j++)
			cout << " ";
		cout << read.sequence;
		for (int j = 0; j<ROA_LENGTH; j++)
			cout << " ";
		cout << " (f " <<read.forward_count <<"-" <<read.matches_ref_on_left()<< " + " <<read.reverse_count << "-"<<read.matches_ref_on_right()<<" / " << ((above > 0) ?  total : ((above == 0) ? frag :  ppm))  << ")" << (( above > 9 )? " > 10% ": (( above > 0) ? " > 1% " : (( above == 0) ? " > ppm " :  " :  NOT PPM "))) << read.is_left_verified() << ":"<<read.is_right_verified()<< endl ;
	}

	//	 Check the middle
	if( position == i+ROA_LENGTH/2
			&& read.left_sequence_half == stringToBitset<Read::half_length>(seq.substr(ROA_LENGTH/2,ROA_LENGTH/2))
			&& read.right_sequence_half == stringToBitset<Read::half_length>(seq.substr(ROA_LENGTH,ROA_LENGTH/2))
	){
		for (int j = 0; j<ROA_LENGTH; j++)
			cout << " ";
		cout << read.sequence;
		for (int j = 0; j<ROA_LENGTH/2; j++)
			cout << " ";
		cout << " (f " <<read.forward_count <<"-" <<read.matches_ref_on_left()<< " + " <<read.reverse_count << "-"<<read.matches_ref_on_right()<<" / " << total  << ")" << (( above > 9 )? " > 10% ": (( above > 0) ? " > 1%" : (( above == 0) ? " > ppm" : "  NOT PPM"))) << endl ;
	}

	// Check the right half
	if( position == i+ROA_LENGTH
			&& read.left_sequence_half == stringToBitset<Read::half_length>(seq.substr(ROA_LENGTH,ROA_LENGTH/2))
			&& read.right_sequence_half == stringToBitset<Read::half_length>(seq.substr(ROA_LENGTH/2+ROA_LENGTH,ROA_LENGTH/2))
	){
		for (int j = 0; j<ROA_LENGTH/2+ROA_LENGTH; j++)
			cout << " ";
		cout << read.sequence;
		cout << " (f " <<read.forward_count <<"-" <<read.matches_ref_on_left()<< " + " <<read.reverse_count << "-"<<read.matches_ref_on_right()<<" / " << ((above > 0) ?  total : ((above == 0) ? frag :  ppm))  << ")" << (( above > 9 )? " > 10% ": (( above > 0) ? " > 1% " : (( above == 0) ? " > ppm " :  " :  NOT PPM "))) << read.is_left_verified() << ":"<<read.is_right_verified()<< endl ;
		check_verify(read, true, gene, (position+ROA_LENGTH/2));

	}

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
									for(int i = 0; i <2*ROA_LENGTH; i++)
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

			for(int i = 0; i < size; i ++)
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

	if( position % (ROA_LENGTH/2) == 0 ){
		verified  =  &verified_histogram_0;
		ppm = &ppm_histogram_0;
	} else {
		verified  =  &verified_histogram_1;
		ppm = &ppm_histogram_1;
	}

	if( read.is_right_left_verified() || read.is_above_ppm() || read.matches_reference())
	{

		// Create a histogram for the left half
		bitset<Read::half_length> s = read.left_sequence_half;
		int i = 0;
		while(s.count())
		{

			if( read.is_right_left_verified() || read.matches_reference())
				(*verified)[gene][position + i][ (s & mask).to_ullong()] += read.total_count();
			if( read.is_above_ppm() || read.matches_reference())
				(*ppm)[gene][position + i][ (s & mask).to_ulong()] += read.total_count();
			s = s>>3;i++;
		}

		// Then the right
		s = read.right_sequence_half;
		while(s.count())
		{
			if( read.is_right_left_verified() || read.matches_reference())
				(*verified)[gene][position + i][ (s & mask).to_ulong()] += read.total_count();
			if( read.is_above_ppm() || read.matches_reference())
				(*ppm)[gene][position + i][ (s & mask).to_ulong()] += read.total_count();
			s=s>>3;i++;
		}

	}
	return 0;
}

ofstream snv_file;
ofstream coverage_file;
ofstream dictionary_file;

int callSNV(int reason, string gene, int pos, int i, std::bitset<Read::half_length> ref_left, std::bitset<Read::half_length> ref,std::bitset<Read::half_length> ref_right, double freq, double cov )
{
	bitset<Read::half_length> convert(i);

	snv_file << setw(10) << gene <<" , ";
	snv_file << reason <<" ,";
	snv_file << setw(5) << (pos+1) <<" , ";
	snv_file << BitsetToString<Read::half_length>(ref & mask) << " , ";
	snv_file << BitsetToString<Read::half_length>(convert) << " , ";
	snv_file << setw(10) << freq << " , ";
	snv_file << BitsetToStringCompare<Read::half_length>(ref_left, 0);
	snv_file << BitsetToString<Read::half_length>(convert);
	snv_file << BitsetToStringCompare<Read::half_length>(ref_right, 0) << " , ";
	snv_file << cov;
	snv_file << endl;
	return 1;
}

void callSNVs(double snv_verified_threshold, double snv_total_threshold, string output )
{
	int snvs = 0;

	snv_file.open (output+"snv.csv");
	snv_file << "RefSeqID, CallReason, Loc, RefBase, CallBase, Freq, LocalSeq, Coverage"<< endl;

	coverage_file.open (output+"coverage.csv");
	coverage_file << "RefSeqID, Loc, Coverage "<< endl;

	auto genes = ppm_histogram_1.begin();
	auto genes1= ppm_histogram_0.begin();

	for(; genes != ppm_histogram_1.end(); ++genes)
	{
		auto positions  =(*genes).second.begin();
		auto positions2 =ppm_histogram_0[(*genes).first].begin();

		for (; positions != (*genes).second.end() && positions2 != ppm_histogram_0[(*genes).first].end();)
		{

			// Set the position : use the lower of the two iterators
			int position = ((*positions).first > (*positions2).first) ? (*positions2).first : (*positions).first;

			map<int, int> verified_counts = verified_histogram_0[(*genes).first][position];
			map<int, int> verified_counts2 = verified_histogram_1[(*genes).first][position];

			double verified_total  	= 0;
			double verified_total2 	= 0;
			double ppm_total 		= 0;
			double ppm_total2		= 0;

			for(int i = 1; i < 6; i++)
			{
				verified_total 	+= verified_counts[i];
				verified_total2 += verified_counts2[i];
				ppm_total  		+= ppm_histogram_0[(*genes).first][position][i];
				ppm_total2  	+= ppm_histogram_1[(*genes).first][position][i];
			}

			double verified = verified_total + verified_total2;
			double ppm = ppm_total + ppm_total2;

			coverage_file << (*genes).first <<","<< position+1 << ","<< (ppm / 2) << endl;

			bitset<Read::half_length> ref       = references[(*genes).first][position];
			bitset<Read::half_length> ref_left  = references[(*genes).first][position-ROA_LENGTH/2];
			bitset<Read::half_length> ref_right = references[(*genes).first][position+1];

			double freq;
			for (int i = 1; i < 6; i++)
				if((references[(*genes).first][position] & mask) != i)
					if(ref.count())
					{
						// --- Call type #1
						// If the reads are in both verified histograms.
						if( (freq = ((double) (verified_counts[i]+ verified_counts2[i]) / verified)) and verified_counts[i] > 0 and verified_counts2[i] > 0)
							snvs += callSNV(1, (*genes).first,position, i, ref_left, ref,ref_right, freq, (verified/2) );

						// --- Call type #2
						// If the reads are only in one histogram
						else if( (freq = ((double) verified_counts[i]) / verified_total) > snv_verified_threshold )
							snvs += callSNV(2, (*genes).first,position, i, ref_left, ref,ref_right, freq, verified_total );

						else if( (freq = ((double) verified_counts2[i]) / verified_total2) > snv_verified_threshold )
							snvs += callSNV(2, (*genes).first,position, i, ref_left, ref,ref_right, freq, verified_total2 );

						// --- Call type #3
						// If the reads exceed a 3rd threshold in either track
						else if( (freq = ((double) ppm_histogram_0[(*genes).first][position][i] + ppm_histogram_1[(*genes).first][position][i]) / ppm ) > snv_total_threshold)
							snvs += callSNV(3, (*genes).first,position, i, ref_left, ref,ref_right, freq, (ppm/2) );
					}


			// Incrementing the indexes.
			// This handles error conditions when traversing two iterators that may put us in an inf loop.
			if((*positions).first == (*positions2).first)
			{
				if(positions != (*genes).second.end())
					positions++;
				if(positions != ppm_histogram_0[(*genes).first].end())
					positions2++;
			}
			else if((*positions).first > (*positions2).first)
			{
				if(positions != ppm_histogram_0[(*genes).first].end())
					positions2++;
				else if(positions != (*genes).second.end())
					positions++;
			}
			else
			{
				if(positions != (*genes).second.end())
					positions++;
				else if(positions != ppm_histogram_0[(*genes).first].end())
					positions2++;
			}

		}
	}
	snv_file.close();
	coverage_file.close();
	cout << " I read " << snvs << " SNVs \n";
}

map<string, ofstream> dictionaries;

int printDictionaries (string gene, int position, string seq, Read& read)
{
	ofstream * dict;

	if(read.is_above_ppm() || read.matches_reference())
	{

		if(DICTIONARY_LEVELS > 1)
			dict = &dictionaries[gene];
		else
			dict = &dictionary_file;

		int roa_coverage = 0; // PPM
		for (auto sequences = reads[gene][position].begin(); sequences != reads[gene][position].end(); ++sequences)
			if((*sequences).second.is_above_ppm())
				roa_coverage+=(*sequences).second.total_count();

		int n_diffs = countDifferences<Read::half_length>(read.left_sequence_half, references[gene][position]) + countDifferences<Read::half_length>(read.right_sequence_half, references[gene][position+ROA_LENGTH/2]);

		(*dict) << gene << ","<< (position+1) <<", ";
		(*dict) << BitsetToStringCompare<Read::half_length>(read.left_sequence_half, references[gene][position]);
		(*dict) << BitsetToStringCompare<Read::half_length>(read.right_sequence_half, references[gene][position+ROA_LENGTH/2]);
		(*dict) << ", " << roa_coverage <<", ";
		(*dict) << n_diffs << ", "<< read.is_left_verified_at_frag();
		(*dict) <<"," << read.is_right_verified_at_frag() << ", ";
		(*dict) << read.is_left_verified() << ", " << read.is_right_verified();
		(*dict) << ", " << read.total_count() << ", " << read.forward_count;
		(*dict) << ", " << read.reverse_count;
		if(DICTIONARY_LEVELS > 0)
		{
			(*dict) << ", " << read.cigar_counts[0] << ", "<< read.cigar_counts[2];
			(*dict) << ", " << read.cigar_counts[1] << ", "<< read.cigar_counts[3];
		}
		if(DICTIONARY_LEVELS > 1)
		{
			for(int i = 0; i < FASTA_ENTRIES; i++)
				if(genes[i] == gene)
					(*dict) << ", " << read.RefID[i] ;
		}
		(*dict) << endl;
		return 1;
	}
	return 0;
}

void printDicitonaries(string output, int level)
{
	DICTIONARY_LEVELS = level;
	int words = 0;
	dictionary_file.open(output+"dictionary.csv");

	// Include seperate files for each gene.
	if(DICTIONARY_LEVELS > 1)
	{
		map<string, int> files;
		for(int i = 0; i < FASTA_ENTRIES; i++)
			if(files[genes[i]]==0){
				files[genes[i]]++;
				dictionaries[genes[i]].open(output+"dictionary-"+genes[i]+".csv");
			}

		for(auto files = dictionaries.begin(); files != dictionaries.end(); ++files)
			(*files).second << "RefSeqID, ROAstart, Sequence, ROAcover, EditDist, LVerFrag, RVerFrag, LVerSNV, RVerSNV, WordCover, FwdCover, RevCover ";
	}
	else
		dictionary_file << "RefSeqID, ROAstart, Sequence, ROAcover, EditDist, LVerFrag, RVerFrag, LVerSNV, RVerSNV, WordCover, FwdCover, RevCover ";

	if(DICTIONARY_LEVELS > 0)
	{
		if(DICTIONARY_LEVELS > 1)
			for(auto files = dictionaries.begin(); files != dictionaries.end(); ++files)
				(*files).second << ", #NoIndel, #DelOnly, #InsOnly, #InsAndDel";
		else
		dictionary_file << ", #NoIndel, #DelOnly, #InsOnly, #InsAndDel";
	}

	if(DICTIONARY_LEVELS > 1)
		for(int i = 0; i < FASTA_ENTRIES; i++)
			dictionaries[genes[i]] << ", " << genes_names[i] ;

	if(DICTIONARY_LEVELS > 1)
		for(auto files = dictionaries.begin(); files != dictionaries.end(); ++files)
			(*files).second << endl;
	else
		dictionary_file<< endl;

	words += iterate(printDictionaries);

	if(DICTIONARY_LEVELS > 1)
		for(auto files = dictionaries.begin(); files != dictionaries.end(); ++files)
			(*files).second.close();

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

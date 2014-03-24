#include "DDiMAP-test.h"

KSEQ_INIT(gzFile, gzread)

void check(int total, double ppm, double frag, string gene,  int position,int i, string name, string seq, Read read, int above)
{
	//	 Check the left half.
	if( position == i &&
			read.left_sequence_half == stringToUINT64(seq.substr(0,17)) &&
			read.right_sequence_half == stringToUINT64(seq.substr(17,17))){
		check_verify(read, false, gene, (position-17));
		for (int j = 0; j<17; j++)
			cout << " ";
		cout << read.sequence;
		for (int j = 0; j<34; j++)
			cout << " ";
		cout << " (f " <<read.forward_count <<"-" <<read.matches_ref_on_left()<< " + " <<read.reverse_count << "-"<<read.matches_ref_on_right()<<" / " << ((above > 0) ?  total : ((above == 0) ? frag :  ppm))  << ")" << (( above > 9 )? " > 10% ": (( above > 0) ? " > 1% " : (( above == 0) ? " > ppm " :  " :  NOT PPM "))) << read.is_left_verified() << ":"<<read.is_right_verified()<< endl ;
	}

	//	 Check the middle
	if( position == i+17
			&& read.left_sequence_half == stringToUINT64(seq.substr(17,17))
			&& read.right_sequence_half == stringToUINT64(seq.substr(34,17))
	){
		for (int j = 0; j<34; j++)
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
		for (int j = 0; j<51; j++)
			cout << " ";
		cout << read.sequence;
		cout << " (f " <<read.forward_count <<"-" <<read.matches_ref_on_left()<< " + " <<read.reverse_count << "-"<<read.matches_ref_on_right()<<" / " << ((above > 0) ?  total : ((above == 0) ? frag :  ppm))  << ")" << (( above > 9 )? " > 10% ": (( above > 0) ? " > 1% " : (( above == 0) ? " > ppm " :  " :  NOT PPM "))) << read.is_left_verified() << ":"<<read.is_right_verified()<< endl ;
		check_verify(read, true, gene, (position+17));

	}

}

void check_frequency_filter(string gene, int position, string name, string sequence)
{
	int VERIFY_THRESHOLD  =  2;
	double PPM = 0.00075;
	double FRAG_THRESHOLD = .01;
	double NON_VERIFIED_THRESHOLD = .1;

	for(int i = 0; i < 3; i++)
		frequency_filter(gene, position, VERIFY_THRESHOLD, PPM, FRAG_THRESHOLD, NON_VERIFIED_THRESHOLD, true, name, sequence, (position + i*17));
}


void test()
{

	map <string, string> johns, johns_alph;
	map <string, string> mine, mine_alph;

	char  *fasta = "/Users/androwis/Downloads/TestCaseFixedHigh-2.fa";
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
		string replacement="_";
		string clean = std::regex_replace (seq_name,e,replacement,std::regex_constants::basic);

		string s = seq->seq.s;
		s.erase( std::remove_if( s.begin(), s.end(), ::isspace ), s.end() );

		std::smatch m;
		std::regex_search(clean, m, frag);

		if(clean.find("Frag") && clean.find("NCBI") == -1 && clean.find("Junction") == -1)
			johns[s] = clean;
	}

	FILE *fast2 = fopen(fasta2,"r");
	fp2 = gzdopen(fileno(fast2), "r");
	seq2 = kseq_init(fp2);
	while (kseq_read(seq2) >= 0){

		string seq_name = seq2->name.s;
		string replacement="_";
		string clean = std::regex_replace (seq_name,e,replacement,std::regex_constants::basic);

		string s = seq2->seq.s;
		s.erase( std::remove_if( s.begin(), s.end(), ::isspace ), s.end() );

		std::smatch m;
		std::regex_search(clean, m, frag);

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
		cout <<  endl;
		for(int i = 0; i < 17; i++)
			cout << " ";
		cout << (*my).second<< " "<< (*my).first << " <-- me " << gene << " : "<< (a++)<< endl;
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
		cout << endl;
		for(int i = 0; i < 17; i++)
			cout << " ";
		cout << (*his).second << " " << (*his).first << " <-- John " << (b++)<< endl;
		check_frequency_filter(gene, position - 1, (*his).first , (*his).second);
	}

	cout << " I had " << a <<" / " << mine.size() << " unique | John had "<< b <<" / " << johns.size()<< " unique"<< endl;



}

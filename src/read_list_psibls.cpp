#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
using namespace std;

int main(int argc,char *argv[])
{
	string line,f,sign,s,uniprot_db;
	ifstream infile(argv[1]);
	uniprot_db = argv[2];
	while(getline(infile,line))
	{
		f="src/psiblast -query in_data/sequences/"+line+" -db database/uniprot_sequence/"+uniprot_db+" -out tmp_data/alignments_psibls/"+line+" -max_target_seqs 20000 -evalue 0.01 -num_threads 3 -outfmt 5 -num_iterations 3";
		system(f.c_str());
		f="src/get_blast_aln_psibls "+line;
		system(f.c_str());
		f="rm -f tmp_data/alignments_psibls/"+line;
		system(f.c_str());
	}
	infile.close();
	return 0;
}


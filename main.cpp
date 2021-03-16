#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "CrisprGroup.h"
#include "Read.h"
#include "pameval.h"
#include "Write.h"
#include "sqlite3.h"

using namespace std;

//int main(int argc, char* argv[])
int main()
{
	clock_t tStart = clock();
	//vector<string> argv = { "Executable","spCas9","NGG","TRUE","FALSE","4","16","0","gal","C:/Users/Tfry/Desktop/Recommended_CSPR_Files/","C:/Users/Tfry/Desktop/CASPERinfo","C:/Users/Tfry/Desktop/Recommended_CSPR_Files/FNA Files/gallus_gallus.fna", "gallus gallus", "notes_go_here" };
	//vector<string> argv = { "Executable","spCas9","NGG", "FALSE","FALSE","TRUE","4","16","0","baccoa","C:/Users/Tfry/Desktop/","C:/Users/Tfry/Desktop/CASPERinfo","C:/Users/Tfry/Desktop/baccoa.fna", "bacillus coagulans", "notes_go_here" };
	//vector<string> argv = { "Executable","spCas9","NGG", "TRUE","FALSE","TRUE","4","16","0","baccoa","C:/Users/Tfry/Desktop/","C:/Users/Tfry/Desktop/CASPERinfo","C:/Users/Tfry/Desktop/Recommended_CSPR_Files/FNA FIles/bacillus_coagulans.fna", "bacillus coagulans", "notes_go_here" };

	//variables
	CrisprGroup genome;
	Read read;
	pameval PamEval;
	Write write;
	string endo, pam, org_code, output_path, score_file, input_file, org_name, notes, cspr_filename, db_filename, pam_regex;
	int five_length, seed_length, three_length, seq_length, pam_length, total = 0;
	bool directionality, reps;
	bool strand = true;
	bool mt = true;
	vector<string> chroms;
	vector<string> sequences;
	vector<unsigned long> compressed_seeds, kstats;
	vector<long> seed_locs;
	vector<int> seed_cnts;
	vector<int> uniques, repeats;

	//read in argv
	endo = string(argv[1]);
	pam = string(argv[2]);
	if (string(argv[3]) == "TRUE")
	{
		mt = true;
	}
	else
	{
		mt = false;
	}
	if (string(argv[4]) == "TRUE")
	{
		directionality = true;
	}
	else
	{
		directionality = false;
	}
	if (string(argv[5]) == "TRUE")
	{
		reps = true;
	}
	else
	{
		reps = false;
	}
	five_length = stoi(string(argv[6]));
	seed_length = stoi(string(argv[7]));
	three_length = stoi(string(argv[8]));
	org_code = string(argv[9]);
	output_path = string(argv[10]);
	score_file = string(argv[11]);
	input_file = string(argv[12]);
	org_name = string(argv[13]);
	notes = string(argv[14]);
	seq_length = five_length + seed_length + three_length;
	pam_length = pam.size();
	total = seq_length + pam_length;
	pam_regex = PamEval.regexPAM(pam);

	//filenames
	cspr_filename = output_path + org_code + "_" + endo + ".cspr";
	db_filename = output_path + org_code + "_" + endo + "_repeats.db";
	
	//read input file
	cout << "Reading input file." << endl;
	read.read_file(input_file, sequences, chroms, kstats);
	seed_cnts.resize(sequences.size());
	
	cout << "Number of Chromosomes/Scaffolds: " << sequences.size() << endl;

	//find pams
	cout << "Finding Targets." << endl;
	genome.findPAMs(directionality, mt, sequences, compressed_seeds, seed_locs, seed_cnts, strand, pam_regex, pam_length, seq_length, five_length, seed_length);

	//process targets
	cout << "Processing Targets." << endl;
	genome.process_targets(uniques, repeats, compressed_seeds, seed_locs);

	//writing data
	if (directionality)
	{
		cout << "Writing out uniques." << endl;
		write.write_uniques_dir(uniques, sequences, seed_locs, seed_cnts, kstats, org_name, cspr_filename, score_file, chroms, notes, pam_length, seq_length);
		cout << "Writing out repeats." << endl;
		write.write_repeats_dir(db_filename, repeats, sequences, seed_locs, compressed_seeds, seed_cnts, score_file, five_length, three_length, seed_length, pam_length, seq_length);
	}
	else
	{
		cout << "Writing out uniques." << endl;
		write.write_uniques(uniques, sequences, seed_locs, seed_cnts, kstats, org_name, cspr_filename, score_file, chroms, notes, pam_length, seq_length);
		cout << "Writing out repeats." << endl;
		write.write_repeats(db_filename, repeats, sequences, seed_locs, compressed_seeds, seed_cnts, score_file, five_length, three_length, seed_length, pam_length, seq_length);
	}
	
	cout << "Finished." << endl;
	printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	//system("pause");
	return 0;

}

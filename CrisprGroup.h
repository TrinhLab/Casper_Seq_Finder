
//
//  CrisprGroup.h
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#include <iostream>
#include <vector>
#include <string>
#include <regex>
#include <thread>
#include <numeric>
#include <math.h>
using namespace std;

class CrisprGroup
{
public:

	CrisprGroup();
	~CrisprGroup();

	void findPAMs(bool dir, bool mt, vector<string> &sequences, vector<unsigned long> &compressed_seeds, vector<long> &seed_locs, vector<int> &chroms, bool &strand, string &pam_regex, int &pam_length, int &seq_length, int &five_length, int &seed_length);
	void process_targets(vector<int> &uniques, vector<int> &repeats, vector<unsigned long> &compressed_seeds, vector<long> &seed_locs);
	void reverseComplement(string &str);
	unsigned long compressSeq(string &s);
	int convertCharBase4(char &c);
	void find_seeds(vector<long> &l, int &c, vector<unsigned long> &comp_seeds, string &seq_pointer, regex &pam, int chrom, int &seq_length, int &five_length, int &seed_length);
	void find_seeds_dir(vector<long> &l, int &c, vector<unsigned long> &comp_seeds, string &seq_pointer, regex &pam, int chrom, int &pam_length, int &seq_length, int &five_length, int &seed_length);

};
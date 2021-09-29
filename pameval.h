
//
//  pameval.h
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#pragma once

#include <vector>
#include <string>
#include <map>
#include <iostream>

using namespace std;

class pameval
{
public:
	pameval(string &pam)
	{
		
		valid_pam_characters.insert(pair<char, vector<char>>('A', vector<char> {'A'}));
		valid_pam_characters.insert(pair<char, vector<char>>('C', vector<char> {'C'}));
		valid_pam_characters.insert(pair<char, vector<char>>('G', vector<char> {'G'}));
		valid_pam_characters.insert(pair<char, vector<char>>('T', vector<char> {'T'}));
		valid_pam_characters.insert(pair<char, vector<char>>('R', vector<char> {'A', 'G'}));
		valid_pam_characters.insert(pair<char, vector<char>>('Y', vector<char> {'C', 'T'}));
		valid_pam_characters.insert(pair<char, vector<char>>('S', vector<char> {'G', 'C'}));
		valid_pam_characters.insert(pair<char, vector<char>>('W', vector<char> {'A', 'T'}));
		valid_pam_characters.insert(pair<char, vector<char>>('K', vector<char> {'G', 'T'}));
		valid_pam_characters.insert(pair<char, vector<char>>('M', vector<char> {'A', 'C'}));
		valid_pam_characters.insert(pair<char, vector<char>>('B', vector<char> {'C', 'G', 'T'}));
		valid_pam_characters.insert(pair<char, vector<char>>('D', vector<char> {'A', 'G', 'T'}));
		valid_pam_characters.insert(pair<char, vector<char>>('H', vector<char> {'A', 'C', 'T'}));
		valid_pam_characters.insert(pair<char, vector<char>>('V', vector<char> {'A', 'C', 'G'}));
		valid_pam_characters.insert(pair<char, vector<char>>('N', vector<char> {'A', 'C', 'G', 'T'}));
		
		global_pam = pam;
		pam_length = pam.size();
	}

	string regexPAM(string &pam);

	// Function for assigning a regex code for a degenerate nucleotide code
	string degenerateRegex(char &c);

	void generatePamsWrapper()
	{
		generatePams("", 0);
	}

	void generatePams(string curr_pam, int i)
	{
		if (i == pam_length)
		{
			pam_list.push_back(curr_pam);
			return;
		}

		for (int j = 0; j < valid_pam_characters[global_pam[i]].size(); j++)
		{
			generatePams(curr_pam + valid_pam_characters[global_pam[i]][j], i + 1);
		}
	}

	vector<string> pam_list;

private:
	map<char, vector<char>> valid_pam_characters;
	string global_pam;
	int pam_length;
};

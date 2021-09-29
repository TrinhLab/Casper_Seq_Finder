//
//  Scoring.h
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#pragma once

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <assert.h>
#include <map>
#include <regex>
#include "pameval.h"

using namespace std;

class Scoring
{
public:
	Scoring(string &file, string &on_score_data, bool &directionality, string &endo, int &gRNA_length, string &pam)
	{ 
		fillScoringAlgorithm(file, on_score_data);

		endo_private = endo;
		gRNA_len_private = gRNA_length;
		directionality_private = directionality;

		pam_scores.push_back(vector<float> {1, 1, 1, 1, 1.66, 2.5, 5});
		pam_scores.push_back(vector<float> {1, 1, 1, 1, 1.66, 2.5, 5});
		pam_scores.push_back(vector<float> {1, 1, 1, 1.66, 3.33, 3.33, 0});
		pam_scores.push_back(vector<float> {2, 2, 2.5, 2.5, 3.33, 0, 0});
		pam_scores.push_back(vector<float> {2, 2, 2.5, 2.5, 0, 0, 0, 0});
		pam_scores.push_back(vector<float> {8, 8, 10, 0, 0, 0, 0});
		pam_scores.push_back(vector<float> {10, 10, 0, 0, 0, 0, 0});

		pam_length = pam.size();
	}
	void fillScoringAlgorithm(string &file, string &on_score_data);
	vector<string> Msplit(const string &text, char sep);
	float scoreSequence(string gRNA, string full_sequence, pameval &PamEval);

private:
	//functions
	float get_sc(string &sequence);
	float get_sij(string &sequence, pameval &PamEval);
	float get_sg(string &sequence);
	float get_p(float &sij_score, float &sg_score);
	float ggg_penalty(string &sequence);
	string reverseComplement(string &str);

	//variables
	bool directionality_private = true;
	string endo_private = "";
	int gRNA_len_private = 0;
	map<int, map<std::string, double>> CRISPRSCAN_data;
	vector<vector<float>> pam_scores;
	int pam_length;
};

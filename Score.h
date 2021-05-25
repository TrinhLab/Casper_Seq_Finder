//
//  Scoring.h
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <assert.h>

using namespace std;

class Scoring
{
public:
	Scoring(string &file, string &on_score_data) { fillScoringAlgorithm(file, on_score_data); };
	void fillScoringAlgorithm(string &file, string &on_score_data);
	vector<string> Msplit(const string &text, char sep);
	void scanScore();
	double calcScore(string &s);

private:
	struct iden
	{
		char nt1;
		char nt2;
		int position;
		double odds_score;
	};
	vector<iden> Idens;
	string sequence;
	double totalScore;
	int returnScore;
};

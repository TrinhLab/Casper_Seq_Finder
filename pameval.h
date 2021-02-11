
//
//  pameval.h
//  Casper_Seq_Finder
//
//  Created by Brian Mendoza on 3/26/16.
//  Copyright (c) 2016 Brian Mendoza. All rights reserved.
//

#include <vector>
#include <string>

using namespace std;

class pameval
{
public:
	string regexPAM(string &pam);

	// Function for assigning a regex code for a degenerate nucleotide code
	string degenerateRegex(char &c);
};

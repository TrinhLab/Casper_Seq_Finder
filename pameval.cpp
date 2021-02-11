#include "pameval.h"

string pameval::regexPAM(string &pam)
{
	string retpam = "(?=(";  //initializes lookahead structure
			// Iterate across the pam and generate regex characters for degenerates
	for (int i = 0; i < pam.size(); i++) {
		if (pam[i] == 'N') {
			retpam += '.';
		}
		else if (pam[i] == 'A' || pam[i] == 'T' || pam[i] == 'C' || pam[i] == 'G') {
			retpam += pam[i];
		}
		else {
			retpam += degenerateRegex(pam[i]);
		}
	}
	// finish the lookahead regex structure:
	retpam += ")).";
	return retpam;
}

string pameval::degenerateRegex(char &c)
{
	string reg_s;
	switch (c)
	{
	case 'W':
		reg_s = "[AT]";
		break;
	case 'S':
		reg_s = "[CG]";
		break;
	case 'M':
		reg_s = "[AC]";
		break;
	case 'K':
		reg_s = "[GT]";
		break;
	case 'R':
		reg_s = "[AG]";
		break;
	case 'Y':
		reg_s = "[CT]";
		break;
	case 'B':
		reg_s = "[TCG]";
		break;
	case 'D':
		reg_s = "[ATG]";
		break;
	case 'H':
		reg_s = "[ATC]";
		break;
	case 'V':
		reg_s = "[AGC]";
		break;
	}
	return reg_s;
}
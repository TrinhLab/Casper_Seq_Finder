#include "Score.h"

using namespace std;
void Scoring::fillScoringAlgorithm(string &file, string &on_score_data)
{
	//open CASPERinfo file
	ifstream fin(file);
	string line;

	//make sure file is open before reading
	if (fin.is_open())
	{
		//loop through each line in CASPERinfo file
		while (getline(fin, line))
		{
			if (line.find(on_score_data) != string::npos)
			{
				//loop through each line in the CRISPRSCAN section
				while (getline(fin, line))
				{
					if (line.find("--") != string::npos)
					{
						break;
					}
					else
					{
						//parse line of CRISPRSCAN data
						vector<string> string_split = Msplit(line, '\t');
						string chars = string_split[0];
						int location = stoi(string_split[1]) - 1;
						double score = stod(string_split[2]);

						//check if primary key exists
						if (CRISPRSCAN_data.find(location) != CRISPRSCAN_data.end())
						{
							CRISPRSCAN_data[location].insert(pair<std::string, double>(chars, score));
						}
						//if key not found, create new entry
						else
						{
							map<std::string, double> temp;
							temp.insert(pair<std::string, double>(chars, score));
							CRISPRSCAN_data.insert(pair<int, map<std::string, double>>(location, temp));
						}
					}
				}
				break;
			}
		}
		//close CASPERinfo file
		fin.close();
	}
	else
	{
		cerr << "Unable to open CASPERinfo file: " << file << endl;
		exit(-1);
	}
};

float Scoring::scoreSequence(string gRNA, string full_sequence, pameval &PamEval)
{
	float score = 0.0;
	float sij_score = get_sij(gRNA, PamEval);
	float sc_score = get_sc(full_sequence);
	float sg_score = get_sg(gRNA);
	float p_score = get_p(sij_score, sg_score);

	if (endo_private != "spCas9" && directionality_private == false)
	{
		p_score += ggg_penalty(gRNA);
	}

	if (p_score == 0)
	{
		score = sc_score;
	}
	else
	{
		score = sc_score / p_score;
	}

	score = score * 100;
	
	if (score <= 100 && score >= 0)
	{
		return score;
	}
	else if (score < 0)
	{
		return 0;
	}
	else
	{
		return 100;
	}
	
}

float Scoring::get_sc(string &sequence)
{
	float score = 0;
	string dnt = "";
	string nt = "";
	if (directionality_private)
	{
		reverse(sequence.begin(), sequence.end());
	}
		
	for (int i = 0; i < sequence.size(); i++)
	{
		nt = sequence[i];
		if (i == sequence.size() - 1)
		{
			dnt = "";
		}
		else
		{
			dnt = sequence.substr(i, 2);
		}
		if (CRISPRSCAN_data.find(i) != CRISPRSCAN_data.end())
		{
			if (CRISPRSCAN_data[i].find(nt + "x") != CRISPRSCAN_data[i].end())
			{
				score += CRISPRSCAN_data[i][nt + "x"];
			}
			if (CRISPRSCAN_data[i].find(dnt) != CRISPRSCAN_data[i].end())
			{
				score += CRISPRSCAN_data[i][dnt];
			}
		}
	}

	if (directionality_private)
	{
		reverse(sequence.begin(), sequence.end());
	}

	score += 0.183930944;
	if (score <= 0)
	{
		return 0;
	}
	
	return score;
	
}

float Scoring::get_sij(string &sequence, pameval &PamEval)
{
	//vars
	string reverse_comp = reverseComplement(sequence);
	string seq = sequence;
	string temp;
	string rev_temp;
	int count = 0;
	int rev_count = 0;
	float pam_penalty = 0;
	long pos = 0;
	long pos_rev = 0;
	int cnt = 0;

	for (int i = 0; i <= seq.size() - pam_length; i++)
	{
		//get subtrings of pam length
		temp = seq.substr(i, pam_length);
		rev_temp = reverse_comp.substr(i, pam_length);
		if (find(PamEval.pam_list.begin(), PamEval.pam_list.end(), temp) != PamEval.pam_list.end())
		{
			count += 1;
		}
		if (find(PamEval.pam_list.begin(), PamEval.pam_list.end(), rev_temp) != PamEval.pam_list.end())
		{
			rev_count += 1;
		}
	}

	if (rev_count <= 6 && count <= 6)
	{
		pam_penalty = pam_scores[count][rev_count];
	}

	return pam_penalty;
}

float Scoring::get_sg(string &sequence)
{
	int score = 0;
	string nt = "";
	for (int i = 0; i < sequence.size(); i++)
	{
		nt = sequence[i];
		if (nt == "G")
		{
			score += 10;
		}
		else if (nt == "C")
		{
			score += 5;
		}
		else if (nt == "A")
		{
			score -= 1;
		}
	}

	return float(score/10) / float(gRNA_len_private);
}

float Scoring::get_p(float &sij_score, float &sg_score)
{
	float p = 0;
	if (sij_score > 1)
	{
		p = sij_score * sg_score;
	}
	else if (sij_score == 1)
	{
		p = 1 - (sg_score / 5);
	}
	else
	{
		p = 0;
	}
	return p;
}

float Scoring::ggg_penalty(string &sequence)
{
	int ggg_count = 0;
	for (int i = 0; i < sequence.size() - 2; i++)
	{
		if (sequence.substr(i, 3) == "GGG")
		{
			ggg_count += 1;
		}
	}

	if (ggg_count < 2)
	{
		return 1.0;
	}
	else if (ggg_count == 2)
	{
		return 0.85;
	}
	else if (ggg_count == 3)
	{
		return 0.7;
	}
	else
	{
		return 0.5;
	}
}

string Scoring::reverseComplement(string &str)
{
	string rc = "";
	for (long i = str.size() - 1; i >= 0; i--)
	{
		char n = str[i];
		char reverse;
		switch (n) {
		case 'A': reverse = 'T'; break;
		case 'T': reverse = 'A'; break;
		case 'G': reverse = 'C'; break;
		case 'C': reverse = 'G'; break;
		default: reverse = 'N';
		}
		rc += reverse;
	}
	return rc;
}

vector<string> Scoring::Msplit(const string &text, char sep)
{
	vector<string> tokens;
	size_t start = 0, end = 0;
	while ((end = text.find(sep, start)) != string::npos)
	{
		tokens.push_back(text.substr(start, end - start));
		start = end + 1;
	}
	tokens.push_back(text.substr(start));
	return tokens;
}
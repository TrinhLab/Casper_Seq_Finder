#include "Score.h"

using namespace std;

void Scoring::fillScoringAlgorithm(string &file)
{
	//Establish file stream
	//ifstream* stream = new ifstream;
	//stream->open(file.c_str(), ifstream::in);

	ifstream stream;
	stream.open(file.c_str());

	//Get to the CRISPRSCAN data section of the file:
	std::string myline = "myline";
	while (getline(stream, myline))
	{
		if(myline.find("RISPRSCAN_DATA") != string::npos)
		{
			break;
		}
		//getline(stream, myline);
	}

	//Load Information
	string nts;
	while (getline(stream, nts))
	{
		if(nts.find("---") != string::npos)
		{
			break;
		}
		iden nid;
		vector<string> mytoke = Msplit(nts, '\t');
		nid.nt1 = mytoke[0][0];
		nid.nt2 = mytoke[0][1];
		nid.position = stoi(mytoke[1]);
		nid.odds_score = stod(mytoke[2]);
		Idens.push_back(nid);
	}
	stream.close();
};

double Scoring::calcScore(string &s)
{
	totalScore = 0;
	returnScore = 0;
	sequence = s;
	scanScore();
	//following line normalizes to best possible score
	totalScore = 1 - ((1.29401 - totalScore) / 1.94947);
	returnScore = (totalScore * 100) + 0.5; //0.5 for proper rounding
	return returnScore;
}

void Scoring::scanScore()
{
	for (int i = 0; i < Idens.size(); i++)
	{
		char nucleo1 = Idens.at(i).nt1;  //may have to change this to pointers b/c of multiple creations of object
		char nucleo2 = Idens.at(i).nt2;
		int pos = Idens.at(i).position;
		if (pos < sequence.size())
		{
			if (nucleo2 != 'x')
			{
				string dinucleo = string() + nucleo1 + nucleo2;
				if (sequence.substr(pos - 1, 2) == dinucleo)
				{
					totalScore += Idens.at(i).odds_score;
				}
			}
			else
			{
				if (sequence.at(pos - 1) == nucleo1)
				{
					totalScore += Idens.at(i).odds_score;
				}
			}
		}
	}
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
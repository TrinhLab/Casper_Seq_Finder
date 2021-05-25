#include "Read.h"

using namespace std;

void Read::read_file(string &filename, vector<string> &sequences, vector<string> &chroms, vector<unsigned long> &kstats)
{
	//variables
	ifstream fin;
	fin.open(filename);
	string line;
	string seq = "";
	int cnt = 1;
	while (getline(fin, line))
	{
		if (line[0] == '>')
		{
			chroms.push_back(line + " (" + to_string(cnt) + ")");
			if (cnt != 1)
			{
				transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
				sequences.push_back(seq);
				seq = "";
			}
			cnt++;
		}
		else
		{
			seq += line;
		}
	}
	transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
	//cout << seq << endl;
	sequences.push_back(seq);
	fin.clear();
	fin.close();
	
	//fill in kstats
	for (int i = 0; i < sequences.size(); i++)
	{
		kstats.push_back(sequences[i].size());
	}
}

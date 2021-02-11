#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

class Read
{
public:
	Read() { };
	void read_file(string &filename, vector<string> &sequences, vector<string> &chroms, vector<unsigned long> &kstats);
};

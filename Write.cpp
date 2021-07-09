#include "Write.h"
using namespace std;

//abs comparison used for sorting vector of longs
bool abs_cmp(long i1, long i2)
{
	return (abs(i1) < abs(i2));
}

void Write::write_uniques(vector<int> &uniques, vector<string> &sequences, vector<long> &seed_locs, vector<int> &seed_cnts, vector<unsigned long> &kstats, string &org_name, string &filename, string &score_file, vector<string> &chroms, string &notes, int &pam_length, int &seq_length, string &on_target_data)
{
	//variables
	Scoring score(score_file, on_target_data);
	string comp, seq, genome, kstat, misc;
	long pos = 0;
	int i = 0;
	vector<long> temp;
	int j = 0;
	int running_cnt = 0;
	ofstream outputfile;
	outputfile.open(filename, ios_base::out | ios_base::binary);
	boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf;
	outbuf.push(boost::iostreams::gzip_compressor());
	outbuf.push(outputfile);
	//Convert streambuf to ostream
	ostream out(&outbuf);

	//output first details - org name, kstats, misc
	genome = "GENOME: " + org_name;
	kstat = "KARYSTATS: ";
	for (int i = 0; i < kstats.size(); i++)
	{
		kstat += to_string(kstats[i]) + ",";
	}
	misc = "MISCELLANEOUS: " + notes;

	out << genome << endl;
	out << kstat << endl;
	out << misc << endl;
	
	//loop through all chromosomes, get locations of pams, write out pams
	for (int curr_chrom = 0; curr_chrom < sequences.size(); curr_chrom++)
	{
		out << chroms[curr_chrom] << endl;
		
		//get locations pertaining to pams found in the current chromosome
		running_cnt += seed_cnts[curr_chrom];

		while (true)
		{
			if (j >= uniques.size())
			{
				break;
			}
			if (uniques[j] >= running_cnt)
			{
				break;
			}
			temp.push_back(seed_locs[uniques[j]]);
			j++;
		}
		
		//sort locations based on absolute value
		sort(temp.begin(), temp.end(), abs_cmp);
		
		//store reverse complement of current chromosome
		comp = reverseComplement(sequences[curr_chrom]);

		//loop through locations in current chromosomes, extract sequence, calculate score, write out
		for (int i = 0; i < temp.size(); i++)
		{
			if (temp[i] > 0)
			{
				seq = sequences[curr_chrom].substr(temp[i] - seq_length, seq_length + pam_length);
				out << temp[i] << ',' << seq.substr(0, seq_length) << ',' << seq.substr(seq_length, pam_length) << ',' << score.calcScore(seq) << endl;
			}
			else
			{
				pos = comp.size() + temp[i] + 1;
				seq = comp.substr(pos - seq_length, seq_length + pam_length);
				out << temp[i] << ',' << seq.substr(0, seq_length) << ',' << seq.substr(seq_length, pam_length) << ',' << score.calcScore(seq) << endl;
			}
		}

		//clear temporary locations vector
		temp.clear();
		temp.shrink_to_fit();
	}

	//close file
	boost::iostreams::close(outbuf);

	//cleanup - clear out uniques vector
	uniques.clear();
	uniques.shrink_to_fit();
}

void Write::write_repeats(string& filename, vector<int> &repeats, vector<string> &sequences, vector<long> &seed_locs, vector<unsigned long> &compressed_seeds, vector<int> &seed_cnts, string &score_file, int &five_length, int &three_length, int &seed_length, int &pam_length, int &seq_length, string &on_target_data)
{
	//variables
	Scoring score(score_file, on_target_data);
	sqlite3 *db;
	char *zErrMsg = 0;
	string sc, sql, seq, locs, seed, scores, fives, threes, pams, cs;
	int rc;
	int i = 0;
	int running_cnt = 0;
	int curr_chrom = 0;
	int cnt = 0;
	long pos = 0;
	vector<string> comps(sequences.size());

	//open and setup DB
	rc = sqlite3_open(filename.c_str(), &db);
	rc = sqlite3_exec(db, "PRAGMA synchronous = OFF;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "PRAGMA journal_mode = OFF;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "PRAGMA locking_mode = EXCLUSIVE;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "DROP TABLE IF EXISTS repeats;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "VACUUM;", NULL, 0, &zErrMsg);
	sql = "CREATE TABLE repeats (seed TEXT PRIMARY KEY, chromosome TEXT, location TEXT, three TEXT, five TEXT, pam TEXT, score TEXT, count INT);";
	rc = sqlite3_exec(db, sql.c_str(), NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, 0, &zErrMsg);
	
	//get all reverse comps of chromosome sequences
	for (i = 0; i < sequences.size(); i++)
	{
		comps[i] = reverseComplement(sequences[i]);
	}
	i = 0;
	//loop through repeats indices (they should be grouped together, so loop until compressed seed value isnt equal anymore)
	while (i < repeats.size())
	{
		sql = "INSERT INTO repeats ('seed', 'chromosome', 'location', 'three', 'five', 'pam', 'score', 'count') VALUES (";
		locs = to_string(seed_locs[repeats[i]]);
		running_cnt = seed_cnts[0];
		for (int j = 0; j < seed_cnts.size(); j++)
		{
			if (repeats[i] < running_cnt)
			{
				curr_chrom = j;
				break;
			}
			running_cnt += seed_cnts[j + 1];
		}

		if (seed_locs[repeats[i]] > 0)
		{
			
			seq = sequences[curr_chrom].substr(seed_locs[repeats[i]] - seq_length, seq_length + pam_length);
		}
		else
		{
			pos = comps[curr_chrom].size() + seed_locs[repeats[i]] + 1;
			seq = comps[curr_chrom].substr(pos - seq_length, seq_length + pam_length);

		}
		cs = to_string(curr_chrom + 1);
		threes = seq.substr(five_length + seed_length, three_length);
		fives = seq.substr(0, five_length);
		pams = seq.substr(seq_length, pam_length);
		scores = to_string(int(score.calcScore(seq)));
		cnt = 1;
		seed = seq.substr(five_length, seed_length);

		while (true)
		{
			if (i >= repeats.size() - 1)
			{
				i++;
				break;
			}
			if (compressed_seeds[repeats[i]] != compressed_seeds[repeats[i + 1]])
			{
				i++;
				break;
			}
			cnt++;
			
			locs += "," + to_string(seed_locs[repeats[i + 1]]);
			
			running_cnt = seed_cnts[0];
			for (int j = 0; j < seed_cnts.size(); j++)
			{
				if (repeats[i + 1] < running_cnt)
				{
					curr_chrom = j;
					break;
				}
				running_cnt += seed_cnts[j + 1];
			}
			if (seed_locs[repeats[i + 1]] > 0)
			{
				seq = sequences[curr_chrom].substr(seed_locs[repeats[i + 1]] - seq_length, seq_length + pam_length);
			}
			else
			{
				pos = comps[curr_chrom].size() + seed_locs[repeats[i + 1]] + 1;
				seq = comps[curr_chrom].substr(pos - seq_length, seq_length + pam_length);
			}
			cs += "," + to_string(curr_chrom + 1);
			threes += "," + seq.substr(five_length + seed_length, three_length);
			fives += "," + seq.substr(0, five_length);
			pams += "," + seq.substr(seq_length, pam_length);
			scores += "," + to_string(int(score.calcScore(seq)));
			i++;
		}

		//build sql insert statement, execute to db file
		seed = "'" + seq.substr(five_length, seed_length) + "'";
		sql += seed + ",'" + cs + "','" + locs + "','" + threes + "','" + fives + "','" + pams + "','" + scores + "'," + to_string(cnt);
		sql += ");";
		rc = sqlite3_exec(db, sql.c_str(), NULL, 0, 0);
	}

	//end db transaction, close db file
	rc = sqlite3_exec(db, "END TRANSACTION;", NULL, 0, &zErrMsg);
	sqlite3_close(db);
}

void Write::write_uniques_dir(vector<int> &uniques, vector<string> &sequences, vector<long> &seed_locs, vector<int> &seed_cnts, vector<unsigned long> &kstats, string &org_name, string &filename, string &score_file, vector<string> &chroms, string &notes, int &pam_length, int &seq_length, string &on_target_data)
{
	//variables
	Scoring score(score_file, on_target_data);
	string comp, seq, genome, kstat, misc;
	int i = 0;
	vector<long> temp;
	int j = 0;
	long pos = 0;
	int running_cnt = 0;
	ofstream outputfile;
	outputfile.open(filename, ios_base::out | ios_base::binary);
	boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf;
	outbuf.push(boost::iostreams::gzip_compressor());
	outbuf.push(outputfile);
	//Convert streambuf to ostream
	ostream out(&outbuf);

	//output first details - org name, kstats, misc
	genome = "GENOME: " + org_name;
	kstat = "KARYSTATS: ";
	for (int i = 0; i < kstats.size(); i++)
	{
		kstat += to_string(kstats[i]) + ",";
	}
	misc = "MISCELLANEOUS: " + notes;

	out << genome << endl;
	out << kstat << endl;
	out << misc << endl;

	//loop through all chromosomes, get locations of pams, write out pams
	for (int curr_chrom = 0; curr_chrom < sequences.size(); curr_chrom++)
	{
		out << chroms[curr_chrom] << endl;

		//get locations pertaining to pams found in the current chromosome
		running_cnt += seed_cnts[curr_chrom];
		while (uniques[j] < running_cnt && j < uniques.size())
		{
			temp.push_back(seed_locs[uniques[j]]);
			j++;
		}

		//sort locations based on absolute value
		sort(temp.begin(), temp.end(), abs_cmp);

		//store reverse complement of current chromosome
		comp = reverseComplement(sequences[curr_chrom]);

		//loop through locations in current chromosomes, extract sequence, calculate score, write out
		for (int i = 0; i < temp.size(); i++)
		{
			if (temp[i] > 0)
			{
				seq = sequences[curr_chrom].substr(temp[i] - 1, seq_length + pam_length);

				out << temp[i] + pam_length << ',' << seq.substr(pam_length, seq_length) << ',' << seq.substr(0, pam_length) << ',' << score.calcScore(seq) << endl;
			}
			else
			{
				pos = comp.size() + temp[i];
				seq = comp.substr(pos, seq_length + pam_length);
				out << temp[i] + pam_length << ',' << seq.substr(pam_length, seq_length) << ',' << seq.substr(0, pam_length) << ',' << score.calcScore(seq) << endl;
			}
		}

		//clear temporary locations vector
		temp.clear();
		temp.shrink_to_fit();
	}

	//close file
	boost::iostreams::close(outbuf);

	//cleanup - clear out uniques vector
	uniques.clear();
	uniques.shrink_to_fit();
}

void Write::write_repeats_dir(string& filename, vector<int> &repeats, vector<string> &sequences, vector<long> &seed_locs, vector<unsigned long> &compressed_seeds, vector<int> &seed_cnts, string &score_file, int &five_length, int &three_length, int &seed_length, int &pam_length, int &seq_length, string &on_target_data)
{
	//variables
	Scoring score(score_file, on_target_data);
	sqlite3 *db;
	char *zErrMsg = 0;
	string sc, sql, seq, locs, seed, scores, fives, threes, pams, cs;
	int rc;
	int i = 0;
	int running_cnt = 0;
	int curr_chrom = 0;
	long pos = 0;
	int cnt = 0;
	vector<string> comps(sequences.size());

	//open and setup DB
	rc = sqlite3_open(filename.c_str(), &db);
	rc = sqlite3_exec(db, "PRAGMA synchronous = OFF;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "PRAGMA journal_mode = OFF;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "PRAGMA locking_mode = EXCLUSIVE;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "DROP TABLE IF EXISTS repeats;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "VACUUM;", NULL, 0, &zErrMsg);
	sql = "CREATE TABLE repeats (seed TEXT PRIMARY KEY, chromosome TEXT, location TEXT, three TEXT, five TEXT, pam TEXT, score TEXT, count INT);";
	rc = sqlite3_exec(db, sql.c_str(), NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, 0, &zErrMsg);

	//get all reverse comps of chromosome sequences
	for (i = 0; i < sequences.size(); i++)
	{
		comps[i] = reverseComplement(sequences[i]);
	}
	i = 0;

	//loop through repeats indices (they should be grouped together, so loop until compressed seed value isnt equal anymore)
	while (i < repeats.size())
	{
		sql = "INSERT INTO repeats ('seed', 'chromosome', 'location', 'three', 'five', 'pam', 'score', 'count') VALUES (";
		locs = to_string(seed_locs[repeats[i]] + pam_length);
		running_cnt = seed_cnts[0];
		for (int j = 0; j < seed_cnts.size(); j++)
		{
			if (repeats[i] < running_cnt)
			{
				curr_chrom = j;
				break;
			}
			running_cnt += seed_cnts[j + 1];
		}
		if (seed_locs[repeats[i]] > 0)
		{
			seq = sequences[curr_chrom].substr(seed_locs[repeats[i]] - 1, seq_length + pam_length);
		}
		else
		{
			pos = comps[curr_chrom].size() + seed_locs[repeats[i]];
			seq = comps[curr_chrom].substr(pos, seq_length + pam_length);

		}
		cs = to_string(curr_chrom + 1);
		threes = seq.substr(pam_length + five_length + seed_length, three_length);
		fives = seq.substr(pam_length, five_length);
		pams = seq.substr(0, pam_length);
		scores = to_string(int(score.calcScore(seq)));
		cnt = 1;
		while (true)
		{
			if (compressed_seeds[repeats[i]] != compressed_seeds[repeats[i + 1]] || i >= repeats.size())
			{
				i++;
				break;
			}
			cnt++;
			locs += "," + to_string(seed_locs[repeats[i + 1]] + pam_length);
			running_cnt = seed_cnts[0];
			for (int j = 0; j < seed_cnts.size(); j++)
			{
				if (repeats[i + 1] < running_cnt)
				{
					curr_chrom = j;
					break;
				}
				running_cnt += seed_cnts[j + 1];
			}
			if (seed_locs[repeats[i + 1]] > 0)
			{
				seq = sequences[curr_chrom].substr(seed_locs[repeats[i + 1]] - 1, seq_length + pam_length);
			}
			else
			{
				pos = comps[curr_chrom].size() + seed_locs[repeats[i + 1]];
				seq = comps[curr_chrom].substr(pos, seq_length + pam_length);
			}
			cs += "," + to_string(curr_chrom + 1);
			threes += "," + seq.substr(pam_length + five_length + seed_length, three_length);
			fives += "," + seq.substr(pam_length, five_length);
			pams += "," + seq.substr(0, pam_length);
			scores += "," + to_string(int(score.calcScore(seq)));
			i++;
		}

		//build sql insert statement, execute to db file
		seed = "'" + seq.substr(pam_length + five_length, seed_length) + "'";
		sql += seed + ",'" + cs + "','" + locs + "','" + threes + "','" + fives + "','" + pams + "','" + scores + "'," + to_string(cnt);
		sql += ");";
		
		rc = sqlite3_exec(db, sql.c_str(), NULL, 0, 0);
	}

	//end db transaction, close db file
	rc = sqlite3_exec(db, "END TRANSACTION;", NULL, 0, &zErrMsg);
	sqlite3_close(db);
}

unsigned long Write::compressSeq(string &s) {
	unsigned long compseq = 0;
	for (int i = 0; i < s.size(); i++)
	{
		compseq += convertCharBase4(s[i])*pow(4, i); //multiplying by power-4 converts to base10
	}
	return compseq; //base10 version of sequence string
}

int Write::convertCharBase4(char &c) {
	switch (c) {
	case 'A': return 0;
	case 'T': return 1;
	case 'C': return 2;
	case 'G': return 3;
	default: return 0;
	}
}

string Write::reverseComplement(string &str)
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

/*
 * SeqFilter.cpp
 *
 *  Created on: 12 Sep 2017
 *      Author: simon
 *
 *      ---
 *
 *      Main functions for performing filtering
 *
 *      //
 *
 *      Current idea list:
 *       * Run together removed regions, for example XXXARGNDEXXX -> XXXXXXXXXXXX so frameshifts are better caught
 *       * Beginning and end regions are too generous -- user should define run_inside
 *       * Need a full set of options and the beginning of documentation
 *       * Store bounds for particular sequences in the pairHMM so custom bounds are internalised and multiple calls are not needed.
 */


#include <algorithm>
#include "Sequence.h"
#include <iomanip>

using namespace::std;

// Global variables. Ugly, but easiest quick solution
vector <CSequence> *testData = NULL;
vector <CSequence> *refData = NULL;
vector <vector <string> > seqBlock;		// The blocks of divvied sequences to sample from
string gapChar = "*-X?";

struct SScore {
	int TP = 0;				// Number of true pairs
	int FN = 0;				// Number of pairs not captured in true alignment
	int totalRef = 0;		// Number pairs in reference alignment
	int totalTest = 0;		// Total number of pairs in test alignment
	SScore operator+=(const SScore &S) {
		TP += S.TP;
		FN += S.FN;
		totalRef += S.totalRef;
		totalTest += S.totalTest;
		return *this;
	}
};

SScore ComparePairs(tuple<string, string> s1, tuple<string, string> s2); // Compare sequences s1/s2 with <test,ref> for each
string RemoveGaps(string &seq);
tuple<vector<int>,vector<int>> MapPositions(string &x, string &y);						// Maps x to y, with -1 for cases where x doesn't occur in y
vector <tuple<int,int>> MakePairs(vector<int> &x,vector<int> &y);
int Intersection(vector <tuple<int,int> > &x, vector <tuple <int,int> > &y);

int main(int argc, char * argv[]) {
	SScore score;

	if(argc < 3) {
		cout << "\nUsage: msascorer testMSA referenceMSA \n\n";
		exit(-1);
	}
	// Options
	string testFile = argv[1];
	string refFile = argv[2];

	// Read data and do some checking
	testData = FASTAReader(testFile); // Reads the sequences
	refData = FASTAReader(refFile);
	if(testData->size() != refData->size()) {  cout << "\nError: test and reference MSAs have different number of sequences"; }
	// Check the test
	assert(!testData->empty());
	for(auto & s : *testData) {
		if(s.length() != testData->at(0).length()) { cout << "\nSequences of uneven length in test MSA file\n"; exit(-1); }
	}
	// Check the reference
	assert(!refData->empty());
	for(auto & s : *refData) {
		if(s.length() != refData->at(0).length()) { cout << "\nSequences of uneven length in reference MSA file\n"; exit(-1); }
	}
	if(refData->at(0).length() > testData->at(0).length()) {
		cout << "\nError: reference MSA is long than test MSA\n"; exit(-1);
	}
	cout << "#Comparing " << testFile << " (seq:" << testData->size() << ";l=" << testData->at(0).length();
	cout << ") => REF " << refFile << "(seq:" << refData->size() << ";l=" << refData->at(0).length() << ")";

	// Do all against all comparison
	for(int i = 0; i < testData->size(); i++) {
		for(int j = i+1; j < testData->size(); j++) {
			score += ComparePairs(tuple<string,string>(testData->at(i).Seq(),refData->at(i).Seq()),tuple<string,string>(testData->at(j).Seq(),refData->at(j).Seq()) );

		}
	}
	int width = 15;
	cout << left << "\n" << setw(width)<< "#TruePos"<< setw(width) <<"FalseNeg"<< setw(width) <<"totalRef";
	cout << "\n" << setw(width) << score.TP << setw(width) << score.FN << setw(width) << score.totalRef;
	cout << "\n";
	return 0;
}

// Compare the pairs x and y <0: test, 1: ref> and return the score
SScore ComparePairs(tuple<string, string> seq1, tuple<string, string> seq2) {
	SScore retScore;
	// The pair of sequences referenced by position in y ; -1 used to denote gap or something not counted in x
	tuple<vector<int>,vector<int>> s1_int = MapPositions(get<0>(seq1),get<1>(seq1));
	tuple<vector<int>,vector<int>> s2_int = MapPositions(get<0>(seq2),get<1>(seq2));
	vector <tuple<int,int> > testPairs = MakePairs(get<0>(s1_int),get<0>(s2_int));
	vector <tuple<int,int> > refPairs = MakePairs(get<1>(s1_int),get<1>(s2_int));
	// Transfer information to the score structure
	retScore.totalTest = testPairs.size();
	retScore.totalRef = refPairs.size();
	retScore.TP = Intersection(testPairs,refPairs);
	retScore.FN = retScore.totalRef - retScore.TP;
//	cout << retScore << "\n" << flush;
	assert(retScore.FN >= 0);
	return retScore;
}

tuple<vector<int>,vector<int>> MapPositions(string &x, string &y) {
	string x_clean = RemoveGaps(x);
	string y_clean = RemoveGaps(y);
	vector <int> x_clean_int(x_clean.size(),-1), y_clean_int(y_clean.size(),-1); // The raw sequences unaligned
	vector <int> x_int(x.size(),-1), y_int(y.size(),-1);		// The return
	int start_x = (int) x_clean.find(y_clean);	// Get the position recast as int
	if(start_x == string::npos) { cout << "\nError: reference sequence is not a valid subset of the test sequence\ntest: " << x << "\nref:  " << y; exit(-1); }
	// Create the map
	for(int i = 0; i < y_clean_int.size(); i++) {
		x_clean_int[i + start_x] = y_clean_int[i] = i;
	}
	// Now put it on the alignments
	int pos = 0;
	for(int i = 0; i < x_int.size(); i++) {
		if(IsGap(x[i])) { continue; }
		x_int[i] = x_clean_int[pos++];
	}
	assert(pos == x_clean_int.size());
	pos = 0;
	for(int i = 0; i < y_int.size(); i++) {
		if(IsGap(y[i])) { continue; }
		y_int[i] = y_clean_int[pos++];
	}
	assert(pos == y_clean_int.size());
	return tuple<vector<int>,vector<int>>(x_int,y_int);
}

// Returns the pairs in strict ordering
vector <tuple<int,int>> MakePairs(vector<int> &x,vector<int> &y) {
	vector <tuple<int,int> > retPairs;
	assert(x.size() == y.size());
	for(int i = 0 ; i < x.size(); i++) {
		if(x[i] == -1 || y[i] == -1) { continue; }
		retPairs.push_back(tuple<int,int>(x[i],y[i]));
	}
	return retPairs;
}

int Intersection(vector <tuple<int,int> > &x, vector <tuple <int,int> > &y) {
	int lastMatch = 0, intersect = 0;
	for(int i = 0; i < x.size(); i++) {
		for(int j = lastMatch; j < y.size(); j++) {
			if(get<0>(y[j]) > get<0>(x[i])) { break; }
			if(get<1>(y[j]) > get<1>(y[j])) { break; }
			if(x[i] == y[i]) {
				intersect ++; lastMatch = j; break;
			}
		}
	}
	return intersect;
}

string RemoveGaps(string &seq) {
	stringstream retSeq;
	for(int i = 0 ; i < seq.size(); i++) {
		if(!IsGap(seq[i])) {
			retSeq << seq[i];
		}
	}
	return retSeq.str();
}

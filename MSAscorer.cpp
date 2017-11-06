/*
 * MSAscorer
 *
 *  Created on: 2 Nov 2017
 *      Author: Simon Whelan
 *      ---
 *	Method for comparing a test MSA to a reference MSA, that may contain only a subset of the testMSA sequence for each OTU
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
	int FP = 0;
	int FN = 0;				// Number of pairs not captured in true alignment
	int totalRef = 0;		// Number pairs in reference alignment
	int totalTest = 0;		// Total number of pairs in test alignment
	SScore operator+=(const SScore &S) {
		TP += S.TP;
		FP += S.FP;
		FN += S.FN;
		totalRef += S.totalRef;
		totalTest += S.totalTest;
		return *this;
	}
};

SScore ComparePairs(tuple<string, string> s1, tuple<string, string> s2); // Compare sequences s1/s2 with <test,ref> for each
string RemoveGaps(string &seq);
tuple<vector<int>,vector<int>> MapPositions(string &x, string &y, string &z); // Maps x to y, with -1 for cases where x doesn't occur in y; uses the second reference sequence z to ignore gaps
vector <tuple<int,int>> MakePairs(vector<int> &x,vector<int> &y);
int CountTP(vector <tuple<int,int> > &x);

int main(int argc, char * argv[]) {
	SScore score;

	if(argc > 1) {
		if(strcmp(argv[1], "-h") == 0) {
			cout << "\n===================================================================";
			cout << "\n\tMSAscorer : written by Simon Whelan";
			cout << "\n===================================================================";
			cout << "\n\nUsage: msascorer TestMSA RefMSA";
			cout << "\n\nTestMSA.fas / RefMSA.fas can be in FASTA/MSF/Phylip/Interleaved format";
			cout << "\n\nResults will look like this:\n";
			cout << "\n#Comparing TestMSA.fas (seq:4;l=112) => REF RefMSA.fas(seq:4;l=78)";
			cout << "\n#TruePos       FalsePos       FalseNeg       totalRef ";
			cout << "\n346            59             110            456";
			cout << "\n\nThe first line (commented with #) is a header for confirming input";
			cout << "\nThe second line references what returned values. They are defined as follows:";
			cout << "\n\tTruePos  :  Number of homologous pairs in both the RefMSA found in TestMSA";
			cout << "\n\tFalsePos :  Number of incorrect identified pairs of characters from the RefMSA found in TestMSA";
			cout << "\n\t(TruePos + FalsePos = Positives = Number of pairs in TestMSA with characters from only the RefMSA.)";
			cout << "\n\tFalseNeg :  Number of homologous pairs in RefMSA not found in TestMSA;";
			cout << "\n\ttotalRef :  Total number of homologous pairs in RefMSA (TruePos + FalseNeg)";
			cout << "\nNOTE: TrueNeg not included because it's the set of all pairs that could be obtained from RefMSA characters";
			cout << "\n\nThe final line are the values you need to use for computing available confusion matrix measures\n\n";
			exit(-1);
		}
	}
	if(argc < 3) {
		cout << "\n===================================================================";
		cout << "\n\tMSAscorer : written by Simon Whelan";
		cout << "\n===================================================================";
		cout << "\n\nStandard usage: msascorer TestMSA RefMSA\n\nType \"./msascorer -h\" for help\n\n";
		exit(-1);
	}
	// Options
	string testFile = argv[1];
	string refFile = argv[2];

	// Read data and do some checking
	testData = ReadSequences(testFile); // Reads the sequences
	refData = ReadSequences(refFile);
	sort(testData->begin(), testData->end(),[](CSequence &a, CSequence&b) { return a.Name() < b.Name(); });
	sort(refData->begin(), refData->end(),[](CSequence &a, CSequence&b) { return a.Name() < b.Name(); });
	if(testData->size() != refData->size()) {  cout << "\nError: test and reference MSAs have different number of sequences"; exit(-1); }
	for(int i = 0 ; i < testData->size(); i++) {
		if(testData->at(i).Name() != refData->at(i).Name()) { cout << "\nError: test ("<< testData->at(i).Name()<< ")and ref ("<< refData->at(i).Name() << ") have different names?\n"; exit(-1); }
	}
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
	cout << left << "\n" << setw(width)<< "#TruePos"<< setw(width) <<"FalsePos"<< setw(width) <<"FalseNeg"<< setw(width) <<"totalRef";
	cout << "\n" << setw(width) << score.TP << setw(width) << score.FP << setw(width) << score.FN << setw(width) << score.totalRef;
	cout << "\n";
	return 0;
}

// Compare the pairs x and y <0: test, 1: ref> and return the score
SScore ComparePairs(tuple<string, string> seq1, tuple<string, string> seq2) {
	SScore retScore;
	// Get labels for each of the sequences so we can identify and compare homology pairs; Detailed in the function
	tuple<vector<int>,vector<int>> s1_int = MapPositions(get<0>(seq1),get<1>(seq1),get<1>(seq2));
	tuple<vector<int>,vector<int>> s2_int = MapPositions(get<0>(seq2),get<1>(seq2),get<1>(seq1));
	// From these labels construct all of the pairs where both characters are present in the reference (defined as both labels >= 0)
	vector <tuple<int,int> > testPairs = MakePairs(get<0>(s1_int),get<0>(s2_int));
	vector <tuple<int,int> > refPairs = MakePairs(get<1>(s1_int),get<1>(s2_int));
	// Transfer information to the score structure
	retScore.totalTest = testPairs.size();
	retScore.totalRef = refPairs.size();
	retScore.TP = CountTP(testPairs);					// Intersection of totalTest and totalRef, but can be computed more simply
	retScore.FN = retScore.totalRef - retScore.TP;		// FNs are the
	retScore.FP = testPairs.size() - retScore.TP;		// FPs are the size of testPairs that are not TP
	assert(retScore.FN >= 0);
	return retScore;
}

// Map the reference sequence y onto the test sequence x. Uses the second reference sequence z to ensure gaps are handled correctly
// ---
// Each character in the pair of sequence x,y (returned in tuple<vector<int>,vector<int>>) is labelled as follows:
// -BIG_NUMBER : Not in reference MSA (only x; all of y will be labelled)
// -int : In the reference, but aligned to a gap
// int : In the reference and aligned to a real character
// The logic is commented into the code
tuple<vector<int>,vector<int>> MapPositions(string &x, string &y, string &z) {
	string x_clean = RemoveGaps(x);	// The raw sequences unaligned
	string y_clean = RemoveGaps(y);
	assert(y.size() == z.size());	// Check the aligned reference is correct
	vector <int> y_clean_int(y_clean.size(),-BIG_NUMBER);	// Storage of character mapping
	vector <int> x_int(x.size(),-BIG_NUMBER), y_int(y.size(),-BIG_NUMBER);		// The return vectors
	int start_x = (int) x_clean.find(y_clean);	// Get the real character position that y starts in clean_x, recast as int
	int end_x = start_x + y_clean.size();		// Get the positon in clean_x where y is expected to end
	if(start_x == string::npos) { cout << "\nError: reference sequence is not a valid subset of the test sequence\ntest: " << x << "\nref:  " << y; exit(-1); }
	// Create the map for reference (y) based on the other reference (z)
	int pos = 0;										// Sequence position number (observable character)
	for(int i = 0; i < y.size(); i++) {
		if(IsGap(z[i])) { y_int[i] = -i; }				// If aligned to a gap => label as -int
		else { y_int[i] = i; }							// If aligned to a character => label as int
		if(!IsGap(y[i])) { y_clean_int[pos++] = i; }    // Transfer that number to y_clean_int for mapping to x
	}
	// Now build the x_int
	pos = 0;											// Sequence position number (observable character)
	for(int i = 0; i < x.size(); i++) {
		if(IsGap(x[i])) { continue; }							// Don't care about gaps
		if(pos < start_x || pos >= end_x) { pos++; continue; }	// If the observable character in x is not also in y, just skip
		x_int[i] = y_clean_int[pos++ - start_x];				// Otherwise add the label from y_clean_int to x_int
	}
	return tuple<vector<int>,vector<int>>(x_int,y_int);
}

// Returns the pairs in strict ordering
vector <tuple<int,int>> MakePairs(vector<int> &x,vector<int> &y) {
	vector <tuple<int,int> > retPairs;
	assert(x.size() == y.size());
	for(int i = 0 ; i < x.size(); i++) {
		if(x[i] < 0 || y[i] < 0) { continue; }
		retPairs.push_back(tuple<int,int>(x[i],y[i]));
	}
	return retPairs;
}

// Simply count the number where the tuples match
int CountTP(vector <tuple<int,int> > &x) {
	int count = 0;
	for(auto &p : x) { if(get<0>(p) == get<1>(p)) { count++; } }
	return count;
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

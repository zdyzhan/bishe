#include<iostream>
#include<fstream>
#include"PCMBuilder.h"
#include"TimeUtility.h"
#include"MQO.h"
#include"SQP.h"
#include"InputCommandLineParser.h"
#include"QueryGraphIsomorphismTest.h"
#include"GlobalConstant.h"
#include"FrequentPatternGroup.h"

using namespace std;


std::vector<AdjacenceListsGRAPH> dataGraphVector, queryGraphVector;


void loadDataGraphs(string dataGraphFileName);
void loadQueryGraphs(string queryGraphFileName);

int main(int argc, char* argv[]) {

	std::ofstream resultFile = std::ofstream("C:/Users/s2813995/Dropbox/HectorResearch/PaperProjects/MQP_Center/MQOTestData/test.result");

	loadQueryGraphs("C:/Users/s2813995/Dropbox/HectorResearch/PaperProjects/MQP_Center/MQOTestData/test5.igraph");


	/*
	 * PCM building approaches tests
	 */
	FrequentPatternGroup grouper;
	std::vector<std::vector<int>> groups;

	grouper.frequentPatternGroup(queryGraphVector, groups);


	system("pause");

	return 0;
}

void loadDataGraphs(string dataGraphFileName) {
	std::ifstream dataGraphFile = std::ifstream(dataGraphFileName);
	if (!dataGraphFile) {
		cout << "The data file '" << dataGraphFileName << "' doesn't exist." << endl;
	}
	AdjacenceListsGRAPH_IO::loadGraphFromFile(dataGraphFile, dataGraphVector);

	for (vector<AdjacenceListsGRAPH>::iterator dataGraphIterator = dataGraphVector.begin(); dataGraphIterator != dataGraphVector.end(); dataGraphIterator++) {
		dataGraphIterator->buildLabelVertexList();
		dataGraphIterator->buildVertexLabelVertexList();
	}

	cout << "0. Load data graphs and build data graph index done. " << dataGraphFileName << endl;
}

void loadQueryGraphs(string queryGraphFileName) {

	std::ifstream queryGraphFile = std::ifstream(queryGraphFileName);
	if (!queryGraphFile) {
		cout << "The query file '" << queryGraphFileName << "' doesn't exist." << endl;
	}
	AdjacenceListsGRAPH_IO::loadGraphFromFile(queryGraphFile, queryGraphVector);


	int queryGraphIndex = 0;
	for (vector<AdjacenceListsGRAPH>::iterator queryGraphIterator = queryGraphVector.begin(); queryGraphIterator != queryGraphVector.end(); queryGraphIterator++) {
		queryGraphIterator->graphId = queryGraphIndex++;
		/*
		* build inverted indexes for each query graph
		*/
		queryGraphIterator->buildDFSTraversalOrder();
		queryGraphIterator->buildVertexLabelVertexList();
		queryGraphIterator->buildLabelVertexList();

		queryGraphIterator->buildVertexLabelEdgeList(); // used by mcs computation
		queryGraphIterator->buildTLSequence(); // used by pcm computation

	}

	cout << "1. Load query graphs and build query graph index done. " << queryGraphFileName << endl;
}


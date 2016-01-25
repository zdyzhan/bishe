#ifndef TLS_GRAPH_H
#define TLS_GRAPH_H

#include"AdjacenceListsGraph.h"
#include"PCMBuilder.h"
#include<vector>
#include<map>

using namespace std;

class TLSGraph{
private:
	int ** tlsGraphMatrix;

	std::vector<AdjacenceListsGRAPH> * queryGraphVector;

	/*
	 * For each TLSequence, we record all the query graphs that have this TLSequence,
	 * First int int the map is the query graph id, second int is how many of this TLSequece inside this query graph
	 */
	std::map<TLSequence, std::map<int, int>> TLSInvertedIndex;

	std::ofstream * resultFile;

public:
	TLSGraph(int ** pTlsGraphMatrix, std::vector<AdjacenceListsGRAPH> * pQueryGraphVector, std::ofstream * pResultFile);
	~TLSGraph();

	/*
	 * Build a TLS graph, each matrix element is the number of common TLS of two graphs
	 */
	void buildTLSGraph();

	/*
	 * Apply the relative similarity ratio, set it as 1 when relatively similar, otherwise 0
	 */
	void applyRelativeSimilarityRatio(std::vector<PCM_Node> * pPatternContainmentMap);
	
	void showTLSGraphMatrix();
	
	void showTLSInvertedIndex();
};


#endif
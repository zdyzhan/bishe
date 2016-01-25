#ifndef ULLMANN_H
#define ULLMANN_H

#include"TraversalAlgorithm.h"
#include"AdjacenceListsGraph.h"
#include"ComboLinkedLists.h"
#include<vector>
#include<string>
#include<map>


class Ullmann {

private:

	AdjacenceListsGRAPH * dataGraph;
	std::vector<int> * numberOfEmbeddings;
	ComboLinkedLists * comboLinkedLists;

private:

	AdjacenceListsGRAPH * queryGraph;
	int * partialEmbedding; /// mapping M: V(q) -> V(g), the arrayId is the query vertex id
	bool needFormatCache;

private:

	std::map<int, vector<int>* > candidates;
	std::map<int,int> inverseEmbedding; ///inverse mapping W: V(g) -> V(q)
	std::vector<AdjacenceListsGRAPH::Vertex *> * DFSTraversalOrder;



public:
	Ullmann();

	Ullmann(AdjacenceListsGRAPH * pDataGraph,  std::vector<int> * pNumberOfEmbeddings, ComboLinkedLists * pComboLinkedLists);
	
	~Ullmann();	

	void setParameters(AdjacenceListsGRAPH * pQueryGraph, int * pPartialEmbedding, bool pNeedToSaveCache);

	void execute();

private:

	void clean();

	void filterCandidates();

	void subgraphSearch();

	AdjacenceListsGRAPH::Vertex nextQueryVertex();

	bool refineCandidates(AdjacenceListsGRAPH::Vertex & u, const int & v);

	bool isJoinable(int u, int v);

	void updateState(int u, int v);

	void restoreState(int u, int v);

	bool degreeFilter(int u, int v);



	//@common
	/* Common Utility Function  */
	void showEmbedding();

};



#endif /* ULLMANN */

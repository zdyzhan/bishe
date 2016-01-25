#ifndef VF_2_H
#define VF_2_H

#include"TraversalAlgorithm.h"
#include<vector>
#include<stack>
#include<string>
#include<map>
#include"ComboLinkedLists.h"

using namespace std;

class VF2 {

private:


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
	std::map<int, int> inverseEmbedding; ///inverse mapping W: V(g) -> V(q)
	std::vector<AdjacenceListsGRAPH::Vertex *> * DFSTraversalOrder;

private:

	//@Unique
	std::map<int, int> Cg;
	std::map<int,std::pair<int,int>> Cq_Mq; // store the |Cq| and |adj(u)/Cq/Mq| 
	int numberOfCgAdj;
	int numberOfCg_Mg_Adj;

public:

	VF2();

	VF2(AdjacenceListsGRAPH * pDataGraph, std::vector<int> * pNumberOfEmbeddings, ComboLinkedLists * pComboLinkedLists);

	~VF2();

	void setParameters(AdjacenceListsGRAPH * pQueryGraph, int * pPartialEmbedding, bool pNeedToSaveCache);

	void execute();

private:

	/* Subgraph Isomorphism Framework Function */
	void clean(); // clean function for each indicidual subgraph isomorphism (1 query , 1 data), will be executed in the beginning of genericQueryProc

	void filterCandidates();

	void subgraphSearch();

	AdjacenceListsGRAPH::Vertex nextQueryVertex();

	bool refineCandidates(AdjacenceListsGRAPH::Vertex & u,const int & v);

	bool isJoinable(int u, int v);

	void updateState(int u, int v);

	void restoreState(int u, int v);

	bool degreeFilter(int u, int v);


private:

	//@Unique VF2
	bool isFlexible(AdjacenceListsGRAPH::Vertex & u, int & v);   // check Cq and Cg requirements

	void preCalculateCq();

	void updateCg(int & v);

	void restoreCg(int & v);

	//@common
	/* Common Utility Function  */

	void showEmbedding();

};



#endif /* VF2 */

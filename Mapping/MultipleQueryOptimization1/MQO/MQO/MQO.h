#ifndef MQO_DEF_H
#define MQO_DEF_H

#include"AdjacenceListsGraph.h"
#include"AdjacenceListsGRAPH_IO.h"
#include"TLSGraph.h"
#include"MQO.h"
#include"PCMBuilder.h"
#include"ComboLinkedLists.h"
#include<vector>
#include<iostream>
#include"RecurParEmbConstruct.h"
#include"GlobalConstant.h"

class MQO{

private:

	/*
	 * service functions
	 */
	TLSGraph * tlsGraph;

	PCMBuilder * pcmBuilder;

	RecurParEmbConstruct * recurParEmbConstruct;

	std::ofstream * resultFile;

private:

	/*
	 * datasets
	 */
	std::vector<AdjacenceListsGRAPH> * queryGraphVector;
	
	AdjacenceListsGRAPH * dataGraph;

	AdjacenceListsGRAPH * queryGraph;

	int queryGraphId;

private: 
	/*
	 * used to contain the containment relationships of the queries and newly generated graphs
	 */
	std::vector<PCM_Node> patternContainmentMap;

	/*
	 * used to store the Triple Label Sequence Matrix
	 */
	int ** tlsGraphMatrix;

	/*
	 * used to store the all the newly generated graphs
	 */
	std::vector<AdjacenceListsGRAPH> newlyGeneratedGraphVector;


	/*
	 * Used to save the number of founded embeddings
	 */
	std::vector<int> numberOfEmbeddings;

	/*
	 * used to save the linked list of cached embeddings
	 * points to the head of each linked list
	 */
	std::vector<std::vector<EmbeddingNode *>> comboLinkedLists;

	/*
	 * store the processing(joining) order of a query's parents
	 */
	vector<int> joiningOrder;

	/*
	* store the uncovered edges by its parents, these edges must be validate seperately
	*/
	std::set<std::pair<int, int>> uncoveredEdges;


	/*
	 * statistic utilities
	 */
	double timeClapse;

public:

	MQO();
	MQO(std::vector<AdjacenceListsGRAPH> * dataGraphVector, std::vector<AdjacenceListsGRAPH> * pQueryGraphVector, std::ofstream * pResultFile);
	~MQO();

	void buildPCM();

	void orederedQueryProcessing();

private:

	int getNextGraphId(vector<int> & candidates);

	void changeParentsWeight(int graphId);

	/*
	* Given a PCM which is a transitive reduction graph, we need an order for processing each PCM node which points to a graph
	*/
	void DFSTopological(int graphId);

	/*
	 * For each PCM node containing multiple parents, we need to compute an order to sequentially join the cached embeddings of its parents
	 * if return false, that means one of its parent doesn't have any embedding. this query graph cannot have any embedding either. We don't 
	 * need to compute it further.
	 */
	bool computeJoiningOrder();

	void showEmbedding(int graphId);

};



#endif 
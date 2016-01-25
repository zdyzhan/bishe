#pragma once
#ifndef RECUR_PARTIAL_EMBEDDING_CONSTUCT_H
#define RECUR_PARTIAL_EMBEDDING_CONSTUCT_H

#include<vector>
#include<map>
#include<set>
#include"AdjacenceListsGraph.h"
#include"ComboLinkedLists.h"
#include"PCMBuilder.h"
#include"Ullmann.h"
#include"VF2.h"
#include"TurboIso.h"
#include"TurboIsoBoosted.h"
#include"TurboIso_MQO.h"
#include"TurboIsoBoosted_MQO.h"

class RecurParEmbConstruct{

private:

	AdjacenceListsGRAPH * dataGraph;

	std::vector<AdjacenceListsGRAPH> * queryGraphVector;

	std::vector<AdjacenceListsGRAPH> * newlyGeneratedGraphVector;

	std::vector<int> * numberOfEmbeddings;

	std::vector<std::vector<EmbeddingNode *>> * comboLinkedLists;

	std::vector<PCM_Node> * patternContainmentMap;

private:
	
	AdjacenceListsGRAPH * queryGraph;

	std::vector<int> * joiningOrder;

	std::set<std::pair<int, int>> * uncoveredEdges;

private :


	/*
	* store the partial embedding generated from the cached parents
	* we use an array to improve the efficiency
	*/
	int * partialEmbedding;

	/*
	* To save the set of mapped data vertices, used to prevent duplication uses
	*/
	std::set<int> mappedDataVertexSet;

private:

	ComboLinkedLists * cachedResultLists;
	
	Ullmann  * ullmann;

	VF2 * vf2;

	TurboIsoBoosted * turboIsoBoosted;

	TurboIso * turboIso;

	TurboIsoBoosted_MQO * turboIsoBoosted_MQO;

	TurboIso_MQO * turboIso_MQO;

public:
	RecurParEmbConstruct();
	RecurParEmbConstruct(AdjacenceListsGRAPH * pDataGraph, std::vector<AdjacenceListsGRAPH> * pQueryGraphVector, std::vector<AdjacenceListsGRAPH> * pNewlyGeneratedGraphVector, std::vector<int> * pNumberOfEmbeddings, std::vector<std::vector<EmbeddingNode *>> * pComboLinkedLists, std::vector<PCM_Node> * pPatternContainmentMap);
	~RecurParEmbConstruct();

	void processing(AdjacenceListsGRAPH * pQueryGraph, std::vector<int> * pJoiningOrder, std::set<std::pair<int, int>> * pUncoveredEdges);

private:
	/*
	* Use the cached embeddings to compute a partial embedding for a given query graph which has a node in the PCM graph
	*/
	void recursiveComputePartialEmbedding(int matchedToOrderIndex, int selfJoinOrderIndex);

	/*
	* Enumerate the embeddings from the linked list cached embeddings
	*/
	void enumberateEmbeddingLinkedList(int comboVertexOrder, int matchedToOrderIndex, int selfJoinOrderIndex, std::vector<ComboNODE> * parentComboGraph, std::vector<int> * parentDFSComboVertexOrder, std::vector<int> * vertexMappingList, std::vector<EmbeddingNode *> * pQarentComboEmbedding);

	/*
	* Final verify this partial embedding by checking uncovered edges by its parents
	*/
	bool checkPartialEmbeddingUncoveredEdge();

	/*
	 * Call one of the subgraph isomorphism algorithms to further process the partial solution
	 */
	void runSubgraphIsomorphismSearch();

};

#endif 
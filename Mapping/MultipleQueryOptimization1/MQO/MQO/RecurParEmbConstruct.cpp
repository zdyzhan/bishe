#include"RecurParEmbConstruct.h"
#include"GlobalConstant.h"
#include"TimeUtility.h"

using namespace std;



RecurParEmbConstruct::RecurParEmbConstruct(){}

RecurParEmbConstruct::RecurParEmbConstruct(AdjacenceListsGRAPH * pDataGraph, std::vector<AdjacenceListsGRAPH>* pQueryGraphVector, std::vector<AdjacenceListsGRAPH>* pNewlyGeneratedGraphVector, std::vector<int>* pNumberOfEmbeddings, std::vector<std::vector<EmbeddingNode*>>* pComboLinkedLists, std::vector<PCM_Node>* pPatternContainmentMap)
{
	dataGraph = pDataGraph;
	queryGraphVector = pQueryGraphVector;
	newlyGeneratedGraphVector = pNewlyGeneratedGraphVector;
	numberOfEmbeddings = pNumberOfEmbeddings;
	comboLinkedLists = pComboLinkedLists;
	patternContainmentMap = pPatternContainmentMap;

	cachedResultLists = new ComboLinkedLists(comboLinkedLists);

	ullmann = new Ullmann(dataGraph, numberOfEmbeddings, cachedResultLists);

	vf2 = new VF2(dataGraph, numberOfEmbeddings, cachedResultLists);

	turboIsoBoosted = new TurboIsoBoosted(dataGraph, numberOfEmbeddings, cachedResultLists);

	turboIso = new TurboIso(dataGraph, numberOfEmbeddings, cachedResultLists);

	turboIsoBoosted_MQO = new TurboIsoBoosted_MQO(dataGraph, numberOfEmbeddings, cachedResultLists);

	turboIso_MQO = new TurboIso_MQO(dataGraph, numberOfEmbeddings, cachedResultLists);
}

RecurParEmbConstruct::~RecurParEmbConstruct()
{
	delete cachedResultLists;
	delete ullmann;

	mappedDataVertexSet.swap(set<int>());
}

void RecurParEmbConstruct::processing(AdjacenceListsGRAPH * pQueryGraph, std::vector<int> * pJoiningOrder, std::set<std::pair<int, int>> * pUncoveredEdges)
{
	queryGraph = pQueryGraph;
	joiningOrder = pJoiningOrder;
	uncoveredEdges = pUncoveredEdges;

	mappedDataVertexSet.swap(set<int>());


	partialEmbedding = new int[queryGraph->getNumberOfVertexes()];

	for (int i = 0; i < queryGraph->getNumberOfVertexes(); i++) {
		partialEmbedding[i] = -1;
	}
	

	recursiveComputePartialEmbedding(0, 0);

	delete[] partialEmbedding;

}


/*
* With the partial embeddings ready for this query graph, we are ready to process the joining process
* The joining process is also a recursively enumeration process
*/
void RecurParEmbConstruct::recursiveComputePartialEmbedding(int matchedToOrderIndex, int selfJoinOrderIndex) {
	if ((*numberOfEmbeddings)[queryGraph->graphId] >= GlobalConstant::G_SUFFICIENT_NUMBER_EMBEDDINGS || (*numberOfEmbeddings)[queryGraph->graphId]==-1) {
		return;
	}
	/*****************************************
	* Terminate conditions *
	*****************************************/
	if (mappedDataVertexSet.size() == queryGraph->getNumberOfVertexes()) {
		if (checkPartialEmbeddingUncoveredEdge()) {
			/*
			* We got a full embedding by only use the caches. Add it.
			*/
			(*numberOfEmbeddings)[queryGraph->graphId] ++;
			/*
			* We cache the results if and only if it has some children
			*/
			if ((*patternContainmentMap)[queryGraph->graphId].children.size() != 0) {
				cachedResultLists->addEmbeddingToCache(queryGraph, partialEmbedding);
			}
		}
		return;
	}
	if (matchedToOrderIndex == joiningOrder->size()) {
		/*
		* we iterate to the end of all its minimum query cover parents
		*/
		if (matchedToOrderIndex > 0 && mappedDataVertexSet.size() == 0) {
			/*
			* it has parents, however, we cannot get any partial embedding from its parent caches, that means, the query graph wouldn't have any embeddings
			*/
			return;
		}
		else {
			/*
			* @subgraph isomorphism. pass it
			* We only save the cache results if and only if this graph has pcm children
			*/
			if (checkPartialEmbeddingUncoveredEdge()) {
				// Call subgraph isomorphism algorithm
				runSubgraphIsomorphismSearch();
			}
			return;
		}
	}


	/*****************************************
	* Set the parent graph in this round
	*****************************************/
	int parentQueryGraphId = joiningOrder->at(matchedToOrderIndex);
	AdjacenceListsGRAPH * parentQueryGraph = NULL;
	if (parentQueryGraphId < queryGraphVector->size()) {
		parentQueryGraph = &(*queryGraphVector)[parentQueryGraphId];
	}
	else {
		parentQueryGraph = &(*newlyGeneratedGraphVector)[parentQueryGraphId - queryGraphVector->size()];
	}
	std::vector<int> * parentDFSComboVertexOrder = parentQueryGraph->getDFSComboVertexOrder();
	std::vector<ComboNODE> * parentComboGraph = parentQueryGraph->getComboGraph();

	/*
	* Find the mappings of this parent graph
	*/
	std::map<int, std::vector<std::vector<int>>>::iterator vertexMappingMapIterator = (*patternContainmentMap)[parentQueryGraphId].containmentRelationshipMappingLists.find(queryGraph->graphId);
	/*
	* The mapping from the parent graph to this graph: (0->? 1->? 2->?)
	*/
	std::vector<int> * vertexMappingList = &vertexMappingMapIterator->second[selfJoinOrderIndex];


	std::vector<EmbeddingNode *> parentComboEmbedding;
	for (int i = 0; i<parentComboGraph->size(); i++) {
		parentComboEmbedding.push_back(NULL);
	}

	/*****************************************
	* Prepare for Next recursieve query round
	*****************************************/
	if (selfJoinOrderIndex == vertexMappingMapIterator->second.size() - 1) {
		matchedToOrderIndex++;
		selfJoinOrderIndex = 0;
	}
	else {
		selfJoinOrderIndex++;
	}


	/*****************************************
	* enumerate and compose the cached embeddings
	*****************************************/
	EmbeddingNode * rootsMatchedVertices = (*comboLinkedLists)[parentQueryGraphId][0];

	int * updatedVertexId = new int[5];
	for (int i = 0; i<5; i++) {
		updatedVertexId[i] = -1;
	}

	while (rootsMatchedVertices != NULL) {
		
		if ((*numberOfEmbeddings)[queryGraph->graphId] >= GlobalConstant::G_SUFFICIENT_NUMBER_EMBEDDINGS || (*numberOfEmbeddings)[queryGraph->graphId]==-1) {
			return;
		}

		//1. Make sure no data vertices can be used multiple times AND the same data vertex to each query vertex
		bool isJoinable = true;
		for (int i = 0; i<(*parentComboGraph)[0].queryVertices.size(); i++) {

			int mappedQueryVertex = (*vertexMappingList)[(*parentComboGraph)[0].queryVertices[i]];

			if (partialEmbedding[mappedQueryVertex] != -1) {
				// same data vertex to each query vertex
				if (partialEmbedding[mappedQueryVertex] != rootsMatchedVertices->comboMappingVertices[i]) {
					isJoinable = false;
					break;
				}
			}
			else {
				if (mappedDataVertexSet.find(rootsMatchedVertices->comboMappingVertices[i]) == mappedDataVertexSet.end()) {
					// no data vertices can be used multiple times 
					partialEmbedding[mappedQueryVertex] = rootsMatchedVertices->comboMappingVertices[i];
					mappedDataVertexSet.insert(rootsMatchedVertices->comboMappingVertices[i]);

					updatedVertexId[i] = mappedQueryVertex;
				}
				else {
					isJoinable = false;
					break;
				}
			}
		}

		if (isJoinable) {
			//2. start iteratiing inner embeddings
			parentComboEmbedding[(*parentDFSComboVertexOrder)[0]] = rootsMatchedVertices;
			enumberateEmbeddingLinkedList(1, matchedToOrderIndex, selfJoinOrderIndex, parentComboGraph, parentDFSComboVertexOrder, vertexMappingList, &parentComboEmbedding);
		}

		//3. restore to initial status
		for (int i = 0; i<5; i++) {
			if (updatedVertexId[i] != -1) {

				mappedDataVertexSet.erase(partialEmbedding[updatedVertexId[i]]);

				partialEmbedding[updatedVertexId[i]] = -1;

				updatedVertexId[i] = -1;
			}
		}

		//4. next root trial
		rootsMatchedVertices = rootsMatchedVertices->adj;

	}
	
	//5. release memory
	delete[] updatedVertexId;
}





void RecurParEmbConstruct::enumberateEmbeddingLinkedList(int comboVertexOrder, int matchedToOrderIndex, int selfJoinOrderIndex, std::vector<ComboNODE> * parentComboGraph, std::vector<int> * parentDFSComboVertexOrder, std::vector<int> * vertexMappingList, std::vector<EmbeddingNode *> * parentComboEmbedding) {
	
	if ((*numberOfEmbeddings)[queryGraph->graphId] >= GlobalConstant::G_SUFFICIENT_NUMBER_EMBEDDINGS || (*numberOfEmbeddings)[queryGraph->graphId]==-1) {
		return;
	}

	if (comboVertexOrder == parentComboGraph->size() || mappedDataVertexSet.size() == queryGraph->getNumberOfVertexes()) {
		/*
		* Push to next parent graph
		*/
		recursiveComputePartialEmbedding(matchedToOrderIndex, selfJoinOrderIndex);
		return;
	}

	int thisComboVertexId = (*parentDFSComboVertexOrder)[comboVertexOrder];
	int parentComboVertexId = (*parentComboGraph)[thisComboVertexId].parent;

	std::map<int, std::set<EmbeddingNode *>>::iterator childrenMapIterator = ((*parentComboEmbedding)[parentComboVertexId])->neighboursList.find(thisComboVertexId);

	/*
	* iteration state helper variables
	*/
	int * updatedVertexId = new int[5];
	for (int i = 0; i<5; i++) {
		updatedVertexId[i] = -1;
	}
	int originalPartialEmbeddingMatchedSize = -1;

	for (std::set<EmbeddingNode *>::iterator childrenSetIterator = childrenMapIterator->second.begin(); childrenSetIterator != childrenMapIterator->second.end(); childrenSetIterator++) {
		
		if ((*numberOfEmbeddings)[queryGraph->graphId] >= GlobalConstant::G_SUFFICIENT_NUMBER_EMBEDDINGS || (*numberOfEmbeddings)[queryGraph->graphId]==-1) {
			return;
		}

		bool isJoinable = true;
		//1. Make sure the non-DFS edge has corresponding mapping
		for (set<int>::iterator neighbourIterator = (*parentComboGraph)[thisComboVertexId].neighbours.begin(); neighbourIterator != (*parentComboGraph)[thisComboVertexId].neighbours.end(); neighbourIterator++) {
			if (*neighbourIterator != parentComboVertexId && (*parentComboEmbedding)[*neighbourIterator] != NULL) {

				std::map<int, std::set<EmbeddingNode *>>::iterator nonDFSEdgeMappingIterator;
				if (thisComboVertexId < *neighbourIterator) {
					nonDFSEdgeMappingIterator = (*childrenSetIterator)->neighboursList.find(*neighbourIterator);
					if (nonDFSEdgeMappingIterator->second.find((*parentComboEmbedding)[*neighbourIterator]) == nonDFSEdgeMappingIterator->second.end()) {
						isJoinable = false;
						break;
					}
				}
				else {
					nonDFSEdgeMappingIterator = (*parentComboEmbedding)[*neighbourIterator]->neighboursList.find(thisComboVertexId);
					if (nonDFSEdgeMappingIterator->second.find(*childrenSetIterator) == nonDFSEdgeMappingIterator->second.end()) {
						isJoinable = false;
						break;
					}
				}
			}
		}

		if (isJoinable) {
			//2. Make sure no data vertices can be used multiple times AND the same data vertex to each query vertex
			for (int i = 0; i<(*parentComboGraph)[thisComboVertexId].queryVertices.size(); i++) {

				int mappedQueryVertex = (*vertexMappingList)[(*parentComboGraph)[thisComboVertexId].queryVertices[i]];

				if (partialEmbedding[mappedQueryVertex] != -1) {
					// same data vertex to each query vertex
					if (partialEmbedding[mappedQueryVertex] != (*childrenSetIterator)->comboMappingVertices[i]) {
						isJoinable = false;
						break;
					}
				}
				else {
					if (mappedDataVertexSet.find((*childrenSetIterator)->comboMappingVertices[i]) == mappedDataVertexSet.end()) {
						// no data vertices can be used multiple times 
						partialEmbedding[mappedQueryVertex] = (*childrenSetIterator)->comboMappingVertices[i];
						mappedDataVertexSet.insert((*childrenSetIterator)->comboMappingVertices[i]);
						updatedVertexId[i] = mappedQueryVertex;
					}
					else {
						isJoinable = false;
						break;
					}
				}
			}

			if (isJoinable) {
				//3. start next iteratiing inner embeddings
				(*parentComboEmbedding)[thisComboVertexId] = (*childrenSetIterator);
				enumberateEmbeddingLinkedList(comboVertexOrder + 1, matchedToOrderIndex, selfJoinOrderIndex, parentComboGraph, parentDFSComboVertexOrder, vertexMappingList, parentComboEmbedding);
			}


			//4. restore to initial status
			for (int i = 0; i<5; i++) {
				if (updatedVertexId[i] != -1) {
					mappedDataVertexSet.erase(partialEmbedding[updatedVertexId[i]]);

					partialEmbedding[updatedVertexId[i]] = -1;

					updatedVertexId[i] = -1;
				}
			}
		}
	}

	//5. reset the combo mapping to NULL after this inner iteration
	(*parentComboEmbedding)[thisComboVertexId] = NULL;
	delete[] updatedVertexId;
}


bool RecurParEmbConstruct::checkPartialEmbeddingUncoveredEdge() {
	for (set<std::pair<int, int>>::iterator uncoveredEdgeIterator = uncoveredEdges->begin(); uncoveredEdgeIterator != uncoveredEdges->end(); uncoveredEdgeIterator++) {
		if (partialEmbedding[uncoveredEdgeIterator->first] != -1 && partialEmbedding[uncoveredEdgeIterator->second] != -1) {
			if (!dataGraph->edge(partialEmbedding[uncoveredEdgeIterator->first], partialEmbedding[uncoveredEdgeIterator->second])) {
				return false;
			}
		}
	}
	return true;
}

void RecurParEmbConstruct::runSubgraphIsomorphismSearch()
{
	if ((*numberOfEmbeddings)[queryGraph->graphId] >= GlobalConstant::G_SUFFICIENT_NUMBER_EMBEDDINGS || (*numberOfEmbeddings)[queryGraph->graphId]==-1) {
		return;
	}

	switch (GlobalConstant::G_RUNNING_OPTION_INDEX) {
	case GlobalConstant::RUN_OP_ULLMANN:
		ullmann->setParameters(queryGraph, partialEmbedding, ((*patternContainmentMap)[queryGraph->graphId].children.size() != 0) ? true : false);
		ullmann->execute();
		break;
	case GlobalConstant::RUN_OP_VF2:
		vf2->setParameters(queryGraph, partialEmbedding, ((*patternContainmentMap)[queryGraph->graphId].children.size() != 0) ? true : false);
		vf2->execute();
		break;
	case GlobalConstant::RUN_OP_TURBOISO:
		if (mappedDataVertexSet.size() > 0) {
			turboIso_MQO->setParameters(queryGraph, partialEmbedding, ((*patternContainmentMap)[queryGraph->graphId].children.size() != 0) ? true : false);
			turboIso_MQO->execute();
		}
		else {	
			turboIso->setParameters(queryGraph, partialEmbedding, ((*patternContainmentMap)[queryGraph->graphId].children.size() != 0) ? true : false);
			turboIso->execute();
		}
		break;
	case GlobalConstant::RUN_OP_TURBO_ISO_BOOSTED:
		if (mappedDataVertexSet.size() > 0) {
			turboIsoBoosted_MQO->setParameters(queryGraph, partialEmbedding, ((*patternContainmentMap)[queryGraph->graphId].children.size() != 0) ? true : false);
			turboIsoBoosted_MQO->execute();
		}
		else {
			turboIsoBoosted->setParameters(queryGraph, partialEmbedding, ((*patternContainmentMap)[queryGraph->graphId].children.size() != 0) ? true : false);
			turboIsoBoosted->execute();
		}
		break;
	}
}

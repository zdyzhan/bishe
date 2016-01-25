#include"PartiLinkedLists.h"
#include<vector>


using namespace std;

/*
 *  1 means first < second
 *  0 means first == second
 *  -1 means first > second
 */
int compareEmbeddingNode(const PartiLinkedLists::EmbeddingNode& first, const PartiLinkedLists::EmbeddingNode& second) {
	for (int i = 0; i < 5; i++) {
		if (first.comboMappingVertices[i] < second.comboMappingVertices[i]) {
			return 1;
		}
		else if (first.comboMappingVertices[i] > second.comboMappingVertices[i]) {
			return -1;
		}
		else if (first.comboMappingVertices[i] == -1) {
			return 0;
		}
	}
	return 0;
}


void PartiLinkedLists::addEmbeddingToCache(AdjacenceListsGRAPH * queryGraph, std::vector<std::vector<EmbeddingNode *>> * cachedResults, int * partialEmbedding) {

	std::vector<ComboNODE> * comboGraph = queryGraph->getComboGraph();
	short int ** comboGraphEdgeMatrix = queryGraph->getComboGraphEdgeMatrix();
	std::vector<int> * DFSComboVertexOrder = queryGraph->getDFSComboVertexOrder();

	PartiLinkedLists::EmbeddingNode ** queryVertexPositions = new PartiLinkedLists::EmbeddingNode *[comboGraph->size()];

	for (int comboVertexIndex = 0; comboVertexIndex<comboGraph->size(); comboVertexIndex++) {

		if ((*cachedResults)[queryGraph->graphId][comboVertexIndex] == NULL) {
			/*
			 * Empty entry
			 */
			(*cachedResults)[queryGraph->graphId][comboVertexIndex] = new PartiLinkedLists::EmbeddingNode();
			for (int j = 0; j < (*comboGraph)[comboVertexIndex].queryVertices.size(); j++) {
				(*cachedResults)[queryGraph->graphId][comboVertexIndex]->comboMappingVertices[j] = partialEmbedding[(*comboGraph)[comboVertexIndex].queryVertices[j]];
			}
			queryVertexPositions[comboVertexIndex] = (*cachedResults)[queryGraph->graphId][comboVertexIndex];
			continue ;
		}


		PartiLinkedLists::EmbeddingNode embeddingNode;
		for (int j = 0; j < (*comboGraph)[comboVertexIndex].queryVertices.size(); j++) {
			embeddingNode.comboMappingVertices[j] = partialEmbedding[(*comboGraph)[comboVertexIndex].queryVertices[j]];
		}
		PartiLinkedLists::EmbeddingNode * insertPosition = (*cachedResults)[queryGraph->graphId][comboVertexIndex];
		PartiLinkedLists::EmbeddingNode * insertPositionParent = NULL;

		int compareValue;
		while (insertPosition != NULL) {
			compareValue = compareEmbeddingNode(embeddingNode, *insertPosition);
			if (compareValue == 1) {
				insertPositionParent = insertPosition;
				insertPosition = insertPosition->leftAdj;
			}
			else if (compareValue == -1) {
				insertPositionParent = insertPosition;
				insertPosition = insertPosition->rightAdj;
			}
			else {
				break;
			}
		}

		if (compareValue != 0) {
			insertPosition = new PartiLinkedLists::EmbeddingNode(embeddingNode);
			if (compareValue == 1) {
				insertPositionParent->leftAdj = insertPosition;
			}
			else {
				insertPositionParent->rightAdj = insertPosition;
			}
		}
		queryVertexPositions[comboVertexIndex] = insertPosition;
	}


	// handle the dfs edges
	for (std::vector<int>::iterator orderComboVertexIterator = DFSComboVertexOrder->begin(); orderComboVertexIterator != DFSComboVertexOrder->end(); orderComboVertexIterator++) {
		for (std::vector<int>::iterator childrenIterator = (*comboGraph)[*orderComboVertexIterator].children.begin(); childrenIterator != (*comboGraph)[*orderComboVertexIterator].children.end(); childrenIterator++) {

			int a = queryVertexPositions[*orderComboVertexIterator]->neighboursList.size();

			std::map<int, std::set<PartiLinkedLists::EmbeddingNode *>>::iterator neighbourLinkIterator = queryVertexPositions[*orderComboVertexIterator]->neighboursList.find(*childrenIterator);

			if (neighbourLinkIterator != queryVertexPositions[*orderComboVertexIterator]->neighboursList.end()) {
				neighbourLinkIterator->second.insert(queryVertexPositions[*childrenIterator]);
			}
			else {
				set<PartiLinkedLists::EmbeddingNode *> newList;
				newList.insert(queryVertexPositions[*childrenIterator]);
				queryVertexPositions[*orderComboVertexIterator]->neighboursList.insert(std::pair<int, set<PartiLinkedLists::EmbeddingNode *>>(*childrenIterator, newList));
			}
		}
	}

	// handle the non-spanning(not in the dfs route) tree edges
	for (int i = 0; i<comboGraph->size(); i++) {
		for (int j = i + 1; j<comboGraph->size(); j++) {
			if (comboGraphEdgeMatrix[i][j] && (*comboGraph)[i].parent != j && (*comboGraph)[j].parent != i) {

				std::map<int, std::set<PartiLinkedLists::EmbeddingNode *>>::iterator neighbourLinkIterator = queryVertexPositions[i]->neighboursList.find(j);
				if (neighbourLinkIterator != queryVertexPositions[i]->neighboursList.end()) {
					neighbourLinkIterator->second.insert(queryVertexPositions[j]);
				}
				else {
					set<PartiLinkedLists::EmbeddingNode *> newList;
					newList.insert(queryVertexPositions[j]);
					queryVertexPositions[i]->neighboursList.insert(std::pair<int, set<PartiLinkedLists::EmbeddingNode *>>(j, newList));
				}
			}
		}
	}

	delete[] queryVertexPositions;
}
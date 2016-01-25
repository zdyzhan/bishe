#include"ComboLinkedLists.h"
#include<vector>


using namespace std;

// Overload 
bool operator<(const EmbeddingNode& first, const EmbeddingNode& second) {
	for (int i = 0; i < 5; i++) {
		if (first.comboMappingVertices[i] < second.comboMappingVertices[i]) {
			return true;
		}
		else if (first.comboMappingVertices[i] > second.comboMappingVertices[i]) {
			return false;
		}
	}
	return false;
}

bool operator==(const EmbeddingNode& first, const EmbeddingNode& second) {
	bool result = true;
	for (int i = 0; i < 5; i++) {
		result == result && (first.comboMappingVertices[i] == second.comboMappingVertices[i]);
	}
	return result;
}


ComboLinkedLists::ComboLinkedLists(std::vector<std::vector<EmbeddingNode*>> * pComboLinkedLists)
{
	comboLinkedLists = pComboLinkedLists;
}

void ComboLinkedLists::addEmbeddingToCache(AdjacenceListsGRAPH * queryGraph, int * partialEmbedding){

	std::vector<ComboNODE> * comboGraph = queryGraph->getComboGraph();
	short int ** comboGraphEdgeMatrix = queryGraph->getComboGraphEdgeMatrix();
	std::vector<int> * DFSComboVertexOrder = queryGraph->getDFSComboVertexOrder();
	std::vector<EmbeddingNode *> queryVertexPositions;

	for (int i = 0; i<comboGraph->size(); i++) {

		EmbeddingNode embeddingNode;
		for (int j = 0; j < (*comboGraph)[i].queryVertices.size(); j++) {
			embeddingNode.comboMappingVertices[j] = partialEmbedding[(*comboGraph)[i].queryVertices[j]];
		}

		/*
		* insert the node in an increasing order.
		* For each embedding docmo, we first find its inserting position
		*/
		EmbeddingNode * insertPosition = (*comboLinkedLists)[queryGraph->graphId][i];
		EmbeddingNode * insertPositionParent = NULL;

		while (insertPosition != NULL && (*insertPosition) < embeddingNode) {
			insertPositionParent = insertPosition;
			insertPosition = insertPosition->adj;
		}

		if (insertPosition != NULL) {
			if ((*insertPosition) == embeddingNode) {
				queryVertexPositions.push_back(insertPosition);

			}
			else if (insertPosition == (*comboLinkedLists)[queryGraph->graphId][i]) {
				(*comboLinkedLists)[queryGraph->graphId][i] = new EmbeddingNode(embeddingNode);
				(*comboLinkedLists)[queryGraph->graphId][i]->adj = insertPosition;
				queryVertexPositions.push_back((*comboLinkedLists)[queryGraph->graphId][i]);
			}
			else {
				insertPositionParent->adj = new EmbeddingNode(embeddingNode);
				insertPositionParent->adj->adj = insertPosition;
				queryVertexPositions.push_back(insertPositionParent->adj);
			}

		}
		else {
			if (insertPositionParent != NULL) {
				insertPositionParent->adj = new EmbeddingNode(embeddingNode);
				insertPositionParent->adj->adj = insertPosition;
				queryVertexPositions.push_back(insertPositionParent->adj);
			}
			else {
				(*comboLinkedLists)[queryGraph->graphId][i] = new EmbeddingNode(embeddingNode);
				(*comboLinkedLists)[queryGraph->graphId][i]->adj = insertPosition;
				queryVertexPositions.push_back((*comboLinkedLists)[queryGraph->graphId][i]);
			}
		}
	}

	// handle the dfs edges
	for (std::vector<int>::iterator orderComboVertexIterator = DFSComboVertexOrder->begin(); orderComboVertexIterator != DFSComboVertexOrder->end(); orderComboVertexIterator++) {
		for (std::vector<int>::iterator childrenIterator = (*comboGraph)[*orderComboVertexIterator].children.begin(); childrenIterator != (*comboGraph)[*orderComboVertexIterator].children.end(); childrenIterator++) {

			int a = queryVertexPositions[*orderComboVertexIterator]->neighboursList.size();

			std::map<int, std::set<EmbeddingNode *>>::iterator neighbourLinkIterator = queryVertexPositions[*orderComboVertexIterator]->neighboursList.find(*childrenIterator);

			if (neighbourLinkIterator != queryVertexPositions[*orderComboVertexIterator]->neighboursList.end()) {
				neighbourLinkIterator->second.insert(queryVertexPositions[*childrenIterator]);
			}
			else {
				set<EmbeddingNode *> newList;
				newList.insert(queryVertexPositions[*childrenIterator]);
				queryVertexPositions[*orderComboVertexIterator]->neighboursList.insert(std::pair<int, set<EmbeddingNode *>>(*childrenIterator, newList));
			}
		}
	}

	// handle the non-spanning(not in the dfs route) tree edges
	for (int i = 0; i<comboGraph->size(); i++) {
		for (int j = i + 1; j<comboGraph->size(); j++) {
			if (comboGraphEdgeMatrix[i][j] && (*comboGraph)[i].parent != j && (*comboGraph)[j].parent != i) {

				std::map<int, std::set<EmbeddingNode *>>::iterator neighbourLinkIterator = queryVertexPositions[i]->neighboursList.find(j);
				if (neighbourLinkIterator != queryVertexPositions[i]->neighboursList.end()) {
					neighbourLinkIterator->second.insert(queryVertexPositions[j]);
				}
				else {
					set<EmbeddingNode *> newList;
					newList.insert(queryVertexPositions[j]);
					queryVertexPositions[i]->neighboursList.insert(std::pair<int, set<EmbeddingNode *>>(j, newList));
				}
			}
		}
	}
}
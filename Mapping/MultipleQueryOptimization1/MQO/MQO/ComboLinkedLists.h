#pragma once
#ifndef COMBO_LINKED_LISTS_H
#define COMBO_LINKED_LISTS_H

#include"AdjacenceListsGraph.h"
#include<vector>
#include<map>
#include<set>

class EmbeddingNode {
public:
	/*
	* At most 5 values
	*/
	int * comboMappingVertices;

	EmbeddingNode * adj;

	std::map<int, std::set<EmbeddingNode *>> neighboursList;

	EmbeddingNode() {
		comboMappingVertices = new int[5];
		for (int i = 0; i<5; i++) {
			comboMappingVertices[i] = -1;
		}
		adj = NULL;
	};

	EmbeddingNode(const EmbeddingNode & pEmbeddingNode) {
		comboMappingVertices = new int[5];
		for (int i = 0; i<5; i++) {
			comboMappingVertices[i] = pEmbeddingNode.comboMappingVertices[i];
		}
		adj = pEmbeddingNode.adj;
		neighboursList = pEmbeddingNode.neighboursList;
	};

	~EmbeddingNode() {
		neighboursList.swap(std::map<int, std::set<EmbeddingNode *>>());
		delete[] comboMappingVertices;
	}
};

class ComboLinkedLists {

private:
	/*
	* used to save the linked list of cached embeddings
	* points to the head of each linked list
	*/
	std::vector<std::vector<EmbeddingNode *>> * comboLinkedLists;

public: 

	ComboLinkedLists(std::vector<std::vector<EmbeddingNode *>> * pComboLinkedLists);

	void ComboLinkedLists::addEmbeddingToCache(AdjacenceListsGRAPH * queryGraph, int * partialEmbedding);
};


#endif
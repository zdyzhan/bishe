#pragma once
#ifndef PARTI_LINKED_LISTS_H
#define PARTI_LINKED_LISTS_H

#include"AdjacenceListsGraph.h"
#include<vector>
#include<map>
#include<set>


class PartiLinkedLists {

public:

class EmbeddingNode {
	public:
		/*
		* At most 5 values
		*/
		int * comboMappingVertices;

		EmbeddingNode * leftAdj;
		EmbeddingNode * rightAdj;

		std::map<int, std::set<EmbeddingNode *>> neighboursList;

		EmbeddingNode() {
			comboMappingVertices = new int[5];
			for (int i = 0; i<5; i++) {
				comboMappingVertices[i] = -1;
			}
			leftAdj = NULL;
			rightAdj = NULL;
		};

		EmbeddingNode(const EmbeddingNode & pEmbeddingNode) {
			comboMappingVertices = new int[5];
			for (int i = 0; i<5; i++) {
				comboMappingVertices[i] = pEmbeddingNode.comboMappingVertices[i];
			}
			leftAdj = pEmbeddingNode.leftAdj;
			rightAdj = pEmbeddingNode.rightAdj;

			// TODO optimized
			neighboursList = pEmbeddingNode.neighboursList;
		};

		~EmbeddingNode() {
			neighboursList.swap(std::map<int, std::set<EmbeddingNode *>>());
			delete[] comboMappingVertices;
		}
	};



public:

	void PartiLinkedLists::addEmbeddingToCache(AdjacenceListsGRAPH * queryGraph, std::vector<std::vector<EmbeddingNode *>> * cachedResults, int * partialEmbedding);
};


#endif
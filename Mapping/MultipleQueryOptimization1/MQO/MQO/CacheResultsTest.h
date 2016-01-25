#pragma once
#ifndef CACHE_RESULTS_TEST_H
#define CACHE_RESULTS_TEST_H


#include"PartiLinkedLists.h"
#include"TurboIso_CacheResults.h"
#include<vector>
#include<set>
#include<iostream>

class CacheResultsTest {

private: 
	

	std::ofstream * resultFile;

	/*
	* Used to save the number of founded embeddings
	*/
	std::vector<int> numberOfEmbeddings;


	std::vector<std::vector<std::vector<int>>> trivialEmbeddingLists;

	std::vector<std::vector<PartiLinkedLists::EmbeddingNode *>> * cachedLinkedLists;
	PartiLinkedLists * partiLinkedLists;

	/*
	 * datasets
	 */
	std::vector<AdjacenceListsGRAPH> * queryGraphVector;


	AdjacenceListsGRAPH * dataGraph;

	TurboIso_CacheResults * turboIso_CacheResults; 

private:

	long long trivialEmbeddingListMSize;
	long long * linkedListsMSize;
	double * linkedListRetriTime;
	double * linkedListSaveTime;

	PartiLinkedLists::EmbeddingNode ** comboVertexMapping;

public :

	CacheResultsTest(AdjacenceListsGRAPH * pDataGraph, std::vector<AdjacenceListsGRAPH> * pQueryGraphVector, std::ofstream * pResultFile);

	void execute();

private:

	void computeMemory(int queryGraphId, int parWidthIndex);

	void inorderTraversalBT(PartiLinkedLists::EmbeddingNode * embeddingNode, int parWidthIndex, int comboVerticeNumber);

	void saveCLLResults(int parWidth, int queryGraphId);

	void computeCLLRetrivingPerfor(int parWidth, int queryGraphId);

	void recursiveEnumLinkedLists(int comboVertexOrder, AdjacenceListsGRAPH * queryGraph);
};



#endif
#include "CacheResultsTest.h"
#include "TimeUtility.h"
#include<stack>

using namespace std;

CacheResultsTest::CacheResultsTest(AdjacenceListsGRAPH * pDataGraph, std::vector<AdjacenceListsGRAPH>* pQueryGraphVector, std::ofstream * pResultFile)
{
	dataGraph = pDataGraph;
	queryGraphVector = pQueryGraphVector;
	resultFile = pResultFile;
	

	linkedListsMSize = new long long[5];
	linkedListRetriTime = new double[5];
	linkedListSaveTime = new double[5];

	cachedLinkedLists = new std::vector<std::vector<PartiLinkedLists::EmbeddingNode *>>[5];
	partiLinkedLists = new PartiLinkedLists();

	for (int parWidthIndex = 0; parWidthIndex < 5; parWidthIndex++) {

		linkedListsMSize[parWidthIndex] = 0;
		linkedListRetriTime[parWidthIndex] = 0;
		linkedListSaveTime[parWidthIndex] = 0;

		for (int queryGraphId = 0; queryGraphId < queryGraphVector->size(); queryGraphId++) {
			cachedLinkedLists[parWidthIndex].push_back(std::vector<PartiLinkedLists::EmbeddingNode *>());
		}
	}

	for (int queryGraphId = 0; queryGraphId < queryGraphVector->size(); queryGraphId++) {

		numberOfEmbeddings.push_back(0);
		trivialEmbeddingLists.push_back(std::vector<std::vector<int>>());
	}

	trivialEmbeddingListMSize = 0;

	turboIso_CacheResults = new TurboIso_CacheResults(dataGraph, &numberOfEmbeddings, &trivialEmbeddingLists);

}

void CacheResultsTest::execute()
{

	for (int queryGraphId = 0; queryGraphId < queryGraphVector->size(); queryGraphId++) {
		
		turboIso_CacheResults->setParameters(&(*queryGraphVector)[queryGraphId]);
		turboIso_CacheResults->execute();
		
		if (trivialEmbeddingLists[queryGraphId].size() == 0) {
			cout << "** Finished One query " << queryGraphId << " Embedding size: " << trivialEmbeddingLists[queryGraphId].size() << endl;
			continue;
		}

		trivialEmbeddingListMSize += trivialEmbeddingLists[queryGraphId].size() * (*queryGraphVector)[queryGraphId].getNumberOfVertexes() * sizeof(int);

		for (int parWidthIndex = 0; parWidthIndex < 5; parWidthIndex++) {
			saveCLLResults(parWidthIndex + 1, queryGraphId);
			computeCLLRetrivingPerfor(parWidthIndex + 1, queryGraphId);
			/*
			* Memory Comparison
			*/
			computeMemory(queryGraphId, parWidthIndex);
		}
		



		trivialEmbeddingLists[queryGraphId].swap(std::vector<std::vector<int>>());
		for (int parWidthIndex = 0; parWidthIndex < 5; parWidthIndex++) {
			cachedLinkedLists[parWidthIndex][queryGraphId].swap(vector<PartiLinkedLists::EmbeddingNode*>());
		}

	}

	(*resultFile) << "**** Cached Results Memory Size: " << endl;
	(*resultFile) << "Trivial Embedding List: " << trivialEmbeddingListMSize << " bytes" << endl;

	for (int parWidthIndex = 0; parWidthIndex < 5; parWidthIndex++) {
		(*resultFile) << parWidthIndex << "-- Width--Trivial Compressed Graph: " << linkedListsMSize[parWidthIndex] << " bytes" << endl;
	}


	(*resultFile) << "****Retriving Time: " << endl;
	for (int parWidthIndex = 0; parWidthIndex < 5; parWidthIndex++) {
		(*resultFile) << parWidthIndex << "-- Width--Retriving Time: " << linkedListRetriTime[parWidthIndex] << " milliseconds" << endl;
	}

	(*resultFile) << "****Saving Time: " << endl;
	for (int parWidthIndex = 0; parWidthIndex < 5; parWidthIndex++) {
		(*resultFile) << parWidthIndex << "-- Width--Saving Time: " << linkedListSaveTime[parWidthIndex] << " milliseconds" << endl;
	}
}


void CacheResultsTest::computeMemory(int queryGraphId, int parWidthIndex)
{
	for (int j = 0; j < cachedLinkedLists[parWidthIndex][queryGraphId].size(); j++) {

		int comboVerticeNumber = (*queryGraphVector)[queryGraphId].getComboGraph()->at(j).queryVertices.size();
		inorderTraversalBT(cachedLinkedLists[parWidthIndex][queryGraphId][j], parWidthIndex, comboVerticeNumber);
	}
}

void CacheResultsTest::inorderTraversalBT(PartiLinkedLists::EmbeddingNode * embeddingNode, int parWidthIndex, int comboVerticeNumber) {
	if (embeddingNode == NULL) {
		return;
	}
	inorderTraversalBT(embeddingNode->leftAdj, parWidthIndex, comboVerticeNumber);
	linkedListsMSize[parWidthIndex] += comboVerticeNumber * sizeof(int);
	inorderTraversalBT(embeddingNode->rightAdj, parWidthIndex, comboVerticeNumber);

}

void CacheResultsTest::saveCLLResults(int parWidth, int queryGraphId) {
	int * embeddingArray = new int[(*queryGraphVector)[queryGraphId].getNumberOfVertexes()];
	(*queryGraphVector)[queryGraphId].cleanComboInfo();
	(*queryGraphVector)[queryGraphId].buildComboGraph(parWidth);

	for (unsigned int comboVertexIndex = 0; comboVertexIndex < (*queryGraphVector)[queryGraphId].getComboGraph()->size(); comboVertexIndex++) {
		cachedLinkedLists[parWidth - 1][queryGraphId].push_back(NULL);
	}

	for (int embeddingIndex = 0; embeddingIndex < trivialEmbeddingLists[queryGraphId].size(); embeddingIndex++) {
		for (int embeddingDataVertexIndex = 0; embeddingDataVertexIndex < (*queryGraphVector)[queryGraphId].getNumberOfVertexes(); embeddingDataVertexIndex++) {
			embeddingArray[embeddingDataVertexIndex] = trivialEmbeddingLists[queryGraphId][embeddingIndex][embeddingDataVertexIndex];
		}

		TimeUtility llTime;
		llTime.StartCounterMill();
		partiLinkedLists->addEmbeddingToCache(&(*queryGraphVector)[queryGraphId], &cachedLinkedLists[parWidth - 1], embeddingArray);
		linkedListSaveTime[parWidth - 1] += llTime.GetCounterMill();
	}
	delete [] embeddingArray;
}

void CacheResultsTest::computeCLLRetrivingPerfor(int parWidth, int queryGraphId) {
	
	AdjacenceListsGRAPH * queryGraph = &(*queryGraphVector)[queryGraphId];

	comboVertexMapping = new PartiLinkedLists::EmbeddingNode*[queryGraph->getComboGraph()->size()];
	for (int i = 0; i < queryGraph->getComboGraph()->size(); i++) {
		comboVertexMapping[i] = NULL;
	}

	TimeUtility llTime;
	llTime.StartCounterMill();
	PartiLinkedLists::EmbeddingNode * comboVertexStarting = cachedLinkedLists[parWidth - 1][queryGraphId][queryGraph->getDFSComboVertexOrder()->at(0)];
	stack<PartiLinkedLists::EmbeddingNode*> inorderStack;
	
	while (!inorderStack.empty() || comboVertexStarting != NULL) {
		if (comboVertexStarting != NULL) {
			inorderStack.push(comboVertexStarting);
			comboVertexStarting = comboVertexStarting->leftAdj;
		}
		else {
			comboVertexStarting = inorderStack.top();
			inorderStack.pop();
			
			comboVertexMapping[queryGraph->getDFSComboVertexOrder()->at(0)] = comboVertexStarting;
			recursiveEnumLinkedLists(1, queryGraph);

			comboVertexStarting = comboVertexStarting->rightAdj;
		}
	}

	linkedListRetriTime[parWidth - 1] += llTime.GetCounterMill();

	delete [] comboVertexMapping;
}


void CacheResultsTest::recursiveEnumLinkedLists(int comboVertexOrder, AdjacenceListsGRAPH * queryGraph) {

	if (comboVertexOrder == queryGraph->getDFSComboVertexOrder()->size()) {
		// TODO anything ??
		return;
	}

	int comboVertexId = queryGraph->getDFSComboVertexOrder()->at(comboVertexOrder);
	ComboNODE & comboVertex = queryGraph->getComboGraph()->at(comboVertexId);

	std::map<int, std::set<PartiLinkedLists::EmbeddingNode *>>::iterator comboMappingIterator = comboVertexMapping[comboVertex.parent]->neighboursList.find(comboVertexId);
	for (std::set<PartiLinkedLists::EmbeddingNode *>::iterator comboMapping = comboMappingIterator->second.begin(); comboMapping != comboMappingIterator->second.end(); comboMapping++) {
		
		for (std::set<int>::iterator comboVertexNeiIterator = comboVertex.neighbours.begin(); comboVertexNeiIterator != comboVertex.neighbours.end(); comboVertexNeiIterator ++) {
			if (*comboVertexNeiIterator != comboVertex.parent) {
				if (comboVertexMapping[*comboVertexNeiIterator] != NULL) {
					if (*comboVertexNeiIterator < comboVertexId) {
						std::map<int, std::set<PartiLinkedLists::EmbeddingNode *>>::iterator nonTreeMappingIter = comboVertexMapping[*comboVertexNeiIterator]->neighboursList.find(comboVertexId);
						if (nonTreeMappingIter->second.find(*comboMapping) == nonTreeMappingIter->second.end()) {
							continue;
						}
					}
					else {
						std::map<int, std::set<PartiLinkedLists::EmbeddingNode *>>::iterator nonTreeMappingIter = (*comboMapping)->neighboursList.find(*comboVertexNeiIterator);
						if (nonTreeMappingIter->second.find(comboVertexMapping[*comboVertexNeiIterator]) == nonTreeMappingIter->second.end()) {
							continue;
						}
					}
				}
			}
		}

		comboVertexMapping[comboVertexId] = *comboMapping;
		recursiveEnumLinkedLists(comboVertexOrder + 1, queryGraph);
		comboVertexMapping[comboVertexId] = NULL;
	}

}
#include"TurboIso_MQO.h"
#include"TimeUtility.h"
#include"MathUtility.h"
#include"AdjacenceListsGRAPH_IO.h"
#include<iostream>
#include<vector>
#include<map>
#include<set>
#include<queue>
#include <limits>
#include"GlobalConstant.h"

using namespace std;


TurboIso_MQO::TurboIso_MQO() {}

TurboIso_MQO::TurboIso_MQO(AdjacenceListsGRAPH * pDataGraph, std::vector<int> * pNumberOfEmbeddings, ComboLinkedLists * pComboLinkedLists) {
	dataGraph = pDataGraph;
	numberOfEmbeddings = pNumberOfEmbeddings;
	comboLinkedLists = pComboLinkedLists;

}

void TurboIso_MQO::setParameters(AdjacenceListsGRAPH * pQueryGraph, int * pPartialEmbedding, bool pNeedToSaveCache) {

	queryGraph = pQueryGraph;
	partialEmbedding = pPartialEmbedding;
	needFormatCache = pNeedToSaveCache;

	inverseEmbedding.swap(std::map<int, int>());

	for (int i = 0; i<queryGraph->getNumberOfVertexes(); i++) {
		if (partialEmbedding[i] != -1) {
			inverseEmbedding.insert(std::pair<int, int>(partialEmbedding[i], i));
		}
	}
	numberOfPassedPartialMappedV = inverseEmbedding.size();


	/*
	* If enoughRecursiveCalls has been called before any embedding founded. We need to abort the process and just return false
	* tune this value based on the data graphs size. Because it can be too long time if the data graph is a large graph
	*/
	enoughRecursiveCalls = 500000;

	/*
	* For TurboIso, we need to set another threshold for the ExploreCR call numbers
	* tune this value based on the data graphs size. Because it can be too long time if the data graph is a large graph
	*/
	exploreCRRecursiveCalls = 50000;

	needReturn = false;
}



void TurboIso_MQO::clean_before() {
	CR.swap(map<std::pair<int, int>, set<int>>());
	necTree.swap(std::vector<NECNode>());
}

void TurboIso_MQO::clean_after() {
	delete[] necTreeVertexIdInverted;
}



void TurboIso_MQO::execute() {

	clean_before();

	rewriteToNecTree();
	
	for (int i = 0; i<numberOfPassedPartialMappedV; i++) {
		vector<int> subRegionCandidates;
		subRegionCandidates.push_back(partialEmbedding[necTree[i].vertexId]);
		if (exploreCR(i, -1, subRegionCandidates) == false) {
			return;
		}
	}

	subgraphSearch();

	clean_after();
}


void TurboIso_MQO::updateState(int u, int v) {
	partialEmbedding[u] = v;
	inverseEmbedding.insert(std::pair<int, int>(v, u));
}

void TurboIso_MQO::restoreState(int u, int v) {
	partialEmbedding[u] = -1;
	inverseEmbedding.erase(v);
}


void TurboIso_MQO::subgraphSearch() {
	if (enoughRecursiveCalls-- <= 0) {
		// break the limit, so the MQO won't iterate here again
		(*numberOfEmbeddings)[queryGraph->graphId] = -1;
		needReturn = true;
		return;
	}

	if (inverseEmbedding.size() == queryGraph->getNumberOfVertexes()) {
		(*numberOfEmbeddings)[queryGraph->graphId] ++;
		// FIND AN EMBDDING : REPORT IT
		//showEmbedding();
		if (needFormatCache) {
			comboLinkedLists->addEmbeddingToCache(queryGraph, partialEmbedding);
		}
		return;
	}

	NECNode  * nectree_u = nextQueryVertex();

	map<std::pair<int, int>, set<int>>::iterator candidateListIterator = CR.find(std::pair<int, int>(nectree_u->id, partialEmbedding[necTree[nectree_u->parent].vertexId]));
	if (candidateListIterator == CR.end()) {
		return;
	}
	// start iteration candidates
	for (set<int>::iterator candidateIterator = candidateListIterator->second.begin(); candidateIterator != candidateListIterator->second.end(); candidateIterator++) {
		if (needReturn) {
			return;
		}

		if ((*numberOfEmbeddings)[queryGraph->graphId] >= GlobalConstant::G_SUFFICIENT_NUMBER_EMBEDDINGS) {
			break; // only calculate 1000 embeddings for each query graph
		}

		// all ready match filter 
		if (inverseEmbedding.find(*candidateIterator) != inverseEmbedding.end()) {
			continue;
		}
		// isJoinable Filtering
		if (!isJoinable(nectree_u->vertexId, *candidateIterator)) {
			continue;
		}
		updateState(nectree_u->vertexId, *candidateIterator);
		subgraphSearch();
		restoreState(nectree_u->vertexId, *candidateIterator);
	}

}

bool  TurboIso_MQO::isJoinable(int u, int v) {

	AdjacenceListsGRAPH::adjIterator adjIterator(queryGraph, u);

	for (AdjacenceListsGRAPH::link t = adjIterator.begin(); !adjIterator.end(); t = adjIterator.next()) {
		if (partialEmbedding[t->v] != -1) {
			// u has an edge with query vertex t->v which has already been matched
			if (necTreeVertexIdInverted[u] < numberOfPassedPartialMappedV) {
				continue;
			}
			if (necTree[necTreeVertexIdInverted[u]].parent == t->v) {
				continue;
			}
			if (dataGraph->edge(v, partialEmbedding[t->v])) {
				continue;
			}
			else {
				return false;
			}
		}
	}
	return true;
}

TurboIso_MQO::NECNode * TurboIso_MQO::nextQueryVertex()
{
	return & necTree[inverseEmbedding.size()];
}


bool TurboIso_MQO::refineCandidates(AdjacenceListsGRAPH::Vertex & u, const int & v) {
	if (inverseEmbedding.find(v) != inverseEmbedding.end()) {
		return false;
	}
	return true;
}

/* u is the query vertex id, v is the data graph vertex id */
bool TurboIso_MQO::degreeFilter(int u, int v) {
	// TODO, no need extra index
	//return true;

	AdjacenceListsGRAPH::Vertex queryVertex = queryGraph->getVertexByVertexId(u);
	AdjacenceListsGRAPH::Vertex dataVertex = dataGraph->getVertexByVertexId(v);


	for (std::map<int, std::vector<int>>::iterator labelVertexListIterator = queryVertex.labelVertexList.begin(); labelVertexListIterator != queryVertex.labelVertexList.end(); labelVertexListIterator++) {
		if (dataVertex.labelVertexList.find(labelVertexListIterator->first) == dataVertex.labelVertexList.end()) {
			return false;
		}
		else if (labelVertexListIterator->second.size() > dataVertex.labelVertexList.find(labelVertexListIterator->first)->second.size()) {
			return false;
		}
	}
	return true;
}



void TurboIso_MQO::computeMatchingOrder() {
	// TODO nothing
}

void TurboIso_MQO::devideNoTreeEdges(int u, int numberOfNoTreeEdges, float * queryVertexScore) {
	// TODO nothing
}


void TurboIso_MQO::getOrderByBFSScore(int u, float * queryVertexScore, vector<int> & nextVertex) {
	// TODO nothing
}



bool TurboIso_MQO::exploreCR(int u, int parentMappedNecTree, vector<int> & subRegionCandidates) {
	if (exploreCRRecursiveCalls-- <= 0) {
		needReturn = true;
		return false;
	}
	for (int vmIndex = 0; vmIndex < subRegionCandidates.size(); vmIndex++) {

		bool matched = true;

		for (vector<int>::iterator childIterator = necTree[u].childList.begin(); childIterator != necTree[u].childList.end(); childIterator++) {

			vector<int> adjChildAdj;
			AdjacenceListsGRAPH::Vertex & dataGraphChild = dataGraph->getVertexByVertexId(subRegionCandidates[vmIndex]);
			map<int, vector<int>>::iterator childLabelList = dataGraphChild.labelVertexList.find(necTree[*childIterator].label);
			
			vector<int> connectedAlreadyMapped;
			AdjacenceListsGRAPH::adjIterator vertexIterator(queryGraph, necTree[*childIterator].vertexId);
			for (AdjacenceListsGRAPH::link t = vertexIterator.begin(); !vertexIterator.end(); t = vertexIterator.next()) {
				if (partialEmbedding[t->v] != -1) {
					connectedAlreadyMapped.push_back(partialEmbedding[t->v]);
				}
			}


			if (childLabelList != dataGraphChild.labelVertexList.end()) {
				for (vector<int>::iterator childLabelItem = childLabelList->second.begin(); childLabelItem != childLabelList->second.end(); childLabelItem++) {
					
					int alreadyMappedIndex = 0;
					for (; alreadyMappedIndex < connectedAlreadyMapped.size(); alreadyMappedIndex++) {
						if (!dataGraph->edge(*childLabelItem, connectedAlreadyMapped[alreadyMappedIndex])) {
							break;
						}
					}
					if (alreadyMappedIndex < connectedAlreadyMapped.size()) {
						continue;
					}

					if (!degreeFilter(necTree[*childIterator].vertexId, *childLabelItem)) {
						continue;
					}

					adjChildAdj.push_back(*childLabelItem);
				}
			}

			/*
			* A important modification for the exploreCR for BoostIso goes here. We need to add itself here to suit for the dynamic candidate loading strategy
			*/
			if (parentMappedNecTree == -1) {
				if (degreeFilter(necTree[*childIterator].vertexId, subRegionCandidates[vmIndex])) {
					adjChildAdj.push_back(subRegionCandidates[vmIndex]);
				}
			}

			if (exploreCR(necTree[*childIterator].id, subRegionCandidates[vmIndex], adjChildAdj) == false) {
				/*
				* If one of its child cannot match under this subregion
				* 1. stop the process, and the elements at "vmIndex" won't be added
				* Two reasons not to delete anything.
				* 1. No need to, we can allow some rubbish data in the memory while no need to bother to clean them before the computing of this query is done.
				* 2. In some special cases, you cannot delete it. A grand-relationship will cause error.
				*/
				matched = false;
				break;
			}
		}
		if (matched == false) {
			continue;
		}
		map<std::pair<int, int>, set<int>>::iterator  CRIterator = CR.find(std::pair<int, int>(u, parentMappedNecTree));
		if (CRIterator == CR.end()) {
			set<int> CRVertices;
			CRVertices.insert(subRegionCandidates[vmIndex]);
			CR.insert(std::pair<std::pair<int, int>, set<int>>(std::pair<int, int>(u, parentMappedNecTree), CRVertices));
		}
		else {
			CRIterator->second.insert(subRegionCandidates[vmIndex]);
		}
	}

	if (CR.find(std::pair<int, int>(u, parentMappedNecTree)) == CR.end()) {
		return false;
	}

	return true;
}




int TurboIso_MQO::chooseStartVertex() {
	// Nothing to do there because of MQO
	return 0;
}


void TurboIso_MQO::rewriteToNecTree() {

	necTreeVertexIdInverted = new int [queryGraph->getNumberOfVertexes()];
	
	bool * flags = new bool[queryGraph->getNumberOfVertexes()];

	queue<int> parentsVertices;

	for (int i = 0; i<queryGraph->getNumberOfVertexes(); i++) {
		if (partialEmbedding[i] != -1) {
			necTree.push_back(NECNode());
			NECNode & necNode = *(necTree.rbegin());
			necNode.vertexId = i;
			necNode.parent = -1;
			necNode.id = necTree.size() - 1;
			necNode.label = queryGraph->getVertexByVertexId(i).label;
			
			parentsVertices.push(necNode.id);
			necTreeVertexIdInverted[i] = necNode.id;

			flags[i] = true;
		}
		else {
			flags[i] = false;
		}
	}
	while (!parentsVertices.empty()) {
		int parentNodeId = parentsVertices.front();
		parentsVertices.pop();

		AdjacenceListsGRAPH::adjIterator vertexIterator(queryGraph, necTree[parentNodeId].vertexId);
		for (AdjacenceListsGRAPH::link t = vertexIterator.begin(); !vertexIterator.end(); t = vertexIterator.next()) {
			if (!flags[t->v]) {

				necTree.push_back(NECNode());
				NECNode & necNode = *(necTree.rbegin());
				necNode.vertexId = t->v;
				necNode.parent = necTree[parentNodeId].id;
				necNode.id = necTree.size() - 1;
				necNode.label = queryGraph->getVertexByVertexId(t->v).label;

				parentsVertices.push(necNode.id);
				necTree[parentNodeId].childList.push_back(necNode.id);

				necTreeVertexIdInverted[t->v] = necNode.id;

				flags[t->v] = true;
			}
		}
	}
	delete[] flags;
}





void TurboIso_MQO::showEmbedding() {
	std::cout << "{";
	for (int i = 0; i < queryGraph->getNumberOfVertexes(); i++) {
		cout << "*: " << i << "->" << partialEmbedding[i] << " , ";
	}
	cout << "}" << endl;
}

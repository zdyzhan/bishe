#include"TurboIso.h"
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


TurboIso::TurboIso(){}

TurboIso::TurboIso(AdjacenceListsGRAPH * pDataGraph, std::vector<int> * pNumberOfEmbeddings, ComboLinkedLists * pComboLinkedLists) {
	dataGraph = pDataGraph;
	numberOfEmbeddings = pNumberOfEmbeddings;
	comboLinkedLists = pComboLinkedLists;


}

void TurboIso::setParameters(AdjacenceListsGRAPH * pQueryGraph, int * pPartialEmbedding, bool pNeedToSaveCache) {

	//cout << "TurboIso Preparing Data -\/-" << endl; //@debug

	queryGraph = pQueryGraph;
	partialEmbedding = pPartialEmbedding;
	needFormatCache = pNeedToSaveCache;

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




void TurboIso::clean_before() {
	necTree.swap(std::vector<NECNode>());
	inverseEmbedding.swap(std::map<int, int>());
	queryMatchingSuquence.swap(std::vector<int>());
}

void TurboIso::clean_after() {
}


void TurboIso::execute() {
	clean_before();

	startQueryVertex = chooseStartVertex();
	if (startQueryVertex == -1) {
		return;
	}
	//@debug
	//startQueryVertex = 1;

	rewriteToNecTree();

	
	vector<int> * startVertexCandidates = &(dataGraph->getLabelVertexList()->find(queryGraph->getVertexByVertexId(startQueryVertex).label)->second);

	for (vector<int>::iterator startVertexCandidateIterator = startVertexCandidates->begin(); startVertexCandidateIterator != startVertexCandidates->end(); startVertexCandidateIterator++) {
		
		if (needReturn) {
			break;
		}

		if ((*numberOfEmbeddings)[queryGraph->graphId] >= GlobalConstant::G_SUFFICIENT_NUMBER_EMBEDDINGS) {
			break; // only calculate 1000 embeddings for each query graph
		}

		vector<int> subRegionCandidates;
		subRegionCandidates.push_back(*startVertexCandidateIterator);

		// explore the ER start from the root of the nectree
		CR.swap(map<std::pair<int, int>, set<int>>());

		if (exploreCR(0, -1, subRegionCandidates) == false) {
			continue;
		}

		/*
		* Computing the vertication order based on the saved CR
		*/
		computeMatchingOrder();

		updateState(startQueryVertex, *startVertexCandidateIterator);
		subgraphSearch();
		restoreState(startQueryVertex, *startVertexCandidateIterator);
	}
	
	clean_after();
}


void TurboIso::updateState(int u, int v) {
	partialEmbedding[u] = v;
	inverseEmbedding.insert(std::pair<int, int>(v, u));
}

void TurboIso::restoreState(int u, int v) {
	partialEmbedding[u] = -1;
	inverseEmbedding.erase(v);
}


void TurboIso::subgraphSearch() {
	if (enoughRecursiveCalls-- <= 0) {
		(*numberOfEmbeddings)[queryGraph->graphId] = -1;
		needReturn = true;
		return;
	}

	if (inverseEmbedding.size() == queryGraph->getNumberOfVertexes()) {
		(*numberOfEmbeddings)[queryGraph->graphId] ++;
		// FIND AN EMBDDING : REPORT IT
		//showEmbedding();
		if (needFormatCache) {
			//@MQO
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

bool  TurboIso::isJoinable(int u, int v) {

	AdjacenceListsGRAPH::adjIterator adjIterator(queryGraph, u);

	for (AdjacenceListsGRAPH::link t = adjIterator.begin(); !adjIterator.end(); t = adjIterator.next()) {
		if (partialEmbedding[t->v] != -1) {
			// u has an edge with query vertex t->v which has already been matched
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

TurboIso::NECNode * TurboIso::nextQueryVertex()
{
	return &necTree[queryMatchingSuquence[inverseEmbedding.size()]];
}


bool TurboIso::refineCandidates(AdjacenceListsGRAPH::Vertex & u, const int & v) {
	if (inverseEmbedding.find(v) != inverseEmbedding.end()) {
		return false;
	}
	return true;
}

/* u is the query vertex id, v is the data graph vertex id */
bool TurboIso::degreeFilter(int u, int v) {
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



void TurboIso::computeMatchingOrder() {

	// set the start vertex
	float * queryVertexScore = new float[necTree.size()];
	for (int queryVertexIndex = 0; queryVertexIndex<necTree.size(); queryVertexIndex++) {
		queryVertexScore[queryVertexIndex] = 0;
	}

	for (map<std::pair<int, int>, set<int>>::iterator crIterator = CR.begin(); crIterator != CR.end(); crIterator++) {
		queryVertexScore[crIterator->first.first] += crIterator->second.size();
	}

	for (int queryVertexIndex = 0; queryVertexIndex<necTree.size(); queryVertexIndex++) {
		// devide the no-tree edges
		int numberOfNoTreeEdges;
		if (necTree[queryVertexIndex].parent != -1) {
			numberOfNoTreeEdges = queryGraph->getVertexByVertexId(necTree[queryVertexIndex].vertexId).inDegree - necTree[queryVertexIndex].childList.size() - 1;
		}
		else {
			numberOfNoTreeEdges = queryGraph->getVertexByVertexId(necTree[queryVertexIndex].vertexId).inDegree - necTree[queryVertexIndex].childList.size();
		}
		devideNoTreeEdges(queryVertexIndex, numberOfNoTreeEdges + 1, queryVertexScore);
	}

	vector<int> nextVertex;
	for (vector<int>::iterator childIterator = necTree[0].childList.begin(); childIterator != necTree[0].childList.end(); childIterator++) {
		nextVertex.push_back(*childIterator);
	}
	getOrderByBFSScore(0, queryVertexScore, nextVertex);


	delete[] queryVertexScore;
}

void TurboIso::devideNoTreeEdges(int u, int numberOfNoTreeEdges, float * queryVertexScore) {
	queryVertexScore[u] /= numberOfNoTreeEdges;
	for (vector<int>::iterator childIterator = necTree[u].childList.begin(); childIterator != necTree[u].childList.end(); childIterator++) {
		devideNoTreeEdges(*childIterator, numberOfNoTreeEdges, queryVertexScore);
	}
}


void TurboIso::getOrderByBFSScore(int u, float * queryVertexScore, vector<int> & nextVertex){

	queryMatchingSuquence.push_back(u);

	if (nextVertex.size() == 0) {
		return;
	}

	float mimimum = FLT_MAX;
	int mimimumVertex = 1;
	int vertexNextId = -1;

	for (int i = 0; i < nextVertex.size(); i++) {
		if (queryVertexScore[nextVertex[i]] == -1) {
			continue;
		}
		else {
			if (queryVertexScore[nextVertex[i]] < mimimum) {
				mimimumVertex = nextVertex[i];
				vertexNextId = i;
				mimimum = queryVertexScore[nextVertex[i]];
			}
		}
	}
	queryVertexScore[mimimumVertex] = -1;
	if (vertexNextId != -1) {
		nextVertex.erase(nextVertex.begin() + vertexNextId);
	}


	for (vector<int>::iterator childIterator = necTree[mimimumVertex].childList.begin(); childIterator != necTree[mimimumVertex].childList.end(); childIterator++) {
		nextVertex.push_back(*childIterator);
	}
	getOrderByBFSScore(mimimumVertex, queryVertexScore, nextVertex);
}



bool TurboIso::exploreCR(int u, int parentMappedNecTree, vector<int> & subRegionCandidates) {
	if (exploreCRRecursiveCalls -- <= 0) {
		needReturn = true;
		return false;
	}

	for (int vmIndex = 0; vmIndex < subRegionCandidates.size(); vmIndex++) {

		bool matched = true;

		for (vector<int>::iterator childIterator = necTree[u].childList.begin(); childIterator != necTree[u].childList.end(); childIterator++) {

			vector<int> adjChildAdj;

			AdjacenceListsGRAPH::Vertex & dataGraphChild = dataGraph->getVertexByVertexId(subRegionCandidates[vmIndex]);
			map<int, vector<int>>::iterator childLabelList = dataGraphChild.labelVertexList.find(necTree[*childIterator].label);

			if (childLabelList != dataGraphChild.labelVertexList.end()) {
				for (vector<int>::iterator childLabelItem = childLabelList->second.begin(); childLabelItem != childLabelList->second.end(); childLabelItem++) {

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




int TurboIso::chooseStartVertex() {

	std::vector<AdjacenceListsGRAPH::Vertex>::iterator queryVertexIterator;

	std::map<int, vector<int>>* labelVertexList = dataGraph->getLabelVertexList();
	std::map<int, vector<int>>::iterator labelVertexIterator;

	vector<float> rankScore;
	for (queryVertexIterator = queryGraph->getVertexList()->begin(); queryVertexIterator != queryGraph->getVertexList()->end(); queryVertexIterator++) {
		if (queryVertexIterator->inDegree != 0) {
			labelVertexIterator = labelVertexList->find(queryVertexIterator->label);
			if (labelVertexIterator != labelVertexList->end()) {
				rankScore.push_back((float)labelVertexIterator->second.size() / queryVertexIterator->inDegree);
			}
			else {
				return -1;
			}
		}
	}
	// minimum query vertex number is 2
	int topThree[3];
	for (int i = 0; i<3; i++) {
		topThree[i] = -1;
	}
	float miniScore = FLT_MAX;
	for (int i = 0; i<rankScore.size(); i++) {
		if (rankScore[i] < miniScore) {
			miniScore = rankScore[i];
			topThree[0] = i;
		}
	}
	miniScore = FLT_MAX;
	for (int i = 0; i<rankScore.size(); i++) {
		if (i != topThree[0]) {
			if (rankScore[i] < miniScore) {
				miniScore = rankScore[i];
				topThree[1] = i;
			}
		}
	}
	miniScore = FLT_MAX;
	if (rankScore.size() >= 3) {
		for (int i = 0; i<rankScore.size(); i++) {
			if (i != topThree[0] && i != topThree[1]) {
				if (rankScore[i] < miniScore) {
					miniScore = rankScore[i];
					topThree[2] = i;
				}
			}
		}
	}
	int topThreeNumberOfCandidates[3];
	int us;
	int miniTopThree = INT_MAX;
	for (int i = 0; i<3; i++) {
		topThreeNumberOfCandidates[i] = 0;
		if (topThree[i] != -1) {
			map<int, vector<int>>::iterator candidateListIterator = dataGraph->getLabelVertexList()->find(queryGraph->getVertexByVertexId(topThree[i]).label);
			if (candidateListIterator == dataGraph->getLabelVertexList()->end()) {
				// the query vertex doesn't have candidates
				return -1;
			}
			for (vector<int>::iterator candidateIterator = candidateListIterator->second.begin();
			candidateIterator != candidateListIterator->second.end(); candidateIterator++) {
				if (degreeFilter(topThree[i], *candidateIterator)) {
					topThreeNumberOfCandidates[i] ++;
				}
			}
			if (miniTopThree > topThreeNumberOfCandidates[i]) {
				us = topThree[i];
				miniTopThree = topThreeNumberOfCandidates[i];
			}
		}
	}
	return us;
}


void TurboIso::rewriteToNecTree() {
	bool *flags = new bool[queryGraph->getNumberOfVertexes()];
	for (int i = 0; i<queryGraph->getNumberOfVertexes(); i++) {
		flags[i] = false;
	}

	queue<std::pair<int, int>> S; // id, parentId
	int vertexId, parentId;

	S.push(std::pair<int, int>(startQueryVertex, -1));
	flags[startQueryVertex] = true;

	while (!S.empty()) {
		vertexId = S.front().first;
		parentId = S.front().second;
		S.pop();
		necTree.push_back(NECNode());
		NECNode & necNode = *(necTree.rbegin());
		necNode.vertexId = vertexId;
		necNode.parent = parentId;
		necNode.id = necTree.size() - 1;
		necNode.label = queryGraph->getVertexByVertexId(vertexId).label;

		if (parentId != -1)
			necTree[parentId].childList.push_back(necNode.id);

		AdjacenceListsGRAPH::adjIterator vertexIterator(queryGraph, vertexId);
		for (AdjacenceListsGRAPH::link t = vertexIterator.begin(); !vertexIterator.end(); t = vertexIterator.next()) {
			if (!flags[t->v]) {
				flags[t->v] = true;
				S.push(std::pair<int, int>(t->v, necNode.id));
			}
		}
	}
	delete[] flags;
}





void TurboIso::showEmbedding() {
	std::cout<<"{";
	for(int i = 0; i < queryGraph -> getNumberOfVertexes(); i++){
		cout<<"*: "<<i<<"->"<< partialEmbedding[i]<<" , ";
	}
	cout<<"}"<<endl;
}

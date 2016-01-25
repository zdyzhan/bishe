#include"MathUtility.h"
#include<iostream>

#include <limits>
#include"TurboIsoBoosted_MQO.h"
#include"AdjacenceListsGRAPH_BOOST.h"
#include"GlobalConstant.h"

using namespace std;



TurboIsoBoosted_MQO::TurboIsoBoosted_MQO() {}

TurboIsoBoosted_MQO::TurboIsoBoosted_MQO(AdjacenceListsGRAPH * pHyperGraph, std::vector<int> * pNumberOfEmbeddings, ComboLinkedLists * pComboLinkedLists) {
	hyperGraph = pHyperGraph;
	numberOfEmbeddings = pNumberOfEmbeddings;
	comboLinkedLists = pComboLinkedLists;


}

void TurboIsoBoosted_MQO::setParameters(AdjacenceListsGRAPH * pQueryGraph, int * pPartialEmbedding, bool pNeedToSaveCache) {

	queryGraph = pQueryGraph;
	partialEmbedding = pPartialEmbedding;
	needFormatCache = pNeedToSaveCache;

	inverseEmbedding.swap(std::map<int, std::stack<int>>());
	mappedQueryVertexSize = 0;


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


	std::map<int, std::stack<int>>::iterator iEIterator;
	for (int i = 0; i<queryGraph->getNumberOfVertexes(); i++) {
		if (partialEmbedding[i] != -1) {
			iEIterator = inverseEmbedding.find(partialEmbedding[i]);
			if (iEIterator == inverseEmbedding.end()) {

				std::stack<int> newMappedDataVertex;
				newMappedDataVertex.push(partialEmbedding[i]);

				inverseEmbedding.insert(std::pair<int, std::stack<int>>(partialEmbedding[i], newMappedDataVertex));

			}
			else {
				iEIterator->second.push(partialEmbedding[i]);
			}

			mappedQueryVertexSize++;
		}
	}
	numberOfPassedPartialMappedV = mappedQueryVertexSize;
}


void TurboIsoBoosted_MQO::clean_before() {

	CR.swap(map<std::pair<int, int>, set<int>>());
	necTree.swap(std::vector<NECNode>());

	scRootCandidates.swap(std::map<int, vector<int>* >());

}

void TurboIsoBoosted_MQO::clean_after() {
	delete[] necTreeVertexIdInverted;
}


void TurboIsoBoosted_MQO::execute() {

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



/**
* The candidate filter for turboIso only compute the candidates for the start vertex
*/
void TurboIsoBoosted_MQO::filterCandidates() {
	// TODO nothing
}

void TurboIsoBoosted_MQO::updateState(AdjacenceListsGRAPH::Vertex * u, int v) {
	partialEmbedding[u->id] = v;

	mappedQueryVertexSize++;
	if (inverseEmbedding.find(v) == inverseEmbedding.end()) {
		stack<int> vertexList;
		vertexList.push(u->id);
		inverseEmbedding.insert(std::pair<int, stack<int>>(v, vertexList));
	}
	else {
		inverseEmbedding.find(v)->second.push(u->id);
	}

}

void TurboIsoBoosted_MQO::restoreState(AdjacenceListsGRAPH::Vertex * u, int v) {
	partialEmbedding[u->id] = -1;
	mappedQueryVertexSize--;

	std::map<int, std::stack<int>>::iterator Global_stack_Iterator = inverseEmbedding.find(v);
	Global_stack_Iterator->second.pop();

	if (Global_stack_Iterator->second.size() == 0) {
		inverseEmbedding.erase(v);
	}
}



void TurboIsoBoosted_MQO::subgraphSearch() {
	if (enoughRecursiveCalls-- <= 0) {
		needReturn = true;
		return;
	}

	if (mappedQueryVertexSize == queryGraph->getNumberOfVertexes()) {
		generateEquivalentEmbedding();
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

		// isJoinable Filtering
		if (!isJoinable(queryGraph->getVertexAddressByVertexId(nectree_u->vertexId), hyperGraph->getVertexAddressByVertexId(*candidateIterator))) {
			continue;
		}

		updateState(queryGraph->getVertexAddressByVertexId(nectree_u->vertexId), *candidateIterator);
		subgraphSearch();
		restoreState(queryGraph->getVertexAddressByVertexId(nectree_u->vertexId), *candidateIterator);
	}

}

bool TurboIsoBoosted_MQO::isJoinable(AdjacenceListsGRAPH::Vertex * u, AdjacenceListsGRAPH::Vertex * v) {

	for (AdjacenceListsGRAPH::link t = u->adj; t != 0; t = t->next) {

		if (partialEmbedding[t->v] != -1) {
			if (partialEmbedding[t->v] == v->id) {
				/*
				* If two connected query vertices are trying to match to the same hypervertex, we need to make sure this hypervertex is a clique
				*/
				if (v->isClique == 2) {
					return false;
				}
			}
			else {
				if (necTreeVertexIdInverted[u->id] < numberOfPassedPartialMappedV) {
					continue;
				}
				if (necTree[necTreeVertexIdInverted[u->id]].parent == t->v) {
					continue;
				}
				if (hyperGraph->edge(v->id, partialEmbedding[t->v])) {
					continue;
				}
				return false;
			}
		}
		/*
		* Check whether the data vertex has already been matched, If it is, we need to make sure it cannot be matched more than the data vertices it has
		*/
		std::map<int, std::stack<int>>::iterator inverseIterator = inverseEmbedding.find(v->id);
		if (inverseIterator != inverseEmbedding.end()) {
			if (inverseIterator->second.size() + 1 > v->sEquivalentVertexList.size()) {
				return false;
			}
		}
	}
	return true;
}


TurboIsoBoosted_MQO::NECNode * TurboIsoBoosted_MQO::nextQueryVertex()
{
	return &necTree[mappedQueryVertexSize];
}

/*
*
*/
void TurboIsoBoosted_MQO::computeMatchingOrder() {
	// TODO nothing
}


void TurboIsoBoosted_MQO::devideNoTreeEdges(int u, int numberOfNoTreeEdges, float * queryVertexScore) {
	// TODO nothing
}

void TurboIsoBoosted_MQO::getOrderByBFSScore(int u, float * queryVertexScore, vector<int> & nextVertex) {
	// TODO nothing
}


bool TurboIsoBoosted_MQO::exploreCR(int u, int parentMappedNecTree, vector<int> & subRegionCandidates) {

	if (exploreCRRecursiveCalls-- <= 0) {
		needReturn = true;
		return false;
	}
	for (int vmIndex = 0; vmIndex < subRegionCandidates.size(); vmIndex++) {

		bool matched = true;

		for (vector<int>::iterator childIterator = necTree[u].childList.begin(); childIterator != necTree[u].childList.end(); childIterator++) {

			vector<int> adjChildAdj;
			AdjacenceListsGRAPH::Vertex & dataGraphChild = hyperGraph->getVertexByVertexId(subRegionCandidates[vmIndex]);
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
						if (!hyperGraph->edge(*childLabelItem, connectedAlreadyMapped[alreadyMappedIndex])) {
							break;
						}
					}
					if (alreadyMappedIndex < connectedAlreadyMapped.size()) {
						continue;
					}
					// TODO, no need extra index
					//if (!AdjacenceListsGRAPH_BOOST::degreeFilter(queryGraph, hyperGraph, necTree[*childIterator].vertexId, *childLabelItem)) {
					//	continue;
					//}

					adjChildAdj.push_back(*childLabelItem);
				}
			}

			/*
			* A important modification for the exploreCR for BoostIso goes here. We need to add itself here to suit for the dynamic candidate loading strategy
			*/
			if (parentMappedNecTree == -1) {
				if (AdjacenceListsGRAPH_BOOST::degreeFilter(queryGraph, hyperGraph, necTree[*childIterator].vertexId, subRegionCandidates[vmIndex])) {
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




/*
* Dynamically append candidates based on the (s and qd) containment relationships
*/
void TurboIsoBoosted_MQO::dynamicLoadCandidates(int dataVertexId) {
	//TODO nothing
}


/*
* Choose the starting query vertex. This process is based on the frequence information for each label
*/
int TurboIsoBoosted_MQO::chooseStartVertex() {

	// TODO nothing
	return 0;
}


void TurboIsoBoosted_MQO::rewriteToNecTree() {

	necTreeVertexIdInverted = new int[queryGraph->getNumberOfVertexes()];

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
				flags[t->v] = true;

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

void TurboIsoBoosted_MQO::generateEquivalentEmbedding() {
	//@MQO
	if (needFormatCache) {
		comboLinkedLists->addEmbeddingToCache(queryGraph, partialEmbedding);
	}

	//showEmbedding();
	int newFinalEmbeddings = 1;
	/*
	* the key represents the mapped hyperVertexId
	* the first int in the pair is the number of data vertices that the hyper vertex has
	* the second int in the pair indicate how many times(how many query vertex) this hypervertex is matched.
	*/
	map<int, std::pair<int, int>> distinctHyperVertexId;

	for (int pEindex = 0; pEindex < mappedQueryVertexSize; pEindex++) {
		map<int, std::pair<int, int>>::iterator mapTwoIntEndMapIterator = distinctHyperVertexId.find(partialEmbedding[pEindex]);
		if (mapTwoIntEndMapIterator == distinctHyperVertexId.end()) {
			distinctHyperVertexId.insert(
				std::pair<int, std::pair<int, int>>(
					partialEmbedding[pEindex],
					std::pair<int, int>(hyperGraph->getVertexByVertexId(partialEmbedding[pEindex]).sEquivalentVertexList.size(), 1)
					)
				);
		}
		else {
			mapTwoIntEndMapIterator->second.second += 1;
		}
	}

	for (map<int, std::pair<int, int>>::iterator mapTwoIntEndMapIterator = distinctHyperVertexId.begin(); mapTwoIntEndMapIterator != distinctHyperVertexId.end(); mapTwoIntEndMapIterator++) {
		if (mapTwoIntEndMapIterator->second.second > 1) {
			newFinalEmbeddings *= Math_Utility::combinations(mapTwoIntEndMapIterator->second.first, mapTwoIntEndMapIterator->second.second);
		}
		else {
			newFinalEmbeddings *= mapTwoIntEndMapIterator->second.first;
		}
	}


	(*numberOfEmbeddings)[queryGraph->graphId] += newFinalEmbeddings;
}

void TurboIsoBoosted_MQO::generateQDEquivalentEmbedding() {
	// TODO nothing
}




void TurboIsoBoosted_MQO::showEmbedding() {
	std::cout << "{";
	for (int i = 0; i < queryGraph->getNumberOfVertexes(); i++) {
		cout << "*: " << i << "->" << partialEmbedding[i] << " , ";
	}
	cout << "}" << endl;
}

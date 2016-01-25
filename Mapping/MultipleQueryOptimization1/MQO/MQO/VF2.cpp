#include"VF2.h"
#include"TimeUtility.h"
#include"AdjacenceListsGRAPH_IO.h"
#include<iostream>
#include<vector>
#include<map>
#include<set>
#include"GlobalConstant.h"

using namespace std;
/*Helper Functions, Use With Care*/

//@unique for VF2 to be used in the update Cg
int vd_parent_level;

//@common
map<int,int>::iterator Global_Map_Int_Iterator; // one use it with in the isJoinable()
std::map<int,std::stack<int>>::iterator Global_stack_Iterator; // used in refineCandidate


VF2::VF2() {}

VF2::VF2(AdjacenceListsGRAPH * pDataGraph, std::vector<int> * pNumberOfEmbeddings, ComboLinkedLists * pComboLinkedLists) {
	
	dataGraph = pDataGraph;
	numberOfEmbeddings = pNumberOfEmbeddings;
	comboLinkedLists = pComboLinkedLists;
}

VF2::~VF2() {
	candidates.swap(std::map<int, vector<int>* >());
	inverseEmbedding.swap(std::map<int, int>());
	Cg.swap(std::map<int, int>());
	Cq_Mq.swap(std::map<int, std::pair<int, int>>());
}

void VF2::setParameters(AdjacenceListsGRAPH * pQueryGraph, int * pPartialEmbedding, bool pNeedToSaveCache) {

	queryGraph = pQueryGraph;
	DFSTraversalOrder = queryGraph->getDFSTraversalOrder();
	partialEmbedding = pPartialEmbedding;
	needFormatCache = pNeedToSaveCache;

	inverseEmbedding.swap(std::map<int, int>());

	for (int i = 0; i<queryGraph->getNumberOfVertexes(); i++) {
		if (partialEmbedding[i] != -1) {
			inverseEmbedding.insert(std::pair<int, int>(partialEmbedding[i], i));
		}
	}
}


void VF2::execute() {

	clean();

	preCalculateCq();

	filterCandidates();

	if (candidates.size() != queryGraph->getNumberOfVertexes() - inverseEmbedding.size()) {
		return;
	}

	/*
	* computing order, starting from the most less candidates
	*/
	int smallestId;
	int leastNumberCandidates = INT_MAX;
	for (std::map<int, vector<int>* >::iterator iterator = candidates.begin(); iterator != candidates.end(); iterator++) {
		if (iterator->second->size() < leastNumberCandidates) {
			leastNumberCandidates = iterator->second->size();
			smallestId = iterator->first;
		}
	}


	subgraphSearch();
}


void VF2::preCalculateCq() {

	std::pair<int,int> Cq_Mq_pair;
	set<int> Cq;
	bool* flags = new bool[DFSTraversalOrder->size()]();

	for(std::vector<AdjacenceListsGRAPH::Vertex *>::iterator queryVertexIterator = DFSTraversalOrder->begin(); queryVertexIterator != DFSTraversalOrder->end(); queryVertexIterator++){
		Cq_Mq_pair.first = 0;
		Cq_Mq_pair.second = 0;
		// update Cq
		flags[(*queryVertexIterator)->id] = true;

		AdjacenceListsGRAPH::adjIterator queryVertexAdj(queryGraph, (*queryVertexIterator)->id);
		for(AdjacenceListsGRAPH::link t = queryVertexAdj.begin(); !queryVertexAdj.end() ; t = queryVertexAdj.next()){
			if(flags[t->v] == true){
				continue;
			}
			if(Cq.insert(t -> v).second == false) {
				Cq_Mq_pair.first ++;
			}
			else {
				Cq_Mq_pair.second ++;
			}
		}
		Cq_Mq.insert(std::pair<int,std::pair<int,int>>((*queryVertexIterator)->id, Cq_Mq_pair));
	}

	delete[] flags;
}



void VF2::clean(){

	candidates.swap(std::map<int, vector<int>* >());

	Cg.swap(std::map<int, int>());
	Cq_Mq.swap(std::map<int, std::pair<int, int>>());

	numberOfCgAdj = -1;
	numberOfCg_Mg_Adj = -1;

}

//@Common
void VF2::filterCandidates(){

	map<int, vector<int>>::iterator candidateSetsIterator;

	vector<AdjacenceListsGRAPH::Vertex> * queryVertexLists = queryGraph->getVertexList();
	map<int, vector<int>> * labelDataVertexList = dataGraph->getLabelVertexList();

	for (vector<AdjacenceListsGRAPH::Vertex>::iterator queryVertexIterator = queryVertexLists->begin(); queryVertexIterator != queryVertexLists->end(); queryVertexIterator++) {
		/*
		* Only consider not matched vertices
		*/
		if (partialEmbedding[queryVertexIterator->id] != -1) {
			continue;
		}
		candidateSetsIterator = labelDataVertexList->find(queryVertexIterator->label);
		if (candidateSetsIterator != labelDataVertexList->end()) {
			candidates.insert(std::pair<int, std::vector<int> *>(queryVertexIterator->id, &candidateSetsIterator->second));
		}
		else {
			// to do ? the query has no embeddings
			return;
		}
	}
}

void VF2::subgraphSearch() {

	if (inverseEmbedding.size() == queryGraph->getNumberOfVertexes()) {
		(*numberOfEmbeddings)[queryGraph->graphId] ++;
		// FIND AN EMBDDING : REPORT IT
		//showEmbedding();
		if (needFormatCache) {
			comboLinkedLists->addEmbeddingToCache(queryGraph, partialEmbedding);
		}
		return;
	}

	AdjacenceListsGRAPH::Vertex u = nextQueryVertex();
	if (partialEmbedding[u.id] != -1) {
		//@MQO already matched
		updateCg(partialEmbedding[u.id]);
		subgraphSearch();
		restoreCg(partialEmbedding[u.id]);
	}
	else {
		vector<int> * candidates_u = candidates.find(u.id)->second;
		// For each v in C(u)
		for (vector<int>::iterator v = candidates_u->begin(); v != candidates_u->end(); v++) {
			if ((*numberOfEmbeddings)[queryGraph->graphId] >= GlobalConstant::G_SUFFICIENT_NUMBER_EMBEDDINGS) {
				break; // only calculate 1000 embeddings for each query graph
			}
			//refine candidate
			if (!refineCandidates(u, *v)) {
				continue;
			}
			// check whether the corresponding edges exists
			if (!isJoinable(u.id, *v)) {
				continue;
			}
			if (!isFlexible(u, *v)) {
				continue;
			}

			updateCg(*v);
			// Next
			updateState(u.id, *v);
			subgraphSearch();
			restoreState(u.id, *v);

			restoreCg(*v);
		}
	}
}


bool VF2::isFlexible(AdjacenceListsGRAPH::Vertex & u, int & v){

	numberOfCgAdj = 0;
	numberOfCg_Mg_Adj = 0;

	AdjacenceListsGRAPH::adjIterator dataVertexAdj(dataGraph, v);
	for(AdjacenceListsGRAPH::link t = dataVertexAdj.begin(); !dataVertexAdj.end() ; t = dataVertexAdj.next()){
		if(inverseEmbedding.find(t->v) == inverseEmbedding.end()) {
			if(Cg[t->v] == 0){
				numberOfCg_Mg_Adj ++;  // if t -> v is not in Embedding and not in Cg
			}
			else {
				numberOfCgAdj ++; // if v in Cg
			}
		}
	}

	std::pair<int,int> Cq_Mq_pair = Cq_Mq.find(u.id) -> second;

	/* According to the description of An inembedding.size() comparison, the folowing pruning rule may be wrong */
	if( Cq_Mq_pair.first > numberOfCgAdj) {
		return false;
	}
	if( Cq_Mq_pair.second > numberOfCg_Mg_Adj){
		return false;
	}

	return true;
}



AdjacenceListsGRAPH::Vertex VF2::nextQueryVertex() {
	for (int i = 0; i<queryGraph->getNumberOfVertexes(); i++) {
		if (partialEmbedding[DFSTraversalOrder->at(i)->id] == -1) {
			return DFSTraversalOrder->at(i)->id;
		}
	}
}


bool VF2::refineCandidates(AdjacenceListsGRAPH::Vertex & u, const int & v){

	if (inverseEmbedding.find(v) != inverseEmbedding.end()) {
		return false;
	}
	// The candidates whose degree is smaller than the query vertex's will be filtered
	if (!degreeFilter(u.id, v)) {
		return false;
	}
	return true;
}

bool VF2::isJoinable(int u, int v) {

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





void VF2::updateCg(int & v) {

	vd_parent_level = Cg[v];
	Cg[v] = 0;

	AdjacenceListsGRAPH::adjIterator dataVertexAdj(dataGraph, v);
	for(AdjacenceListsGRAPH::link t = dataVertexAdj.begin(); !dataVertexAdj.end() ; t = dataVertexAdj.next()){
		if (inverseEmbedding.find(t->v) == inverseEmbedding.end() && Cg[t->v] == 0) // if va is not in Mg and not in Cg
			Cg[t->v] = inverseEmbedding.size() + 1; // add it to Cg and mark it with the current embedding.size()
	}
}

void VF2::restoreCg(int & v){
	AdjacenceListsGRAPH::adjIterator dataVertexAdj(dataGraph, v);
	for(AdjacenceListsGRAPH::link t = dataVertexAdj.begin(); !dataVertexAdj.end() ; t = dataVertexAdj.next()){
		int va = t->v;
		if (Cg[va] == inverseEmbedding.size() + 1){
			Cg[va] = 0;
		}
	}

	Cg[v] = vd_parent_level;
}


/* u is the query vertex id, v is the data graph vertex id */
bool VF2::degreeFilter(int u, int v) {


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


void VF2::updateState(int u, int v) {
	partialEmbedding[u] = v;
	inverseEmbedding.insert(std::pair<int, int>(v, u));
}

void VF2::restoreState(int u, int v) {
	partialEmbedding[u] = -1;
	inverseEmbedding.erase(v);
}


void VF2::showEmbedding() {
	std::cout<<"{";
	for(int i = 0; i < queryGraph -> getNumberOfVertexes(); i++){
		cout<<"Embedding: "<<i<<"->"<<partialEmbedding[i]<<" , ";
	}
	cout<<"}"<<endl;
}
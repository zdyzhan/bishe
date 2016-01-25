#include"Ullmann.h"
#include"TimeUtility.h"
#include"AdjacenceListsGRAPH_IO.h"
#include<iostream>
#include<vector>
#include<map>
#include<queue>
#include<set>
#include<limits.h>
#include"GlobalConstant.h"

using namespace std;

Ullmann::Ullmann(){}

Ullmann::Ullmann(AdjacenceListsGRAPH * pDataGraph, std::vector<int> * pNumberOfEmbeddings, ComboLinkedLists * pComboLinkedLists){
	dataGraph = pDataGraph;
	numberOfEmbeddings = pNumberOfEmbeddings;
	comboLinkedLists = pComboLinkedLists;
}

Ullmann::~Ullmann() {
	candidates.swap(std::map<int, vector<int>* >());
	inverseEmbedding.swap(std::map<int, int>());
}

void Ullmann::setParameters(AdjacenceListsGRAPH * pQueryGraph, int * pPartialEmbedding, bool pNeedToSaveCache){
	
	queryGraph = pQueryGraph;
	DFSTraversalOrder = queryGraph->getDFSTraversalOrder();
	partialEmbedding = pPartialEmbedding;
	needFormatCache = pNeedToSaveCache;

	inverseEmbedding.swap(std::map<int,int>());

	for(int i=0; i<queryGraph->getNumberOfVertexes(); i++){
		if(partialEmbedding[i] != -1){
			inverseEmbedding.insert(std::pair<int, int>(partialEmbedding[i],i));
		}
	}
}

void Ullmann::execute() {

	clean();

	filterCandidates();

	if(candidates.size() != queryGraph -> getNumberOfVertexes() - inverseEmbedding.size()) {
		return;
	}

	/*
	 * computing order, starting from the most less candidates
	 */
	int smallestId;
	int leastNumberCandidates=INT_MAX;
	for(std::map<int, vector<int>* >::iterator iterator = candidates.begin(); iterator != candidates.end(); iterator++){
		if(iterator->second->size() < leastNumberCandidates){
			leastNumberCandidates = iterator->second->size();
			smallestId = iterator->first;
		}
	}


	subgraphSearch();
}


void Ullmann::clean() {
	// clear the candidates
	candidates.swap(std::map<int, vector<int>* > ()); 
}


void Ullmann::filterCandidates() {

	map<int,vector<int>>::iterator candidateSetsIterator;

	vector<AdjacenceListsGRAPH::Vertex> * queryVertexLists = queryGraph -> getVertexList();
	map<int,vector<int>> * labelDataVertexList = dataGraph -> getLabelVertexList();

	for(vector<AdjacenceListsGRAPH::Vertex>::iterator queryVertexIterator = queryVertexLists->begin(); queryVertexIterator != queryVertexLists->end(); queryVertexIterator++){
		/*
		 * Only consider not matched vertices
		 */
		if(partialEmbedding[queryVertexIterator->id] != -1){
			continue;
		}
		candidateSetsIterator = labelDataVertexList -> find(queryVertexIterator->label);
		if(candidateSetsIterator != labelDataVertexList -> end()) {
			candidates.insert(std::pair<int,std::vector<int> *>(queryVertexIterator->id, & candidateSetsIterator->second));
		}
		else{
			// to do ? the query has no embeddings
			return;
		}
	}
}

void Ullmann::subgraphSearch() {
	if(inverseEmbedding.size() == queryGraph -> getNumberOfVertexes()){
		(*numberOfEmbeddings)[queryGraph->graphId] ++;
		// FIND AN EMBDDING : REPORT IT
		//showEmbedding();
		if(needFormatCache){
			comboLinkedLists->addEmbeddingToCache(queryGraph, partialEmbedding);
		}
		return ;
	}

	AdjacenceListsGRAPH::Vertex u = nextQueryVertex();

	if(partialEmbedding[u.id] != -1){
		//@MQO already matched
		subgraphSearch();
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
			if (isJoinable(u.id, *v)) {
				updateState(u.id, *v);
				subgraphSearch();
				restoreState(u.id, *v);
			}
		}
	}
}

AdjacenceListsGRAPH::Vertex Ullmann::nextQueryVertex() {
	for(int i=0; i<queryGraph->getNumberOfVertexes(); i++){
		if(partialEmbedding[DFSTraversalOrder->at(i)->id] == -1){
			return DFSTraversalOrder->at(i)->id;
		}
	}
}

bool Ullmann::refineCandidates(AdjacenceListsGRAPH::Vertex & u, const int & v) {

	if(inverseEmbedding.find(v) != inverseEmbedding.end()) {
		return false;
	}
	// The candidates whose degree is smaller than the query vertex's will be filtered
	if(!degreeFilter(u.id,v)) {
		return false;
	}
	return true;
}

bool Ullmann::isJoinable(int u, int v) {

	AdjacenceListsGRAPH::adjIterator adjIterator(queryGraph, u);

	for(AdjacenceListsGRAPH::link t = adjIterator.begin();  !adjIterator.end(); t=adjIterator.next()) {
		if(partialEmbedding[t->v] != -1){
			// u has an edge with query vertex t->v which has already been matched
			if( dataGraph -> edge( v, partialEmbedding[t->v])){
				continue;
			}else{
				return false;
			}
		}
	}
	return true;
}

void Ullmann::updateState(int u, int v) {
	partialEmbedding[u] = v;
	inverseEmbedding.insert(std::pair<int,int>(v, u));
}

void Ullmann::restoreState(int u, int v){
	partialEmbedding[u] = -1;
	inverseEmbedding.erase(v);
}


/* u is the query vertex id, v is the data graph vertex id */
bool Ullmann::degreeFilter(int u, int v) {


	AdjacenceListsGRAPH::Vertex queryVertex = queryGraph->getVertexByVertexId(u);
	AdjacenceListsGRAPH::Vertex dataVertex = dataGraph->getVertexByVertexId(v);

	
	for(std::map<int,std::vector<int>>::iterator labelVertexListIterator =  queryVertex.labelVertexList.begin();labelVertexListIterator != queryVertex.labelVertexList.end(); labelVertexListIterator++){
		if(dataVertex.labelVertexList.find(labelVertexListIterator->first) == dataVertex.labelVertexList.end()){
			return false;
		}
		else if(labelVertexListIterator->second.size() > dataVertex.labelVertexList.find(labelVertexListIterator->first)->second.size()){
			return false;
		}
	}
	return true;
}

void Ullmann::showEmbedding() {
	std::cout << "{";
	for (int i = 0; i < queryGraph->getNumberOfVertexes(); i++) {
		cout << "Embedding: " << i << "->" << partialEmbedding[i] << " , ";
	}
	cout << "}" << endl;
}
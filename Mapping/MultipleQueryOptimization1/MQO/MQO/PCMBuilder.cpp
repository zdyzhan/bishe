#include"PCMBuilder.h"
#include"AdjacenceListsGRAPH_IO.h"
#include<iostream>
#include<algorithm>
#include<set>
#include"GlobalConstant.h"
#include"TimeUtility.h"


using namespace std;


PCMBuilder::PCMBuilder(int ** pTlsGraphMatrix, std::vector<PCM_Node> * pPatternContainmentMap, std::vector<AdjacenceListsGRAPH> * pQueryGraphVector, std::vector<AdjacenceListsGRAPH> * pNewlyGeneratedGraphVector, std::ofstream * pResultFile){
	
	tlsGraphMatrix = pTlsGraphMatrix;
	
	queryGraphVector = pQueryGraphVector;
	
	newlyGeneratedGraphVector = pNewlyGeneratedGraphVector;

	patternContainmentMap = pPatternContainmentMap;

	graphIsomorphism = QueryGraphIsomorphismTest();

	subgraphIsomorphism = QuerySubgraphIsomorphismSearch();

	findAllMaximalCliques = new FindAllMaximalCliques(pQueryGraphVector);

	resultFile = pResultFile;
}

PCMBuilder::~PCMBuilder(){

}

void PCMBuilder::hideIsomorphicQueries()
{
	for (int i = 0; i < queryGraphVector->size(); i++) {
		for (int j = i + 1; j < queryGraphVector->size(); j++) {
			if ((*patternContainmentMap)[j].isHidden == true 
				|| tlsGraphMatrix[i][j] != (*queryGraphVector)[i].getTotalTLSNumber()
				|| tlsGraphMatrix[i][j] != (*queryGraphVector)[j].getTotalTLSNumber()) {
				continue;
			}
			if (graphIsomorphism.isGraphIsomorphic(&(*queryGraphVector)[i], &(*queryGraphVector)[j], resultFile)) {
				/*
				* Those two query graphs are isomorphic
				*/
				(*patternContainmentMap)[i].equivalentGraph.push_back(j);
				(*patternContainmentMap)[j].isHidden = true;

				for (int z = 0; z < queryGraphVector->size(); z++) {
					tlsGraphMatrix[z][j] = 0;
					tlsGraphMatrix[j][z] = 0;
				}

			}
		}
	}
}

void PCMBuilder::buildContainmentRelations()
{

	for (int i = 0; i < queryGraphVector->size(); i++) {
		if ((*patternContainmentMap)[i].isHidden == true) {
			continue;
		}
		for (int j = i + 1; j < queryGraphVector->size(); j++) {
			/*
			 * Equivalent Filter
			 */
			if ((*patternContainmentMap)[j].isHidden == true) {
				continue;
			}

			AdjacenceListsGRAPH * dataGraph, *queryGraph;
			if ((*queryGraphVector)[i].getNumberOfEdges() >= (*queryGraphVector)[j].getNumberOfEdges()) {
				dataGraph = &(*queryGraphVector)[i];
				queryGraph = &(*queryGraphVector)[j];
			}
			else {
				dataGraph = &(*queryGraphVector)[j];
				queryGraph = &(*queryGraphVector)[i];
			}
			/*
			 * Use TLS filter
			 */
			if (tlsGraphMatrix[i][j] != queryGraph->getTotalTLSNumber()) {
				continue;
			}

			if (queryGraph->getNumberOfVertexes() < GlobalConstant::G_MINIMUM_NUMBER_MCS_VERTEX_RATIO * dataGraph -> getNumberOfVertexes()) {
				continue;
			}

			/*
			 * Already relationship filter
			 */
			if ((*patternContainmentMap)[queryGraph->graphId].descendent.find(dataGraph->graphId) != (*patternContainmentMap)[queryGraph->graphId].descendent.end()) {
				continue;
			}

			vector<vector<int>> mappings;
			if (subgraphIsomorphism.isSubgraphIsomorphic(dataGraph, queryGraph, &mappings, resultFile)) {

				(*patternContainmentMap)[queryGraph->graphId].descendent.insert(dataGraph->graphId);
				(*patternContainmentMap)[queryGraph->graphId].descendent.insert((*patternContainmentMap)[dataGraph->graphId].descendent.begin(), (*patternContainmentMap)[dataGraph->graphId].descendent.end());
				
				(*patternContainmentMap)[queryGraph->graphId].containmentRelationshipMappingLists.insert(std::pair<int, vector<vector<int>>>(dataGraph->graphId, mappings));
			}
		}
	}
}

void PCMBuilder::detectCommonSubgraphs() {
	/*
	 * Call maximal clique funciton to get groups of queries
	 */
	findAllMaximalCliques -> findMaxclique(tlsGraphMatrix, &similiarQueryGroups, patternContainmentMap);
	
	// printSimilarGroups(similiarQueryGroups); //@debug
	/*
	 * Compute a MCS for each group
	 */
	for (std::vector<vector<int>>::iterator similarQueryGroupIterator = similiarQueryGroups.begin(); similarQueryGroupIterator != similiarQueryGroups.end(); similarQueryGroupIterator++) {

		//cout << "Group: "<< similarQueryGroupIterator ->size() << endl; // @debug
		for (int pairQuery = 0; pairQuery < similarQueryGroupIterator->size();) {
			AdjacenceListsGRAPH mcsGraph;
			if (pairQuery +1 < similarQueryGroupIterator->size()) {
				computeMaximumCommonSubgraph(&(*queryGraphVector)[(*similarQueryGroupIterator)[pairQuery]], &(*queryGraphVector)[(*similarQueryGroupIterator)[pairQuery + 1]], &mcsGraph);
				if (mcsGraph.getNumberOfVertexes() > 13) {
					bool qualifiedToBeAdded = true;
					mcsGraph.buildLabelVertexList();
					mcsGraph.buildVertexLabelEdgeList();
					for (int i = 0; i < queryGraphVector->size(); i++) {
						if (graphIsomorphism.isGraphIsomorphic(&mcsGraph, &(*queryGraphVector)[i], resultFile)) {
							qualifiedToBeAdded = false;
							break;
						}
					}
					for (int i = 0; i < newlyGeneratedGraphVector->size(); i++) {
						if (graphIsomorphism.isGraphIsomorphic(&mcsGraph, &(*newlyGeneratedGraphVector)[i], resultFile)) {
							qualifiedToBeAdded = false;
							break;
						}
					}

					if (!qualifiedToBeAdded) {
						continue;
					}
					mcsGraph.graphId = queryGraphVector->size() + newlyGeneratedGraphVector->size();
					mcsGraph.buildComboGraph(GlobalConstant::G_MAXIMUM_WIDTH_COMBO);
					patternContainmentMap->push_back(PCM_Node(mcsGraph.graphId, true));

					newlyGeneratedGraphVector->push_back(mcsGraph);

					for (int i = 0; i < 2; i++) {
						vector<vector<int>> mappings;
						subgraphIsomorphism.isSubgraphIsomorphic(&(*queryGraphVector)[(*similarQueryGroupIterator)[pairQuery + i]], &mcsGraph, &mappings, resultFile);
						if (mappings.size() > 0) {
							// TODO, this is a bug. The mappings should always have values, however, it sometimes empty. fix it later
							(*patternContainmentMap)[mcsGraph.graphId].descendent.insert((*similarQueryGroupIterator)[pairQuery + i]);
							(*patternContainmentMap)[mcsGraph.graphId].containmentRelationshipMappingLists.insert(std::pair<int, vector<vector<int>>>((*similarQueryGroupIterator)[pairQuery + i], mappings));
						}
					}
				}
			}
			pairQuery += 2;
		}
	}
}

/*
 * graph1 act like the relay mcs graph, which can be discarded
 */
void PCMBuilder::computeMaximumCommonSubgraph(AdjacenceListsGRAPH * graph1, AdjacenceListsGRAPH * graph2, AdjacenceListsGRAPH * mcsGraph) {

	/*
	* Each vertex in the product graph represents a edge-pair of graph2
	*/
	map<int, std::pair<std::pair<int, int>, std::pair<int, int>>> productGraphVertexList;

	int numberOfVertex = 0;
	std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> * labelEdgeList1 = graph1->getVertexLabelsEdgeList();
	std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> * labelEdgeList2 = graph2->getVertexLabelsEdgeList();

	/*****************************************************
	* Build the compativity graph
	****************************************************/
	for (std::map<std::pair<int, int>, std::vector<std::pair<int, int>>>::iterator labelEdgeIterator1 = labelEdgeList1->begin(); labelEdgeIterator1 != labelEdgeList1->end(); labelEdgeIterator1++) {
		std::map<std::pair<int, int>, std::vector<std::pair<int, int>>>::iterator labelEdgeIterator2 = labelEdgeList2->find(labelEdgeIterator1->first);
		if (labelEdgeIterator2 != labelEdgeList2->end()) {
			for (unsigned int i = 0; i<labelEdgeIterator1->second.size(); i++) {
				for (unsigned int j = 0; j<labelEdgeIterator2->second.size(); j++) {

					productGraphVertexList.insert(std::pair<int, std::pair<std::pair<int, int>, std::pair<int, int>>>(numberOfVertex, std::pair<std::pair<int, int>, std::pair<int, int>>(labelEdgeIterator1->second[i], labelEdgeIterator2->second[j])));

					numberOfVertex++;
				}
			}
		}
	}

	bool connected = false;
	int ** productGraph = new int*[numberOfVertex];
	for (int i = 0; i < numberOfVertex; i++) {
		productGraph[i] = new int[numberOfVertex];
		for (int j = 0; j < numberOfVertex; j++) {
			if (i == j) {
				productGraph[i][j] = 1;
			}else{
				productGraph[i][j] = 0;
			}
		}
	}


	for (map<int, std::pair<std::pair<int, int>, std::pair<int, int>>>::iterator vertexIterator = productGraphVertexList.begin(); vertexIterator != productGraphVertexList.end(); vertexIterator++) {

		const std::pair<int, int> & e1 = vertexIterator->second.first;
		const std::pair<int, int> & e2 = vertexIterator->second.second;

		map<int, std::pair<std::pair<int, int>, std::pair<int, int>>>::iterator innerVertexIterator = vertexIterator;
		innerVertexIterator++;

		for (; innerVertexIterator != productGraphVertexList.end(); innerVertexIterator++) {
			
			const std::pair<int, int> & e3 = innerVertexIterator->second.first;
			const std::pair<int, int> & e4 = innerVertexIterator->second.second;

			if (e1 == e3 || e2 == e4) {
				continue;
			}

			/*
			* Check connectivity
			*/
			int connect_1 = 0, connect_2 = 0;
			if (e1.first == e3.first) {
				connect_1 = 1;
			}
			else if (e1.second == e3.first) {
				connect_1 = 2;
			}
			else if (e1.second == e3.second) {
				connect_1 = 3;
			}
			if (e2.first == e4.first) {
				connect_2 = 1;
			}
			else if (e2.second == e4.first) {
				connect_2 = 2;
			}
			else if (e2.second == e4.second) {
				connect_2 = 3;
			}

			if (connect_1 == connect_2) {
				productGraph[vertexIterator->first][innerVertexIterator->first] = 1;
				productGraph[innerVertexIterator->first][vertexIterator->first] = 1;
			}
		}
	}

	// printProductGraph(productGraph, numberOfVertex); 	//@debug
	/*
	 * call maximal clique detection algorithm
	 */
	FindMaximumClique findMaximumClique(productGraph, numberOfVertex);
	int * maximumClique = new int[numberOfVertex];
	int maximumCliqueSize = 0;
	findMaximumClique.findMaxclique(maximumClique, &maximumCliqueSize);

	//printMaximalClique(maximumClique, maximumCliqueSize); 	//@debug
	
	/*
	* construct a graph based on the clique result
	* iterate each node in the product graph == iterate each edge of the common subgraph of graph1 (graph2)
	*/
	mcsGraph->clear();
	map<int, int> subIdToMCSGraphVertex;
	map<int, int>::iterator idIteratorSource, idIteratorDesti;
	int newIdSouce, newIdDesti;
	int numberOfCliqueVertex = 0;

	for (int j = 0; j < maximumCliqueSize; j++) {
		map<int, std::pair<std::pair<int, int>, std::pair<int, int>>>::iterator edgePair = productGraphVertexList.find(maximumClique[j]);

		idIteratorSource = subIdToMCSGraphVertex.find(edgePair->second.second.first);
		if (idIteratorSource == subIdToMCSGraphVertex.end()) {
			subIdToMCSGraphVertex.insert(std::pair<int, int>(edgePair->second.second.first, numberOfCliqueVertex));

			AdjacenceListsGRAPH::Vertex vertex(numberOfCliqueVertex, graph2->getVertexByVertexId(edgePair->second.second.first).label);
			mcsGraph->insert(vertex);
			newIdSouce = numberOfCliqueVertex;
			numberOfCliqueVertex++;
		}
		else {
			newIdSouce = idIteratorSource->second;
		}

		idIteratorDesti = subIdToMCSGraphVertex.find(edgePair->second.second.second);
		if (idIteratorDesti == subIdToMCSGraphVertex.end()) {
			subIdToMCSGraphVertex.insert(std::pair<int, int>(edgePair->second.second.second, numberOfCliqueVertex));

			AdjacenceListsGRAPH::Vertex vertex(numberOfCliqueVertex, graph2->getVertexByVertexId(edgePair->second.second.second).label);
			mcsGraph->insert(vertex);
			newIdDesti = numberOfCliqueVertex;
			numberOfCliqueVertex++;
		}
		else {
			newIdDesti = idIteratorDesti->second;
		}

		AdjacenceListsGRAPH::Edge edge(newIdSouce, newIdDesti, 0);
		mcsGraph->insert(edge);
	}
	delete[] maximumClique;
	/*
	 * Make sure the MCS is a connected graph
	 */
	mcsGraph->buildDFSTraversalOrder();
	if (mcsGraph->getDFSTraversalOrder()->size() < mcsGraph->getNumberOfVertexes()) {
		std::map<int, int> mcsVertexIndexMap;
		int * mcsVertexLabels = new int[mcsGraph->getDFSTraversalOrder()->size()];
		int mcsVertexId;
		for (int i = 0; i < mcsGraph->getDFSTraversalOrder()->size(); i++) {
			mcsVertexId = mcsGraph->getDFSTraversalOrder()->at(i)->id;
			mcsVertexIndexMap.insert(std::pair<int, int>(mcsVertexId, i));
			mcsVertexLabels[i] = mcsGraph->getVertexAddressByVertexId(mcsVertexId)->label;
		}

		std::vector<AdjacenceListsGRAPH::Edge> largerConnectedMcsEdges;
		for (std::vector<AdjacenceListsGRAPH::Edge>::iterator mcsEdgeIterator = mcsGraph->getEdgeList()->begin(); mcsEdgeIterator != mcsGraph->getEdgeList()->end(); mcsEdgeIterator++) {
			std::map<int, int>::iterator mcsSourceVertexIndexIterator = mcsVertexIndexMap.find(mcsEdgeIterator->source);
			std::map<int, int>::iterator mcsDestVertexIndexIterator = mcsVertexIndexMap.find(mcsEdgeIterator->destination);
			if (mcsSourceVertexIndexIterator != mcsVertexIndexMap.end() && mcsDestVertexIndexIterator != mcsVertexIndexMap.end()) {
				largerConnectedMcsEdges.push_back(AdjacenceListsGRAPH::Edge(mcsSourceVertexIndexIterator->second, mcsDestVertexIndexIterator->second));
			}
		}
		mcsGraph->clear();
		for (int i = 0; i < mcsVertexIndexMap.size(); i++) {
			AdjacenceListsGRAPH::Vertex vertex(i, mcsVertexLabels[i]);
			mcsGraph->insert(vertex);
		}
		for (std::vector<AdjacenceListsGRAPH::Edge>::iterator largerMCSEdgeIter = largerConnectedMcsEdges.begin(); largerMCSEdgeIter != largerConnectedMcsEdges.end(); largerMCSEdgeIter++) {
			mcsGraph->insert(AdjacenceListsGRAPH::Edge(largerMCSEdgeIter->source, largerMCSEdgeIter->destination, 0));
		}
		mcsGraph->buildDFSTraversalOrder();

		delete[] mcsVertexLabels;
	}


	/*
	 * release memoery
	 */
	for (int i = 0; i < numberOfVertex; i++) {
		delete[] productGraph[i];
	}
	delete [] productGraph;
}

void PCMBuilder::printProductGraph(int ** productGraph, int size) {
	(*resultFile) << "The product graph is:" << endl;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			(*resultFile) << productGraph[i][j] << " ";
		}
		(*resultFile) << endl;
	}
}
void PCMBuilder::printMaximalClique(int * maximumClique, int maximumCliqueSize) {
	(*resultFile) << "The maximal clique size is:" << maximumCliqueSize  << endl;
	for (int i = 0; i < maximumCliqueSize; i++) {
		(*resultFile) << maximumClique[i] << " ";
	}
	(*resultFile) << endl;
}
void PCMBuilder::printSimilarGroups(std::vector<vector<int>> & similiarQueryGroups) {
	for (std::vector<vector<int>>::iterator similarQueryGroupIterator = similiarQueryGroups.begin(); similarQueryGroupIterator != similiarQueryGroups.end(); similarQueryGroupIterator++) {
		(*resultFile) << "Clique size: " << similarQueryGroupIterator->size() << endl;
		for (int i = 0; i < similarQueryGroupIterator->size(); i++) {
			(*resultFile) << " " << (*similarQueryGroupIterator)[i];
		}
		(*resultFile) << endl;
	}
	(*resultFile) << "TLS: Total cliques: " << similiarQueryGroups.size() << endl;
}

/*
 * 1. Do the transitive reduction
 * 2. Set up mapping vertices set helper variables
 */
void PCMBuilder::formatPCM(){
	//1. Do the transitive reduction
	bool onlyDirectChildren;
	for(unsigned int x=0; x<(*patternContainmentMap).size();x++){
		if((*patternContainmentMap)[x].isHidden){				
			continue;
		}
		for(std::set<int>::iterator zIterator = (*patternContainmentMap)[x].descendent.begin(); zIterator != (*patternContainmentMap)[x].descendent.end(); zIterator++){
			if((*patternContainmentMap)[*zIterator].isHidden){				
				continue;
			}

			onlyDirectChildren = true;
			for(std::set<int>::iterator yIterator = (*patternContainmentMap)[x].descendent.begin(); yIterator != (*patternContainmentMap)[x].descendent.end(); yIterator++){
				if(*zIterator == *yIterator || (*patternContainmentMap)[*yIterator].isHidden){
					continue;
				}
				if((*patternContainmentMap)[*yIterator].descendent.find(*zIterator) != (*patternContainmentMap)[*yIterator].descendent.end()){
					onlyDirectChildren = false;
					break;
				}
			}
			if(onlyDirectChildren){
				(*patternContainmentMap)[x].children.push_back(*zIterator);
				(*patternContainmentMap)[*zIterator].parent.push_back(x);
			}
		}
	}

	//2. Set up mapping vertices set helper variables
	for(std::vector<PCM_Node>::iterator pcmNodeIterator = patternContainmentMap->begin(); pcmNodeIterator != patternContainmentMap->end(); pcmNodeIterator++){
		for(std::map<int, std::vector<std::vector<int>>>::iterator mappingListsIterator = pcmNodeIterator->containmentRelationshipMappingLists.begin(); mappingListsIterator != pcmNodeIterator->containmentRelationshipMappingLists.end(); mappingListsIterator++){
			std::set<int> mappingQueryVerticesSet;
			for(vector<vector<int>>::iterator mappingListIterator = mappingListsIterator->second.begin(); mappingListIterator != mappingListsIterator->second.end(); mappingListIterator++){
				for(vector<int>::iterator mappingVertexIterator = mappingListIterator->begin(); mappingVertexIterator != mappingListIterator->end(); mappingVertexIterator++){
					mappingQueryVerticesSet.insert(*mappingVertexIterator);	
				}
			}
			pcmNodeIterator->containmentRelationshipMappingSet.insert(std::pair<int, std::set<int>>(mappingListsIterator->first, mappingQueryVerticesSet));	
		}
	}
}


void printDetail(PCM_Node & pcmNode, std::vector<PCM_Node> * patternContainmentMap, std::ofstream * resultFile){
	if(pcmNode.equivalentGraph.size() > 0){
		(*resultFile) <<"Equivalent graph: ";
		for(std::vector<int>::iterator iterator = pcmNode.equivalentGraph.begin(); iterator != pcmNode.equivalentGraph.end();iterator++){
			(*resultFile) <<(*iterator)<<((*patternContainmentMap)[*iterator].isGeneratedQueryGraph == false ? "(*)" : "(N)")<<",";
		}
		(*resultFile) <<endl;
	}
	if(pcmNode.descendent.size() > 0){
		(*resultFile) <<"Descendent graph: ";
		for(std::set<int>::iterator iterator = pcmNode.descendent.begin(); iterator != pcmNode.descendent.end();iterator++){
			(*resultFile) <<(*iterator)<<((*patternContainmentMap)[*iterator].isGeneratedQueryGraph == false ? "(*)" : "(N)")<<",";
		}
		(*resultFile) <<endl;
	}
	if(pcmNode.children.size() > 0){
		(*resultFile) <<"Children graph: ";
		for(std::vector<int>::iterator iterator = pcmNode.children.begin(); iterator != pcmNode.children.end();iterator++){
			(*resultFile) <<(*iterator)<<((*patternContainmentMap)[*iterator].isGeneratedQueryGraph == false ? "(*)" : "(N)")<<",";
		}
		(*resultFile) <<endl;
	}
	if(pcmNode.parent.size() > 0){
		(*resultFile) <<"Parent graph: ";
		for(std::vector<int>::iterator iterator = pcmNode.parent.begin(); iterator != pcmNode.parent.end();iterator++){
			(*resultFile) <<(*iterator)<<((*patternContainmentMap)[*iterator].isGeneratedQueryGraph == false ? "(*)" : "(N)")<<",";
		}
		(*resultFile) <<endl;
	}
}

void PCMBuilder::showPCM(){

	int numberOfPcmEdges = 0;
	int numberOfEquivalentNodes = 0;
	
	for (std::vector<PCM_Node>::iterator pcmNodeIterator = patternContainmentMap->begin(); pcmNodeIterator != patternContainmentMap->end(); pcmNodeIterator++) {
		if (!pcmNodeIterator -> isHidden) {
			numberOfPcmEdges += pcmNodeIterator->children.size();
			numberOfEquivalentNodes += pcmNodeIterator->equivalentGraph.size();
			/*(*resultFile) << endl << "Graph" << pcmNodeIterator->graphId << endl;
			printDetail(*pcmNodeIterator, patternContainmentMap, resultFile);*/
		}
	}

	int newlyCommonGraphs = patternContainmentMap->size() - queryGraphVector->size();

	(*resultFile) << endl << "****PCM Information: " << endl;
	(*resultFile) << "PCM edges: "<< numberOfPcmEdges << endl;
	(*resultFile) << "PCM newly nodes : " << newlyCommonGraphs << endl;
	(*resultFile) << "PCM query nodes : " << queryGraphVector->size() << endl;
	(*resultFile) << "PCM Euivalent nodes : " << numberOfEquivalentNodes << endl;
}
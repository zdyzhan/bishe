#include"MQO.h"
#include"TimeUtility.h"
#include<queue>

using namespace std;


MQO::MQO(){}

MQO::MQO(std::vector<AdjacenceListsGRAPH> * dataGraphVector, std::vector<AdjacenceListsGRAPH> * pQueryGraphVector, std::ofstream * pResultFile) {
	
	resultFile = pResultFile;

	if (dataGraphVector->size() > 0) {
		dataGraph = &(*dataGraphVector)[0];
	}

	queryGraphVector = pQueryGraphVector;

	/*
	* use a matrix to save the tls information between queries
	*/
	tlsGraphMatrix = new int*[queryGraphVector->size()];
	for (unsigned int i = 0; i<queryGraphVector->size(); i++) {
		tlsGraphMatrix[i] = new int[queryGraphVector->size()];
		for (unsigned int j = 0; j<queryGraphVector->size(); j++) {
			tlsGraphMatrix[i][j] = 0;
		}
		tlsGraphMatrix[i][i] = 1;
	}


	/*
	* create pcm node for each query graph
	*/
	for (unsigned int i = 0; i<queryGraphVector->size(); i++) {
		patternContainmentMap.push_back(PCM_Node(i, false));
	}


	tlsGraph = new TLSGraph(tlsGraphMatrix, queryGraphVector, resultFile);
	pcmBuilder = new PCMBuilder(tlsGraphMatrix, &patternContainmentMap, queryGraphVector, &newlyGeneratedGraphVector, resultFile);
	recurParEmbConstruct = new RecurParEmbConstruct(dataGraph, queryGraphVector, &newlyGeneratedGraphVector, &numberOfEmbeddings, &comboLinkedLists, &patternContainmentMap);
}

MQO::~MQO(){

	for (unsigned int i = 0; i < queryGraphVector->size(); ++i) {
		delete[] tlsGraphMatrix[i];
	}
	delete[] tlsGraphMatrix;



	delete tlsGraph;
	delete pcmBuilder;
	delete recurParEmbConstruct;

	patternContainmentMap.swap(std::vector<PCM_Node>());
	newlyGeneratedGraphVector.swap(std::vector<AdjacenceListsGRAPH>());
	numberOfEmbeddings.swap(std::vector<int>());
	comboLinkedLists.swap(std::vector<std::vector<EmbeddingNode *>>());
	joiningOrder.swap(vector<int>());
	uncoveredEdges.swap(set<std::pair<int, int>>());


}



void MQO::buildPCM(){

	tlsGraph->buildTLSGraph();

	//tlsGraph->showTLSInvertedIndex(); //@debug

	//tlsGraph->showTLSGraphMatrix(); //@debug

	pcmBuilder->hideIsomorphicQueries();

	pcmBuilder->buildContainmentRelations();

	tlsGraph->applyRelativeSimilarityRatio(&patternContainmentMap);

	// tlsGraph -> showTLSGraphMatrix(); //@debug

	//@debug
	pcmBuilder->detectCommonSubgraphs();

	pcmBuilder->formatPCM();

	pcmBuilder->showPCM();
}



/*
 * 1. Starting from the smallest number of embeddings
 * 2. Next parent should has common data vertex with the previous one.
 * 3. The less embeddings the better
 * @return false if any one its parent doesn't have any embedding. 
 */
bool MQO::computeJoiningOrder(){

	joiningOrder.swap(vector<int>());
	uncoveredEdges.swap(set<std::pair<int, int>>());

	if(patternContainmentMap[queryGraphId].parent.size() == 0){
		return true;
	}

	set<int> queryVertexCover;
	/*
	 * Get the uncovered edge
	 */
	for(std::vector<AdjacenceListsGRAPH::Edge>::iterator edgeIterator = queryGraph->getEdgeList()->begin(); edgeIterator != queryGraph->getEdgeList()->end(); edgeIterator++){
		uncoveredEdges.insert(std::pair<int, int>(edgeIterator->source, edgeIterator->destination));
	}

	int thisParentQueryGraphId;
	AdjacenceListsGRAPH * thisParentQueryGraph = NULL;

	/*
	 * 1. Use a flags to indicate whether this parent has been added or not
	 * 2. Retrive the parent vertex mappings to prevent retrieving them in iteration
	 */
	int* flags = new int[patternContainmentMap[queryGraphId].parent.size()]();
	set<int> **  parentVertexMappingLists = new set<int>* [patternContainmentMap[queryGraphId].parent.size()];
	for(int i=0; i<patternContainmentMap[queryGraphId].parent.size(); i++){
		flags[i] = false;
		parentVertexMappingLists[i] = &(patternContainmentMap[patternContainmentMap[queryGraphId].parent[i]].containmentRelationshipMappingSet.find(queryGraphId)->second);
	}

	/*
	 * Starting from the smallest number of embeddings
	 */
	int startParentIndex = -1;
	int MIN_NUMBER = INT_MAX;
	for(int i=0; i<patternContainmentMap[queryGraphId].parent.size(); i++){
		if(numberOfEmbeddings[patternContainmentMap[queryGraphId].parent[i]] <= 0){
			delete [] flags;
			delete [] parentVertexMappingLists;
			return false;
		}
		if(numberOfEmbeddings[patternContainmentMap[queryGraphId].parent[i]] < MIN_NUMBER){
			MIN_NUMBER = numberOfEmbeddings[patternContainmentMap[queryGraphId].parent[i]];
			startParentIndex = i;
		}
	}

	flags[startParentIndex] = true;
	/*
	 * set mapping information for this parent graph. update the query edge cover
	 */
	thisParentQueryGraphId = patternContainmentMap[queryGraphId].parent[startParentIndex];
	if(thisParentQueryGraphId < queryGraphVector->size()){
		thisParentQueryGraph = &(*queryGraphVector)[thisParentQueryGraphId];	
	}else{
		thisParentQueryGraph = &newlyGeneratedGraphVector[thisParentQueryGraphId - queryGraphVector->size()];
	}

	joiningOrder.push_back(thisParentQueryGraphId);

	/*
	 * Update the query vertex and edge covers
	 */
	PCM_Node & parentNode = patternContainmentMap[thisParentQueryGraphId];
	std::map<int, std::vector<std::vector<int>>>::iterator mappingListsIterator = parentNode.containmentRelationshipMappingLists.find(queryGraphId);
	for(std::vector<AdjacenceListsGRAPH::Edge>::iterator edgeIterator = thisParentQueryGraph->getEdgeList()->begin(); edgeIterator != thisParentQueryGraph->getEdgeList()->end(); edgeIterator++){
		for(vector<vector<int>>::iterator mappingListIterator = mappingListsIterator->second.begin(); mappingListIterator != mappingListsIterator->second.end(); mappingListIterator++){
			
			uncoveredEdges.erase(std::pair<int, int>((*mappingListIterator)[edgeIterator->source], (*mappingListIterator)[edgeIterator->destination]));
			uncoveredEdges.erase(std::pair<int, int>((*mappingListIterator)[edgeIterator->destination], (*mappingListIterator)[edgeIterator->source]));

			queryVertexCover.insert((*mappingListIterator)[edgeIterator->source]);
			queryVertexCover.insert((*mappingListIterator)[edgeIterator->destination]);
		}
	}

	/*
	 * Handle other parents
	 */
	int nextParentIndex = -1, nextParentIndexConnected = -1;
	int MINI_NUMBER_CONNECTED = -1;
	while(queryVertexCover.size() < queryGraph->getNumberOfVertexes() && joiningOrder.size() < patternContainmentMap[queryGraphId].parent.size()){

		MIN_NUMBER = INT_MAX;
		MINI_NUMBER_CONNECTED = INT_MAX;
		
		for(int i=0; i<patternContainmentMap[queryGraphId].parent.size(); i++){
			if(!flags[i]){

				if(numberOfEmbeddings[patternContainmentMap[queryGraphId].parent[i]] < MIN_NUMBER){
					MIN_NUMBER = numberOfEmbeddings[patternContainmentMap[queryGraphId].parent[i]];
					nextParentIndex = i;
				}

				bool connected = false;
				for(std::set<int>::iterator singleMappingIterator = (*parentVertexMappingLists[i]).begin(); singleMappingIterator != (*parentVertexMappingLists[i]).end();singleMappingIterator++){
					if(queryVertexCover.find(*singleMappingIterator) != queryVertexCover.end()){
						connected = true;
						break;
					}
				}

				if(connected && MIN_NUMBER < MINI_NUMBER_CONNECTED){
					MINI_NUMBER_CONNECTED = MIN_NUMBER;
					nextParentIndexConnected = i;	
				}
			}
		}

		if(nextParentIndexConnected != -1){
			nextParentIndex = nextParentIndexConnected;
		}
		flags[nextParentIndex] = true;

		/*
		 * Same as to the start query graph, update the query edge and vertex cover
		 */
		thisParentQueryGraphId = patternContainmentMap[queryGraphId].parent[nextParentIndex];
		if (thisParentQueryGraphId < queryGraphVector->size()) {
			thisParentQueryGraph = &(*queryGraphVector)[thisParentQueryGraphId];
		}
		else {
			thisParentQueryGraph = &newlyGeneratedGraphVector[thisParentQueryGraphId - queryGraphVector->size()];
		}

		joiningOrder.push_back(thisParentQueryGraphId);
		PCM_Node & parentNode = patternContainmentMap[thisParentQueryGraphId];

		std::map<int, std::vector<std::vector<int>>>::iterator mappingListsIterator = parentNode.containmentRelationshipMappingLists.find(queryGraphId);
		for(std::vector<AdjacenceListsGRAPH::Edge>::iterator edgeIterator = thisParentQueryGraph->getEdgeList()->begin(); edgeIterator != thisParentQueryGraph->getEdgeList()->end(); edgeIterator++){
			for(vector<vector<int>>::iterator mappingListIterator = mappingListsIterator->second.begin(); mappingListIterator != mappingListsIterator->second.end(); mappingListIterator++){
				
				uncoveredEdges.erase(std::pair<int, int>((*mappingListIterator)[edgeIterator->source], (*mappingListIterator)[edgeIterator->destination]));
				uncoveredEdges.erase(std::pair<int, int>((*mappingListIterator)[edgeIterator->destination], (*mappingListIterator)[edgeIterator->source]));

				queryVertexCover.insert((*mappingListIterator)[edgeIterator->source]);
				queryVertexCover.insert((*mappingListIterator)[edgeIterator->destination]);
			}
		}
	}

	
	delete [] flags;
	delete [] parentVertexMappingLists;
	return true;
}



void MQO::DFSTopological(int graphId) {

	//(*resultFile) << "Ready Compute->" << graphId << endl; //@debug

	/*******************************************************
	@core: this graph is ready to compute Subgraph isomorphism mapping
	*******************************************************/
	queryGraphId = graphId;
	if (graphId < queryGraphVector->size()) {
		queryGraph = &(*queryGraphVector)[graphId];
	}
	else {
		queryGraph = &newlyGeneratedGraphVector[graphId - queryGraphVector->size()];
	}

	/*
	* If there is no partial embedding of its parents, we no need to compute further
	*/
	if (computeJoiningOrder()) {
		recurParEmbConstruct->processing(queryGraph, &joiningOrder, &uncoveredEdges);
	}
	else {
		// (*resultFile) <<graphId<<" : No need to compute: no  partial embeddings from parents"<<endl; //@debug
	}

	/*******************************************************
	* @core: computing the subgraph isomorphism for this query finished,
	* then mark it as visited after computing
	*******************************************************/
	patternContainmentMap[graphId].isVisited = true;


	if (patternContainmentMap[graphId].children.size() == 0) {
		/*
		* We won't save the cache for the query who doesn't have any children. However, we can see the embedding here
		*/
		showEmbedding(graphId); 	//@debug
	}

	/*
	* After computing, increase computed children of its parents
	*/
	for (vector<int>::iterator parentIterator = patternContainmentMap[graphId].parent.begin(); parentIterator != patternContainmentMap[graphId].parent.end(); parentIterator++) {
		patternContainmentMap[*parentIterator].numberOfComputedChildren++;
		if (patternContainmentMap[*parentIterator].numberOfComputedChildren == patternContainmentMap[*parentIterator].children.size()) {
			// release the memory
			//(*resultFile) << "Release Memory->" << graphId  << " " << *parentIterator << " "<< numberOfEmbeddings[*parentIterator] << " "<< numberOfEmbeddings[graphId]  << " "<< (*queryGraphVector)[*parentIterator].getNumberOfVertexes() << endl; // @debug
			showEmbedding(*parentIterator);
			comboLinkedLists[*parentIterator].swap(std::vector<EmbeddingNode *>()); // release the momery of the caches of this query
		}
	}


	/*
	* for each of its children, we increase its number of computed parents.
	*/
	if (patternContainmentMap[graphId].children.size() > 0) {
		for (vector<int>::iterator childIterator = patternContainmentMap[graphId].children.begin(); childIterator != patternContainmentMap[graphId].children.end(); childIterator++) {
			if (!patternContainmentMap[*childIterator].isVisited) {
				patternContainmentMap[*childIterator].numberOfComputedParents++;
				if (patternContainmentMap[*childIterator].numberOfComputedParents != patternContainmentMap[*childIterator].parent.size()) {
					/*
					* Not all of its parents have been computed. It is meaningless still to update if all its parents have already been computed.
					*/
					changeParentsWeight(*childIterator);
				}
			}
		}

		int nextChildQuery = getNextGraphId(patternContainmentMap[graphId].children);
		while (nextChildQuery != -1) {
			DFSTopological(nextChildQuery);
			nextChildQuery = getNextGraphId(patternContainmentMap[graphId].children);
		}
	}

}


int MQO::getNextGraphId(vector<int> & candidates) {
	int index = -1;
	int maxWeight = -1;
	/*
	* given a list of candidates, get the one with the highest score.
	*/
	for (unsigned int i = 0; i<candidates.size(); i++) {
		if (!patternContainmentMap[candidates[i]].isVisited && patternContainmentMap[candidates[i]].numberOfComputedParents == patternContainmentMap[candidates[i]].parent.size()) {
			if (patternContainmentMap[candidates[i]].sameHeightPriority > maxWeight) {
				index = i;
				maxWeight = patternContainmentMap[candidates[i]].sameHeightPriority;
			}
		}
	}
	if (index == -1) {
		return -1;
	}
	else {
		return candidates[index];
	}
}


void MQO::changeParentsWeight(int graphId) {
	for (vector<int>::iterator parentIterator = patternContainmentMap[graphId].parent.begin(); parentIterator != patternContainmentMap[graphId].parent.end(); parentIterator++) {
		if (!patternContainmentMap[*parentIterator].isVisited) {
			patternContainmentMap[*parentIterator].sameHeightPriority++;
			changeParentsWeight(*parentIterator);
		}
	}
}

void MQO::orederedQueryProcessing(){
	/*
	* initialize cache embedding for all the PCM node
	*/
	for (unsigned int i = 0; i < patternContainmentMap.size(); i++) {

		numberOfEmbeddings.push_back(0);

		comboLinkedLists.push_back(vector<EmbeddingNode*>());

		AdjacenceListsGRAPH * tempQueryGraph = NULL;
		if (i < queryGraphVector->size()) {
			tempQueryGraph = &(*queryGraphVector)[i];
		}
		else {
			tempQueryGraph = &newlyGeneratedGraphVector[i - queryGraphVector->size()];
		}
		for (unsigned int j = 0; j < tempQueryGraph->getComboGraph()->size(); j++) {
			comboLinkedLists[i].push_back(NULL);
		}
	}

	(*resultFile) << endl << "******QUERY PROCESSING:" << endl; //@debug																																	 /*																															 */
	vector<int> roots;
	int nextGraphId;
	for (std::vector<PCM_Node>::iterator iterator = patternContainmentMap.begin(); iterator != patternContainmentMap.end(); iterator++) {
		/*
		* As the PCM contains all the graphs no matter whether the query is isomorphic to another or not.
		* the isHidden is a flag to indicate the isomorphic relationship.
		*/
		if (iterator->isHidden) {
			continue;
		}
		if (iterator->parent.size() == 0) {
			roots.push_back(iterator->graphId);
		}
	}
	while ((nextGraphId = getNextGraphId(roots)) != -1) {
		DFSTopological(nextGraphId);
	}
}


void MQO::showEmbedding(int graphId) {
	if(!patternContainmentMap[graphId].isGeneratedQueryGraph){
		cout << graphId << ": Real Query: Total number of embeddings found :" << numberOfEmbeddings[graphId] << " " << endl;
		//(*resultFile) <<graphId<<": Real Query: Total number of embeddings found :"<<numberOfEmbeddings[graphId]<<" "<<endl;
	}else{
		// cout << graphId << ": Generated Query: Total number of embeddings found :" << numberOfEmbeddings[graphId] << " " << endl;
		// (*resultFile) <<graphId<<": Generated Query: Total number of embeddings found :"<<numberOfEmbeddings[graphId]<<" "<<endl;
	}
	/*for(int i=0; i<embeddingCache[graphId].size(); i++){
		std::cout<<"{";
		for(int j = 0; j < embeddingCache[graphId][i].size(); j++){
			cout<<""<<j<<"->"<<embeddingCache[graphId][i][j]<<" , ";
		}
		cout<<"}"<<endl;	
	}*/
}
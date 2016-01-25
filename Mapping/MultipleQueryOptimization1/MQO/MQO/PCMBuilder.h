#ifndef PCM_BUILDER_H
#define PCM_BUILDER_H

#include<vector>
#include<string>
#include<map>
#include<list>
#include"QueryGraphIsomorphismTest.h"
#include"QuerySubgraphIsomorphismSearch.h"
#include"AdjacenceListsGraph.h"
#include"FindAllMaximalCliques.h"
#include"FindMaximumClique.h"

class FindAllMaximalCliques;

class PCM_Node {
public:
	std::set<int> descendent;
	// for each of its children, we keep a mapping list from the vertex-id of this graph to the vertex-ids of its children
	std::map<int, std::vector<std::vector<int>>> containmentRelationshipMappingLists;
	std::map<int, std::set<int>> containmentRelationshipMappingSet;

	std::vector<int> equivalentGraph;
	std::vector<int> children;
	std::vector<int> parent;

	int graphId;
	bool isGeneratedQueryGraph;
	bool isHidden;
	bool isVisited;

	int numberOfComputedParents;
	int numberOfComputedChildren;
	int sameHeightPriority;

	PCM_Node() {
	};

	PCM_Node(int pGraphId, bool pIsGeneratedQueryGraph) {
		graphId = pGraphId;
		isGeneratedQueryGraph = pIsGeneratedQueryGraph;
		isHidden = false;
		isVisited = false;
		numberOfComputedParents = 0;
		numberOfComputedChildren = 0;
		sameHeightPriority = 0;
	};

	~PCM_Node() {
		descendent.swap(std::set<int>());
		containmentRelationshipMappingLists.swap(std::map<int, std::vector<std::vector<int>>>());
		containmentRelationshipMappingSet.swap(std::map<int, std::set<int>>());
		equivalentGraph.swap(std::vector<int>());
		children.swap(std::vector<int>());
		parent.swap(std::vector<int>());
	}

};



class PCMBuilder{
private:

	QueryGraphIsomorphismTest graphIsomorphism;

	QuerySubgraphIsomorphismSearch subgraphIsomorphism;

	FindAllMaximalCliques * findAllMaximalCliques;

	std::vector<AdjacenceListsGRAPH> * queryGraphVector;

	std::vector<PCM_Node> * patternContainmentMap;

	std::vector<AdjacenceListsGRAPH> * newlyGeneratedGraphVector;

	int ** tlsGraphMatrix;
	
	std::vector<vector<int>> similiarQueryGroups;

	std::ofstream * resultFile;

public:

	PCMBuilder(int ** pTlsGraphMatrix, std::vector<PCM_Node> * pPatternContainmentMap, std::vector<AdjacenceListsGRAPH> * pQueryGraphVector, std::vector<AdjacenceListsGRAPH> * pNewlyGeneratedGraphVector, std::ofstream * pResultFile);

	~PCMBuilder();

	void detectCommonSubgraphs();

	void hideIsomorphicQueries();

	void buildContainmentRelations();

	/*
	* 1. Do the transitive reduction
	* 2. Build mappingVerticesSet helper variables from mappingVerticesList. We used a list and a set to save the mapping informations
	*/
	void formatPCM();

	void showPCM();

private:

	void computeMaximumCommonSubgraph(AdjacenceListsGRAPH * graph1, AdjacenceListsGRAPH * graph2, AdjacenceListsGRAPH * mcsGraph);
	
	void printProductGraph(int ** productGraph, int size);

	void printMaximalClique(int * maximumClique, int maximumCliqueSize);

	void printSimilarGroups(std::vector<vector<int>> & similiarQueryGroups);
};

#endif
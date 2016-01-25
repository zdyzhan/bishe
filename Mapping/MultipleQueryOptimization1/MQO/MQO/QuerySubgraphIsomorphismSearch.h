#ifndef QUERY_SUBGRAPH_ISOMORPHISM_H
#define QUERY_SUBGRAPH_ISOMORPHISM_H

#include"TraversalAlgorithm.h"
#include"ComboLinkedLists.h"
#include<vector>
#include<map>


using namespace std;



class QuerySubgraphIsomorphismSearch {
private:

	class NECNode {
	public:
		int vertexId;
		std::vector<int> childList;
		int parent;
		int id;
		int label;
		~NECNode() {
			childList.swap(std::vector<int>());
		}
	};

private: 

	/*
	* If enoughRecursiveCalls has been called before any embedding founded. We need to abort the process and just return false
	*/
	int enoughRecursiveCalls;
	/*
	 * For TurboIso, we need to set another threshold for the ExploreCR call numbers
	 */
	int exploreCRRecursiveCalls;

private:

	std::ofstream * resultFile;
	AdjacenceListsGRAPH * dataGraph;
	AdjacenceListsGRAPH * queryGraph;
	vector<vector<int>> * storageMappings;

private:

	bool enoughMappingsFounded;

	int * embedding; 
	std::map<int, int> inverseEmbedding; 

	int startQueryVertex;
	AdjacenceListsGRAPH::Vertex * startQueryVertexAddress;

	/*
	* The key idea of the CR(candidate subregions) is "locality".
	* We have to use set<int> to storage each tuple in the CR, otherwise there can be some duplications causing error
	*/
	map<std::pair<int, int>, set<int>> CR;
	std::vector<NECNode> necTree;
	std::vector<int> queryMatchingSuquence; //vertexId, parentId; 

public:

	bool isSubgraphIsomorphic(AdjacenceListsGRAPH * pDataGraph, AdjacenceListsGRAPH * pQueryGraph, vector<vector<int>> * pStorageMappings, std::ofstream * pResultFile);

private:
	void isomorphismSearch();

	bool QuerySubgraphIsomorphismSearch::globalStatisticFilter();


	/* QuerySubgraphIsomorphismSearchBoosted query processing algorithm */
	int chooseStartVertex();
	void rewriteToNecTree();




	/*
	* The key idea for exploring the CR is locality, which is to save the subregion into the CR table and then do further verification
	* @param parentMappedNecTree is the data vertex mapped to u's parent in nec tree
	* @param subRegionCandidates is the set of candidates of u within this subRegion
	*/
	bool exploreCR(int u, int parentMappedNecTree, vector<int> & subRegionCandidates);


	void computeMatchingOrder();
	void devideNoTreeEdges(int u, int numberOfNoTreeEdges, float * queryVertexScore);
	void getOrderByBFSScore(int u, float * queryVertexScore, vector<int> & nextVertex); // after get the scroe for each query vertex, get the matching order by BFS

																						/* QuerySubgraphIsomorphismSearchBoosted subgraph isomorphism search */
	void updateState(int u, int v);
	void subgraphSearch();
	bool refineCandidates(AdjacenceListsGRAPH::Vertex & u, const int & v);
	bool degreeFilter(int u, int v);
	bool isJoinable(int u, int v);
	NECNode * nextQueryVertex();
	void restoreState(int u, int v);



	//@common
	/* Common Utility Function  */
	void showEmbedding();

};









#endif
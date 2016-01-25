#ifndef TURBO_ISO_BOOSTED_H
#define TURBO_ISO_BOOSTED_H


#include<vector>
#include<stack>
#include<queue>
#include<string>
#include<map>
#include<set>

#include"ComboLinkedLists.h"
#include"AdjacenceListsGraph.h"
#include"TraversalAlgorithm.h"


using namespace std;



class TurboIsoBoosted {

private:
	/*
	* As we don't consider the neighbourhood equivalent class for the query graphs here.
	* Therefore, the NECNode is only a common node in the BFS tree.
	*/
	class NECNode {
	public:
		int vertexId;
		vector<int> childList;
		int parent;
		int id;
		int label;
	};

private:

	AdjacenceListsGRAPH * hyperGraph;
	std::vector<int> * numberOfEmbeddings;
	ComboLinkedLists * comboLinkedLists;

private:
	AdjacenceListsGRAPH * queryGraph;

	int * partialEmbedding; /// mapping M: V(q) -> V(g), the arrayId is the query vertex id
	bool needFormatCache;

private:

	std::map<int, vector<int>* > scRootCandidates;
	std::map<int,std::stack<int>> inverseEmbedding; ///inverse mapping W: V(g) -> queryList{u1,u2...}
	int mappedQueryVertexSize;



	//@TurboIso
	int startQueryVertex;
	AdjacenceListsGRAPH::Vertex * startQueryVertexAddress;

	/*
	 * The key idea of the CR(candidate subregions) is "locality".
	 * We have to use set<int> to storage each tuple in the CR, otherwise there can be some duplications causing error
	 */
	map<std::pair<int,int>, set<int>> CR;
	std::vector<NECNode> necTree;
	std::vector<int> queryMatchingSuquence; //vertexId, parentId; 

	/*
	 * For turboIso, we only build the DRT table for the starting query vertex.
	 * the d_table contains all the DRT units built over the candidates of this query vertex
	 * the superStartingandidates are the DRT units in the d_table which are not qdc-contained by any other DRT unit.
	 */
	std::queue<int> startingQueryVertexCandidates;

private:
	/*
	* If enoughRecursiveCalls has been called before any embedding founded. We need to abort the process and just return false
	*/
	int enoughRecursiveCalls;
	/*
	* For TurboIso, we need to set another threshold for the ExploreCR call numbers
	*/
	int exploreCRRecursiveCalls;

	bool needReturn;

public:

	TurboIsoBoosted();

	TurboIsoBoosted(AdjacenceListsGRAPH * pHyperGraph, std::vector<int> * pNumberOfEmbeddings, ComboLinkedLists * pComboLinkedLists);

	void setParameters(AdjacenceListsGRAPH * pQueryGraph, int * pPartialEmbedding, bool pNeedToSaveCache);

	void execute();

private:

	/**
	* Subgraph isomorphism
	*/
	void clean_before();
	void clean_after();

	/* TurboIsoBoosted query processing algorithm */
	int chooseStartVertex();
	void rewriteToNecTree();
	void filterCandidates();


	/*
	 * The key idea for exploring the CR is locality, which is to save the subregion into the CR table and then do further verification
	 * @param parentMappedNecTree is the data vertex mapped to u's parent in nec tree
	 * @param subRegionCandidates is the set of candidates of u within this subRegion
	 */
	bool exploreCR(int u, int parentMappedNecTree, vector<int> & subRegionCandidates);


	void computeMatchingOrder();
	void devideNoTreeEdges(int u, int numberOfNoTreeEdges, float * queryVertexScore);
	void getOrderByBFSScore(int u, float * queryVertexScore, vector<int> & nextVertex); // after get the scroe for each query vertex, get the matching order by BFS

	/* TurboIsoBoosted subgraph isomorphism search */
	void updateState(AdjacenceListsGRAPH::Vertex * u, int v);
	void subgraphSearch();
		bool isJoinable(AdjacenceListsGRAPH::Vertex * u, AdjacenceListsGRAPH::Vertex * v);
	void restoreState(AdjacenceListsGRAPH::Vertex * u, int v);
	NECNode * nextQueryVertex();
	
	
	/**
	*	@unique for Boost
	*/
	void generateEquivalentEmbedding();
	void generateQDEquivalentEmbedding();
	void dynamicLoadCandidates(int dataVertexId);



	//@common
	/* Common Utility Function  */
	void showEmbedding();

};







#endif
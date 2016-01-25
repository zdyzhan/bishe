#ifndef ADJACENCE_LISTS_GRAPH
#define ADJACENCE_LISTS_GRAPH
#include<vector>
#include<set>
#include<map>
#include<iostream>
#include<fstream>

class AdjacenceListsGRAPH;



/*
 * The start is always smaller than the end
 * This is to make sure no duplications saved
 */
struct TLSequence{
public: 
	int start;
	int pivot;
	int end;

	TLSequence(){
		start = -1;
		pivot = -1;
		end = -1;
	}
};

bool operator<(const TLSequence& a, const TLSequence& b);

/*
 *  The node for combo vertices, each combo vertices may have up to 3 vertices. 
 */
struct ComboNODE{
public:
	/*
	 * Original query vertices
	 */
	std::vector<int> queryVertices;
	/*
	 * Combo vertex neighours, each neighbour is the id/index of the combo group.
	 */
	std::set<int> neighbours;
	/*
	 * DFS order parent. the id/index of the combo group.
	 */
	int parent;
	/*
	 * DFS order children. the id/index of the combo group.
	 */
	std::vector<int> children;
};


/*
 * Main graph definition
 */
class AdjacenceListsGRAPH{

public:
	struct node{
		int v;
		node * next;
		int edgeLabel;
		node(int x, node *t, int e){v = x; next = t;edgeLabel = e;};
	};
	struct Vertex {
	public:
		int id;
		int label;
		node * adj;
		int inDegree;   // if the graph is indirected graph, inDegree = outdegree
		int outDegree;  // if the graph is indirected graph, both outdegree and indegree are are the degree of the vertex

		Vertex(int id = -1, int label = -1) :id(id), label(label) {
			adj = 0; 
			inDegree = 0; 
			outDegree = 0; 
			boostInitialize();
		};
		
		/*
		 * The following two member will only be built when necessary by calling "buildVertexLabelVertexList"
	     */
		std::map<int,std::vector<int>> labelVertexList; // save a list of vertices for each neighbour labels
		std::set<int> labelSet;

		~Vertex(){
			labelVertexList.swap(std::map<int, std::vector<int>>());
			labelSet.swap(std::set<int>());
		}

	public:
		/*
		* @boostd
		* Used by S-relationships
		*/
		int isClique; // 0. Only one vertex 1. IsClique 2. is not a clique

		int sIndegree;
		int spongeIndegree; // is the same as the sInDegree but is used as a modifiable variable. 

		std::vector<int> sEquivalentVertexList; // only used by hypergraph, it is empty for the containment graph and the query graphs
		std::vector<int> sContainedVertexList;

		/*
		* @boostd
		* Used by QD-relationships
		*/
		int qdIndegree;
		int qdSpongeIndegree;
		bool qdeEquivalentHidden;
		std::vector<int> qdEquivalentVertexList;
		std::vector<int> qdContainedVertexList;

	private:
		void boostInitialize() {
			isClique = 0;
			sIndegree = 0;
			spongeIndegree = 0;
			qdeEquivalentHidden = false;
			qdIndegree = 0;
			qdSpongeIndegree = 0;
		}
	};
	struct Edge{
		int source,destination;
		int label;
		Edge(int sour=-1,int dest=-1, int label=-1):source(sour),destination(dest),label(label){};
	};


	typedef node* link;

	int graphId;

private:
	int numberOfVertex, numberOfEdges;
	bool digraph;
	std::vector<Vertex> vertexList;
	std::vector<Edge> edgeList;


	std::vector<Vertex *> DFSTraversalOrder;
	/*
	 * store all the distinct labels of this graph
	 * The labelset of the graph will be initialized when inserting all the vertices
	 */
	std::set<int> labelSet;
	/*
	 * invert list of label -> vertexes of a specific vertex
	 */
	std::map<int,std::vector<int>> labelVertexList;   


private:

	/*
	 * invert list of label pair -> edges of a specific pair of vertex labels
	 * each edge was represented by from-to int pair
	 */	  
	std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> vertexLabelsEdgeList; 

	/*
     * the triple connection set for this graph
	 * each TLS lable can has multiple sets of three vertice-combinations
	 */
	std::map<TLSequence, std::vector<std::vector<int>>> TLSequenceMap;
	int totalTLSNumber;


private:
	/*
	 * save the spanning tree of this graph
	 */
	std::vector<ComboNODE> comboGraph;
	short int ** comboGraphEdgeMatrix;
	std::vector<int> DFSComboVertexOrder;

private:

	/*
	 * @boostd
	 * inverted list of label -> containment roots list
	 */
	std::map<int, std::vector<int>> labelContainmentRootList; 

public:

	AdjacenceListsGRAPH(int numberOfVertex, bool digraph);
	AdjacenceListsGRAPH(bool digraph);
	AdjacenceListsGRAPH();


	bool directed() const;

	void insert(Vertex v);
	void insert(Edge e);
	void insertRandom(Edge e);


	std::vector<Vertex> * getVertexList();
	std::vector<Edge>* getEdgeList();


	/*
	 * Build the DFS traversal order
	 */
	void buildDFSTraversalOrder();

	/*
	 * build the inverted list of label->vertexes
	 */
	void buildLabelVertexList();

	/*
	* build the inverted list of vertex-label pairs -> edges
	* (A,B)-> (q1,q2), (q3, q4)
	 */
	void buildVertexLabelEdgeList();
	/* 
	* for each vertex build its labelset and for each vertex build vertex list according to the each label
	* For vertex u:
	 * A-> u1, u2
	 * B-> u4, u5
	 */
	void buildVertexLabelVertexList();


public:
	int degree(int);
	bool edge(int, int);

	int getNumberOfVertexes() const;
	int getNumberOfEdges() const;

	std::vector<Vertex *>  * getDFSTraversalOrder();
	std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> * getVertexLabelsEdgeList();
	std::map<int, std::vector<int>> * getLabelVertexList();
	std::set<int>* getLabelSet();

public:
	/*
	 * Build triple connection set for this graph
	 */
	void buildTLSequence();


	/*
	 * Build the spanning tree for this graph
	 */
	void buildComboGraph(int maximumWidth);

	void cleanComboInfo();

	std::map<TLSequence, std::vector<std::vector<int>>> * getTLSequenceMap();

	int getTotalTLSNumber();

	std::vector<ComboNODE> * getComboGraph();
	short int ** getComboGraphEdgeMatrix();
	std::vector<int> * getDFSComboVertexOrder();



	Vertex getVertexByVertexId(int vertexId);
	Vertex* getVertexAddressByVertexId(int vertexId);

public:
	
	std::map<int, std::vector<int>> * getLabelContainmentRootList();

public: 

	void printLabelEdges(std::ofstream * resultFile);

public:

	class adjIterator{
	protected:
		const AdjacenceListsGRAPH *G;
			int v;
			AdjacenceListsGRAPH::link t;
	public:
		adjIterator(const AdjacenceListsGRAPH* G, int v);
			AdjacenceListsGRAPH::link begin();
			AdjacenceListsGRAPH::link next();
			bool end();
	};

	void clear();

private:
	void insert_order_helper(AdjacenceListsGRAPH::Vertex* vertex, int w, int edgeLabel);
};





#endif
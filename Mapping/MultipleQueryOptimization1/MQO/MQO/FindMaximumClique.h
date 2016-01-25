#pragma once
#pragma once
#ifndef FIND_MAXIMUM_CLIQUE_H
#define FIND_MAXIMUM_CLIQUE_H

/*
* Given the MLS matrix, we find all the maximal cliques from the matrix
* Each clique is a group of queries grouped together.
* Notice that, in each group. we don't allow containment-relationship vertices.
*/
#include<vector>
#include"AdjacenceListsGraph.h"
#include"PCMBuilder.h"


class PCM_Node;





class FindMaximumClique
{
private:

	class nodeSet {

	public:
		int N;
		int *node;
		int size;

		nodeSet(void);
		~nodeSet(void);
		void nodeInit(int NInput);
		void add(int nodeA);
		void remove(void);
		void print(void);
	};

public:
	FindMaximumClique();
	FindMaximumClique(int ** pGraph, int graphSize);
	~FindMaximumClique();

	void findMaxclique(int * pMaximumClique, int * pMaximumCliqueSize);

private:
	int ** graph;
	int N;
	nodeSet compsub;
private:
	int * maximumClique;
	int * maximumCliqueSize;

	bool cliqueFounded;
private:
	void bkv2(int* pOld, int ne, int ce);
	void setMaximumClique();
};





#endif
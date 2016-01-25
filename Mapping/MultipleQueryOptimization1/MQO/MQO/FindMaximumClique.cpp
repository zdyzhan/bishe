#include <stdlib.h>
#include <stdio.h>
#include "FindMaximumClique.h"
#include"PCMBuilder.h"
#include"GlobalConstant.h"


// The algorithm for finding all the maximal cliques
// This algorithm is writtedn using C for efficience considerations
// An implementation of Bron-Kerbosch algorithm
// From Algorithm 457 of the Collected Algorithms from CACM
// http://www.netlib.org/tomspdf/457.pdf
// According to the paper: "Finding All Cliques of an Undirected Graph"

/*
 * Important note:  we treat the first found clique as the maximum clique as this algorithm will find the maximum clique as early as possible
 */

FindMaximumClique::FindMaximumClique() {}

FindMaximumClique::FindMaximumClique(int ** pGraph, int graphSize) {
	/*
	 * Make sure the Diagonal are all "1", otherwise you will encounter "s" has been used without initized
	 */
	graph = pGraph;
	N = graphSize;

	compsub.nodeInit(N);
}

FindMaximumClique::~FindMaximumClique() {}



void FindMaximumClique::findMaxclique(int * pMaximumClique, int * pMaximumCliqueSize) {
	maximumClique = pMaximumClique;
	maximumCliqueSize = pMaximumCliqueSize;
	cliqueFounded = false;
	/*
	* Start finding maximal cliques within give graph
	* will use a heuristic method to find the maximum maximal clique first 
	*/
	int i;
	int *all = (int *)malloc(N*sizeof(int));
	for (i = 0; i<N; i++)
		all[i] = i;
	bkv2(all, 0, N);

	free(all);
}

// recursive function version 2 of Bron-Kerbosch algorithm
void FindMaximumClique::bkv2(int* pOld, int ne, int ce) {
	if (cliqueFounded) {
		return;
	}

	int *newSet = (int *)malloc(ce*sizeof(int));
	int nod, fixp;
	int newne, newce, i, j, count, pos, p, s, sel, minnod;

	minnod = ce;
	nod = 0;
	// Determine each counter value and look for minimum
	for (i = 0; i <ce && minnod != 0; i++) {
		p = pOld[i];
		count = 0;
		pos = 0;
		// Count disconnections
		for (j = ne; j < ce && count < minnod; j++)
			if (!graph[p][pOld[j]]) {
				count++;
				// Save position of potential candidate
				pos = j;
			}
		// Test new minimum
		if (count < minnod) {
			fixp = p;
			minnod = count;
			if (i<ne)
				s = pos;
			else {
				s = i;
				// pre-increment
				nod = 1;
			}
		}
	}
	// If fixed point initially chosen from candidates then
	// number of diconnections will be preincreased by one
	// Backtrackcycle
	for (nod = minnod + nod; nod >= 1; nod--) {
		// Interchange
		p = pOld[s];
		pOld[s] = pOld[ne];
		sel = pOld[ne] = p;
		// Fill new set "not"
		newne = 0;
		for (i = 0; i < ne; i++)
			if (graph[sel][pOld[i]])
				newSet[newne++] = pOld[i];

		// Fill new set "cand"
		newce = newne;
		for (i = ne + 1; i<ce; i++)
			if (graph[sel][pOld[i]])
				newSet[newce++] = pOld[i];

		// Add to compsub
		compsub.add(sel);
		if (newce == 0) {
			/*****************************************
			* found a clique: RECORD OR PRINT * We report the first found clique, as it is a heuristic algorithm
			*****************************************/
			setMaximumClique();
			return;
			//compsub.print();
		}
		else if (newne < newce)
			bkv2(newSet, newne, newce);

		// Remove from compsub
		compsub.remove();

		// Add to "not"
		ne++;
		if (nod > 1)
			// Select a candidate disconnected to the fixed point
			for (s = ne; graph[fixp][pOld[s]]; s++)
				;
		// nothing but finding s

	} /* Backtrackcycle */
	free(newSet);
}


void FindMaximumClique::setMaximumClique() {
	cliqueFounded = true;
	*maximumCliqueSize = compsub.size;
	for (int i = 0; i < compsub.size; i++) {
		maximumClique[i] = compsub.node[i];
	}
}



FindMaximumClique::nodeSet::nodeSet(void) {
	N = 0;
	size = 0;
}

FindMaximumClique::nodeSet::~nodeSet(void) {
	if (N != 0)
		free(node);
}

void FindMaximumClique::nodeSet::nodeInit(int NInput) {
	N = NInput;
	node = (int *)malloc(N*sizeof(int));
	size = 0;
}

void FindMaximumClique::nodeSet::add(int nodeA) {
	node[size++] = nodeA;
}

void FindMaximumClique::nodeSet::remove(void) {
	size--;
}

void FindMaximumClique::nodeSet::print(void) {
	int i;
	printf(" the size of clique : %d - [", size);
	for (i = 0; i < size; i++) {
		if (i < size - 1)
			printf("%d,", node[i]);
		else printf("%d ]\n", node[i]);
	}
}



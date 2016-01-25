#ifndef SQO_DEF_H
#define SQO_DEF_H


#include"AdjacenceListsGraph.h"
#include"AdjacenceListsGRAPH_IO.h"
#include"TLSGraph.h"
#include"MQO.h"
#include"PCMBuilder.h"
#include"Ullmann.h"
#include"VF2.h"
#include"TurboIso.h"
#include"TurboIsoBoosted.h"
#include"TurboIso_MQO.h"
#include"TurboIsoBoosted_MQO.h"
#include<vector>
#include<iostream>
#include<list>
#include<array>

class SQP{

private:

	/*
	 * subgraph isomorphism algorithms
	 */

	Ullmann  * ullmann;

	VF2 * vf2;

	TurboIsoBoosted * turboIsoBoosted;

	TurboIso * turboIso;

private:

	/*
	 * datasets
	 */
	std::vector<AdjacenceListsGRAPH> * queryGraphVector;
	
	AdjacenceListsGRAPH * dataGraph;

	AdjacenceListsGRAPH * queryGraph;

	std::ofstream * resultFile;

	/*
	 * Used to save the number of founded embeddings
	 */
	std::vector<int> numberOfEmbeddings;

	/*
	 * store the partial embedding generated from the cached parents
	 * we use an array to improve the efficiency
	 */
	int * partialEmbedding;

	double timeClapse;


public:
	SQP();
	SQP(std::vector<AdjacenceListsGRAPH> * dataGraphVector, std::vector<AdjacenceListsGRAPH> * pQueryGraphVector, std::ofstream * pResultFile);
	~SQP();

	void queryProcessing();


};

#endif 
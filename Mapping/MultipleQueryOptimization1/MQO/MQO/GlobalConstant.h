//
// This file encapsutaltes all the global parameters
//
#pragma once
#ifndef GLOBAL_CONSTANT_H
#define GLOBAL_CONSTANT_H


class GlobalConstant {

public:
	/*
	* The cliques whose size less then the G_GOURP_QUERY_CLIQUE_MINI_SIZE won't be reported
	* The number of vertices in the clique is the number of query graphs in this group
	*/
	static int G_GOURP_QUERY_CLIQUE_MINI_SIZE;

	/*
	* The ratio is the threshold when to compute MCS between to graphs
	*/
	static double G_COMMON_TLS_RATIO;

	/*
	* Given a query group, we compute the MCS within this group. However, if the MCS is less than G_MINIMUM_NUMBER_MCS_VERTEX_RATIO * biggerGraph, this mcs will be ignored
	* TODO: we can consider to use a dynamic value for this parameter
	*/
	static double G_MINIMUM_NUMBER_MCS_VERTEX_RATIO;

	/*
	* Testing options:
	* 0: Only test the scalability of the PCM builder
	* 1: Comparasion test based on Ullmann
	* 2. Comparasion test based on VF2
	* 3. Comparasion test based on TurboIso
	* 4. Comparasion test based on TurboIsoBoosted
	*/
	enum RUN_OPTION {
		RUN_OP_PCM,
		RUN_OP_ULLMANN,
		RUN_OP_VF2,
		RUN_OP_TURBOISO,
		RUN_OP_TURBO_ISO_BOOSTED
	};

	static RUN_OPTION G_RUNNING_OPTION_INDEX;


	static int G_SUFFICIENT_NUMBER_EMBEDDINGS;


	static int G_MAXIMUM_WIDTH_COMBO;

	static int G_GLOBAL_DEBUG_COUNTER;

};




#endif
#include"SQP.h"
#include"TimeUtility.h"
#include"GlobalConstant.h"
#include<queue>


SQP::SQP(){}

SQP::SQP(std::vector<AdjacenceListsGRAPH>* dataGraphVector, std::vector<AdjacenceListsGRAPH>* pQueryGraphVector, std::ofstream * pResultFile)
{

	if (dataGraphVector->size() > 0) {
		dataGraph = &(*dataGraphVector)[0];
	}


	queryGraphVector = pQueryGraphVector;

	resultFile = pResultFile;

	/*
	 * Subgraph Isomorphism algorithms
	 */
	ullmann = new Ullmann(dataGraph, &numberOfEmbeddings, NULL);

	vf2 = new VF2(dataGraph, &numberOfEmbeddings, NULL);

	turboIsoBoosted = new TurboIsoBoosted(dataGraph, &numberOfEmbeddings, NULL);

	turboIso = new TurboIso(dataGraph, &numberOfEmbeddings, NULL);
}

SQP::~SQP()
{
	numberOfEmbeddings.swap(std::vector<int>());

	delete ullmann;

	delete vf2;

	delete turboIsoBoosted;

	delete turboIso;
}


void SQP::queryProcessing(){
	/*
	 * initialize cache embedding for all the query graphs
	 */
	for(unsigned int i=0; i<queryGraphVector->size(); i++){
		numberOfEmbeddings.push_back(0);
	}

	for(unsigned int i=0; i<queryGraphVector->size(); i++){
		//TimeUtility ttt;
		//ttt.StartCounterMill();

		queryGraph =  &(*queryGraphVector)[i];
		
		partialEmbedding = new int[queryGraph->getNumberOfVertexes()];
		for(int j=0; j < queryGraph->getNumberOfVertexes(); j++){
			partialEmbedding[j] = -1;
		}

		switch (GlobalConstant::G_RUNNING_OPTION_INDEX) {
		case GlobalConstant::RUN_OP_ULLMANN:
			ullmann->setParameters(queryGraph, partialEmbedding, false);
			ullmann->execute();
			break;
		case GlobalConstant::RUN_OP_VF2:
			vf2->setParameters(queryGraph, partialEmbedding, false);
			vf2->execute();
			break;
		case GlobalConstant::RUN_OP_TURBOISO:
			turboIso->setParameters(queryGraph, partialEmbedding, false);
			turboIso->execute();
			break;
		case GlobalConstant::RUN_OP_TURBO_ISO_BOOSTED:
			turboIsoBoosted->setParameters(queryGraph, partialEmbedding, false);
			turboIsoBoosted->execute();
			break;
		}

		//cout << " Time Cost " << ttt.GetCounterMill() << endl;
		cout << "Finish One Query: " << i << " Embedding founded " << numberOfEmbeddings[i] << endl;
		//(*resultFile)<<"Finish One Query: "<<i<<" Embedding founded "<< numberOfEmbeddings[i]<<endl;

	}
}
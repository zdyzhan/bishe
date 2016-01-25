#include"TLSGraph.h"
#include"AdjacenceListsGRAPH_IO.h"
#include"GlobalConstant.h"
#include<vector>

using namespace std;




TLSGraph::TLSGraph(int ** pTlsGraphMatrix, std::vector<AdjacenceListsGRAPH> * pQueryGraphVector, std::ofstream * pResultFile) {
	tlsGraphMatrix = pTlsGraphMatrix;
	queryGraphVector = pQueryGraphVector;
	resultFile = pResultFile;
}

TLSGraph::~TLSGraph() {
	TLSInvertedIndex.swap(std::map<TLSequence, std::map<int, int>>());
}


void TLSGraph::buildTLSGraph() {

	for (vector<AdjacenceListsGRAPH>::iterator queryGraphIterator = (*queryGraphVector).begin(); queryGraphIterator != (*queryGraphVector).end(); queryGraphIterator++) {
		for (std::map<TLSequence, std::vector<std::vector<int>>>::iterator tlsIterator = queryGraphIterator->getTLSequenceMap()->begin(); tlsIterator != queryGraphIterator->getTLSequenceMap()->end(); tlsIterator++) {

			std::map<TLSequence, std::map<int, int>>::iterator tcInvertedIndexIterator = TLSInvertedIndex.find(tlsIterator->first);

			if (tcInvertedIndexIterator == TLSInvertedIndex.end()) {
				// no this tc connection
				map<int, int> tcIndexContent;
				tcIndexContent.insert(std::pair<int, int>(queryGraphIterator->graphId, tlsIterator->second.size()));

				TLSInvertedIndex.insert(std::pair<TLSequence, map<int, int>>(tlsIterator->first, tcIndexContent));
			}
			else {
				/*
				* All the query graphs in this triple connection index has a common with this query graph
				*/
				for (map<int, int>::iterator sharedTLSQueryIterator = tcInvertedIndexIterator->second.begin(); sharedTLSQueryIterator != tcInvertedIndexIterator->second.end(); sharedTLSQueryIterator++) {
					if (sharedTLSQueryIterator->second >= tlsIterator->second.size()) {
						tlsGraphMatrix[sharedTLSQueryIterator->first][queryGraphIterator->graphId] += tlsIterator->second.size();
						tlsGraphMatrix[queryGraphIterator->graphId][sharedTLSQueryIterator->first] += tlsIterator->second.size();
					}
					else {
						tlsGraphMatrix[sharedTLSQueryIterator->first][queryGraphIterator->graphId] += sharedTLSQueryIterator->second;
						tlsGraphMatrix[queryGraphIterator->graphId][sharedTLSQueryIterator->first] += sharedTLSQueryIterator->second;
					}
				}
				/*
				* add this graph id to the end of this index entry
				*/
				tcInvertedIndexIterator->second.insert(std::pair<int, int>(queryGraphIterator->graphId, tlsIterator->second.size()));
			}
		}
	}
}

void TLSGraph::applyRelativeSimilarityRatio(std::vector<PCM_Node> * patternContainmentMap)
{
	/*
	* Apply the relative similar ratio here
	*/
	for (int i = 0; i < (*queryGraphVector).size(); i++) {
		if ((*patternContainmentMap)[i].isHidden) {
			continue;
		}

		int iTLSSize = (*queryGraphVector)[i].getTotalTLSNumber();

		for (int j = 0; j < i; j++) {
			if ((*patternContainmentMap)[j].isHidden) {
				continue;
			}

			int jTLSSize = (*queryGraphVector)[j].getTotalTLSNumber();

			if (iTLSSize > jTLSSize) {
				tlsGraphMatrix[i][j] = ((tlsGraphMatrix[i][j] / (double)iTLSSize) > GlobalConstant::G_COMMON_TLS_RATIO) ? 1 : 0;
			}
			else {
				tlsGraphMatrix[i][j] = ((tlsGraphMatrix[i][j] / (double)jTLSSize) > GlobalConstant::G_COMMON_TLS_RATIO) ? 1 : 0;
			}

			tlsGraphMatrix[j][i] = tlsGraphMatrix[i][j];
		}
	}
}



void TLSGraph::showTLSGraphMatrix() {
	
	int numberOfEdges = 0;

	(*resultFile) << "The TLS graph matrix is: " << endl;

	(*resultFile) << "    ";
	for (int i = 0; i < queryGraphVector->size(); i++) {
		(*resultFile) << i << " ";
	}
	(*resultFile) << endl;
	for (int i = 0; i < queryGraphVector->size(); i++) {
		(*resultFile) << "(" << i << ") ";
		for (int j = 0; j < queryGraphVector->size(); j++) {
			(*resultFile) << tlsGraphMatrix[i][j] << " ";
		}
		(*resultFile) << endl;
	}


	for (int i = 0; i < queryGraphVector->size(); i++) {
		for (int j = 0; j < queryGraphVector->size(); j++) {
			if (tlsGraphMatrix[i][j] >= 1) {
				numberOfEdges++;
			}
		}
	}
	(*resultFile) << "TLS: total number of TLS common-edges : "<< (numberOfEdges - queryGraphVector->size()) / 2 << endl;
}

void TLSGraph::showTLSInvertedIndex()
{
	cout << "TLSInvertedIndex: " << endl;
	for (std::map<TLSequence, std::map<int, int>>::iterator tlsMapIterator = TLSInvertedIndex.begin(); tlsMapIterator != TLSInvertedIndex.end(); tlsMapIterator++) {
		cout << tlsMapIterator->first.start << " " << tlsMapIterator->first.pivot << " " << tlsMapIterator->first.end<<":";
		for (std::map<int, int>::iterator tlsQueryIterator = tlsMapIterator->second.begin(); tlsQueryIterator != tlsMapIterator->second.end(); tlsQueryIterator++) {
			cout << tlsQueryIterator->first << "(" << tlsQueryIterator->second << ")" << " " ;
		}
		cout << endl;
	}

}


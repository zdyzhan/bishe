//#include<iostream>
//#include<fstream>
//#include"PCMBuilder.h"
//#include"TimeUtility.h"
//#include"MQO.h"
//#include"SQP.h"
//#include"InputCommandLineParser.h"
//#include"AdjacenceListsGRAPH_BOOST.h"
//#include"GlobalConstant.h"
//#include"CoreSubgraphSpanOverlapQG.h"
//#include"CoreVertexSpanOverlapQG.h"
//#include"OverlapQueryGenerator.h"
//#include"RandomQG.h"
//#include"CacheResultsTest.h"
//
//
//using namespace std;
//
//
//string dataGraphFileName, queryGraphFileName, hyperGraphFileName, containmentGraphFileName, resultFilename;
//std::vector<AdjacenceListsGRAPH> dataGraphVector, queryGraphVector;
//std::ofstream resultFile;
//
//
//int queryType = 0;
//int numberOfCores = 0;
//int groupSize = 0;
//
//void loadDataGraphs();
//void loadQueryGraphs(int argc, char* argv[]);
//void loadHyperGraphs();
//void loadResultFile(int argc, char* argv[]);
//void setTheConfigurationSettings(int argc, char* argv[]);
//void help();
//
//
//bool scalabilityTest;
//
//int maindd(int argc, char* argv[]) {
//	
//	loadResultFile(argc, argv);
//	setTheConfigurationSettings(argc, argv);
//
//	/*******************************************
//	 * Query Graph Generator
//	 *******************************************/
//	if (InputCommandLineParser::cmdOptionExists(argc, argv, "-gen")) {
//		dataGraphFileName = InputCommandLineParser::getCmdOption(argc, argv, "-dg");
//		queryType = std::stoi(InputCommandLineParser::getCmdOption(argc, argv, "-t"));
//		numberOfCores = std::stoi(InputCommandLineParser::getCmdOption(argc, argv, "-nc"));
//		groupSize = std::stoi(InputCommandLineParser::getCmdOption(argc, argv, "-gs"));
//		
//
//		std::ifstream datagraphfile = std::ifstream(dataGraphFileName);
//		if (!datagraphfile) {
//			resultFile << "the data file '" << dataGraphFileName << "' doesn't exist." << endl;
//		}
//		AdjacenceListsGRAPH_IO::loadGraphFromFile(datagraphfile, dataGraphVector);
//
//		//OverlapQueryGenerator overlapQueryGenerator(&dataGraphVector[0], queryType, &resultFile, numberOfCores);
//		//overlapQueryGenerator.generateQueries();
//
//		//CoreSubgraphSpanOverlapQG coreSubgraphSpanOverlapQG(&dataGraphVector[0], queryType, &resultFile, numberOfCores, groupSize);
//		//coreSubgraphSpanOverlapQG.generateQueries();
//
//		//RandomQG randomQG(&dataGraphVector[0], queryType, &resultFile, numberOfCores);
//		//randomQG.generateQueries();
//
//		CoreVertexSpanOverlapQG coreVertexSpanOverlapQG(&dataGraphVector[0], numberOfCores * groupSize, 1 / (double)groupSize, 0, &resultFile);
//		coreVertexSpanOverlapQG.generateQueries();
//
//		resultFile.close();
//		
//		resultFile << "Queries generated. |\-/|" << endl;
//		//system("pause"); //@debug
//		return 0;
//	}
//
//
//
//	if (InputCommandLineParser::cmdOptionExists(argc, argv, "-cacheTest")) {
//		
//		queryGraphFileName = InputCommandLineParser::getCmdOption(argc, argv, "-qg");
//		loadQueryGraphs(argc, argv);
//
//		dataGraphFileName = InputCommandLineParser::getCmdOption(argc, argv, "-dg");
//		loadDataGraphs();
//
//		CacheResultsTest cacheResultsTest(&dataGraphVector[0], &queryGraphVector, &resultFile);
//		cacheResultsTest.execute();
//
//		return 0;
//	}
//
//
//
//	/*******************************************
//	*   MQO Testers *
//	*******************************************/
//	loadQueryGraphs(argc, argv);
//	if (InputCommandLineParser::cmdOptionExists(argc, argv, "-s")) {
//		scalabilityTest = true;
//	}
//	else {
//		if (InputCommandLineParser::cmdOptionExists(argc, argv, "-dg") && InputCommandLineParser::cmdOptionExists(argc, argv, "-hg")) {
//			resultFile << "Wrong Parameters! You cannot use -dg and -hg at the same time. hg is for BoostedIso " << endl;
//			help();
//			exit(1);
//		}
//		else if (InputCommandLineParser::cmdOptionExists(argc, argv, "-dg")) {
//			dataGraphFileName = InputCommandLineParser::getCmdOption(argc, argv, "-dg");
//		}
//		else if (InputCommandLineParser::cmdOptionExists(argc, argv, "-hg")) {
//			hyperGraphFileName = InputCommandLineParser::getCmdOption(argc, argv, "-hg");
//			if (InputCommandLineParser::cmdOptionExists(argc, argv, "-cg")) {
//				containmentGraphFileName = InputCommandLineParser::getCmdOption(argc, argv, "-cg");
//			}
//			else {
//				resultFile << "Wrong Parameters! No containment file specified. " << endl;
//				help();
//				exit(1);
//			}
//		}
//		else {
//			resultFile << "Wrong Parameters! No data file specified. Use -dg for normal subgraph isomorphism and -hg for BoostIso subgraph isomorphism " << endl;
//			help();
//			exit(1);
//		}
//	}
//
//
//
//	if (InputCommandLineParser::cmdOptionExists(argc, argv, "-dg")) {
//		loadDataGraphs();
//	}
//	else if(InputCommandLineParser::cmdOptionExists(argc, argv, "-hg")){
//		loadHyperGraphs();
//	}
//
//
//	//@debug
//	/*dataGraphFileName = "C:/Users/s2813995/Dropbox/HectorResearch/PaperProjects/MQP_Center/MQOTestData/yeast.graph";
//	queryGraphFileName = "C:/Users/s2813995/Dropbox/HectorResearch/PaperProjects/MQP_Center/MQOExperiment/yeast_q9.graph";
//	hyperGraphFileName = "C:/Users/s2813995/Dropbox/HectorResearch/PaperProjects/MQP_Center/MQOTestData/yeast.hgraph";
//	containmentGraphFileName = "C:/Users/s2813995/Dropbox/HectorResearch/PaperProjects/MQP_Center/MQOTestData/yeast.cgraph";
//	resultFilename = "C:/Users/s2813995/Dropbox/HectorResearch/PaperProjects/MQP_Center/MQOTestData/scalability.result";
//	scalabilityTest = true; 
//
//	//loadHyperGraphs();
//	//loadDataGraphs();
//	loadQueryGraphs();
//	setTheConfigurationSettings(argc, argv);*/
//
//	
//	if(!scalabilityTest){
//		SQP sqo(&dataGraphVector, &queryGraphVector, &resultFile);
//
//		TimeUtility tSQP;
//		tSQP.StartCounterMill();
//		sqo.queryProcessing();
//		resultFile << endl << "1. SQO Average Time: " << tSQP.GetCounterMill() << "(milliseconds)" << endl;
//	}
//
//
//
//	MQO mqo(&dataGraphVector, &queryGraphVector, &resultFile);
//	
//	TimeUtility tMQO;
//	tMQO.StartCounterMill();
//	mqo.buildPCM();
//	if(scalabilityTest){
//		resultFile <<"PCM Build Time: "<< tMQO.GetCounterMill()<<"(milliseconds)"<<endl;
//	}else{
//		mqo.orederedQueryProcessing();
//		resultFile << "2. MQO Time: " << tMQO.GetCounterMill() << "(milliseconds)" << endl;
//	}
//	resultFile.close();
//
//	//system("pause"); //@debug
//
//	return 0;
//}
//
//
//void loadDataGraphs() {
//
//	std::ifstream dataGraphFile = std::ifstream(dataGraphFileName);
//	if (!dataGraphFile) {
//		resultFile << "The data file '" << dataGraphFileName << "' doesn't exist." << endl;
//	}
//	AdjacenceListsGRAPH_IO::loadGraphFromFile(dataGraphFile, dataGraphVector);
//
//	for (vector<AdjacenceListsGRAPH>::iterator dataGraphIterator = dataGraphVector.begin(); dataGraphIterator != dataGraphVector.end(); dataGraphIterator++) {
//		dataGraphIterator->buildLabelVertexList();
//		dataGraphIterator->buildVertexLabelVertexList();
//	}
//
//	resultFile << "2. Load data graphs and build data graph index done. " << dataGraphFileName << endl;
//}
//
//void loadQueryGraphs(int argc, char* argv[]) {
//
//	if (InputCommandLineParser::cmdOptionExists(argc, argv, "-qg")) {
//		queryGraphFileName = InputCommandLineParser::getCmdOption(argc, argv, "-qg");
//	}
//	else {
//		resultFile << "Wrong Parameters! Missing query graph [-qg]. " << endl;
//		help();
//		exit(1);
//	}
//
//
//	std::ifstream queryGraphFile = std::ifstream(queryGraphFileName);
//	if (!queryGraphFile) {
//		resultFile << "The query file '" << queryGraphFileName << "' doesn't exist." << endl;
//	}
//	AdjacenceListsGRAPH_IO::loadGraphFromFile(queryGraphFile, queryGraphVector);
//
//
//	int queryGraphIndex = 0;
//	for (vector<AdjacenceListsGRAPH>::iterator queryGraphIterator = queryGraphVector.begin(); queryGraphIterator != queryGraphVector.end(); queryGraphIterator++) {
//
//		queryGraphIterator->graphId = queryGraphIndex++;
//		/*
//		* build inverted indexes for each query graph
//		*/
//		queryGraphIterator->buildDFSTraversalOrder();
//		queryGraphIterator->buildVertexLabelVertexList();
//		queryGraphIterator->buildLabelVertexList();
//
//		queryGraphIterator->buildVertexLabelEdgeList(); // used by mcs computation
//		queryGraphIterator->buildTLSequence(); // used by pcm computation
//
//		queryGraphIterator->buildComboGraph(GlobalConstant::G_MAXIMUM_WIDTH_COMBO);
//	}
//
//	resultFile << "1. Load query graphs and build query graph index done. " << queryGraphFileName << endl;
//}
//
//void loadHyperGraphs() {
//	
//	std::ifstream hyperGraphFile = std::ifstream(hyperGraphFileName);
//	std::ifstream containmentGraphFile = std::ifstream(containmentGraphFileName);
//
//	AdjacenceListsGRAPH_IO::loadGraphFromFile(hyperGraphFile, dataGraphVector);
//	AdjacenceListsGRAPH_BOOST::loadContainmentGraph(dataGraphVector, containmentGraphFile);
//
//	/*
//	* Build indexes for the hypergraph
//	* 1. For each label within label set of the graph, we attach all data vertices with this lable to it
//	* 2. For each label within the label set of each vertex's neighbours, we attach all its neighbours with this label to it
//	* 3. For each label within label set of the graph, we attach the s-roots vertices with this lable to it. s-roots have no s-contained parents
//	*/
//	for (int hyperGraphIndex = 0; hyperGraphIndex < dataGraphVector.size(); hyperGraphIndex++) {
//		dataGraphVector[hyperGraphIndex].buildLabelVertexList();
//		dataGraphVector[hyperGraphIndex].buildVertexLabelVertexList();
//		AdjacenceListsGRAPH_BOOST::buildLabelRootMap(dataGraphVector[hyperGraphIndex]);
//	}
//
//	resultFile << "2. Load hyper graphs and build hyper graph index done. " << hyperGraphFileName << endl;
//}
//
//void loadResultFile(int argc, char* argv[]) {
//	if (InputCommandLineParser::cmdOptionExists(argc, argv, "-out")) {
//		resultFilename = InputCommandLineParser::getCmdOption(argc, argv, "-out");
//		resultFile = std::ofstream(resultFilename, std::ios_base::app);
//	}
//	else {
//		resultFile << "Wrong Parameters! No output file specified. Use [-out] " << endl;
//		help();
//		exit(1);
//	}
//}
//
//void setTheConfigurationSettings(int argc, char* argv[]) {
//	
//	GlobalConstant::G_GOURP_QUERY_CLIQUE_MINI_SIZE = 2;
//
//	if (InputCommandLineParser::cmdOptionExists(argc, argv, "-tlsRatio")) {
//		GlobalConstant::G_COMMON_TLS_RATIO = std::stod(InputCommandLineParser::getCmdOption(argc, argv, "-tlsRatio"));
//	}
//	else {
//		GlobalConstant::G_COMMON_TLS_RATIO = 0.5;
//	}
//	if (InputCommandLineParser::cmdOptionExists(argc, argv, "-minMCSRatio")) {
//		GlobalConstant::G_MINIMUM_NUMBER_MCS_VERTEX_RATIO = std::stod(InputCommandLineParser::getCmdOption(argc, argv, "-minMCSRatio"));
//	}
//	else {
//		GlobalConstant::G_MINIMUM_NUMBER_MCS_VERTEX_RATIO = 0.5;
//	}
//	if (InputCommandLineParser::cmdOptionExists(argc, argv, "-subIso")) {
//		
//		char * algorithmIndex = InputCommandLineParser::getCmdOption(argc, argv, "-subIso");
//
//		if (strcmp(algorithmIndex, "ULLMANN") == 0) {
//			GlobalConstant::G_RUNNING_OPTION_INDEX = GlobalConstant::RUN_OP_ULLMANN;
//		}
//		else if (strcmp(algorithmIndex, "VF2") == 0) {
//			GlobalConstant::G_RUNNING_OPTION_INDEX = GlobalConstant::RUN_OP_VF2;
//		}
//		else if (strcmp(algorithmIndex, "TURBOISO") == 0) {
//			GlobalConstant::G_RUNNING_OPTION_INDEX = GlobalConstant::RUN_OP_TURBOISO;
//		}
//		else if (strcmp(algorithmIndex, "TURBO_ISO_BOOSTED") == 0) {
//			GlobalConstant::G_RUNNING_OPTION_INDEX = GlobalConstant::RUN_OP_TURBO_ISO_BOOSTED;
//		}
//
//	}
//	else {
//		GlobalConstant::G_RUNNING_OPTION_INDEX = GlobalConstant::RUN_OP_TURBOISO;
//	}
//
//	if (InputCommandLineParser::cmdOptionExists(argc, argv, "-suffiEmb")) {
//		GlobalConstant::G_SUFFICIENT_NUMBER_EMBEDDINGS = std::stoi(InputCommandLineParser::getCmdOption(argc, argv, "-suffiEmb"));
//	}
//	else {
//		GlobalConstant::G_SUFFICIENT_NUMBER_EMBEDDINGS = 1000;
//	}
//
//	GlobalConstant::G_MAXIMUM_WIDTH_COMBO = 3;
//}
//
//
//void help() {
//
//	resultFile << "|\-/| Testing Options: " << endl;
//	resultFile << "-s  Only test the scalability of PCM building. If set, only the query graph file will be optional" << endl;
//	resultFile << "-dg  [datagraphFilename] The datagraph file" << endl;
//	resultFile << "-qg  [querygraphFilename] The querygraph file" << endl;
//	resultFile << "-hg  [hypergraphFilename] The hyerpergraph file" << endl;
//	resultFile << "-cg  [containmentgraphFilename] The boosted containment file" << endl;
//	resultFile << "-out  [outputResultFile] Append the result to the this file" << endl;
//
//	resultFile << "|\-/| Optional Testing Configurations: " << endl;
//	resultFile << "-subIso  [ULLMANN | VF2 | TURBOISO | TURBO_ISO_BOOSTED] the subgraph isomorphism algorithm choosed" << endl;
//	resultFile << "-tlsRatio  The tls ratio for comparasing TLS between two graphs" << endl;
//	resultFile << "-minMCSRatio  The min MCS ratio when build the MCS within a query clique " << endl;
//	resultFile << "-suffiEmb  sufficient number of embedddings for each query" << endl;
//
//	//system("pause"); //@debug
//}

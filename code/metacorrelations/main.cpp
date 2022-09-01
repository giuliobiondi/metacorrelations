#include "headers.h"
#include "cxxopts.hpp"
#include "graphtools.h"
#include "measures.h"
#include "Classifier.h"
#include "de.h"
#include "Timer.h"


std::vector<std::vector<double>> parse2DCsvFile(std::string inputFileName) {
    std::cout<<"Loading file: "<<inputFileName<<std::endl;
    std::vector<std::vector<double>> data;
    std::ifstream inputFile(inputFileName);
    int l = 0;
    if (inputFile){
        while (inputFile) {
            l++;
            std::string s;
            if (!getline(inputFile, s)) break;
            if (s[0] != '#') {
                std::istringstream ss(s);
                std::vector<double> record;

                while (ss) {
                    std::string line;
                    if (!getline(ss, line, ','))
                        break;
                    try {
                        record.push_back(stod(line));
                    }
                    catch (const std::invalid_argument e) {
                        std::cout << "NaN found in file " << inputFileName << " line " << l
                            << std::endl;
                        e.what();
                    }
                }

                data.push_back(record);
            }
        }
        const size_t D = data[0].size();
        const size_t NP = data.size();
        double** matrix = new double*[NP];            
        for (size_t i =0; i < NP; i++){
            matrix[i] = new double[D];
            for (size_t j = 0; j<D; j++){
                matrix[i][j] = data[i][j];
            }
        }

        if (!inputFile.eof()) {
            std::cerr << "Could not read file " << inputFileName << "\n";
            throw std::invalid_argument("File not found.");
        }

        return data;
    }
    else{
        std::cout<<"Error: correlations file not found"<<std::endl;
        exit(1);
    }
}

void printConfigInfo(Config a){
    std::cout<<"Data set name: "<<a.datasetName<<std::endl;
    std::cout<<"Data set file name: "<<a.datasetFileName<<std::endl;
    std::cout<<"Data set path: "<<a.pathToFile<<std::endl;
    std::cout<<"Measures: ";
    for (auto i: a.measures){
        cout<<i<<" ";
    }
    std::cout<<std::endl<<"Number of folds: "<<a.nfold<<std::endl;
    std::cout<<"K to use in top-k evaluation: "<<a.k<<std::endl;
}
Config readConfig(std::string configPath){
    PathManager &pathmng = PathManager::getInstance();
    std::ifstream t(configPath);
    std::stringstream buffer;
    buffer << t.rdbuf();
    json j3 = json::parse(buffer.str());
	std::cout << j3.dump(4) << std::endl;
    Config newConfig;
    newConfig.datasetName = j3["name"].get<std::string>();
    newConfig.datasetFileName = j3["filename"].get<std::string>();
    newConfig.pathToFile = j3["path"].get<std::string>();
    newConfig.nodesPath = j3["nodesPath"].get<std::string>();
    newConfig.trainingSetPath = j3["trainingPath"].get<std::string>();
    newConfig.validationSetPath = j3["validationPath"].get<std::string>();
    newConfig.testSetPath = j3["testPath"].get<std::string>();
    newConfig.nfold = j3["nfold"];
    newConfig.k = j3["k"];
	newConfig.aucSamples = j3["aucSamples"];
    newConfig.dump = j3["dump"];
    std::vector<std::string> metrics = j3["metrics"];
    newConfig.metrics = metrics;
    pathmng.setNodesPath(j3["nodesPath"].get<std::string>());
    pathmng.setTrainingSetPath(j3["trainingPath"].get<std::string>());
    pathmng.setValidationSetPath(j3["validationPath"].get<std::string>());
    if(!(j3["trainingValidationPath"].is_null())){
        pathmng.setTrainingValidationSetPath(j3["trainingValidationPath"].get<std::string>());
    }
    pathmng.setTestSetPath(j3["testPath"].get<std::string>());
    if(j3["measures"].is_null()){
        std::cout<<"Error! No measures indicated in config file"<<std::endl;
        std::cout<<"Exiting"<<std::endl;
        exit(1);
    }
    else{
        pathmng.setMeasures(j3["measures"]);
    }
     if(j3["metrics"].is_null()){
        std::cout<<"Error! No metrics indicated in config file"<<std::endl;
        std::cout<<"Exiting"<<std::endl;
        exit(1);
    }
    else{
        pathmng.setMetrics(j3["metrics"]);
    }
    return newConfig;
}

rankedGraph calculate_ranking(const PUNGraph& training_set,const PUNGraph& test_set){
    PathManager &pathmng = PathManager::getInstance();
    rankedGraph ret;
    size_t epot = (static_cast<size_t>(training_set->GetNodes())*static_cast<size_t>((training_set->GetNodes()-1))) / 2-static_cast<size_t>(training_set->GetEdges());
    Time::Timer timer(true);
    std::cout <<"capacity: "<<ret.capacity()<<std::endl;
    std::cout <<": "<<ret.size()<<std::endl;
    ret.reserve(epot);
    for(TUNGraph::TNodeI i = training_set->BegNI(); i < training_set->EndNI();i++){
        int node1 = i.GetId();
        TUNGraph::TNodeI j(i);
        for(j++; j < training_set->EndNI(); j++){
            int node2 = j.GetId();
            if(!training_set->IsEdge(node1,node2)){
                if(test_set->IsEdge(node1, node2)){
                    ret.push_back(std::make_tuple(node1, node2, 0.0, 1));
                }
                else{
                    ret.push_back(std::make_tuple(node1, node2, 0.0, 0));
                }
            }
        }
    }
    std::cout << "Ranking creation time(ms): "<< timer.stop() <<std::endl;
    std::cout << "Elements in ranking: "<<ret.size()<<std::endl;

    timer.start();
    std::string currentMeasure = pathmng.getMeasureName();
    std::cout<<currentMeasure<<std::endl;
    measuresFunction chosenMeasure = std::get<0>(Measures::measuresFunctions[currentMeasure]);
    bool isProximity = std::get<1>(Measures::measuresFunctions[currentMeasure]);
#pragma omp parallel for
    for(size_t i = 0; i < ret.size(); i++){
        TIntV nbrs;
        auto& node1 = std::get<0>(ret[i]);
        auto& node2 = std::get<1>(ret[i]);
        int cmn_nbrs = TSnap::GetCmnNbrs(training_set, node1 ,node2, nbrs);
        std::get<2>(ret[i]) = chosenMeasure(training_set, cmn_nbrs, node1, node2, nbrs);
    }
    std::cout<<"End calculation"<<std::endl;
    std::cout << "Ranking calculation time (s): "<<timer.stop()<<std::endl;
    timer.start();
    if(isProximity){
        #ifdef _MSC_VER
            size_t epot = (static_cast<size_t>(training_set->GetNodes())*static_cast<size_t>((training_set->GetNodes()-1)))/2-static_cast<size_t>(training_set->GetEdges());
            std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
            std::cout <<"capacity: "<<ret.capacity()<<std::endl;
            std::cout <<": "<<ret.size()<<std::endl;
            ret.reserve(epot);
            concurrency::parallel_sort(ret.begin(), ret.end(), [](rankingEdge a, rankingEdge b) {
                    return std::get<2>(a) > std::get<2>(b); });
        #else
                __gnu_parallel::stable_sort(ret.begin(), ret.end(), [](rankingEdge a, rankingEdge b) {
                    return std::get<2>(a) > std::get<2>(b);});
        #endif
    }
    else{
        #ifdef _MSC_VER
                concurrency::parallel_sort(ret.begin(), ret.end(), [](rankingEdge a, rankingEdge b) {
                    return std::get<2>(a) < std::get<2>(b); });
        #else
                __gnu_parallel::stable_sort(ret.begin(), ret.end(), [](rankingEdge a, rankingEdge b) {
                    return std::get<2>(a) < std::get<2>(b); });
        #endif
    }
    std::cout << "Ranking sort time (s): "<<timer.stop()<<std::endl;
    return ret;
}

void test_evolution(int k=-1){
    std::vector<std::string> metrics{"precision","auc","approxAuc","aupr","avgPrec"};
    PathManager &pathmng = PathManager::getInstance();
    PUNGraph training_set = GraphLoader::loadGraphNodesEdgesList<PUNGraph>(TUNGraph::New(), pathmng.getNodesPath(), pathmng.getTrainingSetPath());
    PUNGraph test_set = GraphLoader::loadGraphNodesEdgesList<PUNGraph>(TUNGraph::New(), pathmng.getNodesPath(), pathmng.getTestSetPath());
    std::cout<<"Training Set Edges: "<< training_set->GetEdges()<<std::endl;
    std::cout<<"Test Set Edges: "<< test_set->GetEdges()<<std::endl;
    std::vector<resultsRow> validationResults;
    for(auto it : pathmng.getMeasures())
    {
        pathmng.setMeasureName(it);
        std::cout<<"Measure: "<<it<<std::endl;
        auto ranking = calculate_ranking(training_set, test_set);
        //resultsRow results = Classifier::classify(training_set, test_set, ranking, 0, -1, metrics);
        validationResults.push_back(Classifier::classify(training_set, test_set, ranking, 0, -1,metrics));
        if(pathmng.getDumpRanking()){
            GraphTools::dump_ranking_to_file(ranking, pathmng.getRankingFilePath());
        }
    }
    Random* myRand=Random::getInstance();
    std::ofstream shuffleSeedFile;
    shuffleSeedFile.open (pathmng.getShuffleSeedFilePath(), std::ofstream::app);
    auto shuffleSeeds = Random::getInstance()->getShuffleVector();
    for (unsigned int i = 0; i < shuffleSeeds.size(); i++){
        shuffleSeedFile<<shuffleSeeds[i]<<",";
    }
    shuffleSeedFile<<std::endl;
    std::ofstream aucSeedFile;
    aucSeedFile.open (pathmng.getAucSeedFilePath(), std::ofstream::app);
    auto aucSeeds = Random::getInstance()->getAucVector();
    for (unsigned int i = 0; i < aucSeeds.size(); i++){
        aucSeedFile<<aucSeeds[i]<<",";
    }
    aucSeedFile<<std::endl;
    std::ofstream randomSeedFile;
    randomSeedFile.open (pathmng.getRandomMeasureFilePath(), std::ofstream::app);
    auto randomGenSeed = Random::getInstance()->getRandomGenSeed();
    randomSeedFile<<randomGenSeed<<std::endl;
    Random::getInstance()->reSeed();
    std::cout<<"";
    std::cout<<pathmng.getResultsFilePath()<<std::endl;
    Classifier::dumpMeasuresEvaluationResults(pathmng.getMeasuresResultsFilePath(),validationResults);
}

int main(int argc, char *argv[])
{
    //Random* ran = Random::getInstance(1234,1234);
        char ch;
        cxxopts::Options options("LPDB", "LPDB");
        options.add_options()
          ("mod", "Task", cxxopts::value<std::string>()->default_value("n"))
          ("preprocess", "Perform preprocessing", cxxopts::value<std::string>()->default_value("n"))
          ("run", "Perform run", cxxopts::value<std::string>()->default_value("n"))
          ("correlations_folder","Correlations folder", cxxopts::value<std::string>()->default_value("./"))
          ("metacorr","Metacorrelation", cxxopts::value<std::string>()->default_value("1"))
          ("genmax","DE generations", cxxopts::value<std::string>()->default_value("10"))
          ("CR","CrossoverProbability", cxxopts::value<std::string>()->default_value("0.9"))
          ("FP","F parameter", cxxopts::value<std::string>()->default_value("0.5"))
          ("config","config file", cxxopts::value<std::string>()->default_value(""))
          ("basePath","basePath", cxxopts::value<std::string>()->default_value(""))
          ("dataset","datasetName", cxxopts::value<std::string>()->default_value(""))
          ("validation_fold","validation fold", cxxopts::value<std::string>()->default_value(""))
          ("test_fold","validation fold", cxxopts::value<std::string>()->default_value(""))
          ("strategy","DE strategy", cxxopts::value<std::string>()->default_value("2"))
          ("fitnessf","Evolution metric", cxxopts::value<std::string>()->default_value("precision"))
          ("population","Population", cxxopts::value<std::string>()->default_value("standard"))
          ("aucSeed","Seed for approximated AUC calculation", cxxopts::value<std::string>()->default_value("0"))
          ("DESeed","Seed for DE", cxxopts::value<std::string>()->default_value("1234"))
          ("shuffleSeed","Seed for shuffling operation", cxxopts::value<std::string>()->default_value("0"))
          ("silent","Suppress output", cxxopts::value<std::string>()->default_value("y"));
          ;
        auto commandLineOptions = options.parse(argc, argv);
        std::string mod = commandLineOptions["mod"].as<std::string>();
        PathManager &pathmng = PathManager::getInstance();
        if(mod == "evolve"){
            std::cout<<"Starting evolutionary search"<<std::endl;
            std::string correlationsFolder = commandLineOptions["correlations_folder"].as<std::string>();
            int metacorr=std::stoi(commandLineOptions["metacorr"].as<std::string>());
            int genmax=std::stoi(commandLineOptions["genmax"].as<std::string>());
            double CR=std::stod(commandLineOptions["CR"].as<std::string>());
            double F=std::stod(commandLineOptions["FP"].as<std::string>());
            int strategy = std::stoi(commandLineOptions["strategy"].as<std::string>());
            std::string population = commandLineOptions["population"].as<std::string>();
            unsigned int shuffleSeed=(unsigned int)std::stoul(commandLineOptions["shuffleSeed"].as<std::string>());
            unsigned int aucSeed=(unsigned int)std::stoul(commandLineOptions["aucSeed"].as<std::string>());
            unsigned int DESeed=(unsigned int)std::stoul(commandLineOptions["DESeed"].as<std::string>());
            std::cout<<"DE SEED "<<DESeed<<std::endl;
            std::string basePath=commandLineOptions["basePath"].as<std::string>();
            std::cout<<"BASEPATH "<<basePath<<std::endl;
            std::string dataset=commandLineOptions["dataset"].as<std::string>();            
            int validation_fold=std::stoi(commandLineOptions["validation_fold"].as<std::string>());
            int test_fold=std::stoi(commandLineOptions["test_fold"].as<std::string>());
            std::string fitnessf=commandLineOptions["fitnessf"].as<std::string>();
            Random* ran = Random::getInstance(shuffleSeed, aucSeed, DESeed, 0);
            pathmng.setConfig(basePath,correlationsFolder,dataset,validation_fold,test_fold, population, fitnessf, metacorr,strategy,genmax,F,CR);
            std::vector<std::string> metrics{"precision","auc"};
            pathmng.setMetrics(metrics);
            auto correlationsVector = parse2DCsvFile(pathmng.getCorrelationsFilePath()); 
            std::cout<<"Correlations vector size "<<correlationsVector.size()<<std::endl;
            DE DEInstance(correlationsVector, metacorr, genmax, CR, F, strategy, fitnessf);
            DEInstance.optimize();
            DEInstance.test_evolution();
            return 0;
        }
        else if (mod == "test"){
            std::string pathToConfigFile=commandLineOptions["config"].as<std::string>();
            pathmng.readConfig(pathToConfigFile);
            test_evolution();
            return 0;
        }
}

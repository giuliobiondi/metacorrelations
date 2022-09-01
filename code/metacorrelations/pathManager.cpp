#include "pathManager.h"

PathManager& PathManager::getInstance(){
    static PathManager instance;
    return instance;
}

void PathManager::setMeasureName(std::string measure){
    this->currentMeasure = measure;
}

void PathManager::setNodesPath(std::string path){
    this->nodes_path = path;
}

void PathManager::setTrainingSetPath(std::string path){
    this->training_set_path = path;
}

void PathManager::setValidationSetPath(std::string path){
    this->validation_set_path = path;
}

void PathManager::setTrainingValidationSetPath(std::string path){
    this->training_validation_set_path = path;
}

void PathManager::setTestSetPath(std::string path){
    this->test_set_path = path;
}

void PathManager::setValidationFold(int validation_fold){
    this->validationFold = validation_fold;
}

void PathManager::setTestFold(int test_fold){
    this->testFold = test_fold;
}

void PathManager::setDataset(std::string dataset){
    this->dataset = dataset;
}

std::string PathManager::getCorrelationsFilePath(){
    return this->correlations_file;
}

std::string PathManager::getMeasureName(){
    return this->currentMeasure;
}

std::string PathManager::getNodesPath(){
    return this->nodes_path;
}

std::string PathManager::getTrainingSetPath(){
    return this->training_set_path;
}

std::string PathManager::getValidationSetPath(){
    return this->validation_set_path;
}

std::string PathManager::getTrainingValidationSetPath(){
    return this->training_validation_set_path;
}

std::string PathManager::getTestSetPath(){
    return this->test_set_path;
}

std::string PathManager::getBestCorrelationFilePath(std::string file_prefix){
    return Filesystem::get_directory (this->base_path) + Filesystem::separator() + dataset + Filesystem::separator() + dataset + "_best.csv";
}

std::string PathManager::getShuffleSeedFilePath(std::string file_prefix){
    return Filesystem::get_directory (this->training_set_path) + Filesystem::separator() + file_prefix + "shuffleSeed.csv";
}

std::string PathManager::getAucSeedFilePath(std::string file_prefix){
    return Filesystem::get_directory (this->training_set_path) + Filesystem::separator() + file_prefix + "aucSeed.csv";
}

std::string PathManager::getRankingFilePath(){
    return Filesystem::get_directory (this->training_set_path) + Filesystem::separator() + this->currentMeasure + "_ranking.csv";
}

std::string PathManager::getResultsFilePath(std::string file_prefix){
    if(this->customResultsFilePath != ""){
        return this->customResultsFilePath;
    }
    std::cout<<Filesystem::get_directory (this->base_path) + dataset + Filesystem::separator() + dataset + "_results_table.csv";
    return Filesystem::get_directory (this->base_path) + dataset + Filesystem::separator() + dataset + "_results_table.csv";
}

std::string PathManager::getMeasuresResultsFilePath(std::string file_prefix){
    if(this->customResultsFilePath != ""){
        return this->customResultsFilePath;
    }
    std::cout<<Filesystem::get_directory (this->base_path) + dataset + Filesystem::separator() + dataset + "_measures_results_table.csv";
    return Filesystem::get_directory (this->base_path) + dataset + Filesystem::separator() + dataset + "_measures_results_table.csv";
}

std::string PathManager::getAucPath(){
    return Filesystem::get_directory (this->base_path) + Filesystem::separator()  + dataset + Filesystem::separator() + dataset + "_auc_points.csv";
}

std::string PathManager::getAuPRPath(){
    return Filesystem::get_directory (this->base_path) + Filesystem::separator()  + dataset + Filesystem::separator() + dataset + "_aupr_points.csv";
}

std::string PathManager::getAvgPrecPath(){
    return Filesystem::get_directory (this->base_path) + Filesystem::separator()  + dataset + Filesystem::separator() + dataset + "_avgPrec_points.csv";
}

std::string PathManager::getRandomMeasureFilePath(std::string file_prefix){
    return Filesystem::get_directory (this->training_set_path) + Filesystem::separator()  + this->currentMeasure + "_random_measure_seed.csv";
}
std::string PathManager::getDERankingFilePath(){
    return this->base_path + dataset + Filesystem::separator() + dataset + "_ranking.csv";
}
int PathManager::getValidationFold(){
    return this->validationFold;
}

int PathManager::getTestFold(){
    return this->testFold;
}
std::string PathManager::getDataset(){
    return this->dataset;
}

configRow PathManager::getConfig(){
    configRow temp;
    std::get<0>(temp) = this->dataset;
    std::get<1>(temp) = this->validationFold;
    std::get<2>(temp) = this->testFold;
    std::get<3>(temp) = this->metacorr;
    std::get<4>(temp) = this->strategy;
    std::get<5>(temp) = this->genmax;
    std::get<6>(temp) = this->F;
    std::get<7>(temp) = this->CR;
    std::get<8>(temp) = this->population;
    std::get<9>(temp) = this->emetric;
    return temp;
}

//Run Configuration parameters

void PathManager::setApproxAucSamples(size_t samples){
    if(samples <= 0){
        this->approxAucSamples = -1;
    }
    else{
        this->approxAucSamples = samples;
    }
    
}

void PathManager::setDumpAucPoints(bool dumpAucPoints){
    this->dumpAucPoints = dumpAucPoints;
}

void PathManager::setDumpAuPRPoints(bool dumpAuPRPoints){
    this->dumpAuPRPoints = dumpAuPRPoints;
}

void PathManager::setDumpAvgPrecPoints(bool dumpAvgPrecPoints){
    this->dumpAvgPrecPoints = dumpAvgPrecPoints;
}

void PathManager::setDumpRanking(bool dumpRanking){
    this->dumpRanking = dumpRanking;
}

void PathManager::setMeasures(std::vector<std::string> measures){
        this->measures = std::vector<std::string>(measures);
}

void PathManager::setMu1(std::vector<double> mu1){
        this->mu1 = std::vector<double>(mu1);
}

void PathManager::setMu2(std::vector<double> mu2){
        this->mu2 = std::vector<double>(mu2);
}

void PathManager::setMu3(std::vector<double> mu3){
        this->mu3 = std::vector<double>(mu3);
}

void PathManager::setMu4(std::vector<double> mu4){
        this->mu4 = std::vector<double>(mu4);
}

void PathManager::setMetrics(std::vector<std::string> metrics){
    this->metrics = std::vector<std::string>(metrics);
}

void PathManager::setCustomResultsFilePath(std::string customResultsFilePath){
    this->customResultsFilePath = customResultsFilePath;
}

size_t PathManager::getApproxAucSamples(){
    return this->approxAucSamples;
}

bool PathManager::getDumpAucPoints(){
    return this->dumpAucPoints;
}

bool PathManager::getDumpAuPRPoints(){
    return this->dumpAuPRPoints;
}

bool PathManager::getDumpAvgPrecPoints(){
    return this->dumpAvgPrecPoints;
}

bool PathManager::getDumpRanking(){
    return this->dumpRanking;
}

std::vector<std::string> PathManager::getMeasures(){
    return this->measures;
}

std::vector<double> PathManager::getMu1(){
    return this->mu1;
}

std::vector<double> PathManager::getMu2(){
    return this->mu2;
}

std::vector<double> PathManager::getMu3(){
    return this->mu3;
}

std::vector<double> PathManager::getMu4(){
    return this->mu4;
}

std::vector<std::string> PathManager::getMetrics(){
    return this->metrics;
}

void PathManager::setConfig(std::string basepath, std::string correlations_folder, std::string dataset, int validation_fold, int test_fold, std::string population, std::string emetric, int metacorr, int strategy, int genmax, double F, double CR){
    this->base_path=basepath;
    this->validationFold = validation_fold;
    this->testFold = test_fold;
    this->dataset = dataset;
    this->metacorr = metacorr;
    this->strategy = strategy;
    this->genmax = genmax;
    this->F = F;
    this->CR = CR;
    this->population = population;
    this->emetric = emetric;
    this->correlations_file = correlations_folder + "correlations_" + std::to_string(metacorr) + "_" + population+".csv";
    this->nodes_path = base_path + dataset + "/" + dataset +"_nodes";
    this->training_set_path = base_path + dataset+"/"+dataset+"_"+std::to_string(test_fold)+"/"+dataset+"_"+std::to_string(test_fold)+"_"+std::to_string(validation_fold)+"/"+dataset+"_train_"+std::to_string(test_fold)+"_"+std::to_string(validation_fold);
    this->validation_set_path = base_path + dataset+"/"+dataset+"_"+std::to_string(test_fold)+"/"+dataset+"_"+std::to_string(test_fold)+"_"+std::to_string(validation_fold)+"/"+dataset+"_validation_"+std::to_string(test_fold)+"_"+std::to_string(validation_fold);
    this->training_validation_set_path = base_path + dataset+"/"+dataset+"_"+std::to_string(test_fold)+"/"+dataset+"_train_val_"+std::to_string(test_fold);
    this->test_set_path = base_path + dataset+"/"+dataset+"_"+std::to_string(test_fold)+"/"+dataset+"_test_"+std::to_string(test_fold);
    std::cout<<"Correlations file path set to: "<<correlations_file<<std::endl;
    std::cout<<"Nodes file path set to: "<<nodes_path<<std::endl;
    std::cout<<"Training set file path set to: "<<training_set_path<<std::endl;
}

void PathManager::readConfig(std::string configPath){
    std::ifstream t(configPath);
    std::stringstream buffer;
    buffer << t.rdbuf();
    json j3 = json::parse(buffer.str());
	std::cout << j3.dump(4) << std::endl;
    this->dataset = j3["dataset"].get<std::string>();
    this->base_path = j3["basePath"].get<std::string>();
    this->testFold = j3["test_fold"];
    std::vector<std::string> metrics = j3["metrics"];
    this->nodes_path = this->base_path + this->dataset + "/" + this->dataset +"_nodes";
    this->training_set_path = base_path + dataset+"/"+dataset+"_"+std::to_string(testFold)+"/"+dataset+"_train_val_"+std::to_string(testFold);
    this->test_set_path = base_path + dataset+"/"+dataset+"_"+std::to_string(testFold)+"/"+dataset+"_test_"+std::to_string(testFold);
    if(!(j3["randomGenSeed"].is_null())){
        unsigned int randomGenSeed = (unsigned int)std::stoul(j3["randomGenSeed"].get<std::string>());
        Random* ran = Random::getInstance(0, 0, 0, randomGenSeed);
    }
    if(j3["measures"].is_null()){
        std::cout<<"Error! No measures indicated in config file"<<std::endl;
        std::cout<<"Exiting"<<std::endl;
        exit(1);
    }
    else{
        this->setMeasures(j3["measures"]);
    }
    if(std::find(this->measures.begin(), this->measures.end(), "mu1") != this->measures.end()) {
        this->setMu1(j3["mu1"]);
    }
    if(std::find(this->measures.begin(), this->measures.end(), "mu2") != this->measures.end()) {
        this->setMu2(j3["mu2"]);
    }
    if(std::find(this->measures.begin(), this->measures.end(), "mu3") != this->measures.end()) {
        this->setMu3(j3["mu3"]);
    }
    if(std::find(this->measures.begin(), this->measures.end(), "mu4") != this->measures.end()) {
        this->setMu4(j3["mu4"]);
    }
    if(j3["metrics"].is_null()){
        std::cout<<"Error! No metrics indicated in config file"<<std::endl;
        std::cout<<"Exiting"<<std::endl;
        exit(1);
    }
    else{
        this->setMetrics(j3["metrics"]);
    }

    if(j3["dumpAucPoints"].is_null()){
        this->setDumpAucPoints(false);
    }
    else{
        this->setDumpAucPoints(j3["dumpAucPoints"]);
    }
    if(j3["dumpAuPRPoints"].is_null()){
        this->setDumpAuPRPoints(false);
    }
    else{
        this->setDumpAuPRPoints(j3["dumpAuPRPoints"]);
    }
    if(j3["dumpRanking"].is_null()){
        this->setDumpRanking(false);
    }
    else{
        this->setDumpRanking(j3["dumpRanking"]);
    }

    if(j3["approxAucSamples"].is_null()){
        this->setApproxAucSamples(-1);
    }
    else{
        this->setApproxAucSamples(j3["approxAucSamples"]);
    }
}
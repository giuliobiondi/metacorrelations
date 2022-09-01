#include "Classifier.h"
#include <string>
#include <vector>
#include "metrics.h"

namespace Classifier
{

inline double accuracy_score(long tp, long tn, long p, long n){
    return (tp+tn)/(double)(p+n);
}

inline double precision_score(long tp, long fp){
    return tp/(double)(tp+fp);
}

inline double recall_score(long tp, long fn){
    return tp/(double)(tp+fn);
}

inline double F1_score(double precision, double recall){
    return (2*precision*recall)/(double)(precision+recall);
}
resultsRow classify(PUNGraph trainingSet, PUNGraph testSet, rankedGraph &ranking, int fold, int k, std::vector<std::string> metrics){
    PathManager &pathmng = PathManager::getInstance();
    resultsRow results;
    bool empty = true;
    for (size_t i = 0; i<metrics.size(); i++){
        if (metrics[i] == "precision"){
            if(empty){
                results =  Metrics::precision(ranking, testSet->GetEdges());
                empty = false;
            }
            
        }
        else if (metrics[i] == "auc"){
            resultsRow aucResults = Metrics::auc(ranking, testSet->GetEdges(), 1);
            if(empty){
                results =  aucResults;
                empty = false;
            }
            else{
                std::get<12>(results) = std::get<12>(aucResults);
            }
        }
        else if (metrics[i] == "approxAuc"){
            resultsRow approxAucResults = Metrics::approx_auc(ranking, testSet->GetEdges(), 1);
            if(empty){
                results =  approxAucResults;
                empty = false;
            }
            else{
                std::get<13>(results) = std::get<13>(approxAucResults);
            }
        }
        else if (metrics[i] == "aupr"){
            resultsRow AuPRResults = Metrics::aupr(ranking, testSet->GetEdges(), 1);
            if(empty){
                results =  AuPRResults;
                empty = false;
            }
            else{
                std::get<15>(results) = std::get<15>(AuPRResults);
            }
        }
        else if (metrics[i] == "avgPrec"){
            resultsRow AvgPrecResults = Metrics::average_precision(ranking, testSet->GetEdges(), 1);
            if(empty){
                results =  AvgPrecResults;
                empty = false;
            }
            else{
                std::get<16>(results) = std::get<16>(AvgPrecResults);
            }
        }
    }
    std::get<0>(results) = pathmng.getMeasureName();
    std::get<9>(results) = trainingSet->GetEdges();
    std::get<10>(results) = testSet->GetEdges();
    return results;
}

resultsRow classify(const PUNGraph trainingSet,const PUNGraph testSet, DERanking &ranking, int fold, int k, std::vector<std::string> metrics){
    PathManager &pathmng = PathManager::getInstance();
    //std::vector<std::string> metrics = pathmng.getMetrics();
    resultsRow results;
    bool empty = true;
    for (size_t i = 0; i<metrics.size(); i++){
        if (metrics[i] == "precision"){
            if(empty){
                results =  Metrics::precision(ranking, testSet->GetEdges());
                empty = false;
            }
            
        }
        else if (metrics[i] == "auc"){
            resultsRow aucResults = Metrics::auc(ranking, testSet->GetEdges(), 1);
            if(empty){
                results =  aucResults;
                empty = false;
            }
            else{
                std::get<12>(results) = std::get<12>(aucResults);
            }
        }
        else if (metrics[i] == "aupr"){
            resultsRow AuPRResults = Metrics::aupr(ranking, testSet->GetEdges(), 1);
            if(empty){
                results =  AuPRResults;
                empty = false;
            }
            else{
                std::get<15>(results) = std::get<15>(AuPRResults);
            }
        }
        else if (metrics[i] == "avgPrec"){
            resultsRow AvgPrecResults = Metrics::average_precision(ranking, testSet->GetEdges(), 1);
            if(empty){
                results =  AvgPrecResults;
                empty = false;
            }
            else{
                std::get<16>(results) = std::get<16>(AvgPrecResults);
            }
        }
        else if (metrics[i] == "approxAuc"){
            resultsRow approxAucResults = Metrics::approx_auc(ranking, testSet->GetEdges(), 1);
            if(empty){
                results =  approxAucResults;
                empty = false;
            }
            else{
                std::get<13>(results) = std::get<13>(approxAucResults);
                std::get<14>(results) = std::get<14>(approxAucResults);
            }
        }
    }

    std::get<0>(results) = pathmng.getMeasureName();
    std::get<9>(results) = trainingSet->GetEdges();
    std::get<10>(results) = testSet->GetEdges();
    return results;
}

void dumpEvaluationResults(const std::string& path, const std::vector<resultsRow>& results){
    PathManager &pathmng = PathManager::getInstance();
    std::ofstream ofile;
    //Insert header if file does not exist
    if(!Filesystem::exists(path)){
        ofile.open (path, std::ofstream::app);
        ofile<<"dataset,validationFold,testFold,metacorrelation,DE strategy,generations,F,CR,measure,tp,tn,fp,fn,accuracy,precision,recall,F1,training set edges,test set edges,Epot,auc,approxAuc,approxAucseed,aupr,avgPrec,coefficients\n";
        ofile.close();
    }
    ofile.open (path, std::ofstream::app);
    configRow config = pathmng.getConfig();
    for(auto const& value: results) {
        ofile<<std::get<0>(config)<<",";
        ofile<<std::get<1>(config)<<",";
        ofile<<std::get<2>(config)<<",";
        ofile<<std::get<3>(config)<<",";
        ofile<<std::get<4>(config)<<",";
        ofile<<std::get<5>(config)<<",";
        ofile<<std::get<6>(config)<<",";
        ofile<<std::get<7>(config)<<",";
        ofile<<std::setprecision(20)<<std::get<0>(value)<<",";
        ofile<<std::get<1>(value)<<",";
        ofile<<std::get<2>(value)<<",";
        ofile<<std::get<3>(value)<<",";
        ofile<<std::get<4>(value)<<",";
        ofile<<std::get<5>(value)<<",";
        ofile<<std::get<6>(value)<<",";
        ofile<<std::get<7>(value)<<",";
        ofile<<std::get<8>(value)<<",";
        ofile<<std::get<9>(value)<<",";
        ofile<<std::get<10>(value)<<",";
        ofile<<std::get<11>(value)<<",";
        ofile<<std::get<12>(value)<<",";
        ofile<<std::get<13>(value)<<",";
        ofile<<std::get<14>(value)<<",";
        ofile<<std::get<15>(value)<<",";
        ofile<<std::get<16>(value)<<std::endl;
    }
    ofile.close();
}
void dumpMeasuresEvaluationResults(const std::string& path, const std::vector<resultsRow>& results){
    PathManager &pathmng = PathManager::getInstance();
    std::ofstream ofile;
    //Insert header if file does not exist
    if(!Filesystem::exists(path)){
        ofile.open (path, std::ofstream::app);
        ofile<<"dataset,testFold,measure,tp,tn,fp,fn,accuracy,precision,recall,F1,training set edges,test set edges,Epot,auc,approxAuc,approxAucseed,aupr,avgPrec,coefficients\n";
        ofile.close();
    }
    ofile.open(path, std::ofstream::app);
    configRow config = pathmng.getConfig();
    for(auto const& value: results) {
        ofile<<std::get<0>(config)<<",";
        ofile<<std::get<2>(config)<<",";
        ofile<<std::setprecision(20)<<std::get<0>(value)<<",";
        ofile<<std::get<1>(value)<<",";
        ofile<<std::get<2>(value)<<",";
        ofile<<std::get<3>(value)<<",";
        ofile<<std::get<4>(value)<<",";
        ofile<<std::get<5>(value)<<",";
        ofile<<std::get<6>(value)<<",";
        ofile<<std::get<7>(value)<<",";
        ofile<<std::get<8>(value)<<",";
        ofile<<std::get<9>(value)<<",";
        ofile<<std::get<10>(value)<<",";
        ofile<<std::get<11>(value)<<",";
        ofile<<std::get<12>(value)<<",";
        ofile<<std::get<13>(value)<<",";
        ofile<<std::get<14>(value)<<",";
        ofile<<std::get<15>(value)<<",";
        ofile<<std::get<16>(value);
        if(std::get<0>(value)=="mu1"){
            ofile<<",";
            std::vector<double> coefficients = pathmng.getMu1();
            for(auto const& coeff: coefficients){
                ofile<<std::setprecision(20)<<coeff<<";";
            }
        }
        if(std::get<0>(value)=="mu2"){
            ofile<<",";
            std::vector<double> coefficients = pathmng.getMu2();
            for(auto const& coeff: coefficients){
                ofile<<std::setprecision(20)<<coeff<<";";
            }
        }
        ofile<<std::endl;
    }
    ofile.close();
}
    
}
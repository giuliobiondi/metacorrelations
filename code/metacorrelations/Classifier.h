#pragma once
#include "headers.h"

namespace Classifier 
{
    resultsRow classify(const PUNGraph trainingSet,const PUNGraph testSet, rankedGraph &ranking, int fold, int k, std::vector<std::string> metrics);
    resultsRow classify(const PUNGraph trainingSet,const PUNGraph testSet, DERanking &ranking, int fold, int k, std::vector<std::string> metrics);
    void dumpEvaluationResults(const std::string& path, const std::vector<resultsRow>& results);
    void dumpMeasuresEvaluationResults(const std::string& path, const std::vector<resultsRow>& results);
}

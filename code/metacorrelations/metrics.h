#pragma once
#include "headers.h"
namespace Metrics{
    
    resultsRow auc(rankedGraph &ranking, size_t positiveSamples, int positive_class_label = 1);

    resultsRow aupr(rankedGraph &ranking, size_t positiveSamples, int positive_class_label = 1);
    
    resultsRow average_precision(rankedGraph &ranking, size_t positiveSamples, int positive_class_label = 1);

    resultsRow approx_auc(rankedGraph &ranking, size_t positiveSamples, int positive_class_label = 1);
    
    resultsRow precision(rankedGraph &ranking, uint positiveSamples, int positive_class_label=1);

    resultsRow auc(DERanking &ranking, size_t positiveSamples, int positive_class_label = 1);

    resultsRow aupr(DERanking &ranking, size_t positiveSamples, int positive_class_label = 1);
    
    resultsRow average_precision(DERanking &ranking, size_t positiveSamples, int positive_class_label = 1);
    
    resultsRow approx_auc(DERanking &ranking, size_t positiveSamples, int positive_class_label = 1);
    
    resultsRow precision(DERanking &ranking, uint positiveSamples, int positive_class_label= 1 );
}
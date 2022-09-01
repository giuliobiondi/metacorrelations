#include "metrics.h"
#include <exception>

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

resultsRow Metrics::auc(rankedGraph &ranking, size_t positiveSamples, int positive_class_label){
    PathManager &pathmng = PathManager::getInstance();
    /*Calculate the number of positive and negative samples*/
    size_t p = positiveSamples, n = ranking.size() - positiveSamples;
    std::vector<double> tpr,fpr;
    size_t tp = 0, fp = 0, tn = 0, fn = 0;
    size_t tpcounter = 0, fpcounter = 0;
    size_t i;
    std::vector<double> thresholds;
    std::vector<int> tpcounters;
    std::vector<int> fpcounters;
    thresholds.push_back(0.0);
    tpcounters.push_back(0);
    fpcounters.push_back(0);
    tpr.push_back(0);
    fpr.push_back(0);
    double prev_score = std::get<2>(ranking[0]);
    if(std::get<3>(ranking[0]) == positive_class_label){
        tpcounter++;
    }
    else{
        fpcounter++;
    }
    for(i = 1; i < ranking.size(); i++){
        if (std::get<2>(ranking[i]) != prev_score){
            thresholds.push_back(prev_score);
            tpcounters.push_back(tpcounter);
            fpcounters.push_back(fpcounter);
            tpr.push_back(tpcounter / (double)p);
            fpr.push_back(fpcounter / (double)n);
            prev_score = std::get<2>(ranking[i]);
        }
        if(i<positiveSamples){
            if(std::get<3>(ranking[i]) == positive_class_label){
            tp++;
            tpcounter++;
            }
            else{
                fp++;
                fpcounter++;
            }
        }
        else{
            if(std::get<3>(ranking[i]) == positive_class_label){
            fn++;
            tpcounter++;
            }
            else{
                tn++;
                fpcounter++;
            }
        }
    }
    //Point 1,1
    tpr.push_back(tpcounter / (double)p);
    fpr.push_back(fpcounter / (double)n);
    tpcounters.push_back(tpcounter);
    fpcounters.push_back(fpcounter);
    thresholds.push_back(std::get<3>(ranking[ranking.size()-1]));

    /*Calculate AUC*/
    double auc_score = 0;
    size_t size = tpr.size();
    double q1,q2,p1,p2;
    q1 = fpr[0];
    q2 = tpr[0];
    for(size_t i = 1;i < size;++i){
        p1 = fpr[i];
        p2 = tpr[i];
        auc_score += (p1-q1) * (q2+p2) * 0.5;
        q1=p1;
        q2=p2;   
    }
    /*Dump the ROC points*/
    if (pathmng.getDumpAucPoints()){
        try{
            std::ofstream ofile;
            ofile.exceptions (std::ifstream::failbit | std::ifstream::badbit);
            ofile.open (pathmng.getAucPath(), std::ofstream::app);
            for (size_t i = 0; i < tpr.size(); i++){
                ofile<<std::setprecision(20)<<fpr[i]<<","<<tpr[i]<<std::endl;
            }
        }
        catch(std::ifstream::failure &readErr) {
        std::cerr << "\n\nException occured saving the roc points file\n"
              << readErr.what()
              << std::endl;
        }
    }

    double accuracy = accuracy_score(tp, tn, positiveSamples, ranking.size() - positiveSamples);
    double precision = precision_score(tp, fp);
    double recall = recall_score(tp, fn);
    double F1 = F1_score(precision, recall);

    resultsRow result = std::make_tuple(
        "",
        tp,
        tn,
        fp,
        fn,
        accuracy,
        precision,
        recall,
        F1,
        -1,
        -1,
        ranking.size(),
        auc_score,
        -1.0,
        -1.0,
        0.0,
        0.0
    );
    return result;
}

resultsRow Metrics::aupr(rankedGraph &ranking, size_t positiveSamples, int positive_class_label){
    PathManager &pathmng = PathManager::getInstance();
    /*Calculate the number of positive and negative samples*/
    size_t p = positiveSamples, n = ranking.size() - positiveSamples;
    std::vector<double> precisionVec, recallVec;
    size_t tp = 0, fp = 0, tn = 0, fn = 0;
    size_t tpcounter = 0, fpcounter = 0;
    size_t i;
    std::vector<double> thresholds;
    //add first point (0,1)
    precisionVec.push_back(1.0);
    recallVec.push_back(0.0);
    //initialize the threshold to the first score
    double prev_score = std::get<2>(ranking[0]);
    if(std::get<3>(ranking[0]) == positive_class_label){
        tpcounter++;
    }
    else{
        fpcounter++;
    }
    //loop
    for(i = 1; i < ranking.size(); i++){
        if (std::get<2>(ranking[i]) != prev_score){
            thresholds.push_back(prev_score);
            recallVec.push_back(tpcounter / (double)p);
            precisionVec.push_back(tpcounter / (double)i);
            prev_score = std::get<2>(ranking[i]);
        }
        if(i<positiveSamples){
            if(std::get<3>(ranking[i]) == positive_class_label){
            tp++;
            tpcounter++;
            }
            else{
                fp++;
                fpcounter++;
            }
        }
        else{
            if(std::get<3>(ranking[i]) == positive_class_label){
            fn++;
            tpcounter++;
            }
            else{
                tn++;
                fpcounter++;
            }
        }
    }
    //Point 1,0
    precisionVec.push_back(tpcounter/(double)ranking.size());
    recallVec.push_back(1.0);
    /*Calculate AUC*/
    double aupr_score = 0;
    size_t size = precisionVec.size();
    double q1,q2,p1,p2;
    q1 = recallVec[0];
    q2 = precisionVec[0];
    for(size_t i = 1;i < size;++i){
        p1 = recallVec[i];
        p2 = precisionVec[i];
        aupr_score += abs((p1-q1) * (q2+p2) * 0.5);
        q1=p1;
        q2=p2;
    }
    /*Dump the AUPR points*/
    if (pathmng.getDumpAuPRPoints()){
        try{
            std::ofstream ofile;
            ofile.exceptions (std::ifstream::failbit | std::ifstream::badbit);
            ofile.open (pathmng.getAuPRPath(), std::ofstream::app);
            for (size_t i = 0; i < precisionVec.size(); i++){
                ofile<<std::setprecision(20)<<recallVec[i]<<","<<precisionVec[i]<<std::endl;
            }
        }
        catch(std::ifstream::failure &readErr) {
        std::cerr << "\n\nException occured saving the roc points file\n"
              << readErr.what()
              << std::endl;
        }
    }

    double accuracy = accuracy_score(tp, tn, positiveSamples, ranking.size() - positiveSamples);
    double precision = precision_score(tp, fp);
    double recall = recall_score(tp, fn);
    double F1 = F1_score(precision, recall);

    resultsRow result = std::make_tuple(
        "",
        tp,
        tn,
        fp,
        fn,
        accuracy,
        precision,
        recall,
        F1,
        -1,
        -1,
        ranking.size(),
        0.0,
        -1.0,
        -1.0,
        aupr_score,
        0.0
    );
    return result;
}

resultsRow Metrics::average_precision(rankedGraph &ranking, size_t positiveSamples, int positive_class_label){
    PathManager &pathmng = PathManager::getInstance();
    /*Calculate the number of positive and negative samples*/
    size_t p = positiveSamples, n = ranking.size() - positiveSamples;
    std::vector<double> precisionVec, recallVec;
    size_t tp = 0, fp = 0, tn = 0, fn = 0;
    size_t tpcounter = 0, fpcounter = 0;
    size_t i;
    std::vector<double> thresholds;
    //add first point (0,1)
    precisionVec.push_back(1.0);
    recallVec.push_back(0.0);
    //initialize the threshold to the first score
    double prev_score = std::get<2>(ranking[0]);
    if(std::get<3>(ranking[0]) == positive_class_label){
        tpcounter++;
    }
    else{
        fpcounter++;
    }
    //loop
    for(i = 1; i < ranking.size(); i++){
        if (std::get<2>(ranking[i]) != prev_score){
            thresholds.push_back(prev_score);
            recallVec.push_back(tpcounter / (double)p);
            precisionVec.push_back(tpcounter / (double)i);
            prev_score = std::get<2>(ranking[i]);
        }
        if(i<positiveSamples){
            if(std::get<3>(ranking[i]) == positive_class_label){
            tp++;
            tpcounter++;
            }
            else{
                fp++;
                fpcounter++;
            }
        }
        else{
            if(std::get<3>(ranking[i]) == positive_class_label){
            fn++;
            tpcounter++;
            }
            else{
                tn++;
                fpcounter++;
            }
        }
    }
    //Point 1,0
    precisionVec.push_back(tpcounter/(double)ranking.size());
    recallVec.push_back(1.0);
    /*Calculate AUC*/
    double average_precision_score = 0;
    size_t size = precisionVec.size();
    double q1,q2,p1,p2;
    q1 = recallVec[0];
    q2 = precisionVec[0];
    for(size_t i = 1;i < size;++i){
        average_precision_score += (recallVec[i] - recallVec[i-1]) * precisionVec[i];
    }
    /*Dump the ROC points*/
    if (pathmng.getDumpAvgPrecPoints()){
        try{
            std::ofstream ofile;
            ofile.exceptions (std::ifstream::failbit | std::ifstream::badbit);
            ofile.open (pathmng.getAvgPrecPath(), std::ofstream::app);
            for (size_t i = 0; i < precisionVec.size(); i++){
                ofile<<std::setprecision(20)<<recallVec[i]<<","<<precisionVec[i]<<std::endl;
            }
        }
        catch(std::ifstream::failure &readErr) {
        std::cerr << "\n\nException occured saving the roc points file\n"
              << readErr.what()
              << std::endl;
        }
    }

    double accuracy = accuracy_score(tp, tn, positiveSamples, ranking.size() - positiveSamples);
    double precision = precision_score(tp, fp);
    double recall = recall_score(tp, fn);
    double F1 = F1_score(precision, recall);

    resultsRow result = std::make_tuple(
        "",
        tp,
        tn,
        fp,
        fn,
        accuracy,
        precision,
        recall,
        F1,
        -1,
        -1,
        ranking.size(),
        0.0,
        -1.0,
        -1.0,
        0.0,
        average_precision_score
    );
    return result;
}

resultsRow Metrics::approx_auc(rankedGraph &ranking, size_t positiveSamples, int positive_class_label){
    PathManager &pathmng = PathManager::getInstance();
    size_t tp=0, tn=0, fp=0, fn=0, p = positiveSamples, n = ranking.size() - positiveSamples, approxAucSamples = (pathmng.getApproxAucSamples() == -1 ? positiveSamples : pathmng.getApproxAucSamples());
    int id1, id2;
    size_t i, j;
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    t1 = std::chrono::high_resolution_clock::now();
    std::vector<size_t> positives;
#pragma omp parallel for private(id1,id2) lastprivate(i) reduction(+:tp,fp)
    for (i = 0; i < p; i++){
        id1=std::get<0>(ranking[i]);
        id2=std::get<1>(ranking[i]);
        if(std::get<3>(ranking[i])){
            tp++;
#pragma omp critical
            {
				positives.push_back(i);
            }
        }
        else{
            fp++;
        }
    }
#pragma omp parallel for private(id1,id2) reduction(+:fn,tn)
    for (j = i; j < ranking.size(); j++){
        id1=std::get<0>(ranking[j]);
        id2=std::get<1>(ranking[j]);
        if(std::get<3>(ranking[j])){
            fn++;
#pragma omp critical
            {
				positives.push_back(j);
            }
        }
        else{
            tn++;
        }
    }
    //Sort the positives indexes array
    #ifdef _MSC_VER
            concurrency::parallel_sort(positives.begin(), positives.end());
    #else
            __gnu_parallel::stable_sort(positives.begin(), positives.end());
    #endif
    t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count();
    double accuracy = accuracy_score(tp, tn, p, n);
    double precision = precision_score(tp, fp);
    double recall = recall_score(tp, fn);
    double F1 = F1_score(precision, recall);
    double auc = 0;
    int missingIndex, nonExistentIndex;
    double missingScore, nonExistentScore;
    //std::random_device rd;  //Will be used to obtain a seed for the random number engine
    //std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::mt19937& gen = Random::getInstance()->aucGen;
    auto aucSeed = Random::getInstance()->getAucSeed();
    std::uniform_int_distribution<> pos_dis(0, positives.size()-1);
    std::uniform_int_distribution<> ranking_dis(0, ranking.size()-1);
    try{
    for (int n = 0; n < approxAucSamples; ++n){
        //randomly select a sample from the positives and get the score
        missingIndex = pos_dis(gen);
        missingScore = std::get<2>(ranking[positives[missingIndex]]);
        nonExistentIndex = ranking_dis(gen);
        id1=std::get<0>(ranking[nonExistentIndex]);
        id2=std::get<1>(ranking[nonExistentIndex]);
        while(std::get<3>(ranking[nonExistentIndex])){
            nonExistentIndex = ranking_dis(gen);
            id1=std::get<0>(ranking[nonExistentIndex]);
            id2=std::get<1>(ranking[nonExistentIndex]);
        }
        nonExistentScore = std::get<2>(ranking[nonExistentIndex]);
        if(missingScore > nonExistentScore){
            auc += 1.0;
        }
        else if(missingScore == nonExistentScore){
            auc += 0.5;
        }
    }
    }
    catch (std::exception& e){
        std::cout<<e.what()<<std::endl;
    }
    auc /= approxAucSamples;
    resultsRow result = std::make_tuple(
        "",
        tp,
        tn,
        fp,
        fn,
        accuracy,
        precision,
        recall,
        F1,
        -1,
        -1,
        ranking.size(),
        -1.0,
        auc,
        aucSeed,
        0.0,
        0.0
    );
    return result;
}

resultsRow Metrics::precision(rankedGraph &ranking, uint positiveSamples, int positive_class_label){
    size_t tp = 0, tn = 0, fp = 0, fn = 0, i;
    //#pragma omp parallel for lastprivate(i) reduction(+:tp,fp)
    for (i = 0; i < positiveSamples; i++){
        if(std::get<3>(ranking[i])==positive_class_label){
            tp++;
        }
        else{
            fp++;
        }
    }
    /*#pragma omp parallel for reduction(+:fn,tn)
    for (size_t j = i; j < ranking.size(); j++){
        if(std::get<3>(ranking[j])==positive_class_label){
            fn++;
        }
        else{
            tn++;
        }
    }*/
    fn=positiveSamples-tp;
    tn=ranking.size()-(positiveSamples+fn);
    double accuracy = accuracy_score(tp, tn, positiveSamples, ranking.size()-positiveSamples);
    double precision = precision_score(tp, fp);
    double recall = recall_score(tp, fn);
    double F1 = F1_score(precision, recall);
    resultsRow result = std::make_tuple(
        "",
        tp,
        tn,
        fp,
        fn,
        accuracy,
        precision,
        recall,
        F1,
        -1,
        -1,
        ranking.size(),
        -1.0,
        -1.0,
        -1.0,
        0.0,
        0.0
    );
    return result;
}


resultsRow Metrics::auc(DERanking &ranking, size_t positiveSamples, int positive_class_label){
    PathManager &pathmng = PathManager::getInstance();
    /*Calculate the number of positive and negative samples*/
    size_t p = positiveSamples, n = ranking.size() - positiveSamples;
    std::vector<double> tpr,fpr;
    size_t tp = 0, fp = 0, tn = 0, fn = 0;
    size_t tpcounter = 0, fpcounter = 0;
    size_t i = 0;
    std::vector<double> thresholds;
    std::vector<int> tpcounters;
    std::vector<int> fpcounters;
    thresholds.push_back(0.0);
    tpcounters.push_back(0);
    fpcounters.push_back(0);
    tpr.push_back(0);
    fpr.push_back(0);
    double prev_score = std::get<9>(ranking[0]);
    if(std::get<10>(ranking[0]) == positive_class_label){
        tp++;
        tpcounter++;
    }
    else{
        fp++;
        fpcounter++;
    }
    for(i = 1; i < ranking.size(); i++){
        if (std::get<9>(ranking[i]) != prev_score){
                thresholds.push_back(prev_score);
                tpcounters.push_back(tpcounter);
                fpcounters.push_back(fpcounter);
                tpr.push_back(tpcounter / (double)p);
                fpr.push_back(fpcounter / (double)n);
                prev_score = std::get<9>(ranking[i]);
        }
        if(i<positiveSamples){
            if(std::get<10>(ranking[i]) == positive_class_label){
            tp++;
            tpcounter++;
            }
            else{
                fp++;
                fpcounter++;
            }
        }
        else{
            if(std::get<10>(ranking[i]) == positive_class_label){
            fn++;
            tpcounter++;
            }
            else{
                tn++;
                fpcounter++;
            }
        }
    }
    //Point 1,1
    tpr.push_back(tpcounter / (double)p);
    fpr.push_back(fpcounter / (double)n);
    tpcounters.push_back(tpcounter);
    fpcounters.push_back(fpcounter);
    thresholds.push_back(std::get<9>(ranking[ranking.size()-1]));

    /*Calculate AUC*/
    double auc_score = 0;
    size_t size = tpr.size();
    double q1,q2,p1,p2;
    q1 = fpr[0];
    q2 = tpr[0];
    for(size_t i = 1;i < size;++i){
        p1 = fpr[i];
        p2 = tpr[i];
        auc_score += (p1-q1) * (q2+p2) * 0.5;
        q1=p1;
        q2=p2;   
    }
    /*Dump the ROC points*/
    if (pathmng.getDumpAucPoints()){
        try{
            std::ofstream ofile;
            ofile.exceptions (std::ifstream::failbit | std::ifstream::badbit);
            ofile.open (pathmng.getAucPath(), std::ofstream::app);
            for (size_t i = 0; i < tpr.size(); i++){
                ofile<<std::setprecision(20)<<fpr[i]<<","<<tpr[i]<<","<<thresholds[i]<<","<<fpcounters[i]<<","<<tpcounters[i]<<std::endl;
            }
        }
        catch(std::ifstream::failure &readErr) {
        std::cerr << "\n\nException occured saving the roc points file\n"
              << readErr.what()
              << std::endl;
        }
    }

    double accuracy = accuracy_score(tp, tn, positiveSamples, ranking.size() - positiveSamples);
    double precision = precision_score(tp, fp);
    double recall = recall_score(tp, fn);
    double F1 = F1_score(precision, recall);

    resultsRow result = std::make_tuple(
        "",
        tp,
        tn,
        fp,
        fn,
        accuracy,
        precision,
        recall,
        F1,
        -1,
        -1,
        ranking.size(),
        auc_score,
        -1.0,
        -1.0,
        0.0,
        0.0
    );
    return result;
}

resultsRow Metrics::aupr(DERanking &ranking, size_t positiveSamples, int positive_class_label){
    PathManager &pathmng = PathManager::getInstance();
    /*Calculate the number of positive and negative samples*/
    size_t p = positiveSamples, n = ranking.size() - positiveSamples;
    std::vector<double> precisionVec, recallVec;
    size_t tp = 0, fp = 0, tn = 0, fn = 0;
    size_t tpcounter = 0, fpcounter = 0;
    size_t i;
    //add first point (0,1)
    std::vector<double> thresholds;
    precisionVec.push_back(1.0);
    recallVec.push_back(0.0);
    thresholds.push_back(0.0);
    //initialize the threshold to the first score
    double prev_score = std::get<9>(ranking[0]);
    if(std::get<10>(ranking[0]) == positive_class_label){
        tp++;
        tpcounter++;
    }
    else{
        fp++;
        fpcounter++;
    }
    //loop
    for(i = 1; i < ranking.size(); i++){
        if (std::get<9>(ranking[i]) != prev_score){
                thresholds.push_back(prev_score);
                recallVec.push_back(tpcounter / (double)p);
                precisionVec.push_back(tpcounter / (double)i);
                prev_score = std::get<9>(ranking[i]);
        }
        if(i<positiveSamples){
            if(std::get<10>(ranking[i]) == positive_class_label){
            tp++;
            tpcounter++;
            }
            else{
                fp++;
                fpcounter++;
            }
        }
        else{
            if(std::get<10>(ranking[i]) == positive_class_label){
            fn++;
            tpcounter++;
            }
            else{
                tn++;
                fpcounter++;
            }
        }
    }
    //Point 1,0
    precisionVec.push_back(tpcounter/(double)ranking.size());
    recallVec.push_back(1.0);
    /*Calculate AUC*/
    double aupr_score = 0;
    size_t size = precisionVec.size();
    double q1,q2,p1,p2;
    q1 = recallVec[0];
    q2 = precisionVec[0];
    for(size_t i = 1;i < size;++i){
        p1 = recallVec[i];
        p2 = precisionVec[i];
        aupr_score += abs((p1-q1) * (q2+p2) * 0.5);
        q1=p1;
        q2=p2;   
    }
    /*Dump the AUPR points*/
    if (pathmng.getDumpAuPRPoints()){
        try{
            std::ofstream ofile;
            ofile.exceptions (std::ifstream::failbit | std::ifstream::badbit);
            ofile.open (pathmng.getAuPRPath(), std::ofstream::app);
            for (size_t i = 0; i < precisionVec.size(); i++){
                ofile<<std::setprecision(20)<<recallVec[i]<<","<<precisionVec[i]<<","<<thresholds[i]<<std::endl;
            }
        }
        catch(std::ifstream::failure &readErr) {
        std::cerr << "\n\nException occured saving the roc points file\n"
              << readErr.what()
              << std::endl;
        }
    }

    double accuracy = accuracy_score(tp, tn, positiveSamples, ranking.size() - positiveSamples);
    double precision = precision_score(tp, fp);
    double recall = recall_score(tp, fn);
    double F1 = F1_score(precision, recall);

    resultsRow result = std::make_tuple(
        "",
        tp,
        tn,
        fp,
        fn,
        accuracy,
        precision,
        recall,
        F1,
        -1,
        -1,
        ranking.size(),
        0.0,
        -1.0,
        -1.0,
        aupr_score,
        0.0
    );
    return result;
}

resultsRow Metrics::average_precision(DERanking &ranking, size_t positiveSamples, int positive_class_label){
     PathManager &pathmng = PathManager::getInstance();
    /*Calculate the number of positive and negative samples*/
    size_t p = positiveSamples, n = ranking.size() - positiveSamples;
    std::vector<double> precisionVec, recallVec;
    size_t tp = 0, fp = 0, tn = 0, fn = 0;
    size_t tpcounter = 0, fpcounter = 0;
    size_t i;
    //add first point (0,1)
    std::vector<double> thresholds;
    precisionVec.push_back(1.0);
    recallVec.push_back(0.0);
    thresholds.push_back(0.0);
    //initialize the threshold to the first score
    double prev_score = std::get<9>(ranking[0]);
    if(std::get<10>(ranking[0]) == positive_class_label){
        tp++;
        tpcounter++;
    }
    else{
        fp++;
        fpcounter++;
    }
    //loop
    for(i = 1; i < ranking.size(); i++){
        if (std::get<9>(ranking[i]) != prev_score){
                thresholds.push_back(prev_score);
                recallVec.push_back(tpcounter / (double)p);
                precisionVec.push_back(tpcounter / (double)i);
                prev_score = std::get<9>(ranking[i]);
        }
        if(i<positiveSamples){
            if(std::get<10>(ranking[i]) == positive_class_label){
            tp++;
            tpcounter++;
            }
            else{
                fp++;
                fpcounter++;
            }
        }
        else{
            if(std::get<10>(ranking[i]) == positive_class_label){
            fn++;
            tpcounter++;
            }
            else{
                tn++;
                fpcounter++;
            }
        }
    }
    //Point 1,0
    precisionVec.push_back(tpcounter/(double)ranking.size());
    recallVec.push_back(1.0);
    /*Calculate AUC*/
    double average_precision_score = 0;
    size_t size = precisionVec.size();
    double q1,q2,p1,p2;
    q1 = recallVec[0];
    q2 = precisionVec[0];
    for(size_t i = 1;i < size;++i){
        average_precision_score += (recallVec[i] - recallVec[i-1]) * precisionVec[i];
    }
    /*Dump the ROC points*/
    if (pathmng.getDumpAvgPrecPoints()){
        try{
            std::ofstream ofile;
            ofile.exceptions (std::ifstream::failbit | std::ifstream::badbit);
            ofile.open (pathmng.getAvgPrecPath(), std::ofstream::app);
            for (size_t i = 0; i < precisionVec.size(); i++){
                ofile<<std::setprecision(20)<<recallVec[i]<<","<<precisionVec[i]<<","<<thresholds[i]<<std::endl;
            }
        }
        catch(std::ifstream::failure &readErr) {
        std::cerr << "\n\nException occured saving the roc points file\n"
              << readErr.what()
              << std::endl;
        }
    }

    double accuracy = accuracy_score(tp, tn, positiveSamples, ranking.size() - positiveSamples);
    double precision = precision_score(tp, fp);
    double recall = recall_score(tp, fn);
    double F1 = F1_score(precision, recall);

    resultsRow result = std::make_tuple(
        "",
        tp,
        tn,
        fp,
        fn,
        accuracy,
        precision,
        recall,
        F1,
        -1,
        -1,
        ranking.size(),
        0.0,
        -1.0,
        -1.0,
        0.0,
        average_precision_score
    );
    return result;
}

resultsRow Metrics::approx_auc(DERanking &ranking, size_t positiveSamples, int positive_class_label){
    PathManager &pathmng = PathManager::getInstance();
    size_t tp = 0, tn = 0, fp = 0, fn = 0, p = positiveSamples, n = ranking.size() - positiveSamples, approxAucSamples = (pathmng.getApproxAucSamples() == -1 ? positiveSamples : pathmng.getApproxAucSamples());
    int id1, id2;
    size_t i = 0;
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    t1 = std::chrono::high_resolution_clock::now();
    std::vector<size_t> positives;
//#pragma omp parallel for private(id1,id2) lastprivate(i) reduction(+:tp,fp,p,n)
    for (i = 0; i < p; i++){
        id1=std::get<0>(ranking[i]);
        id2=std::get<1>(ranking[i]);
        if(std::get<10>(ranking[i]) == positive_class_label){
            tp++;
//#pragma omp critical
//            {
				positives.push_back(i);
//            }
        }
        else{
            fp++;
        }
    }
#pragma omp parallel for private(id1,id2) reduction(+:fn,tn)
    for (size_t j = i; j < ranking.size(); j++){
        id1=std::get<0>(ranking[j]);
        id2=std::get<1>(ranking[j]);
        if(std::get<10>(ranking[j]) == positive_class_label){
            fn++;
#pragma omp critical
            {
				positives.push_back(j);
            }
        }
        else{
            tn++;
        }
    }
    //Sort the positives indexes array
    #ifdef _MSC_VER
            concurrency::parallel_sort(positives.begin(), positives.end());
    #else
            __gnu_parallel::stable_sort(positives.begin(), positives.end());
    #endif

    t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count();
    double accuracy = accuracy_score(tp, tn, p, n);
    double precision = precision_score(tp, fp);
    double recall = recall_score(tp, fn);
    double F1 = F1_score(precision, recall);
    double auc = 0;
    int missingIndex, nonExistentIndex;
    double missingScore, nonExistentScore;
    //std::random_device rd;  //Will be used to obtain a seed for the random number engine
    //std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::mt19937& gen = Random::getInstance()->aucGen;
    auto aucSeed = Random::getInstance()->getAucSeed();
    std::uniform_int_distribution<> pos_dis(0, positives.size()-1);
    std::uniform_int_distribution<> ranking_dis(0, ranking.size()-1);
    for (int n = 0; n < approxAucSamples; ++n){
        //randomly select a sample from the positives and get the score
        missingIndex = pos_dis(gen);
        missingScore = std::get<9>(ranking[positives[missingIndex]]);
        nonExistentIndex = ranking_dis(gen);
        id1=std::get<0>(ranking[nonExistentIndex]);
        id2=std::get<1>(ranking[nonExistentIndex]);
        while(std::get<10>(ranking[nonExistentIndex])){
            nonExistentIndex = ranking_dis(gen);
            id1=std::get<0>(ranking[nonExistentIndex]);
            id2=std::get<1>(ranking[nonExistentIndex]);
        }
        nonExistentScore=std::get<9>(ranking[nonExistentIndex]);
        if(missingScore>nonExistentScore){
            auc += 1.0;
        }
        else if(missingScore==nonExistentScore){
            auc += 0.5;
        }
    }
    auc /= approxAucSamples;
    resultsRow result = std::make_tuple(
        "",
        tp,
        tn,
        fp,
        fn,
        accuracy,
        precision,
        recall,
        F1,
        -1,
        -1,
        ranking.size(),
        -1.0,
        auc,
        aucSeed,
        0.0,
        0.0
    );
    return result;
}

resultsRow Metrics::precision(DERanking &ranking, uint positiveSamples, int positive_class_label){
    /*Calculate the number of positive and negative samples*/
    size_t tp = 0, tn = 0, fp = 0, fn = 0, i=0;
    #pragma omp parallel for lastprivate(i) reduction(+:tp,fp)
    for (i = 0; i < positiveSamples; i++){
        if(std::get<10>(ranking[i])==positive_class_label){
            tp++;
        }
        else{
            fp++;
        }
    }
    #pragma omp parallel for reduction(+:fn,tn)
    for (size_t j = i; j < ranking.size(); j++){
        if(std::get<10>(ranking[j])==positive_class_label){
            fn++;
        }
        else{
            tn++;
        }
    }
    double accuracy = accuracy_score(tp, tn, positiveSamples, ranking.size() - positiveSamples);
    double precision = precision_score(tp, fp);
    double recall = recall_score(tp, fn);
    double F1 = F1_score(precision, recall);

    resultsRow result = std::make_tuple(
        "",
        tp,
        tn,
        fp,
        fn,
        accuracy,
        precision,
        recall,
        F1,
        -1,
        -1,
        ranking.size(),
        -1.0,
        -1.0,
        -1,
        0.0,
        0.0
    );
    return result;
}
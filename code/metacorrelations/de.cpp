#include <iostream>
#include "stdafx.h"
#include "de.h"
#include <fstream>
#include "headers.h"
#include "Loader.h"
#include "Timer.h"

/**
 *
 * @param paths A vector of paths. Contains the nodes file path, the training set file path, the validation set file path, the training+validation file path, the test path
 * @return
 */

DE::DE(const std::vector<std::vector<double>>& correlations, int corr, int genmax, double CR, double F, int strategy, std::string evolution_metric){
    PathManager &pathmng = PathManager::getInstance();
    correlation = corr;
    this->strategy = strategy;
    std::cout<< "Correlation set to "<<correlation<<std::endl;
    D = static_cast<int>(correlations[0].size());
    NP = static_cast<int>(correlations.size());
    for (size_t i=0;i<correlations.size();i++){
        memcpy(initialValues[i],correlations[i].data(),sizeof (double)*static_cast<size_t>(D));
    }
    this->evolution_metric = evolution_metric;
    this->genmax=genmax;
    this->CR = CR;
    this->F = F;
    this->corr = corr;
    std::cout<<"Created instance"<<std::endl;
    std::cout<<"Population members: "<<correlations.size()<<std::endl;
    std::cout<<"Population members size: "<<correlations[0].size()<<std::endl;
    std::cout<<"Generations: "<<genmax<<std::endl;
    std::cout<<"F: "<<F<<std::endl;
    std::cout<<"CR: "<<CR<<std::endl;
    std::cout<<"Population members size:"<<correlations[0].size()<<std::endl;
}

void DE::initialize_ranking(const PUNGraph training_set,const PUNGraph test_set){
    PathManager &pathmng = PathManager::getInstance();
    ranking.clear();
    TIntV nbrs, nbrs1, nbrs2;
    DEEdge tempEdge;
    int node1, node2, cmn_nbrs, neighbour;
    int nodes = training_set->GetNodes();
    double a1,b1,c1;
    int a,b,c,d;
    //Parameter to reset the same order before each sort 
    int counter = 0;
    int positivesLoaded=0, negativesLoaded=0;
    //Create the ranking
    for(TUNGraph::TNodeI node1it = training_set->BegNI(); node1it != training_set->EndNI(); node1it++){
        node1 = node1it.GetId();
        int deg_1 = node1it.GetDeg();
        TUNGraph::TNodeI node2it(node1it);
        for(node2it++; node2it < training_set->EndNI(); node2it++){
            int deg_2 = node2it.GetDeg();
            node2 = node2it.GetId();
            if(!training_set->IsEdge(node1,node2)){
                cmn_nbrs = TSnap::GetCmnNbrs(training_set, node1, node2, nbrs);
                //Fill node1nbrs with all neighbours of Node1
                for (int i = 0; i < deg_1; i++){
                    neighbour=node1it.GetNbrNId(i);
                    if(neighbour != node2){
                        nbrs1.Add(neighbour);
                    }
                }
                //Fill node2nbrs with all neighbours of Node2
                for (int i = 0; i < deg_2; i++){
                    neighbour=node2it.GetNbrNId(i);
                    if(neighbour!=node1){
                        nbrs2.Add(node2it.GetNbrNId(i));
                    }
                }
                //sort vectors (all vectors must be sorted for diff operation
                nbrs.Sort();
                nbrs1.Sort();
                nbrs2.Sort();
                nbrs1.Diff(nbrs);
                nbrs2.Diff(nbrs);
                a=cmn_nbrs;
                b=deg_1-cmn_nbrs;
                c=deg_2-cmn_nbrs;
                d=nodes-(a+b+c)-2;
                a1 = Measures::adamic_adar_index(training_set, cmn_nbrs, node1, node2, nbrs);
                b1 = Measures::adamic_adar_index(training_set, cmn_nbrs, node1, node2, nbrs1);
                c1 = Measures::adamic_adar_index(training_set, cmn_nbrs, node1, node2, nbrs2);
                std::get<0>(tempEdge)=node1;
                std::get<1>(tempEdge)=node2;
                std::get<2>(tempEdge)=a; //a
                std::get<3>(tempEdge)=b; //b
                std::get<4>(tempEdge)=c; //c
                std::get<5>(tempEdge)=d; //d
                std::get<6>(tempEdge)=a1;
                std::get<7>(tempEdge)=b1;
                std::get<8>(tempEdge)=c1;
                std::get<9>(tempEdge)=0.0; //rank
                //To avoid checking later whether the edge is in the validation set or not
                if(test_set->IsEdge(node1,node2)){
                    std::get<10>(tempEdge) = true;
                    positivesLoaded++;
                }
                else{
                    std::get<10>(tempEdge) = false;
                    negativesLoaded++;
                }
                std::get<11>(tempEdge) = counter; 
                this->ranking.push_back(tempEdge);
                nbrs1.Clr();
                nbrs2.Clr();
                counter++;
            }
        }
    }
    rankingAux.reserve(ranking.size());
    memcpy(rankingAux.data(),ranking.data(),sizeof (DEEdge)*ranking.size());
}

void DE::optimize(){
    PathManager &pathmng = PathManager::getInstance();
    auto config=pathmng.getConfig();
    std::string bestCorrelationFilePath=pathmng.getBestCorrelationFilePath("");
    std::ofstream ofile;
    if(!Filesystem::exists(bestCorrelationFilePath)){
        ofile.open (bestCorrelationFilePath, std::ofstream::app);
        ofile<<"dataset,validationFold,testFold,metacorrelation,DE strategy,generations,F,CR,population,emetric,measure,tp,tn,fp,fn,accuracy,precision,recall,F1,training set edges,test set edges,Epot,auc,approxAuc,approxAucseed,aupr,avgPrec,generation,maxFitness,coefficients\n";
        ofile.close();
    }
    training_set = GraphLoader::loadGraphNodesEdgesList<PUNGraph>(TUNGraph::New(), pathmng.getNodesPath(), pathmng.getTrainingSetPath());
    validation_set = GraphLoader::loadGraphNodesEdgesList<PUNGraph>(TUNGraph::New(), pathmng.getNodesPath(), pathmng.getValidationSetPath());
    initialize_ranking(training_set, validation_set);
    std::mt19937& DEGen = Random::getInstance()->DEGen;
    rnd_uni_init = -static_cast<long>(seed);  // initialization of rnd_uni()
    nfeval       =  0;           // reset number of function evaluations
    double r;
    resultsRow fitness;
//    // Initialization
//    // Right now this part is kept fairly simple and just generates
//    // random numbers in the range [-initfac, +initfac]. You might
//    // want to extend the init part such that you can initialize
//    // each parameter separately.

    // spread initial population members
    /*for (i=0; i<NP; i++) {
        for (j=0; j<D; j++) {
            r = rnd_uni(&rnd_uni_init);
            c[i][j] = inibound_l + r*(inibound_h - inibound_l);
        }
        energy[i] = evaluate(D, c[i], &nfeval,dbInstance);
        // printf("%2d %20.8f %3d\n", i, energy[i], nfeval);
        // cin.get(ch);
    }*/

    for (i=0; i<NP; i++) {
        for (j=0; j<D; j++) {
            c[i][j] = initialValues[i][j];

        }
        fitness = fitness_f(c[i],{evolution_metric},training_set,validation_set);
        if(evolution_metric == "precision"){
            energy[i] = std::get<6>(fitness);
        }
        else if (evolution_metric == "auc"){
            energy[i] = std::get<12>(fitness);
        }
        else if(evolution_metric == "aupr"){
            energy[i] = std::get<15>(fitness);
        }
        else if (evolution_metric == "avgprec"){
            energy[i] = std::get<16>(fitness);
        }
        //energy[i] = evaluate(D, c[i], &nfeval);
        // printf("%2d %20.8f %3d\n", i, energy[i], nfeval);
        // cin.get(ch);
    }

    emax = energy[0];
    imax = 0;
    for (i=1; i<NP; i++) {
        if (energy[i] > emax) {
            emax = energy[i];
            imax = i;
        }
    }

    CopyVector(best, c[imax]);
    CopyVector(bestit, c[imax]);

    // old population (generation G)
    // new population (generation G+1)
    CopyArray(oldarray, c);
    // new population (generation G+1)
    CopyArray(newarray, d);

    // Iteration loop
    gen = 0; // generation counter reset
    while ((gen < genmax)) {
        gen++;
        //std::cout<<"Generation: "<<gen<<std::endl;
        imax = 0;

        for (i=0; i<NP; i++) {
            // Pick a random population member
            do {
                // Endless loop for NP < 2 !!!
                //r = rnd_uni(&rnd_uni_init);
                r = ((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min());
                r1 =static_cast<int>(r*NP);
            } while(r1 == i);

            do {
                // Endless loop for NP < 3 !!!
                r = ((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min());
                r2 = static_cast<int>(r*NP);
            } while((r2 == i) || (r2 == r1));

            do {
                // Endless loop for NP < 4 !!!
                r = ((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min());
                r3 = static_cast<int>(r*NP);
            } while((r3 == i) || (r3 == r1) || (r3 == r2));

            do {
                // Endless loop for NP < 5 !!!
                r = ((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min());
                r4 = static_cast<int>(r*NP);
            } while((r4 == i) || (r4 == r1) || (r4 == r2) || (r4 == r3));

            do {
                // Endless loop for NP < 6 !!!
                r = ((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min());
                r5 = static_cast<int>(r*NP);
            } while((r5==i) || (r5==r1) || (r5==r2) || (r5==r3) || (r5==r4));

            // Choice of strategy
            // We have tried to come up with a sensible naming-convention: DE/x/y/z
            //   DE :  stands for Differential Evolution
            //   x  :  a string which denotes the vector to be perturbed
            //   y  :  number of difference vectors taken for perturbation of x
            //   z  :  crossover method (exp = exponential, bin = binomial)
            //
            // There are some simple rules which are worth following:
            //   1)  F is usually between 0.5 and 1 (in rare cases > 1)
            //   2)  CR is between 0 and 1 with 0., 0.3, 0.7 and 1. being worth to be tried first
            //   3)  To start off NP = 10*D is a reasonable choice. Increase NP if misconvergence=
            //       happens.
            //   4)  If you increase NP, F usually has to be decreased
            //   5)  When the DE/best... schemes fail DE/rand... usually works and vice versa


            // EXPONENTIAL CROSSOVER
            // DE/best/1/exp
            // Our oldest strategy but still not bad. However, we have found several
            // optimization problems where misconvergence occurs.

            // strategy DE0 (not in our paper)
            if (strategy == 1) {
                for (int k=0; k<MAXDIM; k++) {
                    tmp[k] = oldarray[i][k];
                }
                n = static_cast<int>((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min()))*D);
                L = 0;
                do {
                    tmp[n] = bestit[n] + F*(oldarray[r2][n] - oldarray[r3][n]);
                    n = (n+1)%D;
                    L++;
                } while(((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min())) < CR) && (L < D));
            }
            // DE/rand/1/exp
            // This is one of my favourite strategies. It works especially well when the
            // "bestit[]"-schemes experience misconvergence. Try e.g. F=0.7 and CR = 0.5
            // as a first guess.
            // strategy DE1 in the techreport
            else if (strategy == 2) {
                for (int k=0; k<MAXDIM; k++) {
                    tmp[k] = oldarray[i][k];
                }
                n = static_cast<int>((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min()))*D);
                L = 0;
                do {
                    tmp[n] = oldarray[r1][n] + F*(oldarray[r2][n] - oldarray[r3][n]);
                    n = (n+1)%D;
                    L++;
                } while((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min())) < CR && (L < D));
            }
            // DE/rand-to-best/1/exp
            // This strategy seems to be one of the best strategies. Try F=0.85 and CR = 1.0
            // If you get misconvergence try to increase NP. If this doesn't help you
            // should play around with all three control variables.
            // similiar to DE2 but generally better
            else if (strategy == 3) {
                for (int k=0; k<MAXDIM; k++) {
                    tmp[k] = oldarray[i][k];
                }
                n = static_cast<int>((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min()))*D);
                L = 0;
                do {
                    tmp[n] = tmp[n] + F*(bestit[n] - tmp[n]) + F*(oldarray[r1][n] - oldarray[r2][n]);
                    n = (n+1)%D;
                    L++;
                } while(((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min())) < CR) && (L < D));
            }
            // DE/best/2/exp is another powerful strategy worth trying
            else if (strategy == 4) {
                for (int k=0; k<MAXDIM; k++) {
                    tmp[k] = oldarray[i][k];
                }
                n = static_cast<int>((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min()))*D);
                L = 0;
                do {
                    tmp[n] = bestit[n] + (oldarray[r1][n] + oldarray[r2][n] - oldarray[r3][n] - oldarray[r4][n])*F;
                    n = (n+1)%D;
                    L++;
                } while(((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min())) < CR) && (L < D));
            }
            // DE/rand/2/exp seems to be a robust optimizer for many functions
            else if (strategy == 5) {
                for (int k=0; k<MAXDIM; k++) {
                    tmp[k] = oldarray[i][k];
                }
                n = static_cast<int>((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min()))*D);
                L = 0;
                do {
                    tmp[n] = oldarray[r5][n] + (oldarray[r1][n] + oldarray[r2][n] - oldarray[r3][n] - oldarray[r4][n])*F;
                    n = (n+1)%D;
                    L++;
                } while(((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min())) < CR) && (L < D));
            }
            // Essentially same strategies but BINOMIAL CROSSOVER
            // DE/best/1/bin
            else if (strategy == 6) {
                for (int k=0; k<MAXDIM; k++) {
                    tmp[k] = oldarray[i][k];
                }
                n = static_cast<int>((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min()))*D);
                // perform D binomial trials
                for (L=0; L<D; L++) {
                    // change at least one parameter
                    if (((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min())) < CR) || L == (D-1)) {
                        tmp[n] = bestit[n] + F*(oldarray[r2][n] - oldarray[r3][n]);
                    }
                    n = (n+1)%D;
                }
            }
            // DE/rand/1/bin
            else if (strategy == 7) {
                for (int k=0; k<MAXDIM; k++) {
                    tmp[k] = oldarray[i][k];
                }
                n = static_cast<int>((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min()))*D);
                // perform D binomial trials */
                for (L=0; L<D; L++) {
                    // change at least one parameter
                    if (((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min())) < CR) || L == (D-1)) {
                        tmp[n] = oldarray[r1][n] + F*(oldarray[r2][n] - oldarray[r3][n]);
                    }
                    n = (n+1)%D;
                }
            }
            // DE/rand-to-best/1/bin
            else if (strategy == 8) {
                for (int k=0; k<MAXDIM; k++) {
                    tmp[k] = oldarray[i][k];
                }
                n = static_cast<int>((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min()))*D);
                for (L=0; L<D; L++) {
                    if (((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min())) < CR) || L == (D-1)) {
                        tmp[n] = tmp[n] + F*(bestit[n] - tmp[n]) + F*(oldarray[r1][n] - oldarray[r2][n]);
                    }
                    n = (n+1)%D;
                }
            }
            // DE/best/2/bin
            else if (strategy == 9) {
                for (int k=0; k<MAXDIM; k++) {
                    tmp[k] = oldarray[i][k];
                }
                n = static_cast<int>((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min()))*D);
                for (L=0; L<D; L++) {
                    if (((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min())) < CR) || L == (D-1)) {
                        tmp[n] = bestit[n] + (oldarray[r1][n] + oldarray[r2][n] - oldarray[r3][n] - oldarray[r4][n])*F;
                    }
                    n = (n+1)%D;
                }
            }
            // DE/rand/2/bin
            else {
                for (int k=0; k<MAXDIM; k++) {
                    tmp[k] = oldarray[i][k];
                }
                n = static_cast<int>((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min()))*D);
                for (L=0; L<D; L++) {
                    if (((((double)DEGen() -  DEGen.min()) /  (DEGen.max()-  DEGen.min())) < CR) || L == (D-1)) {
                        tmp[n] = oldarray[r5][n] + (oldarray[r1][n] + oldarray[r2][n] - oldarray[r3][n] - oldarray[r4][n])*F;
                    }
                    n = (n+1)%D;
                }
            }

            // Trial mutation now in tmp[]. Test how good this choice really was.
            //std::cout<<"Mutation of individual: "<<i<<": ";
            fitness = fitness_f(tmp,{evolution_metric}, training_set,validation_set);
            if(evolution_metric == "precision"){
                trial_energy = std::get<6>(fitness);
            }
            else if (evolution_metric == "auc"){
                trial_energy = std::get<12>(fitness);
            }
            else if(evolution_metric == "aupr"){
            energy[i] = std::get<15>(fitness);
            }
            else if (evolution_metric == "avgprec"){
                energy[i] = std::get<16>(fitness);
            }
            // improved objective function value?
            //std::cout<<"trial_energy: "<<trial_energy<<std::endl;
            if (trial_energy >= energy[i]) {
                //std::cout<<"Replacing individual "<<i<<" with its mutant"<<std::endl;
                energy[i] = trial_energy;
                for (int k=0; k<MAXDIM; k++) {
                    newarray[i][k] = tmp[k];
                }
                // Was this a new max?
                if (trial_energy>emax) {
                    // reset emax to new low...
                    emax = trial_energy;
                    imax = i;
                    for (int k=0; k<MAXDIM; k++) {
                        best[k] = tmp[k];
                    }
                }
            } else {
                // replace target with old value
                for (int k=0; k<MAXDIM; k++) {
                    newarray[i][k] = oldarray[i][k];
                }
            }
        }
        resultsRow results = fitness_f(best,{"auc","aupr","avgPrec"},training_set,validation_set);
        ofile.open (bestCorrelationFilePath, std::ios::app);
        ofile<<std::get<0>(config)<<",";
        ofile<<std::get<1>(config)<<",";
        ofile<<std::get<2>(config)<<",";
        ofile<<std::get<3>(config)<<",";
        ofile<<std::get<4>(config)<<",";
        ofile<<std::get<5>(config)<<",";
        ofile<<std::get<6>(config)<<",";
        ofile<<std::get<7>(config)<<",";
        ofile<<std::get<8>(config)<<",";
        ofile<<std::get<9>(config)<<",";
        ofile<<std::setprecision(20)<<std::get<0>(results)<<",";
        ofile<<std::get<1>(results)<<",";
        ofile<<std::get<2>(results)<<",";
        ofile<<std::get<3>(results)<<",";
        ofile<<std::get<4>(results)<<",";
        ofile<<std::get<5>(results)<<",";
        ofile<<std::get<6>(results)<<",";
        ofile<<std::get<7>(results)<<",";
        ofile<<std::get<8>(results)<<",";
        ofile<<std::get<9>(results)<<",";
        ofile<<std::get<10>(results)<<",";
        ofile<<std::get<11>(results)<<",";
        ofile<<std::get<12>(results)<<",";
        ofile<<std::get<13>(results)<<",";
        ofile<<std::get<14>(results)<<",";
        ofile<<std::get<15>(results)<<",";
        ofile<<std::get<16>(results)<<",";
        ofile<<gen<<",";
        ofile<<std::setprecision(20)<<emax<<",";
        for(int t=0;t<D;t++){
            ofile<<std::setprecision(20)<<best[t]<<";";
        }
        ofile<<std::endl;
        ofile.close();

        /*std::cout<<"Best: ";
        for(int t=0;t<D;t++){
            std::cout<<best[t]<<" ";
        }
        std::cout<<std::endl<<"With best fitness: "<<emax<<std::endl;*/
        CopyVector(bestit, best);  // Save best population member of current iteration

        // swap population arrays. New generation becomes old one
        CopyArray(swaparray, oldarray);
        CopyArray(oldarray, newarray);
        CopyArray(newarray, swaparray);

        // display after every refresh generations

        // if (gen%refresh == 1) {
        if (false) {
            printf("\n\n Best-so-far obj. funct. value: %15.10f",emax);
            for (j=0; j<D; j++) {
                printf("\n best[%d]: %14.7f", j, best[j]);
            }
            printf("\n Generation: %d  NFEs: %ld", gen, nfeval);
            printf("\n Strategy: %d  NP: %d  F: %f  CR: %f", strategy, NP, F, CR);
        }
    }

    printf( "\n\n Best-so-far objective funct. value: %-15.10g", emax);
    for (j=0; j<D; j++) {
        printf("\n best[%d]: %18.10f", j, best[j]);
    }
    printf("\n No. of Generations: %d", gen);
    printf("\n Number of Function Evaluations: %ld", nfeval);
    printf("\n Strategy: %d, NP: %d, F: %f, CR: %f", strategy, NP, F, CR);
}

resultsRow DE::fitness_f(double params[],std::vector<std::string> metrics, const PUNGraph training_set, PUNGraph test_set){
    memcpy(ranking.data(),rankingAux.data(),sizeof (DEEdge)*ranking.size());
    //std::copy(rankingAux.begin(),rankingAux.end(),ranking);
    switch (corr) {
        /**/
    case 1:
        for(size_t t=0; t < ranking.size(); ++t){
            //Consider only pairs of edges with common neighbours
            if(std::get<2>(ranking[t])>0){
                std::get<9>(ranking[t])=
                    (
                        std::get<2>(ranking[t])*params[0]+ //a
                        std::get<3>(ranking[t])*params[1]+ //b
                        std::get<4>(ranking[t])*params[1]+ //c
                        std::get<5>(ranking[t])*params[2] //d
                    )/(
                        std::get<2>(ranking[t])*params[3]+ //a
                        std::get<3>(ranking[t])*params[4]+ //b
                        std::get<4>(ranking[t])*params[4]+ //c
                        std::get<5>(ranking[t])*params[5]+ //d
                        1*params[6]
                    );
                if(std::isnan(std::get<9>(ranking[t]))){

                    std::cout<<"NAN"<<std::endl;
                }
            }
        }
        break;
        
    case 2:
        for(size_t t=0; t<ranking.size(); ++t){
            if(std::get<2>(ranking[t])>0){
                std::get<9>(ranking[t])=
                    (
                        std::get<2>(ranking[t])*params[0]+ //a
                        std::get<2>(ranking[t])*std::get<2>(ranking[t])*params[1]+ //a^2
                        std::get<2>(ranking[t])*std::get<3>(ranking[t])*params[2]+ //ab
                        std::get<2>(ranking[t])*std::get<4>(ranking[t])*params[2]+ //ac
                        std::get<3>(ranking[t])*std::get<4>(ranking[t])*params[3] //bc
                    )/(
                        std::get<2>(ranking[t])*params[4]+ //a
                        std::get<2>(ranking[t])*std::get<2>(ranking[t])*params[5]+ //a^2
                        std::get<2>(ranking[t])*std::get<3>(ranking[t])*params[6]+ //ab
                        std::get<2>(ranking[t])*std::get<4>(ranking[t])*params[6]+ //ac
                        std::get<3>(ranking[t])*std::get<4>(ranking[t])*params[7]+ //bc
                        1*params[8]
                    );
                if(std::isnan(std::get<9>(ranking[t]))){

                    std::cout<<"NAN"<<std::endl;
                }
            }
        }
        break;
    case 3:
        for(size_t t=0; t < ranking.size(); ++t){
            //Consider only pairs of edges with common neighbours
            if(std::get<2>(ranking[t])>0){
                std::get<9>(ranking[t])=
                    (
                        std::get<2>(ranking[t])*params[0]+ //a
                        std::get<3>(ranking[t])*params[1]+ //b
                        std::get<4>(ranking[t])*params[1]+ //c
                        std::get<5>(ranking[t])*params[2]+ //d
                        std::get<6>(ranking[t])*params[3]+ //a1
                        std::get<7>(ranking[t])*params[4]+ //b1
                        std::get<8>(ranking[t])*params[4] //c1
                    )/(
                        std::get<2>(ranking[t])*params[5]+ //a
                        std::get<3>(ranking[t])*params[6]+ //b
                        std::get<4>(ranking[t])*params[6]+ //c
                        std::get<5>(ranking[t])*params[7]+ //d
                        std::get<6>(ranking[t])*params[8]+ //a1
                        std::get<7>(ranking[t])*params[9]+ //b1
                        std::get<8>(ranking[t])*params[9]+ //c1
                        1*params[10]
                    );
                if(std::isnan(std::get<9>(ranking[t]))){

                    std::cout<<"NAN"<<std::endl;
                }
            }
        }
        break;
    case 4:
        for(size_t t=0; t<ranking.size(); ++t){
            if(std::get<2>(ranking[t])>0){
                std::get<9>(ranking[t])=
                    (
                        std::get<2>(ranking[t])*params[0]+ //a
                        std::get<2>(ranking[t])*std::get<2>(ranking[t])*params[1]+ //a^2
                        std::get<2>(ranking[t])*std::get<3>(ranking[t])*params[2]+ //ab
                        std::get<2>(ranking[t])*std::get<4>(ranking[t])*params[2]+ //ac
                        std::get<3>(ranking[t])*std::get<4>(ranking[t])*params[3]+ //bc
                        std::get<6>(ranking[t])*params[4]+ //a1
                        std::get<7>(ranking[t])*params[5]+ //b1
                        std::get<8>(ranking[t])*params[5] //c1
                    )/(
                        std::get<2>(ranking[t])*params[6]+ //a
                        std::get<2>(ranking[t])*std::get<2>(ranking[t])*params[7]+ //a^2
                        std::get<2>(ranking[t])*std::get<3>(ranking[t])*params[8]+ //ab
                        std::get<2>(ranking[t])*std::get<4>(ranking[t])*params[8]+ //ac
                        std::get<3>(ranking[t])*std::get<4>(ranking[t])*params[9]+ //bc
                        std::get<6>(ranking[t])*params[10]+ //a1
                        std::get<7>(ranking[t])*params[11]+ //b1
                        std::get<8>(ranking[t])*params[11]+ //c1
                        1*params[12]
                    );
                if(std::isnan(std::get<9>(ranking[t]))){

                    std::cout<<"NAN"<<std::endl;
                }
            }
        }
        break;
    }
    /*//first sort according to ids to reset the ranking to its initial state
    #ifdef _MSC_VER
            concurrency::parallel_sort(ranking.begin(), ranking.end(), [](DEEdge a, DEEdge b) {
                return std::get<11>(a) > std::get<11>(b); });
    #else
            __gnu_parallel::sort(ranking.begin(), ranking.end(), [](DEEdge a, DEEdge b) {
                return std::get<11>(a) > std::get<11>(b);});
    #endif*/
    //then sort according to the calculated similarity values
    #ifdef _MSC_VER
            concurrency::parallel_sort(ranking.begin(), ranking.end(), [](DEEdge a, DEEdge b) {
                return std::get<9>(a) > std::get<9>(b);});
    #else
            __gnu_parallel::stable_sort(ranking.begin(), ranking.end(), [](DEEdge a, DEEdge b) {
                return std::get<9>(a) > std::get<9>(b);});
    #endif
    resultsRow results = Classifier::classify(training_set, test_set, ranking, 0, -1, metrics);
    return results;
}

int DE::CopyVector(double a[], double b[]) {
    for (int k=0; k<MAXDIM; k++) {
        a[k] = b[k];
    }
    return 0;
}

int DE::CopyArray(double dest[MAXPOP][MAXDIM], double src[MAXPOP][MAXDIM]) {
    for (int j=0; j<MAXPOP; j++) {
        for (int k=0; k<MAXDIM; k++) {
            dest[j][k] = src[j][k];
        }
    }
    return 0;
}
//Definition for Destructor

//       D I F F E R E N T I A L     E V O L U T I O N
//
// Program: de.c
// Version: 3.6
//
// Authors: Dr. Rainer Storn
//         c/o ICSI, 1947 Center Street, Suite 600
//         Berkeley, CA 94707
//         Tel.:   510-642-4274 (extension 192)
//         Fax.:   510-643-7684
//         E-mail: storn@icsi.berkeley.edu
//         WWW: http://http.icsi.berkeley.edu/~storn/
//         on leave from
//         Siemens AG, ZFE T SN 2, Otto-Hahn Ring 6
//         D-81739 Muenchen, Germany
//         Tel:    636-40502
//         Fax:    636-44577
//         E-mail: rainer.storn@zfe.siemens.de
//
//         Kenneth Price
//         836 Owl Circle
//         Vacaville, CA 95687
//         E-mail: kprice@solano.community.net
//
//  Changes for Visual Sutdio 2010, by Charles Brauer
//
// This program implements some variants of Differential
// Evolution (DE) as described in part in the techreport
// tr-95-012.ps of ICSI. You can get this report either via
// ftp.icsi.berkeley.edu/pub/techreports/1995/tr-95-012.ps.Z
// or via WWW: http://http.icsi.berkeley.edu/~storn/litera.html*
// A more extended version of tr-95-012.ps is submitted for
// publication in the Journal Evolutionary Computation.
//
// ou may use this program for any purpose, give it to any
// person or change it according to your needs as long as you
// are referring to Rainer Storn and Ken Price as the origi-
// nators of the the DE idea.
// If you have questions concerning DE feel free to contact
// us. We also will be happy to know about your experiences
// with DE and your suggestions of improvement.
// SRC-FUNCTION   :main()
// LONG_NAME      :main program)
// FUNCTIONS      :rnd_uni(), evaluate(), printf(), fprintf(),
//                fopen(), fclose(), fscanf().
//
// GLOBALS        :rnd_uni_init    input variable for rnd_uni()
//
// PARAMETERS     :argc            #arguments = 3
//                argv            pointer to argument strings
//
// PRECONDITIONS  :main must be called with three parameters
//                e.g. like de1 <input-file> <output-file>, if
//                the executable file is called de1.
//                The input file must contain valid inputs accor-
//                ding to the fscanf() section of main().
//
// POSTCONDITIONS :main() produces consecutive console outputs and
//                writes the final results in an output file if
//                the program terminates without an error.

/*fscanf(stream, "%d", &strategy);       // choice of strategy
  fscanf(stream, "%d", &genmax);         // maximum number of generations
  fscanf(stream, "%d", &refresh);        // output refresh cycle
  fscanf(stream, "%d", &D);              // number of parameters
  fscanf(stream, "%d", &NP);             // population size.
  fscanf(stream, "%f", &inibound_h);     // upper parameter bound for init
  fscanf(stream, "%f", &inibound_l);     // lower parameter bound for init
  fscanf(stream, "%f", &F);              // weight factor
  fscanf(stream, "%f", &CR);             // crossing over factor
  fscanf(stream, "%d", &seed);           // random seed
  fclose(stream);


  exit(0);

  return(0);*/
void DE::test_evolution(){
    PathManager &pathmng = PathManager::getInstance();
    training_set = GraphLoader::loadGraphNodesEdgesList<PUNGraph>(TUNGraph::New(), pathmng.getNodesPath(), pathmng.getTrainingValidationSetPath());
    test_set = GraphLoader::loadGraphNodesEdgesList<PUNGraph>(TUNGraph::New(), pathmng.getNodesPath(), pathmng.getTestSetPath());
    initialize_ranking(training_set, test_set);
    /*Fitness*/
    /*pathmng.setDumpAucPoints(true);
    pathmng.setDumpAuPRPoints(true);
    pathmng.setDumpAvgPrecPoints(true);*/
    resultsRow results = fitness_f(best,{"auc","approxAuc","aupr","avgPrec"}, training_set,test_set);
    std::cout<<"Best AUC: "<<std::setprecision(20)<<std::get<12>(results)<<std::endl;
    std::ofstream aucSeedFile;
    aucSeedFile.open(pathmng.getAucSeedFilePath(std::to_string(this->correlation) + "_" + std::to_string(this->strategy) + "_"), std::ofstream::app);
    auto aucSeed = Random::getInstance()->getAucSeed();
    auto DESeed = Random::getInstance()->getDESeed();    
    aucSeedFile<<aucSeed<<std::endl<<DESeed;
    std::string resultsFilePath=pathmng.getResultsFilePath();
    std::ofstream ofile;
    if(!Filesystem::exists(resultsFilePath)){
        ofile.open (resultsFilePath, std::ofstream::app);
        ofile<<"dataset,validationFold,testFold,metacorrelation,DE strategy,generations,F,CR,population,emetric,measure,tp,tn,fp,fn,accuracy,precision,recall,F1,training set edges,test set edges,Epot,auc,approxAuc,approxAucseed,aupr,avgPrec,coefficients\n";
        ofile.close();
    }
    auto config=pathmng.getConfig();
    ofile.open (resultsFilePath, std::ios::app);
    ofile<<std::get<0>(config)<<",";
    ofile<<std::get<1>(config)<<",";
    ofile<<std::get<2>(config)<<",";
    ofile<<std::get<3>(config)<<",";
    ofile<<std::get<4>(config)<<",";
    ofile<<std::get<5>(config)<<",";
    ofile<<std::get<6>(config)<<",";
    ofile<<std::get<7>(config)<<",";
    ofile<<std::get<8>(config)<<",";
    ofile<<std::get<9>(config)<<",";
    ofile<<std::setprecision(20)<<std::get<0>(results)<<",";
    ofile<<std::get<1>(results)<<",";
    ofile<<std::get<2>(results)<<",";
    ofile<<std::get<3>(results)<<",";
    ofile<<std::get<4>(results)<<",";
    ofile<<std::get<5>(results)<<",";
    ofile<<std::get<6>(results)<<",";
    ofile<<std::get<7>(results)<<",";
    ofile<<std::get<8>(results)<<",";
    ofile<<std::get<9>(results)<<",";
    ofile<<std::get<10>(results)<<",";
    ofile<<std::get<11>(results)<<",";
    ofile<<std::get<12>(results)<<",";
    ofile<<std::get<13>(results)<<",";
    ofile<<std::get<14>(results)<<",";
    ofile<<std::get<15>(results)<<",";
    ofile<<std::get<16>(results)<<",";
    for(int t=0;t<D;t++){
        ofile<<std::setprecision(20)<<best[t]<<";";
    }
    ofile<<std::endl;
    ofile.close();
    //GraphTools::dump_DEranking_to_file(this->ranking,pathmng.getDERankingFilePath());
}
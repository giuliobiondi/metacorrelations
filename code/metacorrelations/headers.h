#pragma once

#ifdef _MSC_VER
#define NOMINMAX
#endif
#include "Snap.h"
#include <map>
#include <functional>
#include <algorithm>
#include <vector>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <tuple>
#include "json.hpp"
#include <chrono>
#include <cmath> 
#ifdef _MSC_VER
#include <ppl.h>
#else
#include <parallel/algorithm>
#endif
#include <random>
#include "Random.h"
#include <exception>
#include "pathManager.h"
#include "Filesystem.h"
using namespace std;
using json = nlohmann::json;

using edge = std::pair<int,int>;
using graph = std::vector<edge>;
//ID1, ID2, score, IsTestSet
using rankingEdge = std::tuple<uint,uint,double,uint>;
using rankedGraph = std::vector<rankingEdge>;
//Node1, Node2, a, b, c, d, a1 ,b1 ,c1, rank, testset, addOrder
using DEEdge=std::tuple<int,int,int,int,int,int,double,double,double,double,bool,int>;
using DERanking = std::vector<DEEdge>;
using configRow=std::tuple<std::string, int, int, int, int, int, double, double,std::string,std::string>;

using measuresFunction = std::function<double(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs)>;
//measure, tp, tn, fp, fn, accuracy, precision, recall, F1, trainingSetEdges, testSetEdges, ranking size, auc, approxAuc, approxAucSeed,aupr,avg_prec
using resultsRow=std::tuple<std::string, long, long, long, long, double, double, double, double, long, long, long, double, double, long, double, double>;

struct Config{
    std::string datasetName;
    std::string datasetFileName;
    std::string pathToFile;
    std::string nodesPath;
    std::string trainingSetPath;
    std::string validationSetPath;
    std::string testSetPath;
    std::vector<std::string> measures;
    std::vector<std::string> metrics;
    int nfold;
    int k;
	int aucSamples;
    bool dump;
};


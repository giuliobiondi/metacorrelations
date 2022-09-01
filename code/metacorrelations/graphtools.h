#ifndef GRAPHTOOLS_H
#define GRAPHTOOLS_H
#include "headers.h"

using vector_pair_int_int = std::vector<std::pair<int,int>>;
using sper_vector_pair_int_int = std::shared_ptr<std::vector<std::pair<int,int>>>;

namespace GraphTools{
    PUNGraph load_dataset_from_file_to_graph(std::string pathToFile);
    PNGraph load_dataset_from_file_to_directed_graph(std::string pathToFile);
    std::pair<PUNGraph,PUNGraph> load_dataset_from_files(std::string nodesFile, std::string trainingFile, std::string testFile);
    graph load_dataset_from_file(std::string path_to_file);
    #if 0
    std::pair<PUNGraph,PUNGraph> load_dataset_from_db_file(std::string path_to_file);
    #endif
    void dump_ranking_to_file(rankedGraph myGraph, std::string path_to_file);
    void dump_DEranking_to_file(DERanking myGraph, std::string path_to_file);
    void dump_graph_to_file(PUNGraph P, std::string path_to_file);
    void dump_directed_graph_to_file(PNGraph P, std::string path_to_file);
}

#endif // GRAPHTOOLS_H

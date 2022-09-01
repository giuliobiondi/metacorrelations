#include "graphtools.h"

graph GraphTools::load_dataset_from_file(std::string pathToFile){
    std::cout<<"Loading data set from file: "<<pathToFile<<std::endl;
    std::ifstream infile;
    //infile.exceptions(std::ifstream::failbit);
    infile.open(pathToFile);
    std::string line;
    graph myGraph;
    uint node1, node2;
    while (!infile.eof() && std::getline(infile, line))
    {
        std::istringstream iss(line);
        if((iss.str())[0] != '#'&& (iss.str())[0] != '%'){
            iss >> node1 >> node2;
            myGraph.push_back(std::make_pair(node1, node2));
        }
    }
    std::cout<<"Loaded data set from file: "<<pathToFile<<std::endl;
    std::cout<<"Added "<<myGraph.size()<<" edges to table"<<std::endl;
    infile.close();
    return myGraph;
}
PUNGraph GraphTools::load_dataset_from_file_to_graph(std::string pathToFile){
    std::cout<<"Loading data set from file: "<<pathToFile<<std::endl;
    std::ifstream infile;
    //infile.exceptions(std::ifstream::failbit);
    infile.open(pathToFile);
    std::string line;
    PUNGraph myGraph=TUNGraph::New();
    int node1, node2;
    while (!infile.eof() && std::getline(infile, line))
    {
        std::istringstream iss(line);
        if((iss.str())[0] != '#'&& (iss.str())[0] != '%'){
            iss >> node1 >> node2;
            myGraph->AddNodeUnchecked(node1);
            myGraph->AddNodeUnchecked(node2);
            myGraph->AddEdge(node1,node2);
        }
    }
    std::cout<<"Loaded data set from file: "<<pathToFile<<std::endl;
    std::cout<<"Added "<<myGraph->GetNodes()<<" nodes to graph"<<std::endl;
    std::cout<<"Added "<<myGraph->GetEdges()<<" edges to graph"<<std::endl;
    infile.close();
    return myGraph;
}
PNGraph GraphTools::load_dataset_from_file_to_directed_graph(std::string pathToFile){
    std::cout<<"Loading data set from file: "<<pathToFile<<std::endl;
    std::ifstream infile;
    //infile.exceptions(std::ifstream::failbit);
    infile.open(pathToFile);
    std::string line;
    PNGraph myGraph=TNGraph::New();
    int node1, node2;
    while (!infile.eof() && std::getline(infile, line))
    {
        std::istringstream iss(line);
        if((iss.str())[0] != '#'&& (iss.str())[0] != '%'){
            iss >> node1 >> node2;
            myGraph->AddNodeUnchecked(node1);
            myGraph->AddNodeUnchecked(node2);
            myGraph->AddEdge(node1,node2);
        }
    }
    std::cout<<"Loaded data set from file: "<<pathToFile<<std::endl;
    std::cout<<"Added "<<myGraph->GetNodes()<<" nodes to graph"<<std::endl;
    std::cout<<"Added "<<myGraph->GetEdges()<<" edges to graph"<<std::endl;
    infile.close();
    return myGraph;
}
std::pair<PUNGraph,PUNGraph> GraphTools::load_dataset_from_files(std::string nodesFile, std::string trainingFile, std::string testFile){
    std::cout<<"Loading nodes from: "<<nodesFile<<std::endl;
    std::cout<<"Loading training edges from: "<<trainingFile<<std::endl;
    std::cout<<"Loading test edges set from: "<<testFile<<std::endl;
    std::ifstream infile;
    //infile.exceptions(std::ifstream::failbit);
    infile.open(nodesFile);
    std::string line;
    PUNGraph trainingGraph=TUNGraph::New();
    PUNGraph testGraph=TUNGraph::New();
    int node1, node2;
    while (!infile.eof() && std::getline(infile, line))
    {
        std::istringstream iss(line);
        if((iss.str())[0] != '#'&& (iss.str())[0] != '%'){
            iss >> node1;
            trainingGraph->AddNodeUnchecked(node1);
            testGraph->AddNodeUnchecked(node1);
        }
    }
    infile.close();
    infile.open(trainingFile);
    while (!infile.eof() && std::getline(infile, line))
    {
        std::istringstream iss(line);
        if((iss.str())[0] != '#'&& (iss.str())[0] != '%'){
            iss >> node1 >> node2;
            trainingGraph->AddEdge(node1,node2);
        }
    }
    infile.close();
    infile.open(testFile);
    while (!infile.eof() && std::getline(infile, line))
    {
        std::istringstream iss(line);
        if((iss.str())[0] != '#'&& (iss.str())[0] != '%'){
            iss >> node1 >> node2;
            testGraph->AddEdge(node1,node2);
        }
    }
    std::cout<<"Loaded nodes from: "<<nodesFile<<std::endl;
    std::cout<<"Loaded training edges from: "<<trainingFile<<std::endl;
    std::cout<<"Loaded test edges set from: "<<testFile<<std::endl;
    std::cout<<"Added "<<trainingGraph->GetNodes()<<" nodes to training graph"<<std::endl;
    std::cout<<"Added "<<trainingGraph->GetEdges()<<" edges to training graph"<<std::endl;
    std::cout<<"Added "<<testGraph->GetNodes()<<" nodes to test graph"<<std::endl;
    std::cout<<"Added "<<testGraph->GetEdges()<<" edges to test graph"<<std::endl;
    infile.close();
    std::pair<PUNGraph,PUNGraph> ret=std::make_pair(trainingGraph,testGraph);
    return ret;
}

void GraphTools::dump_ranking_to_file(rankedGraph myGraph, std::string path_to_file){
    std::ofstream ofile;
    std::cout<<"Dump ranking: "<<path_to_file<<std::endl;
    ofile.open (path_to_file, std::ofstream::app);
    for(std::vector<rankingEdge>::iterator it = myGraph.begin(); it!=myGraph.end(); ++it){
        ofile<<std::setprecision(20)<<std::get<0>(*it)<<","<<std::get<1>(*it)<<","<<std::get<2>(*it)<<","<<std::get<3>(*it)<<","<<std::endl;
    }
    ofile.close();
}
void GraphTools::dump_DEranking_to_file(DERanking myGraph, std::string path_to_file){
    std::ofstream ofile;
    ofile.open (path_to_file, std::ofstream::app);
    for(std::vector<DEEdge>::iterator it = myGraph.begin(); it!=myGraph.end(); ++it){
        ofile<<std::setprecision(18)<<std::get<0>(*it)<<","<<std::get<1>(*it)<<","<<std::get<2>(*it)<<","<<std::get<3>(*it)<<","<<std::get<4>(*it)<<","<<std::get<5>(*it)<<","<<std::get<6>(*it)<<","<<std::get<7>(*it)<<","<<std::get<8>(*it)<<","<<std::get<9>(*it)<<","<<std::get<10>(*it)<<","<<std::get<11>(*it)<<","<<std::endl;
    }
    ofile.close();
}
void GraphTools::dump_graph_to_file(PUNGraph P, std::string path_to_file){
    std::ofstream ofile;
    ofile.open (path_to_file, std::ofstream::app);
    for (TUNGraph::TEdgeI EI = P->BegEI(); EI < P->EndEI(); EI++) {
        ofile<<EI.GetSrcNId()<<"\t"<<EI.GetDstNId()<<std::endl;
    }
    std::cout<<std::endl<<path_to_file<<std::endl;
    ofile.close();
}
void GraphTools::dump_directed_graph_to_file(PNGraph P, std::string path_to_file){
    std::ofstream ofile;
    ofile.open (path_to_file, std::ofstream::app);
    for (TNGraph::TEdgeI EI = P->BegEI(); EI < P->EndEI(); EI++) {
        ofile<<EI.GetSrcNId()<<"\t"<<EI.GetDstNId()<<std::endl;
    }
    std::cout<<std::endl<<path_to_file<<std::endl;
    ofile.close();
}

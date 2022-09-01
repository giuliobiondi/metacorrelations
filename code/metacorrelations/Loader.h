#pragma once
#include <fstream>
namespace GraphLoader
{
    /*template <class PGraph>  
    PGraph loadGraph(TStr& filePath, TStr graphType = "EdgeList", char& separator=' '){
        if(graphType == "Pajek"){
            return TSnap::LoadPajek(filePath); 
        }
        else if (graphType == "EdgeList"){
            return TSnap::LoadEdgeList(filePath,0,1, separator);
        }
    }*/

    template <class PGraph> 
    PGraph loadGraphNodesEdgesList(PGraph graph,std::string nodesPath, std::string edgesPath, char separator=','){
        int counter = 0;
        int nodeId;
        std::string line;
        std::ifstream infile;
        infile.open(nodesPath);
        std::string tempId;
        int node1, node2;
        while (!infile.eof() && std::getline(infile, line))
        {
            std::stringstream iss(line);
            if((iss.str())[0] != '#'&& (iss.str())[0] != '%'){
                iss >> node1;
                graph->AddNode(node1);
            }
        }
        infile.close();
        std::string temp;
        infile.open(edgesPath);
        int retval;
        while (!infile.eof() && std::getline(infile, line))
        {
            std::stringstream ss(line);
            std::getline(ss, temp, separator);
            node1 = stoi(temp);
            std::getline(ss, temp, separator);
            node2 = stoi(temp);
            if (node1!=node2){
                    retval = graph->AddEdge(node1,node2);
                    if(retval!=-1){
                        std::cout<<"err node1 "<<node1<<" node2 "<<node2<<std::endl;
                    }
                    counter++;
            }
            
        }
        infile.close();
        std::cout<<"Nodes: "<<nodesPath<<" "<<graph->GetNodes()<<" edges:"<<edgesPath<<" "<<graph->GetEdges()<<std::endl;
        return graph;
    }
};

/*template <class PGraph> 
PGraph GraphLoader<PGraph>::loadEdgeList(){

}*/
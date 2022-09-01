#pragma once
#include "Filesystem.h"
#include "headers.h"
using configRow=std::tuple<std::string, int, int, int, int, int, double, double,std::string,std::string>;

class PathManager
{
    private:
        PathManager() {}
        std::string base_path;
        std::string correlations_path;
        std::string nodes_path;
        std::string training_set_path;
        std::string validation_set_path;
        std::string training_validation_set_path;
        std::string test_set_path;
        std::string currentMeasure;
        std::string silent;
        std::string customResultsFilePath;
        std::string dataset;
        std::string population;
        std::string emetric;
        std::string correlations_file;
        int validationFold;
        int testFold;
        int metacorr;
        int strategy;
        int genmax;
        double F;
        double CR;
        
        //Run configuration parameters
        bool dumpAucPoints = false;
        bool dumpAuPRPoints = false;
        bool dumpAvgPrecPoints =false;
        bool dumpRanking = false;
        size_t approxAucSamples = -1;
        std::vector<std::string> measures;
        std::vector<std::string> metrics;
        std::vector<double> mu1;
        std::vector<double> mu2;
        std::vector<double> mu3;
        std::vector<double> mu4;

    public:
        static PathManager& getInstance();
        PathManager(PathManager const&) = delete;
        void operator=(PathManager const&)  = delete;
        void setMeasureName(std::string measure);
        void setConfig(std::string basepath, std::string correlations_folder, std::string dataset, int validation_fold, int test_fold, std::string population, std::string emetric, int metacorr, int strategy, int genmax, double F, double CR);
        void setNodesPath(std::string path);
        void setTrainingSetPath(std::string path);
        void setValidationSetPath(std::string path);
        void setTrainingValidationSetPath(std::string path);
        void setTestSetPath(std::string path);
        std::string getCorrelationsFilePath();
        std::string getMeasureName();
        std::string getNodesPath();
        std::string getTrainingSetPath();
        std::string getValidationSetPath();
        std::string getTrainingValidationSetPath();
        std::string getTestSetPath();
        std::string getBestCorrelationFilePath(std::string file_prefix="");
        std::string getShuffleSeedFilePath(std::string file_prefix="");
        std::string getAucSeedFilePath(std::string file_prefix="");
        std::string getRandomMeasureFilePath(std::string file_prefix="");
        std::string getDERankingFilePath();
        int getValidationFold();
        int getTestFold();
        std::string getDataset();
        std::string getRankingFilePath();
        std::string getResultsFilePath(std::string file_prefix="");
        std::string getMeasuresResultsFilePath(std::string file_prefix="");
        std::string getAucPath();
        std::string getAuPRPath();
        std::string getAvgPrecPath();
        configRow getConfig();
        //Run configuration parameters
        void setApproxAucSamples(size_t samples);
        void setDumpAucPoints(bool dumpAucPoints);
        void setDumpAuPRPoints(bool dumpAuPRPoints);
        void setDumpAvgPrecPoints(bool dumpAvgPrecPoints);
        void setDumpRanking(bool dumpRanking);
        void setMeasures(std::vector<std::string> measures);
        void setMu1(std::vector<double> mu1);
        void setMu2(std::vector<double> mu2);
        void setMu3(std::vector<double> mu3);
        void setMu4(std::vector<double> mu4);
        void setMetrics(std::vector<std::string> metrics);
        void setCustomResultsFilePath(std::string);
        void setValidationFold(int validation_fold);
        void setTestFold(int test_fold);
        void setDataset(std::string dataset);


        size_t getApproxAucSamples();
        bool getDumpAucPoints();
        bool getDumpAuPRPoints();
        bool getDumpAvgPrecPoints();
        bool getDumpRanking();
        std::vector<std::string> getMeasures();
        std::vector<double> getMu1();
        std::vector<double> getMu2();
        std::vector<double> getMu3();
        std::vector<double> getMu4();
        std::vector<std::string> getMetrics();

        void readConfig(std::string configPath);
};
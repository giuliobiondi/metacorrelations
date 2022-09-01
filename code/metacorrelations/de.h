#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "memory.h"
#include "headers.h"
#include <vector>
#include "graphtools.h"
#include "measures.h"
#include "Filesystem.h"
#include "Loader.h"
#include "Classifier.h"

#define MAXPOP 40
#define MAXDIM 20
#define IM1 2147483563
#define IM2 2147483399
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

class DE
{
public:
    DE(const std::vector<std::vector<double>>& correlations,int corr,int genmax, double CR, double F, int strategy, std::string evolution_metric);
    void optimize();
    void test_evolution();
    resultsRow fitness_f(double params[], std::vector<std::string> metrics, const PUNGraph training_set, PUNGraph test_set);
private:
    PUNGraph training_set;
    PUNGraph validation_set;
    PUNGraph test_set;
    std::vector<DEEdge> ranking, rankingAux;
    int correlation;
    double initialValues[MAXPOP][MAXDIM];
    long  rnd_uni_init;
    double c[MAXPOP][MAXDIM], d[MAXPOP][MAXDIM];
    double oldarray[MAXPOP][MAXDIM];
    double newarray[MAXPOP][MAXDIM];
    double swaparray[MAXPOP][MAXDIM];
    int corr;
    int   i, j, L, n;      // counting variables
    int   r1, r2, r3, r4;  // placeholders for random indexes
    int   r5;              // placeholders for random indexes
    int   D;               // Dimension of parameter vector
    int   NP;              // number of population members
    int   imax;            // index to member with highest energy
    int   refresh;         // refresh rate of screen output
    int   strategy = 7;        // choice parameter for screen output
    int   gen, genmax = 100, seed = 1234;
    long  nfeval;          // number of function evaluations
    double trial_energy;    // buffer variable
    double inibound_h=2;      // upper parameter bound
    double inibound_l=0;      // lower parameter bound
    double tmp[MAXDIM], best[MAXDIM], bestit[MAXDIM]; // members
    double energy[MAXPOP];  // obj. funct. values
    double F=0.5, CR=0.9;           // control variables of DE
    double emax;            // help variables
    std::string evolution_metric;
    char ch;
    void initialize_ranking(const PUNGraph training_set, PUNGraph test_set);
    int CopyVector(double a[], double b[]);
    int CopyArray(double dest[MAXPOP][MAXDIM], double src[MAXPOP][MAXDIM]);
};

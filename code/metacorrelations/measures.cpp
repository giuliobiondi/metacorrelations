#include "measures.h"
#include "headers.h"

double Measures::cmn_nbrs(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    return (double)cmn_nbrs;
}

double Measures::jaccard_index(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    if(deg_1 == 0 || deg_2 == 0){
        return 0;
    }
    double jaccard = (double)cmn_nbrs/(deg_1+deg_2-cmn_nbrs);
    return jaccard;
}

double Measures::preferential_attachment(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    return (double)deg_1*deg_2;
}

double Measures::adamic_adar_index(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    double sum = 0, den;
    int nodeDeg;
    for (int i = 0; i < nbrs.Len();i++){
        nodeDeg = T->GetNI(nbrs[i]).GetDeg();
        den = std::log(std::max(2,nodeDeg));
        sum += 1.0/den;
    }
    return sum;
}

double Measures::resource_allocation(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    double sum = 0;
    for (int i = 0; i < nbrs.Len();i++){
        sum += 1.0/T->GetNI(nbrs[i]).GetDeg();
    }
    return sum;
}

double Measures::salton_index(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    if(deg_1==0 || deg_2==0){
        return 0;
    }
    return (double)cmn_nbrs/sqrt((double)deg_1*deg_2);
}

double Measures::sorensen(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    if(deg_1==0 && deg_2==0){
        return 0;
    }
    return 2*cmn_nbrs/((double)deg_1+(double)deg_2);
}

double Measures::hub_promoted(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    int min = std::min(deg_1,deg_2);
    if(min==0){
        return 0;
    }
    return (double)cmn_nbrs/min;
}

double Measures::hub_depressed(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    int max = std::max(deg_1,deg_2);
    if(max==0){
        return 0;
    }
    return (double)cmn_nbrs/max;
}

double Measures::lhn(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    if(deg_1==0 || deg_2==0){
        return 0;
    }
    return (double)cmn_nbrs/(deg_1*deg_2);
}

double Measures::random(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs)
{
    std::mt19937& gen = Random::getInstance()->randomGen;
    double r = static_cast <double> (gen()) / static_cast <double> (gen.max());
    return r;
}

double Measures::avg_confidence(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    size_t deg_1 = T->GetNI(id_1).GetDeg();
    size_t deg_2 = T->GetNI(id_2).GetDeg();
    if(deg_1 == 0 || deg_2 == 0){
        return 0;
    }
    return ((double)cmn_nbrs/deg_1+(double)cmn_nbrs/deg_2)/2;
}

double Measures::max_confidence(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    size_t deg_1 = T->GetNI(id_1).GetDeg();
    size_t deg_2 = T->GetNI(id_2).GetDeg();
    if(deg_1 == 0 || deg_2 == 0){
        return 0;
    }
    return std::max((double)cmn_nbrs/deg_1,(double)cmn_nbrs/deg_2);
}

double Measures::PMI(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs)
{
    size_t l_d1 = T->GetNI(id_1).GetDeg();
    size_t l_d2 = T->GetNI(id_2).GetDeg();
    size_t prod = l_d1 * l_d2;
    double pmi;
    if (l_d1 == 0 || l_d2 == 0 || cmn_nbrs == 0){
        return std::numeric_limits<double>::lowest();
    }
    else{
        pmi = log2((double)T->GetNodes()*cmn_nbrs/(prod));
    }
    return pmi;
}

double Measures::NGD(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs)
{
    size_t deg_1 = T->GetNI(id_1).GetDeg();
    size_t deg_2 = T->GetNI(id_2).GetDeg();
    double ngd;
    if (deg_1 == 0 || deg_2 == 0 || cmn_nbrs == 0){
        ngd = std::numeric_limits<double>::max();
    }
    else{
        ngd = (std::max(log2((double)deg_1),log2((double)deg_2))-log2(cmn_nbrs))/(log2((double)T->GetNodes())-std:: min(log2((double)deg_1),log2((double)deg_2)));

    }
    return ngd;
}
//[0,inf]
double Measures::X2(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs)
{
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (pow(a*d-b*c,2)*n)/((a+b)*(a+c)*(b+d)*(c+d));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Dice (or Czekanowski or Sorensen) similarity measure. [0,1]
double Measures::t_Dice(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    if (a == 0){
        return 0;
    }
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double ret = (2 * a)/(2 * a + b + c); 
    return ret;
}

//3WJaccard [0,1]
double Measures::t_3WJaccard(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    if (a == 0){
        return 0;
    }
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double ret = (3 * a)/(3 * a + b + c);
    return ret;
}

//Sokal and Sneath similarity measure 1. [0,1]
double Measures::t_SokalSneath1(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    if (a == 0){
        return 0;
    }
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b + c + a) - 2;
    double ret = a / (a + 2 * (b + c));
    return ret;
}
//Sokal Michener [0,1]
double Measures::t_SokalMichener(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a + b + c + d;
    return (a + d) / n;
}

//Sokal and Sneath similarity measure 2. [0,1]
double Measures::t_SokalSneath2(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    return (2 * (a + d)) / (2 * (a + d) + b + c);
}

//Rogers and Tanimoto similarity measure. [0,1]
double Measures::t_RogersTanimoto(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    return (a + d) / (a + d + 2 * (b + c));
}

//Faith similarity measure. [0,1]
double Measures::t_Faith(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a + b + c + d;
    return (a + 0.5 * d) / n;
}

//Gower and Legendre similarity measure. [0,1]
double Measures::t_GowerLegendre(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    return (a + d) / (a + 0.5 * (b + c) + d);
}

//Innerproduct similarity measure. [0,n]
double Measures::t_InnerProduct(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    return a + d;
}

//Russell and Rao similarity measure. This measure is the binary dot product. [0,1]
double Measures::t_RussellRao(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a + b + c + d;
    return a / n;
}

//Cosine similarity [0,1]
double Measures::t_Cosine(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b + c + a) - 2;
    double ret = a / (sqrt(( a + b) * (a + c)));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//GilbertWells Similarity [-inf,+inf]
double Measures::t_GilbertWells(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    if (a == 0){
        return std::numeric_limits<double>::lowest();
    }
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    //double ret = log(a) - log(n) - log((a + b) / n) - log((a + c) / n);
    double ret = log((a * n) / ((a + b) * (a + c)));
    unsigned int nodes = T->GetNodes();
    return ret;
}



//Ochiai similarity measure [0,1]
//todeschini 2012
double Measures::t_Ochiai(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double ret = a / sqrt((a + b)*(a + c));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Fossum similarity [0,inf]
double Measures::t_Fossum(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a + b + c + d;
    double ret = (n * pow(a - 0.5,2)) / ((a + b) * (a + c));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Fossum similarity [0,1]
//todeschini 2012
double Measures::t_Fossum_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a + b + c + d;
    double ret = (n * pow(a - 0.5,2)) / ((a + b) * (a + c));
    if (!std::isfinite(ret)){
        return 0;
    }
    ret = ret * n / pow(n-0.5,2);
    return ret;
}

//Sorgenfrei similarity [0,1]
//todeschini 2012
double Measures::t_Sorgenfrei(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a + b + c + d;
    double ret = (pow(a, 2)) / ((a + b) * (a + c));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Mountford similarity [0,1]
//todeschini 2012
double Measures::t_Mountford(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a + b + c + d;
    double ret = a / (0.5 * (a * b + a * c) + b * c);
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}
//Mountford similarity [0,1]
//todeschini 2012
double Measures::t_Mountford_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a + b + c + d;
    double ret = a / (0.5 * (a * b + a * c) + b * c);
    if (!std::isfinite(ret)){
        return a / n;
    }
    return ret / 2;
}

//McConnaughey similarity [-1,1]
double Measures::t_McConnaughey(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes() - (b + c + a) - 2;
    double n = a + b + c + d;
    double ret = (pow(a, 2) - b * c) / ((a + b) * (a + c));
    if (!std::isfinite(ret)){
        return -1;
    }
    return ret;
}

//McConnaughey similarity [-1,1]
//todeschini 2012 [0,1]
double Measures::t_McConnaughey_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    if(a == 0){
        return 0;
    }
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes() - (b + c + a) - 2;
    double n = a + b + c + d;
    double ret = (pow(a, 2) - b * c) / ((a + b) * (a + c));
    return (ret + 1.0) / 2.0;
}

//Tarwid similarity [-1,1]
double Measures::t_Tarwid(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    if(a == 0){
        return -1;
    }
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (n*a -((a + b) * (a + c))) / (n * a +((a + b) * (a + c)));
    return ret;
}


//Kulczynski similarity measure 2.[0,1]
//todeschini 2012
double Measures::t_Kulczynski2(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    if(a == 0){
        return 0;
    }
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double ret = (a / 2.0) * ((2.0 * a + b + c)/((a + b) * (a + c)));
    return ret;
}

//Driver Kroeber similarity measure.[0,1]
double Measures::t_DriverKroeber(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (a / 2.0) * (1 / (a + b) + 1 / (a + c));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Johnson similarity measure.[0,2]
//todeschini 2012
double Measures::t_Johnson(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = a / (a + b) + a / (a + c);
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}
//Johnson similarity measure.[0,2]
//todeschini 2012
double Measures::t_Johnson_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    if (a == 0){
        return 0;
    }
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = a / (a + b) + a / (a + c);
    return ret / 2;
}

//Dennis similarity measure
double Measures::t_Dennis(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (a * d - b * c) / sqrt(n * (a + b) * (a + c));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Dennis similarity measure
//todeschini 2012 [0,1]
double Measures::t_Dennis_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    if (a == n || d==n){
        return 1;
    }
    double ret = (a * d - b * c) / sqrt(n * (a + b) * (a + c));
    if (!std::isfinite(ret)){
        return 0;
    }
    ret = (ret + sqrt(n/2))/sqrt(n); 
    return ret;
}

//Simpson similarity measure.[0,1]
//todeschini 2012
double Measures::t_Simpson(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    if (a == 0){
        return 0;
    }
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = a / std::min(a + b, a + c);
    return ret;
}

//Braun Banquet similarity measure.[0,1]
//todeschini 2012
double Measures::t_BraunBanquet(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    if (a == 0){
        return 0;
    }
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = a / std::max(a + b, a + c);
    return ret;
}

//Fager McGowan similarity measure [-n/2,0.5]
double Measures::t_FagerMcGowan(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (a / sqrt((a + b) * (a + c))) - std::max(a + b, a + c) / 2.0;
    if (!std::isfinite(ret)){
        return (-n) / 2;
    }
    return ret;
}

//Forbes 2 similarity measure 2. [-1,1]
//hubalek
double Measures::t_Forbes2(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (n * a - (a + b) * (a + c))/(n * std::min(a + b, a + c) - ((a + b) * (a + c)));
    if (!std::isfinite(ret)){
        return -1;
    }
    return ret;
}

//Sokal and Sneath similarity measure 4. [0,1]
//todeschini 2012
double Measures::t_SokalSneath4(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (a / (a + b) + a / (a + c) + d / (b + d) + d / (c + d)) / 4.0;
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Sokal and Sneath similarity measure 4. [0,1]
//todeschini 2012
double Measures::t_SokalSneath4_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    if (a == n || d == n){
        return 1;
    }
    if (a == 0 && d == 0){
        return 0;
    }
    double ret = (a / (a + b) + a / (a + c) + d / (b + d) + d / (c + d)) / 4.0;
    return ret;
}


//Gower similarity measure. (Verified)
double Measures::t_Gower(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    if (a == n || d == n){
        return n / (n-1);
    }
    double ret = (a + d)/(sqrt((a + b) * (a + c) * (b + d) * (c + d)));
    return ret;
}

//Pearson-1 similarity measure. [0,+inf] 
double Measures::t_Pearson1(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (n * pow((a * d) - (b * c) , 2)) / ((a + b) * (a + c) * (c + d) * (b + d));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}
//Pearson-1 similarity measure. [0,+inf] 
// var like todeschini
double Measures::t_Pearson1_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    if (a == n || d == n){
        return n;
    }
    double ret = (n * pow((a * d) - (b * c) , 2)) / ((a + b) * (a + c) * (c + d) * (b + d));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Pearson-2 similarity measure. [0,1]
double Measures::t_Pearson2(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double X2 = Measures::t_Pearson1(T, cmn_nbrs, id_1, id_2, nbrs);
    double ret = sqrt(X2 / (n + X2));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Pearson-2 similarity measure. [0,1]
//var
double Measures::t_Pearson2_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double X2 = Measures::t_Pearson1_var(T, cmn_nbrs, id_1, id_2, nbrs);
    double ret = sqrt(X2 / (n + X2));
    return ret;
}

//Pearson-3 similarity measure. (OK)
//var like todeschini
double Measures::t_Pearson3(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    if (a == n || d == n){
        return 1;
    }
    double ro = (a * d - b * c) / (sqrt((a + b) * (a + c) * (c + d) * (b + d)));
    double ret = sqrt((ro) / (n + ro));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Pearson Heron 1 similarity measure. [-1,1]
double Measures::t_PearsonHeron1(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (a * d - b * c) / (sqrt((a + b) * (a + c) * (c + d) * (b + d)));
    if (!std::isfinite(ret)){
        return -1;
    }
    return ret;
}

//Pearson Heron 1 similarity measure. [-1,1]
//todeschini 2012 [0,1]
double Measures::t_PearsonHeron1_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    if (a == n || d == n){
        return 1;
    }
    if (b == n || c == n){
        return 0;
    }
    double ret = (a * d - b * c) / (sqrt((a + b) * (a + c) * (c + d) * (b + d)));
    if (!std::isfinite(ret)){
        return 0;
    }
    ret = (ret + 1.0)/2.0;
    return ret;
}


//Pearson Heron 2 similarity measure. [-1;1]
double Measures::t_PearsonHeron2(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = std::cos((M_PI * sqrt(b * c)) / (sqrt(a * d) + sqrt(b * c)));
    if (!std::isfinite(ret)){
        return -1;
    }
    return ret;
}

//Pearson Heron 2 similarity measure. [-1;1]
//var like todeschini
double Measures::t_PearsonHeron2_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    if (a == n || d == n){
        return 1;
    }
    double ret = std::cos((M_PI * sqrt(b * c)) / (sqrt(a * d) + sqrt(b * c)));
    return ret;
}

//Sokal and Sneath similarity measure 3. Undefined when no nonmatches. [0,+inf]. (Verified)
double Measures::t_SokalSneath3(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double ret = (a + d) / (b + c);
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Sokal and Sneath similarity measure 3. Undefined when no nonmatches. [0,+inf]. (Verified)
double Measures::t_SokalSneath3_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double ret = (a + d) / (b + c);
    if (!std::isfinite(ret)){
        return a+d;
    }
    return ret;
}

//Sokal and Sneath similarity measure 5.[0,1]
//Todeschini2012
double Measures::t_SokalSneath5(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double ret = (a * d) / sqrt((a + b) * (a + c)*(b + d)*(c + d));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Sokal and Sneath similarity measure 5.[0,1]
//Todeschini2012
double Measures::t_SokalSneath5_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    if(a == n || d == n){
        return 1;
    }
    if(a == 0 && d == 0){
        return 0;
    }
    double ret = (a * d) / sqrt((a + b) * (a + c)*(b + d)*(c + d));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

/*//Cole similarity measure 
double Measures::t_Cole(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (sqrt(2) * (a * d - b * c)) / sqrt(pow(a * d - b * c,2) - (a+b) * (a+c) * (b+d) * (c+d));
    if (!std::isfinite(ret)){
        return std::numeric_limits<double>::lowest();
    }
    return ret;
}*/

/*//Stiles similarity measure 
double Measures::t_Stiles(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = std::log10((n * pow((abs(a * d - b * c) - n / 2), 2)) / ((a + b) * (a + c) * (b + d) * (c + d)));
    if (!std::isfinite(ret)){
        return std::numeric_limits<double>::lowest();
    }
    return ret;
}*/

//Ochiai2 similarity measure (Verified)
/*double Measures::t_Ochiai2(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes() - (b + c + a) - 2;
    double n = a+b+c+d;
    double ret = (a * d) / sqrt((a + b) * (a + c) * (b + d) * (c + d));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}*/

//This is the 2×2 version of Goodman and Kruskal’s ordinal measure gamma.
//Like Yule’s Y, Q is a function of the cross-product ratio for a 2×2 table and has a range of –1 to +1. [-1,1]
double Measures::t_YuleQ(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double ret = (a * d - b * c)/(a * d + b * c);
    if (!std::isfinite(ret)){
        return -1;
    }
    return ret;
}

//This is the 2×2 version of Goodman and Kruskal’s ordinal measure gamma.
//Like Yule’s Y, Q is a function of the cross-product ratio for a 2×2 table and has a range of –1 to +1. [-1,1]
//transformed Todeschini2012 [0,1]
double Measures::t_YuleQ_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (a * d - b * c)/(a * d + b * c);
    if (a == n || d == n || b * c == 0){
        return 1;
    }
    return (ret + 1.0) / 2.0;
}

//YuleW Similarity Measure [-1,1]
double Measures::t_YuleW(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes() - (b + c + a) - 2;
    double n = a+b+c+d;
    double ret = (sqrt(a * d) - sqrt(b * c)) / (sqrt(a * d) + sqrt(b * c));
    if (!std::isfinite(ret)){
        return -1;
    }
    return ret;
}

//YuleW Similarity Measure [-1,1]
//Todeschini2012 [0,1]
double Measures::t_YuleW_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes() - (b + c + a) - 2;
    double n = a+b+c+d;
    double ret = (sqrt(a * d) - sqrt(b * c)) / (sqrt(a * d) + sqrt(b * c));
    if (a == n || d == n || b*c == 0){
        return 1;
    }
    return (ret + 1.0) / 2.0;
}

//Kulczynski similarity measure 1. Undefined when no nonmatches. [0,+inf].
double Measures::t_Kulczynski1(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double ret = a/(b+c);
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Dispersion similarity measure. [-1;1]
double Measures::t_Disper(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (a * d - b * c) / pow(a + b + c + d,2);
    return ret;
}

//Dispersion similarity measure. [-1;1]
//todeschini 2012[0,1]
double Measures::t_Disper_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    if (a == n || d == n){
        return 1;
    }
    double ret = (a * d - b * c) / pow(a + b + c + d,2);
    return (ret + 0.25)*2;
}

//Hamann similarity measure. [-1,1]
double Measures::t_Hamann(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    return ((a + d) - (b + c)) / (a + b + c + d);
}

//Hamann similarity measure. [-1,1]
//todeschini 2012[0,1]
double Measures::t_Hamann_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double ret = ((a + d) - (b + c)) / (a + b + c + d);
    return (ret + 1.0) / 2.0;
}

//Michael similarity measure.[-1,1]
double Measures::t_Michael(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (4 * (a * d - b * c)) / (pow(a + d, 2)+ pow(b + c, 2));
    return ret;
}

//Michael similarity measure.[-1,1]
//Todeschini2012 [0,1]
double Measures::t_Michael_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (4 * (a * d - b * c)) / (pow(a + d, 2)+ pow(b + c, 2));
    if (a == n || d == n){
        return 1;
    }
    if (b + c == 0){
        return 1;
    }
    return (ret + 1.0) / 2.0;
}

//Goodman and Kruskal’s lambda (similarity).[0,1]
//https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/55190/versions/1/previews/simbin.m/index.html
double Measures::t_GoodManKruskalLambda(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    int t1 = std::max(a, b) + std::max(c, d) + std::max(a, c) + std::max(b,d);
    int t2 = std::max(a + c, b + d) + std::max(a + d, c + d);
    double ret = (t1 - t2) / (2 * (a + b + c + d) - t2);
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Goodman and Kruskal’s lambda (similarity).[0,1]
//var like Todeschini
double Measures::t_GoodManKruskalLambda_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    if (a == n || d == n){
        return 1;
    }
    int t1 = std::max(a, b) + std::max(c, d) + std::max(a, c) + std::max(b,d);
    int t2 = std::max(a + c, b + d) + std::max(a + d, c + d);
    double ret = (t1 - t2) / (2 * (a + b + c + d) - t2);
    return ret;
}

//Anderberg’s D (similarity)[0,1] (Verified)
double Measures::t_AnderbergD(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    int t1 = std::max(a, b) + std::max(c, d) + std::max(a, c) + std::max(b,d);
    int t2 = std::max(a + c, b + d) + std::max(a + d, c + d);
    return (t1 - t2) / (2 * (a + b + c + d));
}

//BaroniUrbaniBuser1 similarity measure. [0,1]
double Measures::t_BaroniUrbaniBuser1(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a + b + c + d;
    double ret = (sqrt(a * d) + a ) / (sqrt(a * d) + a + b + c);
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//BaroniUrbaniBuser1 similarity measure. [0,1]
//todeschini 2012
double Measures::t_BaroniUrbaniBuser1_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a + b + c + d;
    double ret = (sqrt(a * d) + a ) / (sqrt(a * d) + a + b + c);
    if (d == n){
        return 1;
    }
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//BaroniUrbaniBuser2 similarity measure. [-1,1]
double Measures::t_BaroniUrbaniBuser2(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a + b + c + d;
    double ret = (sqrt(a * d) + a - (b + c)) / (sqrt(a * d) + a + b + c);
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//BaroniUrbaniBuser2 similarity measure. [-1,1]
//Todeschini 2012 [0,1]
double Measures::t_BaroniUrbaniBuser2_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a + b + c + d;
    if (d == n){
        return 1;
    }
    double ret = (sqrt(a * d) + a - (b + c)) / (sqrt(a * d) + a + b + c);
    return (ret + 1.0 ) / 2.0;
}

//Eyraud similarity measure. [-1,0]
double Measures::t_Eyraud(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (a - (a + b) * (a + c)) / ((a + b) * (a + c) * (b + d) * (c + d));
        if (!std::isfinite(ret)){
        return -1;
    }
    return ret;
}

//Rogot_Goldberg [0,1]
double Measures::t_RogotGoldberg(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = a / (2 * a + b +c) + d / (2 * d + b +c);
        if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Rogot_Goldberg [0,1]
//Todeschini 2012 [0,1]
double Measures::t_RogotGoldberg_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    if (a == n || d == n){
        return 1;
    }
    double ret = a / (2 * a + b +c) + d / (2 * d + b +c);
    return ret;
}

//HawkinsDotson[0,1]
double Measures::t_HawkinsDotson(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = 0.5 *((a / (a + b + c)) + (d / (b + c + d)));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//HawkinsDotson[0,1]
//Todeschini 2012 [0,1]
double Measures::t_HawkinsDotson_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    if (a == n || d == n){
        return 1;
    }
    double ret = 0.5 *((a / (a + b + c)) + (d / (b + c + d)));
    return ret;
}

//Cohen [-1,1]
double Measures::t_Cohen(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (2 * (a * d - b * c))/((a + b)*(b + d) + (a + c)*(c + d));
        if (!std::isfinite(ret)){
        return -1;
    }
    return ret;
}

//Cohen [-1,1]
//Todeschini 2012 [0,1]
double Measures::t_Cohen_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    if (a == n || d == n){
        return 1;
    }
    double ret = (2 * (a * d - b * c))/((a + b)*(b + d) + (a + c)*(c + d));
        if (!std::isfinite(ret)){
        return 0;
    }
    return (ret + 1.0)/2.0;
}

//Maxwell-Pilliner [-1,1]
double Measures::t_MaxwellPilliner(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (2 * (a * d - b * c))/((a + b) * (c + d) + (a + c)* (b + d));
    if (!std::isfinite(ret)){
        return -1;
    }
    return ret;
}

//Maxwell-Pilliner
//Todeschini 2012 [0,1]
double Measures::t_MaxwellPilliner_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    if (a == n || d == n){
        return 1;
    }
    double ret = (2 * (a * d - b * c))/((a + b) * (c + d) + (a + c)* (b + d));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Harris Lahey
//Todeschini 2012 [0,n]
double Measures::t_HarrisLahey(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (a * (2 * d + b + c)) / (2 * (a + b + c)) + (d * (2 * a + b + c)) / (2 * (b + c + d));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Harris Lahey
//Todeschini 2012 [0,1]
double Measures::t_HarrisLahey_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    if (a == n || d == n){
        return 1;
    }
    double ret = (a * (2 * d + b + c)) / (2 * (a + b + c)) + (d * (2 * a + b + c)) / (2 * (b + c + d));
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret / n;
}

//
double Measures::t_CT1(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (log(1 + a + d))/(log(1+n));
    return ret;
}

//
double Measures::t_CT2(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (log(1 + n) -log(1 + b +c))/(log(1 + n));
    return ret;
}
//
double Measures::t_CT3(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (log(1 + a)) / (log(1 + n));
    return ret;
}
//
double Measures::t_CT4(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (log(1 + a)) / (log(1 + a + b + c));
    return ret;
}

//
double Measures::t_CT5(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (log(1 + a * d)-log(1 + b * c)) / (log(1 + pow(n , 2) / 4.0));
    return ret;
}

//
double Measures::t_AustinColwell(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (2.0 / M_PI) * asin(sqrt((a + d) / n));
    return ret;
}

//Scott [-1,0]
double Measures::t_Scott(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (4 * a * d - pow(b + c, 2)) / ((2*a + b + c)*(2 * d + b + c));
        if (!std::isfinite(ret)){
        return -1;
    }
    return ret;
}

//Scott [-1,0]
//Todeschini 2012 [0,1]
double Measures::t_Scott_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    if (a == n || d == n){
        return 1;
    }
    double ret = (4 * a * d - pow(b + c, 2)) / ((2*a + b + c)*(2 * d + b + c));
    return ret + 1.0;
}

//
double Measures::t_VanDerMaarel(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (2 * a - b - c) / (2 * a - b - c);
        if (!std::isfinite(ret)){
        return -1;
    }
    return ret;
}

//
double Measures::t_VanDerMaarel_var(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    double ret = (2 * a - b - c) / (2 * a - b - c);
        if (!std::isfinite(ret)){
        return -1;
    }
    return (ret + 1.0) / 2.0;
}

//Peirce similarity measure. (OK) (NON-SYMMETRIC)
double Measures::t_Peirce(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    return (a * b + b * c)/(a * b + 2 * b * c + c * d);
}

//Tarantula similarity measure. (NON-SYMMETRIC)
double Measures::t_Tarantula(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    return (a *(c+d))/(c*(a+b));
}

//Ample similarity measure. (NON-SYMMETRIC)
double Measures::t_Ample(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    return ((a+d)-(b+c))/(a+b+c+d);
}

//Yule’s Y coefficient of colligation (similarity)[-1,1]
/*double Measures::t_YuleY(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double ret = (sqrt(a*d)-sqrt(b*c))/(sqrt(a*d)+sqrt(b*c));
    if (!std::isfinite(ret)){
        return -1;
    }
    return ret;
}*/


//Binary Euclidean distance[0,+inf)
double Measures::t_Beuclid(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    return sqrt(b+c);
}
//Binary squared Euclidean distance.[0,+inf)
double Measures::t_Bseuclid(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    return b+c;
}


//Simple matching similarity measure. This measure is the ratio of the number of matches to the total number of characteristics.
double Measures::t_SimpleMatching(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    return (a+d)/(a+b+c+d);
}
//Jaccard similarity measure. This measure is also known as the similarity ratio.
double Measures::t_Jaccard(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double ret = a/(a+b+c);
    if (!std::isfinite(ret)){
        return 0;
    }
    return ret;
}

//Hamming dissimilarity measure. (Verified)
double Measures::t_Hamming(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    return b + c;
}

//Variance dissimilarity measure. (Verified)
double Measures::t_Variance(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a + b + c + d;
    return (b + c) / (4 * n);
}

//Size dissimilarity measure.[0,+inf) (Verified)
double Measures::t_Size(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    return pow((b + c), 2) / pow(n, 2);
}

//Binary shape difference(-inf,+inf) (Verified)
double Measures::t_Bshape(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a+b+c+d;
    return (n * (b + c) - pow(b - c, 2)) / pow( n ,2);
}

//Pattern difference.[0,1] (Verified)
double Measures::t_Pattern(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double n = a + b + c + d;
    return (4 * b * c) / pow(n, 2);
}

//Binary Lance-and-Williams nonmetric dissimilarity measure [0,1] (Verified)
double Measures::t_Blwmn(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double ret = (b + c) / (2 * a + b + c);
    if (!std::isfinite(ret)){
        return 1.0;
    }
    return ret;
}

//Hellinger dissimilarity [0,2] (Verified)
double Measures::t_Hellinger(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double ret = 2 * sqrt(1 - (a / sqrt((a + b) * (a + c))));
    if (!std::isfinite(ret)){
        return 2.0;
    }
    return ret;
}

//YuleQ Distance (Verified)
double Measures::t_YuleQDistance(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes() - (b + c + a) - 2;
    double n = a+b+c+d;
    double ret = (2 * b * c) / (a * d + b * c);
    if (!std::isfinite(ret)){
        return 2;
    }
    return ret;
}

double Measures::mu4(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    PathManager &pathmng = PathManager::getInstance();
    std::vector<double> coefficients = pathmng.getMu4();
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double a1 = adamic_adar_index(T, cmn_nbrs, id_1, id_2, nbrs);
    TUNGraph::TNodeI node1it=T->GetNI(id_1);
    TUNGraph::TNodeI node2it=T->GetNI(id_2);
    TIntV node1nbrs;
    TIntV node2nbrs;
    int neighbour;
    //Fill node1nbrs with all neighbours of Node1
    for (int i = 0; i < deg_1; i++){
        neighbour=node1it.GetNbrNId(i);
        if(neighbour!=id_2){
            node1nbrs.Add(neighbour);
        }
    }
    //Fill node2nbrs with all neighbours of Node2
    for (int i = 0; i < deg_2; i++){
        neighbour=node2it.GetNbrNId(i);
        if(neighbour!=id_1){
            node2nbrs.Add(neighbour);
        }
    }
    //sort vectors (all vectors must be sorted for diff operation
    nbrs.Sort();
    node1nbrs.Sort();
    node2nbrs.Sort();
    node1nbrs.Diff(nbrs);
    node2nbrs.Diff(nbrs);
    double b1 = adamic_adar_index(T, cmn_nbrs, id_1, id_2, node1nbrs);
    double c1 = adamic_adar_index(T, cmn_nbrs, id_1, id_2, node2nbrs);
    return (a*coefficients[0]+a*a*coefficients[1]+a*b*coefficients[2]+a*c*coefficients[2]+b*c*coefficients[3]+a1*coefficients[4]+b1*coefficients[5]+c1*coefficients[5])/(a*coefficients[6]+a*a*coefficients[7]+a*b*coefficients[8]+a*c*coefficients[8]+b*c*coefficients[9]+a1*coefficients[10]+b1*coefficients[11]+c1*coefficients[11]+coefficients[12]);
}

double Measures::mu3(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    PathManager &pathmng = PathManager::getInstance();
    std::vector<double> coefficients = pathmng.getMu3();
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    double a1 = adamic_adar_index(T, cmn_nbrs, id_1, id_2, nbrs);
    TUNGraph::TNodeI node1it=T->GetNI(id_1);
    TUNGraph::TNodeI node2it=T->GetNI(id_2);
    TIntV node1nbrs;
    TIntV node2nbrs;
    int neighbour;
    //Fill node1nbrs with all neighbours of Node1
    for (int i = 0; i < deg_1; i++){
        neighbour=node1it.GetNbrNId(i);
        if(neighbour!=id_2){
            node1nbrs.Add(neighbour);
        }
    }
    //Fill node2nbrs with all neighbours of Node2
    for (int i = 0; i < deg_2; i++){
        neighbour=node2it.GetNbrNId(i);
        if(neighbour!=id_1){
            node2nbrs.Add(node2it.GetNbrNId(i));
        }
    }
    //sort vectors (all vectors must be sorted for diff operation
    nbrs.Sort();
    node1nbrs.Sort();
    node2nbrs.Sort();
    //Preserve only non-shared neighbours
    node1nbrs.Diff(nbrs);
    node2nbrs.Diff(nbrs);
    //Calculate pseudo-adamic adar indices (note: only the last parameter is used)
    double b1 = adamic_adar_index(T, cmn_nbrs, id_1, id_2, node1nbrs);
    double c1 = adamic_adar_index(T, cmn_nbrs, id_1, id_2, node2nbrs);
    return (a*coefficients[0]+b*coefficients[1]+c*coefficients[1]+d*coefficients[2]+a1*coefficients[3]+b1*coefficients[4]+c1*coefficients[4])/(a*coefficients[5]+b*coefficients[6]+c*coefficients[6]+d*coefficients[7]+a1*coefficients[8]+b1*coefficients[9]+c1*coefficients[9]+coefficients[10]);
}

double Measures::mu1(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    PathManager &pathmng = PathManager::getInstance();
    std::vector<double> coefficients = pathmng.getMu1();
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    if(a == 0){
        return 0;
    }
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    TUNGraph::TNodeI node1it=T->GetNI(id_1);
    TUNGraph::TNodeI node2it=T->GetNI(id_2);
    TIntV node1nbrs;
    TIntV node2nbrs;
    int neighbour;
    //Fill node1nbrs with all neighbours of Node1
    for (int i = 0; i < deg_1; i++){
        neighbour=node1it.GetNbrNId(i);
        if(neighbour!=id_2){
            node1nbrs.Add(neighbour);
        }
    }
    //Fill node2nbrs with all neighbours of Node2
    for (int i = 0; i < deg_2; i++){
        neighbour=node2it.GetNbrNId(i);
        if(neighbour!=id_1){
            node2nbrs.Add(node2it.GetNbrNId(i));
        }
    }
    //sort vectors (all vectors must be sorted for diff operation
    nbrs.Sort();
    node1nbrs.Sort();
    node2nbrs.Sort();
    //Preserve only non-shared neighbours
    node1nbrs.Diff(nbrs);
    node2nbrs.Diff(nbrs);
    return 
    (a*coefficients[0]
    +b*coefficients[1]
    +c*coefficients[1]
    +d*coefficients[2])/
    (a*coefficients[3]
    +b*coefficients[4]
    +c*coefficients[4]
    +d*coefficients[5]
    +coefficients[6]);
}

double Measures::mu2(const PUNGraph& T, int cmn_nbrs, int id_1, int id_2, TIntV nbrs){
    PathManager &pathmng = PathManager::getInstance();
    std::vector<double> coefficients = pathmng.getMu2();
    int deg_1 = T->GetNI(id_1).GetDeg();
    int deg_2 = T->GetNI(id_2).GetDeg();
    double a = cmn_nbrs;
    if(a == 0){
        return 0;
    }
    double b = deg_1-cmn_nbrs;
    double c = deg_2-cmn_nbrs;
    double d = T->GetNodes()-(b+c+a)-2;
    TUNGraph::TNodeI node1it=T->GetNI(id_1);
    TUNGraph::TNodeI node2it=T->GetNI(id_2);
    TIntV node1nbrs;
    TIntV node2nbrs;
    int neighbour;
    //Fill node1nbrs with all neighbours of Node1
    for (int i = 0; i < deg_1; i++){
        neighbour=node1it.GetNbrNId(i);
        if(neighbour!=id_2){
            node1nbrs.Add(neighbour);
        }
    }
    //Fill node2nbrs with all neighbours of Node2
    for (int i = 0; i < deg_2; i++){
        neighbour=node2it.GetNbrNId(i);
        if(neighbour!=id_1){
            node2nbrs.Add(neighbour);
        }
    }
    //sort vectors (all vectors must be sorted for diff operation
    nbrs.Sort();
    node1nbrs.Sort();
    node2nbrs.Sort();
    node1nbrs.Diff(nbrs);
    node2nbrs.Diff(nbrs);
    return 
    (a*coefficients[0]
    +a*a*coefficients[1]
    +a*b*coefficients[2]
    +a*c*coefficients[2]
    +b*c*coefficients[3])/
    (a*coefficients[4]
    +a*a*coefficients[5]
    +a*b*coefficients[6]
    +a*c*coefficients[6]
    +b*c*coefficients[7]
    +coefficients[8]);
}


//Function pointer, boolean isProximity
std::map<std::string,std::pair<measuresFunction,bool>> Measures::measuresFunctions={
    {"Common Neighbours", {&Measures::cmn_nbrs,true}},
    {"Jaccard", {&Measures::jaccard_index,true}},
    {"Preferential Attachment", {&Measures::preferential_attachment,true}},
    {"Adamic Adar", {&Measures::adamic_adar_index,true}},
    {"Resource Allocation", {&Measures::resource_allocation,true}},
    {"Salton", {&Measures::salton_index,true}},
    {"Sorensen", {&Measures::sorensen,true}},
    {"Hub Promoted", {&Measures::hub_promoted,true}},
    {"Hub Depressed", {&Measures::hub_depressed,true}},
    {"LHN", {&Measures::lhn,true}},
    {"Random", {&Measures::random,true}},
    {"Average Confidence", {&Measures::avg_confidence,true}},
    {"Max Confidence", {&Measures::max_confidence,true}},
    {"PMI", {&Measures::PMI,true}},
    {"NGD", {&Measures::NGD,false}},
    {"X2", {&Measures::X2,true}},
    {"T Dice", {&Measures::t_Dice,true}},
    {"T 3WJaccard", {&Measures::t_3WJaccard,true}},
    {"T Sokal Sneath 1", {&Measures::t_SokalSneath1,true}},
    {"T Sokal Michener", {&Measures::t_SokalMichener,true}},
    {"T Sokal Sneath 2", {&Measures::t_SokalSneath2,true}},
    {"T Rogers Tanimoto", {&Measures::t_RogersTanimoto,true}},
    {"T Faith", {&Measures::t_Faith,true}},
    {"T Gower Legendre", {&Measures::t_GowerLegendre,true}},
    {"T InnerProduct", {&Measures::t_InnerProduct,true}},
    {"T Russell Rao", {&Measures::t_RussellRao,true}},
    {"T Cosine", {&Measures::t_Cosine,true}},
    {"T Gilbert Wells", {&Measures::t_GilbertWells,true}},
    {"T Ochiai", {&Measures::t_Ochiai,true}},
    {"T Fossum", {&Measures::t_Fossum,true}},
    {"T Fossum Var", {&Measures::t_Fossum_var,true}},
    {"T Sorgenfrei", {&Measures::t_Sorgenfrei,true}},
    {"T Mountford", {&Measures::t_Mountford,true}},
    {"T Mountford Var", {&Measures::t_Mountford_var,true}},
    {"T McConnaughey", {&Measures::t_McConnaughey,true}},
    {"T McConnaughey Var", {&Measures::t_McConnaughey_var,true}},
    {"T Tarwid", {&Measures::t_Tarwid,true}},
    {"T Kulczynski 2", {&Measures::t_Kulczynski2,true}},
    {"T Driver Kroeber", {&Measures::t_DriverKroeber,true}},
    {"T Johnson", {&Measures::t_Johnson,true}},
    {"T Johnson Var", {&Measures::t_Johnson_var,true}},
    {"T Dennis", {&Measures::t_Dennis,true}},
    {"T Dennis Var", {&Measures::t_Dennis_var,true}},
    {"T Simpson", {&Measures::t_Simpson,true}},
    {"T Braun Banquet", {&Measures::t_BraunBanquet,true}},
    {"T Fager McGowan", {&Measures::t_FagerMcGowan,true}},
    {"T Forbes 2", {&Measures::t_Forbes2,true}},
    {"T Sokal Sneath 4", {&Measures::t_SokalSneath4,true}},
    {"T Sokal Sneath 4 Var", {&Measures::t_SokalSneath4_var,true}},
    {"T Gower", {&Measures::t_Gower,true}},
    {"T Pearson1", {&Measures::t_Pearson1,true}},
    {"T Pearson1 Var", {&Measures::t_Pearson1,true}},
    {"T Pearson2", {&Measures::t_Pearson2,true}},
    {"T Pearson2 Var", {&Measures::t_Pearson1,true}},
    {"T Pearson3", {&Measures::t_Pearson3,true}},
    {"T PearsonHeron1", {&Measures::t_PearsonHeron1,true}},
    {"T PearsonHeron1 Var", {&Measures::t_PearsonHeron1_var,true}},
    {"T PearsonHeron2", {&Measures::t_PearsonHeron2,true}},
    {"T PearsonHeron2 Var", {&Measures::t_PearsonHeron2_var,true}},
    {"T Sokal Sneath 3", {&Measures::t_SokalSneath3,true}},
    {"T Sokal Sneath 3 Var", {&Measures::t_SokalSneath3_var,true}},
    {"T Sokal Sneath 5", {&Measures::t_SokalSneath5,true}},
    {"T Sokal Sneath 5 Var", {&Measures::t_SokalSneath5_var,true}},
    //{"T Cole", {&Measures::t_Cole,true}},
    //{"T Stiles", {&Measures::t_Stiles,true}},
    //{"T Ochiai2", {&Measures::t_Ochiai2,true}},
    {"T YuleQ", {&Measures::t_YuleQ,true}},
    {"T YuleQ Var", {&Measures::t_YuleQ_var,true}},
    {"T YuleW", {&Measures::t_YuleW,true}},
    {"T YuleW Var", {&Measures::t_YuleW_var,true}},
    {"T Disper", {&Measures::t_Disper,true}},
    {"T Disper Var", {&Measures::t_Disper_var,true}},
    {"T Hamann", {&Measures::t_Hamann,true}},
    {"T Hamann Var", {&Measures::t_Hamann_var,true}}, 
    {"T Michael", {&Measures::t_Michael,true}},
    {"T Michael Var", {&Measures::t_Michael_var,true}},
    {"T GoodMan Kruskall Lambda", {&Measures::t_GoodManKruskalLambda,true}},
    {"T GoodMan Kruskall Lambda Var", {&Measures::t_GoodManKruskalLambda_var,true}},
    //{"T Yule Y", {&Measures::t_YuleY,true}},
    {"T Yule Q", {&Measures::t_YuleQ,true}},
    {"T Kulczynski 1", {&Measures::t_Kulczynski1,true}},
    {"T Anderberg", {&Measures::t_AnderbergD,true}},
    {"T Baroni Urbani Buser 1", {&Measures::t_BaroniUrbaniBuser1,true}},
    {"T Baroni Urbani Buser 1 Var", {&Measures::t_BaroniUrbaniBuser1_var,true}},
    {"T Baroni Urbani Buser 2", {&Measures::t_BaroniUrbaniBuser2,true}},
    {"T Baroni Urbani Buser 2 Var", {&Measures::t_BaroniUrbaniBuser2_var,true}},
    {"T Rogot Goldberg", {&Measures::t_RogotGoldberg,true}},
    {"T Rogot Goldberg Var", {&Measures::t_RogotGoldberg_var,true}},
    {"T Hawkins Dotson", {&Measures::t_HawkinsDotson,true}},
    {"T Hawkins Dotson Var", {&Measures::t_HawkinsDotson_var,true}},
    {"T Cohen", {&Measures::t_Cohen,true}},
    {"T Cohen Var", {&Measures::t_Cohen_var,true}},
    {"T Maxwell Pilliner", {&Measures::t_MaxwellPilliner,true}},
    {"T Maxwell Pilliner Var", {&Measures::t_MaxwellPilliner_var,true}},
    {"T Harris Lahey", {&Measures::t_HarrisLahey,true}},
    {"T Harris Lahey Var", {&Measures::t_HarrisLahey_var,true}},
    {"T CT1", {&Measures::t_CT1,true}},
    {"T CT2", {&Measures::t_CT2,true}},
    {"T CT3", {&Measures::t_CT3,true}},
    {"T CT4", {&Measures::t_CT4,true}},
    {"T CT5", {&Measures::t_CT5,true}},
    {"T Austin Colwell", {&Measures::t_AustinColwell,true}},
    {"T Scott", {&Measures::t_Scott,true}},
    {"T Scott Var", {&Measures::t_Scott_var,true}},
    {"T Van Der Maarel", {&Measures::t_VanDerMaarel,true}},
    {"T Van Der Maarel Var", {&Measures::t_VanDerMaarel_var,true}},
    {"T Eyraud", {&Measures::t_Eyraud,true}},
    {"T Peirce", {&Measures::t_Peirce,true}},
    {"T Tarantula", {&Measures::t_Tarantula,true}},
    {"T Ample", {&Measures::t_Ample,true}},
    {"T Hamming", {&Measures::t_Hamming,false}},
    {"T Variance", {&Measures::t_Variance,false}},
    {"T Size", {&Measures::t_Size,false}},
    {"T BShape", {&Measures::t_Bshape,false}},
    {"T Pattern", {&Measures::t_Pattern,false}},
    {"T BLWMN", {&Measures::t_Blwmn,false}},
    {"T Hellinger", {&Measures::t_Hellinger,false}},
    {"T YuleQ Distance", {&Measures::t_YuleQDistance,false}},
    //{"T Phi", {&Measures::t_Phi,true}},
    {"T Beuclid", {&Measures::t_Beuclid,false}},
    {"T Bseuclid", {&Measures::t_Bseuclid,false}},
    {"T Simple Matching", {&Measures::t_SimpleMatching,true}},
    {"T Jaccard", {&Measures::t_Jaccard,true}},
    {"mu2", {&Measures::mu2,true}},
    {"mu1", {&Measures::mu1,true}},
    {"mu3", {&Measures::mu3,true}},
    {"mu4", {&Measures::mu4,true}},
};
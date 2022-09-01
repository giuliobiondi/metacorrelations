#include<iostream>
#include "Random.h"
static Random* uniqueinstance;


unsigned int Random::getShuffleSeed()
{
	return this->shuffleSeed;
}

unsigned int Random::getAucSeed()
{
	return this->aucSeed;
}

unsigned int Random::getDESeed()
{
	return this->DESeed;
}

unsigned int Random::getRandomGenSeed()
{
	return this->randomGenSeed;
}

Random * Random::getInstance()
{
	if (!uniqueinstance) {
		uniqueinstance = new Random();
        std::cout<<"Creating random generator"<<std::endl;
    }
	return uniqueinstance;
}


Random * Random::getInstance(unsigned int shuffleSeed = 0, unsigned int aucSeed = 0, unsigned int DESeed = -1, unsigned int randomGenSeed = 0)
{
	if (!uniqueinstance) {
        std::cout<<"Creating random generator"<<std::endl;
		uniqueinstance = new Random(shuffleSeed, aucSeed, DESeed, randomGenSeed);
	}
	return uniqueinstance;
}

void Random::reSeed(unsigned int shuffleSeed, unsigned int aucSeed, unsigned int randomGenSeed){
    if (!uniqueinstance) {
        std::cout<<"Error, generator not initialized"<<std::endl;
    }
    this->shuffleSeed = shuffleSeed;
    this->aucSeed = aucSeed;
    std::mt19937 shuffGen(shuffleSeed);
    this->shuffleGen= shuffGen;
    std::mt19937 aGen(aucSeed);
    this->aucGen = aGen;
}

void Random::reSeed(){
    if (!uniqueinstance) {
        std::cout<<"Error, generator not initialized"<<std::endl;
    }
    std::random_device rd;
    this->shuffleSeed = rd();
    this->aucSeed = rd();
    this->shuffleSeedVec.push_back(this->shuffleSeed);
    this->aucSeedVec.push_back(this->aucSeed);
    std::mt19937 shuffGen(this->shuffleSeed);
    this->shuffleGen = shuffGen;
    std::mt19937 aGen(this->aucSeed);
    this->aucGen = aGen;
}

std::vector<unsigned int> Random::getShuffleVector(){
    return shuffleSeedVec;
}

std::vector<unsigned int> Random::getAucVector(){
    return aucSeedVec;
}

Random::Random(unsigned int shuffleSeed, unsigned int aucSeed, unsigned int DESeed, unsigned int randomGenSeed)
{
    std::random_device rd;
    this->shuffleSeed = shuffleSeed == 0 ? rd() : shuffleSeed;
    this->aucSeed = aucSeed == 0 ? rd() : aucSeed;
    this->DESeed = DESeed == 0 ? rd() : DESeed;
    this->randomGenSeed = randomGenSeed == 0 ? rd() : randomGenSeed;
    std::mt19937 shuffGen(this->shuffleSeed);
    this->shuffleGen = shuffGen;
    std::mt19937 aGen(this->aucSeed);
    this->aucGen = aGen;
    std::mt19937 DEGen(this->DESeed);
    this->DEGen = DEGen;
    std::mt19937 randomGen(this->randomGenSeed);
    this->randomGen = randomGen;
    double num = (double)randomGen();
    std::cout<<"Initializing random shuffle generator with seed "<<shuffleSeed<<std::endl;
    std::cout<<"Initializing random auc generator with seed "<<aucSeed<<std::endl;
    std::cout<<"Initializing random DE generator with seed "<<this->DESeed<<std::endl;
    std::cout<<"Initializing random measure generator with seed "<<randomGenSeed<<std::endl;
}

Random::Random()
{
    std::random_device rd;
    this->shuffleSeed = rd();
    this->aucSeed = rd();
    this->shuffleSeedVec.push_back(this->shuffleSeed);
    this->aucSeedVec.push_back(this->aucSeed);
    std::mt19937 shuffGen(this->shuffleSeed);
    this->shuffleGen= shuffGen;
    std::mt19937 aGen(this->aucSeed);
    this->aucGen = aGen;
}
Random::~Random()
{

}






#pragma once
#include <random>

class Random
{
public:
	unsigned int getShuffleSeed();
	unsigned int getAucSeed();
	unsigned int getDESeed();
	unsigned int getRandomGenSeed();
	static Random* getInstance(unsigned int seed, unsigned int aucSeed, unsigned int DESeed, unsigned int randomGenSeed);
	static Random* getInstance();
	std::mt19937 shuffleGen;
	std::mt19937 aucGen;
	std::mt19937 DEGen;
	std::mt19937 randomGen;
    void reSeed(unsigned int shuffleSeed, unsigned int aucSeed, unsigned int randomGenSeed);
    void reSeed();
    std::vector<unsigned int> getShuffleVector();
    std::vector<unsigned int> getAucVector();
private:
	unsigned int shuffleSeed;
	unsigned int aucSeed;
	unsigned int DESeed;
	unsigned int randomGenSeed;
    std::vector<unsigned int> shuffleSeedVec;
    std::vector<unsigned int> aucSeedVec;
	Random(Random const&) = delete;             // Copy construct
	Random(Random&&) = delete;                  // Move construct
	Random& operator=(Random const&) = delete;  // Copy assign
	Random& operator=(Random &&) = delete;      // Move assign
	Random(unsigned int shuffleSeed, unsigned int aucSeed, unsigned int DESeed, unsigned int randomGenSeed);
	Random();
	~Random();
};


#pragma once
#include <list>
#include <vector>
#include "Chain.h"
#include "FileData.h"

// class containing static methods for informed search
class SearchStrategy
{
public:
	static std::list<Chain> informedSearch(FileData* f);
private:
	// private methods used by informedSearch method
	static void assignSignals(FileData* f);
	static int assignSignals(FileData* f, double error);
	static void findUniqueCombination(FileData* f,  std::list<Chain>* l);
	static bool isSignal(std::vector<residueType>::iterator i);
	static void makeResidueChains(FileData* f,  std::list<Chain>* l);
	static void minError(FileData*, std::list<Chain>*);
	static double assignmentError(Chain, Chain);
};
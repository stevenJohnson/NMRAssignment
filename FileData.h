#pragma once
#include "Residue.h"
#include <vector>

//serves as a wrapper for all data in the input file
class FileData
{
public:
	//FileData constructors
	FileData();
	FileData(int nOR, int eT, std::vector<Residue> rV, std::vector<residueType> rC);

	//public data members allow for easy access in SearchStrategy class
	int numberOfResidues;
	int experimentTypes;
	std::vector<Residue> residueVector;
	std::vector<residueType> referenceChain;
};
#include "FileData.h"
#include <iostream>

using namespace std;

//FileData default constructor. Used to create an instance of FileData 
//without explicitly setting all data members. Does not take in any 
//parameters or inputs from the user. Does not display or return anything. 
FileData::FileData(void)
{
	//initial all data members to default conditions
	numberOfResidues = 0;
	experimentTypes = 0;
	vector<Residue> residueVector = vector<Residue>();
	vector<residueType> referenceChain = vector<residueType>();
}

//FileData constructor. Used to set all data parameters concisely. Takes 
// the number of residue (nOR), a vector containing the residues (rV), and 
// the reference chain (rC). Does not take in  any inputs from the user. 
//Does not display or return anything. 
FileData::FileData(int nOR, int eT, vector<Residue> rV, vector<residueType> rC)
{
	numberOfResidues = nOR;
	experimentTypes = eT;
	vector<Residue> residueVector = rV;
	vector<residueType> referenceChain = rC;
}
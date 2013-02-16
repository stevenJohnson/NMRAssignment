#pragma once
#include <list>
#include "Residue.h"

//wrapper for a list of consecutive residues 
class Chain
{
private:
	//list of consecutive residues 
    std::list<Residue> residues;
    
public:
	//constructors for the various circumstances of new chain creation
	Chain();
    Chain(Residue, Residue);
    Chain(Residue);
    Chain(Chain, Residue);
    Chain(Residue, Chain);
    Chain(Chain, Chain);

	//accessor functions for the chain and important parts of the chain
    std::list<Residue>& getResidues();
	Residue first();
	Residue last();

	//error of the chain function
	double chainError();
};

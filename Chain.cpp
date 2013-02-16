#include "Chain.h"

using namespace std;//IMPORT THE STANDARD TEMPLATE LIBRARY!!!!

//chain constructor. Used to create a chain of length zero. Takes in
//no parameters or user inputs. Does not return or output anything. 
Chain::Chain()
{
	//nothing. Chain of length zero is created 
}

//chain constructor. Used to create a chain of two individual residues. 
//Takes in no parameters or user inputs. Does not return or output 
//anything. 
Chain::Chain(Residue a, Residue b)
{
    residues.push_back(a);
    residues.push_back(b);
}

//chain constructor. Used to create a chain of one individual residue. 
//Takes in no parameters or user inputs. Does not return or output 
//anything. 
Chain::Chain(Residue c)
{
    residues.push_back(c);
}

//chain constructor. Used to create a chain of one chain and one 
//residue. Takes in no parameters or user inputs. Does not return or output 
//anything. 
Chain::Chain(Chain a,Residue b)
{
    residues.splice(residues.end(), a.getResidues());
    residues.push_back(b);
}

//chain constructor. Used to create a chain of one residue and one 
//chain. Takes in no parameters or user inputs. Does not return or output 
//anything. 
Chain::Chain(Residue a, Chain b)
{
    residues.push_back(a);
    residues.splice(residues.end(), b.getResidues());
}

//chain constructor. Used to create a chain of two chains. Takes
//in no parameters or user inputs. Does not return or output anything. 
Chain::Chain(Chain a, Chain b)
{
    residues.splice(residues.end(), a.getResidues());
    residues.splice(residues.end(), b.getResidues());
}

//residue list accessor function. Used to access the list of residues 
//stored in this chain. Does not take in any parameters or inputs from
//the user. Does not display anything. Returns the list of residues.
list<Residue>& Chain::getResidues()
{
    return residues;
}

//first residue accessor function. Used to access the first element of
//residues stored in this chain. Does not take in any parameters or inputs 
//from the user. Does not display anything. Returns the first residues.
Residue Chain::first()
{
	return residues.front();
}

//last residue accessor function. Used to access the last element of
//residues stored in this chain. Does not take in any parameters or inputs 
//from the user. Does not display anything. Returns the last residues.
Residue Chain::last()
{
	return residues.back();
}

//chain error function. Used to calculate the error within this chain using
//the percent error between each residue. Does not take in any parameters or inputs 
//from the user. Does not display anything. Returns the chain error. 
double Chain::chainError()
{
	if(residues.size() < 2)
	{
		return 0;
	}

	double totalError = 0;
	list<Residue>::iterator iteri = residues.begin();
	advance(iteri, 1);
	for(list<Residue>::iterator iteriM1 = residues.begin(); iteri != residues.end(); advance(iteri, 1))
	{
		iteriM1 = iteri; 
		advance(iteriM1, -1);
		totalError += abs((iteriM1->getCAI()[0] - iteri->getCBIM()[0]) / iteriM1->getCBI()[0]) + 
			abs((iteriM1->getCBI()[0] - iteri->getCBIM()[0]) / iteriM1->getCBI()[0]);
	}
	return totalError;
}

#include "SearchStrategy.h"
#include <vector>
#include <iostream>
#include "Residue.h"
#include <math.h>

using namespace std;

// method that repeatedly calls assignSignals to
// assign a known number of signal values
//
// Filedata* f : filedata object read in from file - contains reference chain,
// residue chain, the number of residues, and the experiment types.
void SearchStrategy::assignSignals(FileData* f)
{
	int assignment = 0; 
	double error = 0; 
	do
	{
		assignment = assignSignals(f, error);
		if(assignment < 0)
		{
			error += 0.001;
		}
		else
		{
			error -= 0.001;
		}
	} while (assignment != 0);
}

// converts the filedata's residue vector into a list of chains l
//
// Filedata* f : filedata object read in from file - contains reference chain,
// residue chain, the number of residues, and the experiment types.
//
// l : pointer to list of chains that is to be used for assignment
void SearchStrategy::makeResidueChains(FileData* f, list<Chain>* l)
{
	for(vector<Residue>::iterator i = f->residueVector.begin(); i < f->residueVector.end(); i++)
	{
		Chain c(*i);
		l->push_front(c);
	}
}

// method that tests for alanine, serine, and threonine, three residues with unique alpha and beta values for initial assignment
//
// Filedata* f : filedata object read in from file - contains reference chain,
// residue chain, the number of residues, and the experiment types.
//
// double error: error tolerance allowed from reference values
//
// returns number of signals not assigned properly (ideally 0)
int SearchStrategy::assignSignals(FileData* f, double error)
{
	int ala = 0, thr = 0 ,ser = 0/*, val = 0, ile = 0*/;//check for possible overassignment
	for(vector<Residue>::iterator i = f->residueVector.begin(); i < f->residueVector.end(); i++)
	{
		double cai = i->getCAI()[0];
		double cbi = i->getCBI()[0];
		//Test for alanine
		if((cai > 51.5 * (1.0 - error) && cai < 54.8 * (1.0 + error)) && (cbi > 18.3 * (1.0 - error) && cbi < 21.1 * (1.0 + error)))
		{
			i->setRT(ALA);
			ala++;
		}
		//Test for Serine
		else if((cai > 57.5 * (1.0 - error) && cai < 60.9 * (1.0 + error)) && (cbi > 63.1 * (1.0 - error) && cbi < 65.2 * (1.0 + error)))
		{
			i->setRT(SER);
			ser++;
		}
		//Test for Threonine
		else if((cai > 61.1 * (1.0 - error) && cai < 65.6 * (1.0 + error)) && (cbi > 68.9 * (1.0 - error) && cbi < 70.8 * (1.0 + error)))
		{
			i->setRT(THR);
			thr++;
		}
		/*
		else if((cai > 60.1 * (1.0 - error) && cai < 64.6 * (1.0 + error)) && (cbi > 37.6 * (1.0 - error) && cbi < 39.9 * (1.0 + error)))
		{
		i->setRT(ILE);
		ile++;
		}
		else if((cai > 60.8 * (1.0 - error) && cai < 66.2 * (1.0 + error)) && (cbi > 31.5 * (1.0 - error) && cbi < 33.9 * (1.0 + error)))
		{
		i->setRT(VAL);
		val++;
		}
		*/
	}
	for(size_t i = 0; i < f->referenceChain.size(); i++)
	{
		if(f->referenceChain[i] == ALA)
		{
			ala--;
		}
		else if(f->referenceChain[i] == SER)
		{
			ser--;
		}
		else if(f->referenceChain[i] == THR)
		{
			thr--;
		}
		/*
		else if(f->referenceChain[i] == VAL)
		{
		val--;
		}
		else if(f->referenceChain[i] == ILE)
		{
		ile--;
		}
		*/
	}
	if((ala + ser + thr /*+ val + ile*/) > 0)//If any of these aren't zero (over/under)assignment has occured
	{
		cerr << "\nOVERassignment has occured\n";
	}
	else if((ala + ser + thr /*+ val + ile*/) < 0)//If any of these aren't zero (over/under)assignment has occured
	{
		cerr << "\nUNDERassignment has occured\n";
	}
	return ala + ser + thr /*+ val + ile*/;
}

// contains calls to the methods to process the filedata object for search purposes
//
// Filedata* f : filedata object read in from file - contains reference chain,
// residue chain, the number of residues, and the experiment types.
//
// returns: list of chain(s) that the search finds possible
list<Chain> SearchStrategy::informedSearch(FileData* f)
{
	list<Chain>* listOfChains = new list<Chain>(); 
	// PROCESS OF ASSIGNMENT
	// [1] assigns the signal residues
	assignSignals(f);

	// [2] fills list of chains from filedata's residueVector
	makeResidueChains(f, listOfChains);

	// [3] assignment unique combinations of residues
	findUniqueCombination(f, listOfChains);

	// [4] creates chain with lowest error
	minError(f, listOfChains);

	return *listOfChains;
}

// takes a list of chains and performs an iterative greedy search on them, returning the list of chains
// containing only the result from that search. search allows each residue to be the first in the chain, and then
// does the greedy search finding the residue that minimizes the error in the chain to be placed next. then takes
// the results and returns the one with the least chain error
//
// Filedata* f : filedata object read in from file - contains reference chain,
// residue chain, the number of residues, and the experiment types.
void SearchStrategy::minError(FileData* f, list<Chain>* l)
{
	list<Chain> results;

	for(size_t i = l->size(); i > 0; i--)
	{
		list<Chain> tmpList(*l);

		list<Chain>::iterator ity = tmpList.begin();
		advance(ity, i-1);
		Chain master = *ity;
		tmpList.erase(ity);
		// ity now points to the chain we want to start with and is removed from the tmpList

		// add chain with lowest error to master repeatedly
		while(tmpList.size() > 0)
		{
			double minError = DBL_MAX;
			list<Chain>::iterator minErrorIter;

			for(list<Chain>::iterator it = tmpList.begin(); it != tmpList.end(); it++)
			{
				double d = assignmentError(master, *it);
				if(d < minError)
				{
					minError = d;
					minErrorIter = it;
				}
			}
			master = Chain(master, *minErrorIter);
			tmpList.erase(minErrorIter);
		}

		results.push_back(master);
	}

	// loop through results, and return the chain with the lowest sum of error
	double minError = DBL_MAX;
	Chain result;

	for(list<Chain>::iterator it = results.begin(); it != results.end(); it++)
	{
		double d = it->chainError();
		if(d < minError)
		{
			minError = d;
			result = *it;
		}
	}

	// erase the list of chains and add the result to it
	l->clear();
	l->push_back(result);
}

// method looks for pairs of signal values and attmepts to form chains with them using the reference chain
//
// Filedata* f : filedata object read in from file - contains reference chain,
// residue chain, the number of residues, and the experiment types.
//
// list<Chain>* l: pointer to chain list that contains the chains of the dataset being processed
void SearchStrategy::findUniqueCombination(FileData* f, list<Chain>* l)
{
	for(vector<residueType>::iterator i = f->referenceChain.begin(); i < f->referenceChain.end() - 1; i++)
	{
		if(isSignal(i) && isSignal(i+1))
		{
			residueType first = *i;
			residueType second = *(i+1);
			list<list<Chain>::iterator> *firstList = new list<list<Chain>::iterator>();
			list<list<Chain>::iterator> *secondList = new list<list<Chain>::iterator>();

			// store the iterators to the possible residues for the two residues being checks
			for(list<Chain>::iterator resit = l->begin(); resit != l->end(); resit++)
			{
				if(resit->last().getRT() == first) 
				{
					firstList->push_back(resit);
				}
				if(resit->first().getRT() == second)
				{
					secondList->push_back(resit);
				}
			}

			double tolerance = .02;
			double firstAlpha = 0.0;
			double firstBeta = 0.0;
			double secondAlphaM1 = 0.0;
			double secondBetaM1 = 0.0;
			list<list<list<Chain>::iterator>::iterator> *matchesA = new list<list<list<Chain>::iterator>::iterator>();
			list<list<list<Chain>::iterator>::iterator> *matchesB = new list<list<list<Chain>::iterator>::iterator>();

			// should be in a tolerance loop later
			for(size_t k = 0; k < firstList->size(); k++)
			{
				list<list<Chain>::iterator>::iterator tmpIter = firstList->begin();
				advance(tmpIter, k);

				firstAlpha = (*tmpIter)->last().getCAI()[0];
				firstBeta = (*tmpIter)->last().getCBI()[0];

				for(size_t j = 0; j < secondList->size(); j++)
				{
					list<list<Chain>::iterator>::iterator tmpBIter = secondList->begin();
					advance(tmpBIter, j);

					secondAlphaM1 = (*tmpBIter)->first().getCAIM()[0];
					secondBetaM1 = (*tmpBIter)->first().getCBIM()[0];

					// if alpha and beta values are within tolerance of eachother
					if(((firstAlpha - secondAlphaM1) / firstAlpha) < tolerance && ((firstAlpha - secondAlphaM1) / firstAlpha) > tolerance*-1)
					{
						if(((firstBeta - secondBetaM1) / firstBeta) < tolerance && ((firstBeta - secondBetaM1) / firstBeta) > tolerance*-1)
						{
							if(*tmpIter != *tmpBIter)
							{
								matchesA->push_back(tmpIter);
								matchesB->push_back(tmpBIter);
							}
						}
					}
				}
			}

			// processing of matches
			if(matchesA->size() > 0 && matchesB->size() > 0)
			{
				list<list<list<Chain>::iterator>::iterator>::iterator bestA = matchesA->begin();
				list<list<list<Chain>::iterator>::iterator>::iterator bestB = matchesB->begin();

				double bestEval = assignmentError(***bestA, ***bestB);

				list<list<list<Chain>::iterator>::iterator>::iterator tmpIter = matchesA->begin();
				list<list<list<Chain>::iterator>::iterator>::iterator tmpBIter = matchesB->begin();

				int count = matchesA->size();
				for(size_t count = matchesA->size(); count > 0; count--)
				{
					if(assignmentError(***tmpIter, ***tmpBIter) < bestEval)
					{
						bestEval = assignmentError(***tmpIter, ***tmpBIter);
						bestA = tmpIter;
						bestB = tmpBIter;
					}

					advance(tmpIter, 1);
					advance(tmpBIter, 1);
				}

				l->push_back(Chain(***bestA, ***bestB));
				l->erase(**bestA);
				l->erase(**bestB);
			}

			delete matchesA;
			delete matchesB;
			delete firstList;
			delete secondList;
		}
	}
}

// method calculates error between two chains, used in finduniquecombination method
//
// Chain upper: the first chain, uses i values
// Chain lower: the second chain, uses i-1 values
//
// returns: double of error amount
double SearchStrategy::assignmentError(Chain upper, Chain lower)
{
	double AiM1 = (lower.first()).getCAIM()[0];
	double BiM1 = (lower.first()).getCBIM()[0];
	double Ai = (upper.last()).getCAI()[0];
	double Bi = (upper.last()).getCBI()[0];

	double squashedError = abs((AiM1 - Ai) / Ai) + abs((BiM1 - Bi)/ Bi);

	return squashedError;
}

// method to return whether or not a residue is of type alanine, thr, or serine
//
// vector<residueType>::iterator i: iterator to the residuetype you want to check for if it's a signal
//
// returns bool: true = is a signal
bool SearchStrategy::isSignal(vector<residueType>::iterator i)
{
	bool signal = false;
	if(*i == ALA || *i == THR || *i == SER /*|| *i == ILE || *i == VAL*/)
	{
		signal = true;
	}
	return signal;
}
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "Residue.h"
#include "SearchStrategy.h"
#include "FileData.h"

using namespace std;

// residue type reader function. Converts the input file string to a vector of reference
// residue types for for help in data analysis. 
//
// string referenceString: string passed in to be parsed into vector of residueTypes
// returns: vector of residueTypes parsed from input string
vector<residueType> rCStringToVector(string referenceString)
{
	vector<residueType> r; 
	for(int i = 0; referenceString[i] != '#'; i += 3)//# symbols serves as break in the loop
	{
		char rS[3];
		rS[0] = referenceString[i];
		rS[1] = referenceString[i + 1];
		rS[2] = referenceString[i + 2];

		//Set the enum residue type for the reference chain
		if(rS[0] == 'U' && rS[1] == 'D' && rS[2] == 'F')
		{
			r.push_back(UDF);
		}
		else if(rS[0] == 'L' && rS[1] == 'E' && rS[2] == 'U')
		{
			r.push_back(LEU);
		}
		else if(rS[0] == 'C' && rS[1] == 'Y' && rS[2] == 'o')
		{
			r.push_back(CYo);
		}
		else if(rS[0] == 'C' && rS[1] == 'Y' && rS[2] == 'r')
		{
			r.push_back(CYr);
		}
		else if(rS[0] == 'A' && rS[1] == 'L' && rS[2] == 'A')
		{
			r.push_back(ALA);
		}
		else if(rS[0] == 'A' && rS[1] == 'S' && rS[2] == 'P')
		{
			r.push_back(ASP);
		}
		else if(rS[0] == 'G' && rS[1] == 'L' && rS[2] == 'U')
		{
			r.push_back(GLU);
		}
		else if(rS[0] == 'P' && rS[1] == 'H' && rS[2] == 'E')
		{
			r.push_back(PHE);
		}
		else if(rS[0] == 'G' && rS[1] == 'L' && rS[2] == 'Y')
		{
			r.push_back(GLY);
		}
		else if(rS[0] == 'H' && rS[1] == 'I' && rS[2] == 'S')
		{
			r.push_back(HIS);
		}
		else if(rS[0] == 'I' && rS[1] == 'L' && rS[2] == 'E')
		{
			r.push_back(ILE);
		}
		else if(rS[0] == 'L' && rS[1] == 'Y' && rS[2] == 'S')
		{
			r.push_back(LYS);
		}
		else if(rS[0] == 'M' && rS[1] == 'E' && rS[2] == 'T')
		{
			r.push_back(MET);
		}
		else if(rS[0] == 'A' && rS[1] == 'S' && rS[2] == 'N')
		{
			r.push_back(ASN);
		}
		else if(rS[0] == 'P' && rS[1] == 'R' && rS[2] == 'O')
		{
			r.push_back(PRO);
		}
		else if(rS[0] == 'G' && rS[1] == 'L' && rS[2] == 'N')
		{
			r.push_back(GLN);
		}
		else if(rS[0] == 'A' && rS[1] == 'R' && rS[2] == 'G')
		{
			r.push_back(ARG);
		}
		else if(rS[0] == 'S' && rS[1] == 'E' && rS[2] == 'R')
		{
			r.push_back(SER);
		}
		else if(rS[0] == 'T' && rS[1] == 'H' && rS[2] == 'R')
		{
			r.push_back(THR);
		}
		else if(rS[0] == 'V' && rS[1] == 'A' && rS[2] == 'L')
		{
			r.push_back(VAL);
		}
		else if(rS[0] == 'T' && rS[1] == 'R' && rS[2] == 'P')
		{
			r.push_back(TRP);
		}
		else if(rS[0] == 'T' && rS[1] == 'Y' && rS[2] == 'R')
		{
			r.push_back(TYR);
		}
		else
		{
			cerr << "\nRESIUDE:rCStringToVector:residue type in reference chain not known/n";
		}
	}
	return r;
}

// file data reader function. Reads in the input text file. Put data in FileData wrapper
// class for easy transfer of data
//
// string inputFileName: name of file containing dataset
//
// returns: filedata object containing information read in from file
FileData readDataFile(string inputFileName)
{
	FileData f = FileData();
	fstream inFile(inputFileName.c_str(), ios::in);//opens input file
	//READ IN THE FILE
	string comment;
	//skip over the comments
	getline(inFile, comment);
	getline(inFile, comment);
	getline(inFile, comment);
	//read in actual data
	inFile >> f.numberOfResidues;//STORE: number of residues
	string referenceChainString;
	getline(inFile, comment);//Advacne to the next line
	getline(inFile, referenceChainString);//reference chain as a string
	f.referenceChain = rCStringToVector(referenceChainString);//STORE: reference chain as a vector

	inFile >> f.experimentTypes;//STORE: experiment type (which experiments are being used)
	vector<Residue> r; 
	for(int i = 0; i < f.numberOfResidues; i++)
	{
		double n, h;
		vector<double> cai = vector<double>();
		vector<double> cbi = vector<double>();
		vector<double> caim1 = vector<double>();
		vector<double> cbim1 = vector<double>(); 
		inFile >> n;
		inFile >> h;
		if(f.experimentTypes == 10)//This experiment type corresponds to only having cai, cbi, cai-1, cbi-1
		{
			double caiD1, cbiD1, caim1D1, cbim1D1;
			inFile >> caiD1;
			inFile >> cbiD1;
			inFile >> caim1D1;
			inFile >> cbim1D1;
			cai.push_back(caiD1);
			cbi.push_back(cbiD1);
			caim1.push_back(caim1D1);
			cbim1.push_back(cbim1D1);
		}
		else if(f.experimentTypes == 11)
		{
			double caiD1, cbiD1, caim1D1, cbim1D1, caim1D2, cbim1D2;
			inFile >> caiD1;
			inFile >> cbiD1;
			inFile >> caim1D1;
			inFile >> cbim1D1;
			inFile >> caim1D2;
			inFile >> cbim1D2;
			cai.push_back(caiD1);
			cbi.push_back(cbiD1);
			caim1.push_back(caim1D1);
			cbim1.push_back(cbim1D1);
			caim1.push_back(caim1D2);
			cbim1.push_back(cbim1D2);
		}
		else
		{
			cerr << "Experiment type not known";
		}
		Residue res = Residue(h, n, cai, cbi, caim1, cbim1); 
		r.push_back(res);// STORE: the chain of residues in a vector
	}
	f.residueVector = r;// STORE: residue change is stores into the data structure
	return f;
}

// output results function. Iterates through the resulting residue chain in order to
// concisely print it to the console for testing purposes.
//
// list<Chain> listy: list of chains to print out
void outputResults(list<Chain> listy)
{
	cout << endl << "Ending Results ::: " << endl; // start line for output
	for(list<Chain>::iterator ity = listy.begin(); ity != listy.end(); ity++)
	{
		list<Residue> residues = ity->getResidues();
		for(list<Residue>::iterator resit = residues.begin(); resit != residues.end(); resit++)
		{
			cout << *resit << endl; // call output function for the residue
		}
		cout << endl;
	}
}

// main function. Serves only as the kick-off point for the data input function, informed 
// search class, and results output.
int main()
{
	FileData test = readDataFile("testFile4.txt"); // store data in FileData wrapper class
	outputResults(SearchStrategy::informedSearch(&test)); // Pass data to informed search for analysis and output
	system("PAUSE");
}
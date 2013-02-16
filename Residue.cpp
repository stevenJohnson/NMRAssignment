#include "Residue.h"
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

// residue default constructor. Should never be used in actual implementation but present
// to assist in debugging effort.
Residue::Residue()
{
	//Should never be used
	cerr << "ERROR:RESIDUE: Using default residue constructor"; 
}

// residue constructor. Used to store data from experiments without explicitly knowing the
// actual residue type of the residue. passes in hydrogen and nitrogen values, and vectors of doubles for alpha, beta, and their minus ones
Residue::Residue(double h, double n, vector<double> cai, vector<double> cbi, vector<double> caiminus1, vector<double> cbiminus1)
{
	rT = UDF;//stores the residue type UNDEFINED (UDF)
	cAi = cai;
	cBi = cbi;
	cAiMINUS1 = caiminus1;
	cBiMINUS1 = cbiminus1;
	hydrogen = h;
	nitrogen = n;
}

// residue constructor. Used to store data from experiments with explicitly knowing the
// actual residue type of the residue. passes in hydrogen and nitrogen values, and vectors of doubles for alpha, beta, and their minus ones
Residue::Residue(double h, double n, vector<double> cai, vector<double> cbi, vector<double> caiminus1, vector<double> cbiminus1, residueType rt)
{
	rT = rt;//stores the explictly defined residue type
	cAi = cai;
	cBi = cbi;
	cAiMINUS1 = caiminus1;
	cBiMINUS1 = cbiminus1;
	hydrogen = h;
	nitrogen = n;
}

//residue type mutator function. Used to change the residue type with the enumerated data type
void Residue::setRT(residueType rt)
{
	rT = rt;
}

//residue type mutator function. Used to change the residue type with a corresponding integer
void Residue::setRT(int rt)
{
	rT = (residueType)rt;
}

//residue type accessor. Used to access to the assigned residue type of this residue
residueType Residue::getRT()
{
	return rT;
}

//carbon alpha of residue i accessor. Used to access all CAI chemcial shift values for this residue
vector<double> Residue::getCAI() 
{
	return 	cAi;
}

//carbon beta of residue i accessor. Used to access all CBI chemcial shift values for this residue
vector<double> Residue::getCBI() 
{
	return 	cBi;
}

//carbon alpha of residue i-1 accessor. Used to access all CAIM chemcial shift values for this residue
vector<double> Residue::getCAIM() 
{
	return 	cAiMINUS1;
}

//carbon beta of residue i-1 accessor. Used to access all CBIM chemcial shift values for this residue
vector<double> Residue::getCBIM() 
{
	return 	cBiMINUS1;
}

//hydrogen of residue (i implicit) accessor. Used to access hydrogen chemical shift values for this residue
double Residue::getH() 
{
	return 	hydrogen;
}

//nitrogen of residue (i implicit) accessor. Used to access nitrogen chemical shift values for this residue
double Residue::getN() 
{
	return 	nitrogen;
}

//toString output function. Used to compactly output all relevant residue data (e.g. residue type, chemical shifts, etc.)
string Residue::toString() const
{
	stringstream rS;

	string residueType = "ERROR:RESIDUE:toString():residue type not defined";//initialized to a error message for debugging
	//find the corresponding string to the enumerated residue type assigned
	switch (rT)
	{
	case 0:
		residueType = "UDF";
		break;
	case 1:
		residueType = "ALA";
		break;
	case 2:
		residueType = "CYr";
		break;
	case 3:
		residueType = "CYo";
		break;
	case 4:
		residueType = "ASP";
		break;
	case 5:
		residueType = "GLU";
		break;
	case 6:
		residueType = "PHe";
		break;
	case 7:
		residueType = "GlY";
		break;
	case 8:
		residueType = "HIS";
		break;
	case 9:
		residueType = "ILE";
		break;
	case 10:
		residueType = "LYS";
		break;
	case 11:
		residueType = "LEU";
		break;
	case 12:
		residueType = "MET";
		break;
	case 13:
		residueType = "ASN";
		break;
	case 14:
		residueType = "PRO";
		break;
	case 15:
		residueType = "GLN";
		break;
	case 16:
		residueType = "ARG";
		break;
	case 17:
		residueType = "SER";
		break;
	case 18:
		residueType = "THR";
		break;
	case 19:
		residueType = "VAL";
		break;
	case 20:
		residueType = "TRP";
		break;
	case 21:
		residueType = "TYR";
		break;
	}
	//only outputs the chemical shift value for the first experiment (modification necesary in the future)
	rS << "RT: " << residueType << " N:"  << fixed << setprecision(3) << nitrogen << " H:" << hydrogen << " Cai:" << cAi[0] << 
		" Cbi:" << cBi[0] << " Ca(i-1):" << cAiMINUS1[0] << " Cb(i-1):" << cBiMINUS1[0];

	return rS.str();
}

//output stream function. Used to display the relevant resiude information to the terminal 
std::ostream& operator<<(std::ostream &strm, const Residue &a)
{
	strm << a.toString();
	return strm;
}

//residue destructor. Deallocates any dynamically created variables 
Residue::~Residue(void)
{
	//no dynamically allocated memory yet!
}
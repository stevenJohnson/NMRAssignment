#pragma once
#include <vector>

//enumerated values corresponding to the residues found in protein chains
enum residueType {UDF, ALA, CYr, CYo, ASP, GLU, PHE, GLY, HIS, ILE, LYS, LEU, MET, ASN, PRO, GLN, ARG, SER, THR, VAL, TRP, TYR};
class Residue
{
public:
	//constructors (only used when reading in data)
	Residue();
	Residue(double h, double n, std::vector<double> cai, std::vector<double> cbi, std::vector<double> caiminus1, std::vector<double> cbiminus1);
	Residue(double h, double n, std::vector<double> cai, std::vector<double> cbi, std::vector<double> caiminus1, std::vector<double> cbiminus1, residueType rT);

	//residue type accessor and mutator functions
	void setRT(residueType rt);
	void setRT(int rt);
	residueType getRT();

	//accessor functions to the residue data
	double getH();
	double getN();
	std::vector<double> getCAI();
	std::vector<double> getCBI();
	std::vector<double> getCAIM();
	std::vector<double> getCBIM();

	//output streams function 
	std::string toString() const;
	friend std::ostream& operator<<(std::ostream &strm, const Residue &a);

	//destructor
	~Residue(void);
private:
	//residue type (emunerated type)
	residueType rT;

	//residue numerical data
	std::vector<double> cAi;
	std::vector<double> cBi;
	std::vector<double> cAiMINUS1;
	std::vector<double> cBiMINUS1;
	double hydrogen;
	double nitrogen;
};
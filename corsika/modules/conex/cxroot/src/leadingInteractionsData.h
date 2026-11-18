#ifndef _include_leadingInteractionsData_h_
#define _include_leadingInteractionsData_h_

#include <iostream>
#include <iomanip>
#include <vector>

class LeadingInteractionsData {
public:
 LeadingInteractionsData() : pId(0), pEnergy(0), energyCM(0), kinel(0), mult(0),
    matg(0), Depth(0), Height(0) {}
 LeadingInteractionsData(int _pid, double _kinel, double _pEnergy, double _energyCM,
                 int _mult, int _matg, double _depth, double _height) :
  pId(_pid), pEnergy(_pEnergy), energyCM(_energyCM), kinel(_kinel), mult(_mult),
    matg(_matg), Depth(_depth),
    Height(_height) {}
  int pId;           // parent id
  double pEnergy;    // parent energy
  double energyCM;   // center-of-mass energy
  double kinel;      // inelasticity
  int mult;          // multiplicity
  int matg;          // target nucleus mass
  double Depth;      // depth
  double Height;     // height
};

class LeadingInteractionsParticleData {
 public:
  LeadingInteractionsParticleData(double e, double m, double px, double py, double pz,
				  int id, int nInt)
  {fEnergy=e; fMass=m; fpx=px; fpy=py; fpz=pz; fID=id; fnInt=nInt;}
  double fEnergy;
  double fMass;
  double fpx;
  double fpy;
  double fpz;
  int fID;
  int fnInt;
  friend bool operator<(const LeadingInteractionsParticleData& lhs,
			const LeadingInteractionsParticleData& rhs);
private:
  LeadingInteractionsParticleData();
};

bool
operator<(const LeadingInteractionsParticleData& lhs,
	  const LeadingInteractionsParticleData& rhs);


extern "C" {
  void outpart1_(double& energy, int& id, double& mass,
		 double& px, double& py, double& pz,
		 int &origin, int& n);
  void outpart2_(double& energyCMS, double &ptlIntIn, int& multxs, int& matargxs,
		 double& depthxs, double& heightxs);
  void outpart3_(double &id, double& generation, double& energy);
}

extern std::vector<LeadingInteractionsParticleData> gInteractionParticleData;
extern std::vector<LeadingInteractionsData> gInteractionData;
extern LeadingInteractionsData gInteractionDataPartial;

#endif

#include <leadingInteractionsData.h>

#include <iostream>
#include <iomanip>

using namespace std;

std::vector<LeadingInteractionsParticleData> gInteractionParticleData;
std::vector<LeadingInteractionsData> gInteractionData;
LeadingInteractionsData gInteractionDataPartial;


bool
operator<(const LeadingInteractionsParticleData& lhs,
	  const LeadingInteractionsParticleData& rhs)
{
  if (lhs.fnInt==rhs.fnInt)
    return (lhs.fEnergy<rhs.fEnergy);
  else
    return (lhs.fnInt>rhs.fnInt);
}


// called from cxafinal (ADD PARTILCE)
void
outpart1_(double &energy, int &id, double &mass,
	  double &px, double &py, double &pz,
	  int &/*origin*/,int &nInt)
/*
 */
{
  gInteractionParticleData.push_back(LeadingInteractionsParticleData(energy,mass,px,py,pz,id,(int)nInt));
}


// called from cxafinal (END INTERACTION)
void
outpart2_(double& energyCMS, double& ptlIntIn, int& multxs, int& matargxs,
          double& depthxs, double& heightxs)
{
  gInteractionDataPartial.energyCM = energyCMS;
  gInteractionDataPartial.kinel = ptlIntIn;
  gInteractionDataPartial.mult = multxs;
  gInteractionDataPartial.matg = matargxs;
  gInteractionDataPartial.Depth = depthxs;
  gInteractionDataPartial.Height = heightxs;
  gInteractionData.push_back(gInteractionDataPartial);
}

// called from s2d (NEW INTERACTION)
void
outpart3_(double& id, double& /*generation*/, double& energy)
{
  gInteractionDataPartial = LeadingInteractionsData((int)id, 0.0, energy, 0, 0, 0, 0.0, 0.0);
}

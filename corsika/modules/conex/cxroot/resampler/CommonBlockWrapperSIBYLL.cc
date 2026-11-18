#include <CommonBlockSIBYLL.h>
#include <CommonBlockWrapperSIBYLL.h>
#include <CommonBlockParticleSIBYLL.h>
#include <Verbosity.h>

#include <iostream>
#include <iomanip>
#include <map>

using namespace std;
using namespace resample;

CommonBlockWrapperSIBYLL::CommonBlockWrapperSIBYLL(CommonBlockSIBYLL* data, bool cms) 
: ParticleBlock(gMaxParticlesSIBYLL, cms),
  fData(data) {
  // IdentifyLeadingParticle( (cms ? eBoth : eForward ) );
}

CommonBlockWrapperSIBYLL::CommonBlockWrapperSIBYLL(CommonBlockSIBYLL& data, bool cms) 
: ParticleBlock(gMaxParticlesSIBYLL, cms),
  fData(&data) {
  // IdentifyLeadingParticle( (cms ? eBoth : eForward ) );
}

ParticleBlockEntry& CommonBlockWrapperSIBYLL::GetEntry(int i) {
  if (fParticles.count(i)==0) {
    fParticles[i] = new CommonBlockParticleSIBYLL(&GetData(), i); 
  }
  return *(fParticles[i]);
}

const ParticleBlockEntry& CommonBlockWrapperSIBYLL::GetEntry(int i) const {
  if (fParticles.count(i)==0) {
    fParticles[i] = new CommonBlockParticleSIBYLL(const_cast<CommonBlockSIBYLL*>(&GetData()), i); //, (i==fLeadingParticleIndex));
  }
  return *(fParticles[i]);
}

void CommonBlockWrapperSIBYLL::Copy(const int i1, const int i2) {
  
  CommonBlockParticleSIBYLL p1(fData, i1);
  CommonBlockParticleSIBYLL p2(fData, i2); // todo: sub-optimal
  
  // p2.fCharge = p1.fCharge; // todo: not needed !!??
  
  p2.SetMomentum(0, p1.GetMomentum(0));
  p2.SetMomentum(1, p1.GetMomentum(1));
  p2.SetMomentum(2, p1.GetMomentum(2));
  p2.SetMomentum(3, p1.GetMomentum(3));
  p2.SetMomentum(4, p1.GetMomentum(4));
  p2.SetId(p1.GetId());
}


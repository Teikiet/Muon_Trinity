#include <CommonBlockCONEX.h>
#include <CommonBlockWrapperCONEX.h>
#include <CommonBlockParticleCONEX.h>
#include <Verbosity.h>

#include <iostream>
#include <iomanip>
#include <map>

using namespace std;
using namespace resample;

CommonBlockWrapperCONEX::CommonBlockWrapperCONEX(CommonBlockCONEX* data, bool cms) 
: ParticleBlock(gMxpltlxs, cms), 
  fData(data) {
  // IdentifyLeadingParticle( (cms ? eBoth : eForward ) );
}
CommonBlockWrapperCONEX::CommonBlockWrapperCONEX(CommonBlockCONEX& data, bool cms) 
: ParticleBlock(gMxpltlxs, cms),
  fData(&data) {
  // IdentifyLeadingParticle( (cms ? eBoth : eForward ) );
}

ParticleBlockEntry& CommonBlockWrapperCONEX::GetEntry(int i) {
  if (fParticles.count(i)==0) {
    fParticles[i] = new CommonBlockParticleCONEX(&GetData(), i); 
  }
  return *(fParticles[i]);
}

const ParticleBlockEntry& CommonBlockWrapperCONEX::GetEntry(int i) const {
  if (fParticles.count(i)==0) {
    fParticles[i] = new CommonBlockParticleCONEX(const_cast<CommonBlockCONEX*>(&GetData()), i);//, (i==fLeadingParticleIndex));
  }
  return *(fParticles[i]);
}

void CommonBlockWrapperCONEX::Copy(int i1, int i2) {
  
  CommonBlockParticleCONEX p1(fData, i1);
  CommonBlockParticleCONEX p2(fData, i2); // todo: sub-optimal
  
  p2.SetMomentum(0, p1.GetMomentum(0));
  p2.SetMomentum(1, p1.GetMomentum(1));
  p2.SetMomentum(2, p1.GetMomentum(2));
  p2.SetMomentum(3, p1.GetMomentum(3));
  p2.SetMomentum(4, p1.GetMomentum(4));
  p2.SetFormationTime(p1.GetFormationTime());
  p2.SetDestructionTime(p1.GetDestructionTime());
  p2.SetFormationPoint(0, p1.GetFormationPoint(0));
  p2.SetFormationPoint(1, p1.GetFormationPoint(1));
  p2.SetFormationPoint(2, p1.GetFormationPoint(2));
  p2.SetFormationPoint(3, p1.GetFormationPoint(3));
  p2.SetIbtlxs(0, p1.GetIbtlxs(0));
  p2.SetIbtlxs(1, p1.GetIbtlxs(1));
  p2.SetIbtlxs(2, p1.GetIbtlxs(2));
  p2.SetIbtlxs(3, p1.GetIbtlxs(3));
  p2.SetIfrptlxs(0, p1.GetIfrptlxs(0));
  p2.SetIfrptlxs(1, p1.GetIfrptlxs(1));
  p2.SetStatus(p1.GetStatus());
  p2.SetMotherId(p1.GetMotherId());
  p2.SetId(p1.GetId());
  p2.SetFatherId(p1.GetFatherId());
  p2.SetOrigin(p1.GetOrigin());
}


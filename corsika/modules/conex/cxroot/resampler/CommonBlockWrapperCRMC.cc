#include <CommonBlockCRMC.h>
#include <CommonBlockWrapperCRMC.h>
#include <CommonBlockParticleCRMC.h>
#include <Verbosity.h>

#include <iostream>
#include <iomanip>
#include <map>

using namespace std;
using namespace resample;

CommonBlockWrapperCRMC::CommonBlockWrapperCRMC(CommonBlockCRMC* data, bool cms) 
: ParticleBlock(gMxpltlxs, cms), 
  fData(data) {
  // IdentifyLeadingParticle( (cms ? eBoth : eForward ) );
}
CommonBlockWrapperCRMC::CommonBlockWrapperCRMC(CommonBlockCRMC& data, bool cms) 
: ParticleBlock(gMxpltlxs, cms),
  fData(&data) {
  // IdentifyLeadingParticle( (cms ? eBoth : eForward ) );
}

ParticleBlockEntry& CommonBlockWrapperCRMC::GetEntry(int i) {
  if (fParticles.count(i)==0) {
    fParticles[i] = new CommonBlockParticleCRMC(&GetData(), i); 
  }
  return *(fParticles[i]);
}

const ParticleBlockEntry& CommonBlockWrapperCRMC::GetEntry(int i) const {
  if (fParticles.count(i)==0) {
    fParticles[i] = new CommonBlockParticleCRMC(const_cast<CommonBlockCRMC*>(&GetData()), i);//, (i==fLeadingParticleIndex));
  }
  return *(fParticles[i]);
}

void CommonBlockWrapperCRMC::Copy(int from, int to) {
  
  //CommonBlockParticleCRMC pFrom(fData, from);
  //CommonBlockParticleCRMC pTo(fData, to); // todo: sub-optimal  
  fData->Copy(from, to); // this is better
}


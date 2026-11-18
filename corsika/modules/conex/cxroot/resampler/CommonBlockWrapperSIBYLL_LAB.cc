#include <CommonBlockSIBYLL_LAB.h>
#include <CommonBlockWrapperSIBYLL_LAB.h>
#include <CommonBlockParticleCONEX.h>
#include <Verbosity.h>
#include <CommonBlockWrapperCONEX.h> // only for debugging

#include <iostream>
#include <iomanip>
#include <map>
#include <cstdlib>
#include <algorithm>

using namespace std;
using namespace resample;

// summary-check>      sumEnergy (w/o leading)=164218978.906 (delta=9.493e+06) sumCharge=12 (delta=1) sumN=372 (delta=-145)
 
CommonBlockWrapperSIBYLL_LAB::CommonBlockWrapperSIBYLL_LAB(CommonBlockCONEX& data, 
							   CommonBlockSIBYLL_LAB& mapping, 
							   const int idInteraction,
							   const bool cms) 
: CommonBlockWrapperCONEX(data, cms), 
  fIdInteraction(idInteraction), 
  fMapping(&mapping) {
  fIndexMap.reserve(fMapping->numSecondariesOfInteraction[idInteraction]);
}

ParticleBlockEntry& CommonBlockWrapperSIBYLL_LAB::GetEntry(int i) {
  if (i<0 || int(fIndexMap.size())<=i) {
   FindIndex(i);
  }
  return CommonBlockWrapperCONEX::GetEntry(fIndexMap[i]);
}

const ParticleBlockEntry& CommonBlockWrapperSIBYLL_LAB::GetEntry(int i) const {
  if (i<0 || int(fIndexMap.size())<=i) {
    FindIndex(i);
  }
  return CommonBlockWrapperCONEX::GetEntry(fIndexMap[i]);
}

int CommonBlockWrapperSIBYLL_LAB::FindIndex(const int i) const {
  if (i<0) {
    cerr << "CommonBlockWrapperSIBYLL_LAB::FindIndex(i="<< i << ") "
	 << " ERROR (FATAL): No particles at negative indices !!!" 
	 << endl;
    exit(2);
  }
  
  if (int(fIndexMap.size())>i) {
    return fIndexMap[i];
  }
  
  fIndexMap.resize(i+1);
  int index = 0;
  int search = 0;
  for( ; index<fData->number; ++index) {
    /*cout << " find-index: " << index << " n_block=" << fData->number 
      << " mapping=" << fMapping->id[index] << endl;
    */
    if (fMapping->id[index] == fIdInteraction+1) { // index: fortan -> C++ 
      fIndexMap[search] = index;
      if ( (search++) == i ) { // count particle of this interaction
	break;
      }
    }
  }
  if (index==fData->number) {
    cerr << "CommonBlockWrapperSIBYLL_LAB::FindIndex(i="<< i << ") "
	 << " ERROR (FATAL): No particle found!!!" 
	 << endl;
    exit(2);
  }
  return index;
}

int CommonBlockWrapperSIBYLL_LAB::FindIndexReverse(const int index) const { 
  vector<int>::iterator it = find(fIndexMap.begin(), fIndexMap.end(), index); 
  if (it==fIndexMap.end()) 
    return -1;
  return it-fIndexMap.begin();
}

// this function uses indices in the full conex-stack
void CommonBlockWrapperSIBYLL_LAB::Copy(const int index1, const int index2) {
  
  if (ruVerbosity>=1) {
    cout << " CommonBlockWrapperSIBYLL_LAB::Copy " << index1 << " to " << index2 << endl;
  }
  
  if (index1==index2) {
    if (ruVerbosity>=1) {
      cerr << " CommonBlockWrapperSIBYLL_LAB::Copy WARNING: source and destiny are identical " 
	   << " index=" << index1 
	   << endl;
    }
    return;
  } 
  
  CommonBlockParticleCONEX p1(fData, index1);
  CommonBlockParticleCONEX p2(fData, index2); // todo: sub-optimal
  
  // -- copy data
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
  
  // -- copy mapping in CommonBlockSIBYLL_LAB
  fMapping->id[index2] = fMapping->id[index1];
  
  // -- check index map + leading particles
  const int i2 = FindIndexReverse(index2);
  if (i2>=0) {
    // check for overriding the leading particle(s)
    if (i2==fLeadingParticleForwardIndex) {
      fLeadingParticleForwardIndex = -1;
      if (ruVerbosity>=0) {
	cerr << " WARNING: leading particle (" << "FORWARD" 
	     << ") deleted : i=" << i2 
	     << endl;
      }
    }
    if (i2==fLeadingParticleBackwardIndex) {
      fLeadingParticleBackwardIndex = -1;
      if (ruVerbosity>=0) {
	cerr << " WARNING: leading particle (" << "BACKWARD" 
	     << ") deleted : i=" << i2 
	     << endl;
      }
    }
    // check for shifting the leading particle indices
    if (fMapping->id[index1] != fIdInteraction+1) { // only if we insert particle from different interaction !
      fIndexMap.erase(fIndexMap.begin()+i2);
      if (fLeadingParticleForwardIndex>i2) {
      if (ruVerbosity>=3) {
	cout << " CommonBlockWrapperSIBYLL_LAB::Copy INFO shift fwd leading particle index to: " << fLeadingParticleForwardIndex-1 << endl;
      }
      --fLeadingParticleForwardIndex;
      }
      if (fLeadingParticleBackwardIndex>i2) {
      if (ruVerbosity>=3) {
	cout << " CommonBlockWrapperSIBYLL_LAB::Copy INFO  shift bwd leading particle index to: " << fLeadingParticleBackwardIndex-1 << endl;
      }
      --fLeadingParticleBackwardIndex;
      }
    }
  }
  
  // check for moving the leading particle
  const int i1 = FindIndexReverse(index1);
  if (i1>=0) {
    // fIndexMap.erase(fIndexMap.begin()+i1);
    if (i1==fLeadingParticleForwardIndex) {
      if (i2<0) {
	cerr << " CommonBlockWrapperSIBYLL_LAB::Copy ERROR: lost index of leading particle (FORWARD) !!!!!!"
	     << endl;
      }
      fLeadingParticleForwardIndex = i2;
      if (ruVerbosity>4) {
	cout << " CommonBlockWrapperSIBYLL_LAB::Copy  INFO: leading particle (" << "FORWARD" 
	     << ") moved : " << i1 << "->" << i2
	     << endl;
      }
    }
    if (i1==fLeadingParticleBackwardIndex) {
      if (i2<0) {
	cerr << " CommonBlockWrapperSIBYLL_LAB::Copy ERROR: lost index of leading particle (Backward) !!!!!!"
	     << endl;
      }
      fLeadingParticleBackwardIndex = i2;
      if (ruVerbosity>4) {
	cout << " CommonBlockWrapperSIBYLL_LAB::Copy INFO: leading particle (" << "BACKWARD" 
	     << ") moved : " << i1 << "->" << i2
	     << endl;
      }
    }
  } /*else {
      for (unsigned int k=0; k<fIndexMap.size(); ++k) {
      cout << " indexmap[" << k << "]=" << fIndexMap[k] << endl;
      }
      }*/
}


int CommonBlockWrapperSIBYLL_LAB::Size() const {
  return fMapping->numSecondariesOfInteraction[fIdInteraction];
}

void CommonBlockWrapperSIBYLL_LAB::IncrementSize() {
  fMapping->id[fData->number] = fIdInteraction+1; // assure proper mapping (fortan index)
  fMapping->numSecondariesOfInteraction[fIdInteraction]++; 
  ++fData->number;
}

void CommonBlockWrapperSIBYLL_LAB::DecrementSize() {
  --fData->number;
  fMapping->numSecondariesOfInteraction[fIdInteraction]--;
  const int lastIndex = fData->number;
  const int iLast = FindIndexReverse(lastIndex);
  if (iLast>=0) {
    fIndexMap.erase(fIndexMap.begin()+iLast);
  }
}

void CommonBlockWrapperSIBYLL_LAB::Delete(int i) {
  
  if (ruVerbosity>4) {
    cout << " CommonBlockWrapperSIBYLL_LAB::Delete i=" << i << endl;
  }
  
  const int index = FindIndex(i);
  
  if (ruVerbosity>5) {
    cout << " CommonBlockWrapperSIBYLL_LAB::Delete, mapped to index=" << index << endl;
  }
  
  const int lastIndex = fData->number-1;
  // const int iLast = FindIndexReverse(lastIndex);

  if (index <= lastIndex && lastIndex>0) {
    
    Copy(lastIndex, index);
    
    if (fParticles.count(lastIndex)!=0) {
      if (ruVerbosity>4) {
	cout << " CommonBlockWrapperSIBYLL_LAB::Delete INFO: remove particle from cache : " << lastIndex << endl;
      }
      delete fParticles[lastIndex];
      fParticles.erase(lastIndex);
    }  
    
    if (fParticles.count(index)!=0) {
      if (ruVerbosity>4) {
	cout << " CommonBlockWrapperSIBYLL_LAB::Delete INFO: remove particle from cache : " << index << endl;
      }
      delete fParticles[index];
      fParticles.erase(index);
    }
  }
  
  DecrementSize();
}

bool CommonBlockWrapperSIBYLL_LAB::Duplicate(int i) {

  const int index = FindIndex(i);
  
  const int lastIndex = fData->number-1;
  
  if (lastIndex>=GetMaxCapacity()-1) {
    ++fParticleBlockSizeOverflowCounter;
    if (ruVerbosity>0) {
      cerr << " CommonBlockWrapperSIBYLL_LAB::Duplicate> WARNING cannot generate particle, no space left!" 
	   << endl;
    }
    return false;
  }
   
  IncrementSize();    

  if (ruVerbosity>4) {
    cout << " CommonBlockWrapperSIBYLL_LAB::Duplicate particle " 
	 << index << " into " << lastIndex+1
	 << endl;
  }
  
  Copy(index, lastIndex+1);
  
  if (ruVerbosity>4) {
    // if index 'i' was leading particle, you duplicated it!
    if (IsLeading(i)) {
      cerr << " CommonBlockWrapperSIBYLL_LAB::Duplicate WARNING: duplicate leading particle : " << index << endl;
      }
  }
  return true;
}

void CommonBlockWrapperSIBYLL_LAB::Dump() const {
  cout << " -------- BEGIN DUMP (FULL) STACK ------- " << endl;
  cout << " *** mapping: " << endl;
  for (int i=0; i<fMapping->numberOfInteractions; ++i) {
    cout << " n_sec = " << fMapping->numSecondariesOfInteraction[i] << endl; 
  }
  cout << " *** particles: " << endl;
  CommonBlockWrapperCONEX wrap(fData, false);
  const int n = wrap.Size();
  for (int iParticle=0; iParticle<n; ++iParticle) {
    const ParticleBlockEntry& p = wrap.GetEntry(iParticle);
    
    const double E = p.GetEnergy();

    cout << " index=" << setw(3) << iParticle
	 << " id=" << setw(5) << p.GetId()
	 << " (" << setw(11) << p.GetName() << ")"
	 << " E=" << setw(18) << E // p.GetMomentum(3)
	 << " m=" << setw(12) << setprecision(6) << p.GetMass()
	 << " nucleon=" << setw(3) << fMapping->id[iParticle]
	 << "/" << fMapping->numberOfInteractions;

    if (!p.IsGood()) { // These are real particles, all the rest is crap
      cout << " BAD PARTICLE"; // , status=" << p.GetStatus();
    } else {
      const int i = FindIndexReverse(iParticle);   
      if (i>=0) {
	cout << " [mapped to " << i << "]";
	if (fLeadingParticleForwardIndex==i) 
	  cout << " [fwd-leading]";
	if (fLeadingParticleBackwardIndex==i) 
	  cout << " [bkw-leading]";	  
      }      
    }
    cout << endl;
  } // loop particles
  cout << " -------- END DUMP (FULL) STACK ------- " << endl;
}

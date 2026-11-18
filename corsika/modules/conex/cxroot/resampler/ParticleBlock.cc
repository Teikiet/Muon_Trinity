#include <ParticleBlock.h>
#include <ParticleBlockEntry.h>
#include <SummaryInfo.h>
#include <Verbosity.h>

#include <map>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
using namespace std;
using namespace resample;

ParticleBlock::~ParticleBlock() 
{
  if (fParticleBlockSizeOverflowCounter>0) {
    cerr << "ParticleBlock> WARNING: " << fParticleBlockSizeOverflowCounter 
	 << " particles could not have been created, since the fortran common block was full! " << endl;
  }
  for(map<int, ParticleBlockEntry*>::iterator iParticle = fParticles.begin();
      iParticle != fParticles.end(); ++iParticle) {
    delete iParticle->second;
  }
}

bool
ParticleBlock::IsLeading(int i, ELeading where) const 
{  
  if (!fIsCMS && where==eBackward) {
    cerr << " ParticleBlock::IsLeading ERROR: Cannot get leading particle in backward direction\n"
	 << "                                 for particles in lab-system!" 
	 << endl;
    exit(1);
  }
  bool result = false;
  switch(where) {
  case eForward: 
    result = fLeadingParticleForwardIndex==i;
    break;
  case eBackward: 
    result = fLeadingParticleBackwardIndex==i;
    break;
  case eBoth: 
    result = fLeadingParticleForwardIndex==i || fLeadingParticleBackwardIndex==i;
    break;
  }
  return result;
}

int
ParticleBlock::GetLeadingParticleIndex(ELeading where) const 
{
  if (!fIsCMS && where!=eForward) {
    cerr << " ParticleBlock::GetLeadingParticleIndex ERROR: Cannot access leading particle in non-forward direction\n"
	 << "                                 for particles in lab-system!" 
	 << endl;
    exit(1);
  }
  int result = -1;
  switch(where) {
  case eForward: 
    result = fLeadingParticleForwardIndex;
    break;
  case eBackward: 
    result = fLeadingParticleBackwardIndex;
    break;
  case eBoth: 
    cerr << " ParticleBlock::GetLeadingParticleIndex ERROR: Cannot access both leading particle indices at the same time!"
	 << endl;
    exit(1);
    break;
  }
  return result;
}

void
ParticleBlock::SetLeadingParticleIndex(const int index, ELeading where) const 
{
  if (!fIsCMS && where!=eForward) {
    cerr << " ParticleBlock::SetLeadingParticleIndex ERROR: Cannot access leading particle in non-forward direction\n"
	 << "                                 for particles in lab-system!" 
	 << endl;
    exit(1);
  }
  switch(where) {
  case eForward: fLeadingParticleForwardIndex=index;
    break;
  case eBackward: fLeadingParticleBackwardIndex=index;
    break;
  case eBoth: 
    cerr << " ParticleBlock::GetLeadingParticleIndex ERROR: Cannot access both leading particle indices at the same time!"
	 << endl;
    exit(1);
    break;
  }
}


void
ParticleBlock::Delete(int i) 
{  
  if (ruVerbosity>4) {
    cout << " ParticleBlock::Delete index=" << i << endl;
  }
  
  const int lastIndex = Size()-1;
  if (i <= lastIndex && lastIndex>0) {
    
    Copy(lastIndex, i);
    
    // iModus==0 : forward leading particle
    // iModus==1 : backward leading particle (only CMS)
    for (int iModus = 0; iModus <= (fIsCMS ? 1 : 0); ++iModus) {
      
      const ELeading where = (iModus==0 ? eForward : eBackward);
      
      if (IsLeading(i, where)) {
	
	// if leading particle was 'i', it is deleted and lost
	SetLeadingParticleIndex(-1, where);
	if (ruVerbosity>=0) {
	  cout << " WARNING: leading particle (" << (iModus ? "FORWARD" : "BACKWARD")
	       << ") deleted : " << i 
	       << endl;
	}
	
      } else if (IsLeading(lastIndex, where)) {
	
	// if leading particle was 'lastIndex' is is now 'i'
	if (ruVerbosity>4) {
	  cout << " INFO: leading particle (" << (iModus ? "FORWARD" : "BACKWARD")
	       << ") moved : " << lastIndex << "->" << i
	       << endl;
	}
	SetLeadingParticleIndex(i, where); 
      }
    } // loop forward/backward leading

	
    if (fParticles.count(lastIndex)!=0) {
      if (ruVerbosity>4) {
	cout << " ParticleBlock::Delete INFO: remove particle from cache : " << lastIndex << endl;
      }
      delete fParticles[lastIndex];
      fParticles.erase(lastIndex);
    }  
    
    if (fParticles.count(i)!=0) {
      if (ruVerbosity>4) {
	cout << " ParticleBlock::Delete INFO: remove particle from cache : " << i << endl;
      }
      delete fParticles[i];
      fParticles.erase(i);
    }
  }
  
  DecrementSize();
}


bool
ParticleBlock::Duplicate(int i) 
{ 
  const int lastIndex = Size()-1;
  
  if (lastIndex>=fMaxParticleCapacity-1) {
    ++fParticleBlockSizeOverflowCounter;
    if (ruVerbosity>0) {
      cerr << " ParticleBlock::Duplicate> WARNING cannot generate particle, no space left!" 
	   << endl;
    }
    
    return false;
  } 
    
  IncrementSize();
  
  if (ruVerbosity>4) {
    cout << " ParticleBlock::Duplicate particle " 
	 << i << " into " << lastIndex+1
	 << endl;
  }
  
  Copy(i, lastIndex+1);
  
  if (ruVerbosity>4) {
    // if index 'i' was leading particle, you duplicated it!
    if (IsLeading(i)) {
      cout << " warning: duplicate leading particle : " << i << endl;
    }
  }    
  return true;
}


void
ParticleBlock::IdentifyLeadingParticle(ELeading where) const 
{
  /*
    if (ruVerbosity>4) {
    cout << " ParticleBlock::IdentifyLeadingParticle " << where << endl;
    }
  */
  
  // iModus==0 : forward leading particle
  // iModus==1 : backward leading particle
  for (int iModus = (where==eForward||where==eBoth ? 0 : 1);
       iModus <= (where==eBackward||where==eBoth ? 1 : 0);
       ++iModus) {
    
    int leadingIndex = 0;
    int leadingId = 0;
    double leadingEnergy = 0;
    string leadingName;
    
    bool first = true;
    for (int iSecondary = 0; 
	 iSecondary<Size();
	 ++iSecondary) {
      
      const ParticleBlockEntry& secondary = GetEntry(iSecondary);
      
      if (!secondary.IsGood())
	continue;
      
      if (iModus==0 && !secondary.IsForward()) 
	continue;
      
      if (iModus==1 && secondary.IsForward()) 
	continue;
      
      const int id = secondary.GetId();
      const double energy = secondary.GetEnergy();
      
      if (first) {
	first = false;
	leadingIndex = iSecondary;
	leadingId = id;
	leadingEnergy = energy;
	leadingName = secondary.GetName();
      } else if (energy>leadingEnergy) {
	leadingIndex = iSecondary;
	leadingId = id;
	leadingEnergy = energy;
	leadingName = secondary.GetName();
      }
      
    } // loop secondaries
    
    if (ruVerbosity>=1) {
      cout << " IdentifyLeadingParticle> "
	   << (iModus==0 ? " FORWARD " : " BACKWARD")
	   << " index=" << setw(5) << leadingIndex
	   << " " << setw(10) << string('\''+leadingName+'\'')
	   << " Id=" << setw(5) << leadingId 
	   << " E=" << setw(12) << setprecision(3) << leadingEnergy
	   << endl;
    }
    
    if (iModus==0)
      fLeadingParticleForwardIndex = leadingIndex;
    else 
      fLeadingParticleBackwardIndex = leadingIndex;

  }
}


map<int, SummaryInfo> 
  ParticleBlock::MakeSummary(const EIdMode idMode,
			     const EExcludeLeading excludeLeadingMode) const 
{
  
  const bool excludeLeading = (excludeLeadingMode == eExcludeLeading);

  map<int, SummaryInfo> sum;
  
  for(int iSecondary = 0; iSecondary<Size();
      ++iSecondary) {
    
    const ParticleBlockEntry& secondary = GetEntry(iSecondary);
    
    if (!secondary.IsGood()) 
      continue;
    
    if (excludeLeading && IsLeading(iSecondary, eBoth))
      continue;
    
    const int id = (idMode == eTypes ? 
		    (int)abs((double)secondary.GetId()) :
		    secondary.GetId());
    const double totEnergy = secondary.GetEnergy();
    const double mass = secondary.GetMass();
    const double kinEnergy = secondary.GetEnergy() - mass;
    const double charge = secondary.GetCharge();
    
    // ## DEBUG ############################################
    if (ruVerbosity>2) {
      cout << setprecision(3);
      cout << " MakeSummary> "
	   << "[" << setw(4) << iSecondary << "]"
	   << " id=" << setw(5) << secondary.GetId() 
	   << " (" << setw(11) << secondary.GetName() << ")"
	   << " q=" << setw(2) << charge 
	   << " Etot=" << setw(14) << totEnergy  
	   << " Ekin=" << setw(14) << kinEnergy  
	   << " mass=" << setw(7) << setprecision(4) << mass << setprecision(3)
	   << " " << (IsLeading(iSecondary, eForward) ?  "[leading (FOWD)]" : "");
      if (fIsCMS)
	cout << " " << (IsLeading(iSecondary, eBackward) ?  "[leading (BKWD)]" : "");
      cout << endl;
    }
    // ## end-DEBUG ############################################
    
    if (!sum.count(id)) {
      
      sum[id].charge    = charge;
      sum[id].totEnergy = totEnergy;
      sum[id].kinEnergy = kinEnergy;
      sum[id].number    = 1;
      sum[id].minKinEnergy = kinEnergy;
      sum[id].maxKinEnergy = kinEnergy;
      /*
      ostringstream hname;
      hname << "hist_pid_" << id;
      sum[id].hE = new TH1D(hname.str().c_str(), hname.str().c_str(), 50,10,-10);
      sum[id].hE->SetDirectory(0);
      sum[id].hE->Fill(log10(energy));
      */
      
    } else {
      
      sum[id].charge    += charge;
      sum[id].totEnergy += totEnergy;
      sum[id].kinEnergy += kinEnergy;
      sum[id].number ++;
      if (sum[id].minKinEnergy>kinEnergy) sum[id].minKinEnergy = kinEnergy;
      if (sum[id].maxKinEnergy<kinEnergy) sum[id].maxKinEnergy = kinEnergy;
      //sum[id].hE->Fill(log10(energy));
      
    }
  }
  
  return sum;
}


void
ParticleBlock::Dump() const 
{
  cout << " -------- BEGIN DUMP STACK ------- " << endl;
  const int n = Size();
  for (int iParticle=0; iParticle<n; ++iParticle) {
    const ParticleBlockEntry& p = GetEntry(iParticle);
    if (p.IsGood()) { // These are real particles, all the rest is crap
      cout << " index=" << setw(3) << iParticle;
      p.Dump();
    }
  } // loop particles
  cout << " -------- END DUMP STACK ------- " << endl;
}

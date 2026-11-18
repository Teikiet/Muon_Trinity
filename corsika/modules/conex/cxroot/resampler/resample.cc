#include <ParticleBlock.h>
#include <ParticleBlockEntry.h>
#include <CommonBlockParticleCONEX.h>
#include <ResamplingMode.h>
#include <Verbosity.h>

#include <TRandom.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <set>
#include <map>
#include <cmath>
#include <vector>
#include <cstdlib>
using namespace std;

#include <resample.h>


namespace resample {

  /* InitResampling 
   * **************
   *
   * call this to set the random number seed
   */
  void
  InitResampling(const int seed) 
  {
    if (ruVerbosity>0)
      cout << "InitResampling: seed=" << seed << endl;
    gRandom->SetSeed(seed);
  }


  /* MultiplicityResampling
   * *********************
   */
  void
  MultiplicityResampling(const double factor, ParticleBlock& particlesList) 
  {  
    const int nPart = particlesList.Size();
  
    // get leading particles in forward/backward direction
    vector<int> leadingParticleIndices;
    particlesList.IdentifyLeadingParticle(ParticleBlock::eForward);
    const int leadingForwardIndex = particlesList.GetLeadingParticleIndex(ParticleBlock::eForward);
    const int leadingForwardId = particlesList.GetEntry(leadingForwardIndex).GetId(); // for checking 
    const double leadingForwardEnergy = particlesList.GetEntry(leadingForwardIndex).GetEnergy(); // for checking 
    leadingParticleIndices.push_back(leadingForwardIndex);
    int leadingBackwardIndex = -1;
    int leadingBackwardId = -1;
    double leadingBackwardEnergy = -1;
    if (particlesList.IsCenterOfMassSystem()) {
      particlesList.IdentifyLeadingParticle(ParticleBlock::eBackward);
      leadingBackwardIndex = particlesList.GetLeadingParticleIndex(ParticleBlock::eBackward);
      leadingBackwardId = particlesList.GetEntry(leadingBackwardIndex).GetId(); // for checking 
      leadingBackwardEnergy = particlesList.GetEntry(leadingBackwardIndex).GetEnergy(); // for checking 
      leadingParticleIndices.push_back(leadingBackwardIndex);
    }
  
    // get overview of the particle block
    double sumEnergy = 0;
    double sumCharge = 0;
    int sumN = 0;
    map<int, SummaryInfo> sum = particlesList.MakeSummary(ParticleBlock::eTypes, ParticleBlock::eExcludeLeading); 
    for (map<int, SummaryInfo>::iterator i = sum.begin(); i != sum.end(); ++i) {
      sumN      +=  i->second.number;
      sumEnergy += i->second.totEnergy;
      sumCharge += i->second.charge;
      if (ruVerbosity>=1) {
	cout << setprecision(3);
	cout << " summary-particles> id=" << setw(5) << i->first 
	     << " N=" << setw(4) << i->second.number 
	     << " z=" << setw(2) << i->second.charge 
	     << " Etot=" << setw(14) << i->second.totEnergy
	     << " Ekin=" << setw(14) << i->second.kinEnergy
	     << " EkinMin=" << setw(14) << i->second.minKinEnergy
	     << " EkinMax=" << setw(14) << i->second.maxKinEnergy
	     << endl;
      }
    }
  
    double EtotLeading = 0;
    for (unsigned int iLeading=0; 
	 iLeading != leadingParticleIndices.size(); 
	 ++iLeading) {
      const ParticleBlockEntry& leading = particlesList.GetEntry(leadingParticleIndices[iLeading]);
      EtotLeading += leading.GetEnergy();
      /*
	cout << " leading particle (" << (iLeading==0 ? "FORWARD ": "BACKWARD" )
	<< ") " << leading.GetName() << " (id: " << leading.GetId() << ") "
	<< " E_tot=" << leading.GetEnergy()
	<< endl;
      */
    }
    if (ruVerbosity>=1) {
      cout << setprecision(3);
      cout << " summary-particles> sumTotEnergy(w/o leading)=" << sumEnergy 
	   << " Etot_leading=" << EtotLeading
	   << " sumTotEnergy=" << sumEnergy+EtotLeading
	   << " sumCharge=" << sumCharge 
	   << " sumN=" << sumN 
	   << endl;
    }
  
    // number of particles to delete/create
    int deltaMultiplicity = (int)abs(double((factor-1.0) * sumN + 0.5));
    const bool deleteMode = (factor<=1);
  
    if (ruVerbosity>0) {
      cout << " Multiplicity-resample mode: " << (deleteMode?"delete":"create") << " " << deltaMultiplicity << " particles."
	   << " (N_part=" << sumN << ")" << endl; 
    }

    map<int, SummaryInfo> sumChanged;  
    int sizeCurrent = nPart;
    int iDelta = 0;
    int iTry = 0;
    while ((iDelta<deltaMultiplicity) &&    // until all particles are deleted/duplicated
	   (iTry<deltaMultiplicity*2)) {   // or the loop-limit is reached
    
      ++iTry;
    
      if (ruVerbosity>=5) {
	particlesList.Dump(); // extensive debugging
      }
    
      if (deleteMode) {
	sizeCurrent = particlesList.Size(); // changing size (particles are disapearing)
      } else {
	// don't adapt to the growing particle block, only duplicate _initial_ particles
      }
    
      const int iPart = int(gRandom->Rndm() * sizeCurrent);
      const ParticleBlockEntry& secondary = particlesList.GetEntry(iPart);
    
      if (!secondary.IsGood() ||             // do not take BAD particles
	  particlesList.IsLeading(iPart, ParticleBlock::eBoth)) {  // do not harm leading particle
	continue;
      }
    
      const int id = (int)abs((double)secondary.GetId());
    
      if (sumChanged.count(id)==0) { // init for each particle id
	sumChanged[id].charge = 0;
	sumChanged[id].totEnergy = 0;
	sumChanged[id].kinEnergy = 0;
	sumChanged[id].number = 0;
      }
    
      if (deleteMode) {
      
	// we don't want to remove particle of a given Id completely ...
	if (sumChanged[id].number+sum[id].number<=1) 
	  continue;
      
	sumChanged[id].charge -= secondary.GetCharge();
	sumChanged[id].totEnergy -= secondary.GetEnergy();
	sumChanged[id].kinEnergy -= secondary.GetEnergy()-secondary.GetMass();
	sumChanged[id].number--;
	particlesList.Delete(iPart);
	++iDelta;
        
      } else { // create mode

	if (particlesList.Duplicate(iPart)) {
	  sumChanged[id].charge += secondary.GetCharge();
	  sumChanged[id].totEnergy += secondary.GetEnergy();
	  sumChanged[id].kinEnergy += secondary.GetEnergy()-secondary.GetMass();
	  sumChanged[id].number++;
	  ++iDelta;
	}
      
      }
    
    } // end of multiplicity changing loop
  
  
    /********************************************
     * 
     * check energy balance (per particle type, therefore also global)
     *
     ********************************************/
  
    double sumChargeDelta = 0;
    map<int,double> modify;
    for (map<int, SummaryInfo>::const_iterator iChange = sumChanged.begin();
	 iChange != sumChanged.end(); ++iChange) {
    
      const unsigned int id = (int)abs((double)iChange->first);
      const int deltaN = iChange->second.number;
      const double deltaCharge = iChange->second.charge;
      const double deltaEtot = iChange->second.totEnergy;
      const double deltaEkin = iChange->second.kinEnergy;
      const double charge = iChange->second.charge;
      sumChargeDelta += charge;
    
      const double sumEkin = sum[id].kinEnergy + deltaEkin;
      const double sumEtot = sum[id].totEnergy + deltaEtot;
      const double sumMass = sumEtot - sumEkin;
    
      if (sumEkin==0) {
	modify[id] = 1;
      } else {
	modify[id] = (sum[id].totEnergy - sumMass) / sumEkin;
      }
    
      if (modify[id]<0) {
	modify[id] = 0;
      }
    
      // ## DEBUG ############################################
      if (ruVerbosity>=1) {
	cout << setprecision(2);
	cout << " summary-modify>    id=" << setw(5) << id 
	     << " N="     << setw(4) << sum[id].number+deltaN
	     << " z="     << setw(4) << sum[id].charge+deltaCharge;
	if (ruVerbosity>3) {
	  cout << " Etot="  << setw(14)<< sumEtot
	       << " EtotIn=" << setw(14) << sum[id].totEnergy
	       << " sumMass=" << setw(14) << sumMass;
	}
	cout << " Ekin="  << setw(14)<< sumEkin
	     << " dN="    << setw(4) << deltaN
	     << " dz="    << setw(3) << deltaCharge
	     << " dEtot=" << setw(14)<< deltaEtot
	     << " dEkin=" << setw(14)<< deltaEkin
	     << " factor="<< setw(10) << setprecision(8) << modify[id] 
	     << endl;
      }
      // ## end-DEBUG ############################################
    }
  
  
    for(int iSecondary = 0; 
	iSecondary < particlesList.Size();
	++iSecondary) {
    
      ParticleBlockEntry& particle = particlesList.GetEntry(iSecondary);
    
      if (!particle.IsGood() ||
	  particlesList.IsLeading(iSecondary, ParticleBlock::eBoth)) {
	continue;
      }
    
      const int id = (int)abs((double)particle.GetId());
      const double mass = particle.GetMass();
      const double energy = particle.GetEnergy();
      const double px = particle.GetPx();
      const double py = particle.GetPy();
      const double pz = particle.GetPz();
      const double Pi = sqrt(px*px + py*py + pz*pz);
    
      const double factor = ((modify.count(id)!=0) ?
			     modify[id] :
			     1.0);
    
      const double Ekin = energy - mass;
      const double energyNew = mass + Ekin*factor;
      const double Pf = sqrt((energyNew+mass) * (energyNew-mass)); // P = sqrt(e^2 - m^2)
    
      double newPx = 0;
      double newPy = 0;
      double newPz = Pf;
      if (Pi<=0) {
	if (ruVerbosity>0)
	  cout << " WARNING: P_f=" << Pf 
	       << " P_i=" << Pi
	       << " Ekin_i=" << Ekin
	       << " factor=" << factor
	       << " E_f=" << energyNew
	       << " id=" << id
	       << " [-> new momentum in z-direction]"
	       << endl;
      } else {
	if (Pf<=0) {
	  if (ruVerbosity>0)
	    cout << " WARNING: P_f=" << Pf 
		 << " P_i=" << Pi
		 << " Ekin_i=" << Ekin
		 << " factor=" << factor
		 << " E_f=" << energyNew
		 << " id=" << id
		 << " [-> new momentum forecd to ZERO]"
		 << endl;
	  newPz = 0;
	} else {
	  const double momentumFactor = Pf/Pi;
	  newPx = px * momentumFactor;
	  newPy = py * momentumFactor;
	  newPz = pz * momentumFactor;
	}
      }
      particle.SetEnergy(energyNew);
      particle.SetPx(newPx);
      particle.SetPy(newPy);
      particle.SetPz(newPz);    
    }
  
    if (ruVerbosity>=1) {
      cout << " charge> total charge balance = " <<  sumChargeDelta 
	   << " (total number of particles on stack = " << particlesList.Size() << ")" 
	   << endl;
    }
  
  
    /******************************************************
     * 
     * check charge balance (global only, with a maximum deviation of +-1e)
     *
     **************************************************/
  
    for (int iTry=0; iTry<2; ++iTry) {
    
      const int maxTest = particlesList.Size() * 10;
      for (int iTest=0; iTest<maxTest; ++iTest) {
	
	if (fabs(sumChargeDelta)<=1) {
	  // ## DEBUG ############################################
	  if (ruVerbosity>2) {
	    cout << " charge> success: delta=" << sumChargeDelta << endl;
	  }
	  // ## end-DEBUG ############################################
	  iTry = 10; // break outer loop
	  break;
	}
      
	// choose random particle 

	int iParticle = int(gRandom->Rndm() * particlesList.Size());
	ParticleBlockEntry& particle = particlesList[iParticle];

	if (!particle.IsGood() || 
	    particlesList.IsLeading(iParticle, ParticleBlock::eBoth)) { // do not harm leading particle
	  continue;
	}
      
	const double particleCharge = particle.GetCharge();
	const int particleId = particle.GetId();
      
	// check whether particle-type will help to reduce global charge change 
      
	SummaryInfo& typeInfo = sumChanged[particleId];
	double& deltaChargeType = typeInfo.charge;
      
	if ((sumChargeDelta*deltaChargeType>0) ||   // check if type is usefull to reduce charge change
	    (iTry==1) ) {                           // don't check in second try
        
	  if (sumChargeDelta*particleCharge>0) {  // check if particle is usefull to reduce charge change
          
	    // ## DEBUG ############################################
	    if (ruVerbosity>2) {
	      cout << " charge> try=" << setw(4) << iTest << "|" << iTry << " i=" << setw(4) << iParticle 
		   << " id=" << setw(5) << particleId << " charge=" << particleCharge
		   << " deltaCharge=" << deltaChargeType
		   << " sumDeltaCharge=" << sumChargeDelta
		   << endl;
	    }
	    // ## end-DEBUG ############################################
	    
	    // particle.SetCharge(particleCharge*-1); // reverse charge
	    particle.SetId(particle.GetIdOfAntiParticle());     // -> anti-particle (reverse charge)
          
	    if (sumChargeDelta<0) {
	      deltaChargeType += 2;
	      sumChargeDelta += 2;
	    } else {
	      deltaChargeType -= 2;
	      sumChargeDelta -= 2;
	    }
          
	  }
        
	} // if type will do
      
      } // loop iTest

    } // loop iTry
    
    if (ruVerbosity>0)
      cout << " charge> final charge balance = " <<  sumChargeDelta << endl;
    
    map<int, SummaryInfo> sumCheck = particlesList.MakeSummary(ParticleBlock::eTypes, ParticleBlock::eExcludeLeading);
    double sumEnergyCheck = 0;
    double sumChargeCheck = 0;
    int sumNCheck = 0;
    for (map<int, SummaryInfo>::iterator i = sumCheck.begin(); i != sumCheck.end(); ++i) {
      sumEnergyCheck += i->second.totEnergy;
      sumChargeCheck += i->second.charge;
      sumNCheck      += i->second.number;
      if (ruVerbosity>=1) {
	cout << setprecision(3);
	cout << " summary-check>     id=" << setw(5) << i->first 
	     << " N=" << setw(4) << i->second.number 
	     << " z=" << setw(2) << i->second.charge 
	     << " Etot=" << setw(14) << i->second.totEnergy
	     << " Ekin=" << setw(14) << i->second.kinEnergy
	     << " EkinMin=" << setw(14) << i->second.minKinEnergy
	     << " EkinMax=" << setw(14) << i->second.maxKinEnergy
	     << endl;
      }
    }
    
    if (ruVerbosity>0) {
      cout << setprecision(3);
      cout << " summary-check>      sumEnergy (w/o leading)=" << sumEnergyCheck << " (delta=";
      cout << scientific;
      cout << sumEnergyCheck-sumEnergy << ")";
      cout << fixed;
      cout << " sumCharge=" << sumChargeCheck << " (delta=" << sumChargeCheck-sumCharge << ")"
	   << " sumN=" << sumNCheck << " (delta=" << sumNCheck-sumN << ")"
	 << endl;
      
      // ------- for checking ---------
      // particlesList.IdentifyLeadingParticle(ParticleBlock::eForward);
      const int leadingForwardIndexCheck = particlesList.GetLeadingParticleIndex(ParticleBlock::eForward);
      const int leadingForwardIdCheck = particlesList.GetEntry(leadingForwardIndexCheck).GetId();
      const double leadingForwardEnergyCheck = particlesList.GetEntry(leadingForwardIndexCheck).GetEnergy();
      if (leadingForwardId != leadingForwardIdCheck ||
	  leadingForwardEnergy != leadingForwardEnergyCheck ) {
	cerr << "MultiplicityResampling: ERROR LEADING-PARTICLE (forward) lost! " << endl;
	cerr << "       was: index=" << leadingForwardIndex << " id=" << leadingForwardId << " E=" << leadingForwardEnergy << endl
	     << "        is: index=" << leadingForwardIndexCheck << " id=" << leadingForwardIdCheck << " E=" << leadingForwardEnergyCheck << endl;
	exit(3);    
      }    
      if (particlesList.IsCenterOfMassSystem()) {
	// particlesList.IdentifyLeadingParticle(ParticleBlock::eBackward);
	const int leadingBackwardIndexCheck = particlesList.GetLeadingParticleIndex(ParticleBlock::eBackward);
	const int leadingBackwardIdCheck = particlesList.GetEntry(leadingBackwardIndexCheck).GetId();
	const double leadingBackwardEnergyCheck = particlesList.GetEntry(leadingBackwardIndexCheck).GetEnergy();
	if (leadingBackwardId != leadingBackwardIdCheck ||
	    leadingBackwardEnergy != leadingBackwardEnergyCheck) {
	  cerr << "MultiplicityResampling: ERROR LEADING-PARTICLE (backward) lost! " << endl;
	  cerr << "       was: index=" << leadingBackwardIndex << " id=" << leadingBackwardId << " E=" << leadingBackwardEnergy << endl
	       << "        is: index=" << leadingBackwardIndexCheck << " id=" << leadingBackwardIdCheck << " E=" << leadingBackwardEnergyCheck << endl;
	  exit(3);
	}
      }      
      // particlesList.IdentifyLeadingParticle(ParticleBlock::eBoth);
    }
  
  } // end of MultiplicityResampling



  /* ElasticityResampling
   * *********************
   */
  void
  ElasticityResampling(const double factor, ParticleBlock& particlesList) 
  {  
    const int nPart = particlesList.Size();

    // get overview ( for checking )
    map<int, SummaryInfo> sum = particlesList.MakeSummary(ParticleBlock::eTypes, ParticleBlock::eIncludeLeading);
    double sumEnergyTotIni = 0;
    double sumEnergyKinIni = 0;
    //bool first = true;
    for (map<int, SummaryInfo>::iterator i = sum.begin(); i != sum.end(); ++i) {
      sumEnergyTotIni += i->second.totEnergy;
      sumEnergyKinIni += i->second.kinEnergy;
      if (ruVerbosity>=1) {
	cout << setprecision(3);
	cout << " summary-particles> id=" << setw(5) << i->first 
	     << " Etot=" << setw(14) << i->second.totEnergy
	     << " Ekin=" << setw(14) << i->second.kinEnergy
	     << " EkinMin=" << setw(14) << i->second.minKinEnergy
	     << " EkinMax=" << setw(14) << i->second.maxKinEnergy
	     << endl;
      }
    }
    if (ruVerbosity>=1) {
      cout << setprecision(3);
      cout << " summary-particles> sumTotEnergy=" << sumEnergyTotIni 
	   << " sumKinEnergy=" << sumEnergyKinIni
	   << endl;
    }
  
    // iMode==0 : forward
    // iMode==1 : backward
    for (unsigned int iMode=0; 
	 iMode <= (particlesList.IsCenterOfMassSystem() ? 1 : 0 );
	 ++iMode) {
    
      if (ruVerbosity>0 && particlesList.IsCenterOfMassSystem()) {
	cout << " ElasticityResampling in " << (iMode==0 ? "FORWARD " : "BACKWARD" ) 
	     << " direction  -------------------- " << endl;
      }
    
      const ParticleBlock::ELeading direction = (iMode==0 ? 
						 ParticleBlock::eForward : 
						 ParticleBlock::eBackward);
    
      particlesList.IdentifyLeadingParticle(direction);
      const int leadingParticleIndex = particlesList.GetLeadingParticleIndex(direction);
    
      // -----------------------------------------------------------------
      // find second to leading particle
      bool secondFirst = true;
      int secondLeadingIndex = 0;
      int secondLeadingId = 0;
      double secondLeadingEnergy = 0;
      string secondLeadingName;
      int secondSecondIndex = 0;
      int secondSecondId = 0;
      double secondSecondEnergy = 0;
      string secondSecondName;
      // get total/kinetic energy in forward/backward direction
      double sumEnergyTot = 0;
      double sumEnergyKin = 0;
    
      double maxMass = 0;
      int nSecondaries = 0;
      for (int iSecondary = 0; 
	   iSecondary<particlesList.Size();
	   ++iSecondary) {
      
	const ParticleBlockEntry& secondary = particlesList.GetEntry(iSecondary);
      
	if (!secondary.IsGood() ||
	    (particlesList.IsCenterOfMassSystem() &&
	     ((direction==ParticleBlock::eForward && !secondary.IsForward()) ||
	      (direction==ParticleBlock::eBackward && secondary.IsForward())))) {
	  continue;
	}
      
	nSecondaries++;
      
	const string& name = secondary.GetName();
	const double energy = secondary.GetEnergy();
	const double mass = secondary.GetMass();
	const int id = secondary.GetId();
      
	// - energy -
	sumEnergyTot += energy;
	sumEnergyKin += energy - secondary.GetMass();
      
	// - 2nd leading particle -
      
	if (secondFirst) {
	  secondFirst = false;
	  secondLeadingIndex = iSecondary;
	  secondLeadingId = id;
	  secondLeadingEnergy = energy;
	  secondLeadingName = name;
	  maxMass = mass;
	} else {
	  maxMass = std::max(maxMass, mass);
	  if (energy>secondLeadingEnergy) {
	    secondSecondIndex = secondLeadingIndex;
	    secondSecondId = secondLeadingId;
	    secondSecondEnergy = secondLeadingEnergy;
	    secondSecondName = secondLeadingName;
	    secondLeadingIndex = iSecondary;
	    secondLeadingId = id;
	    secondLeadingEnergy = energy;
	    secondLeadingName = name;
	  }
	}
      } // end loop secondaries
    
      const double minLeadingParticleEnergy = std::max(maxMass, sumEnergyTot/nSecondaries);
      
      if (ruVerbosity>0)
	cout << " minLeadingParticleEnergy= " << minLeadingParticleEnergy << " (nSecondaries=" << nSecondaries << ")" << endl;
      
      if (ruVerbosity>=1) {
	cout << " Energy flow in " << (iMode==0 ? "forward" : "backward" ) << "-direction"
	     << " total=" << sumEnergyTot << " kinetic=" << sumEnergyKin << endl; 
	
	cout << " FindSecondLeadingParticle> " 
	     << " index=" << setw(5) << secondSecondIndex
	     << " " << setw(10) << string('\''+secondSecondName+'\'')
	     << " id=" << setw(4) << secondSecondId
	     << " energy=" << setw(15) << setprecision(3) << secondSecondEnergy
	     << endl;
      }
      // ----------------------------------------------------------------
    
    
      const ParticleBlockEntry& leading = particlesList.GetEntry(leadingParticleIndex);
    
      const double massLeading = leading.GetMass();
      const double EleadingTot = leading.GetEnergy();
      const double EleadingKin = leading.GetEnergy() - massLeading;
      const double Etot = sumEnergyTot;
    
      const double minimalElasticity = minLeadingParticleEnergy*1.001 / Etot; // absolute minimum
      const double maximalElasticity = (massLeading + sumEnergyKin)*0.999 / Etot; // absolute maximum 
    
      const double elasticity = EleadingTot/Etot; 
      double newElasticity = elasticity*factor;
      //if (newElasticity<=minimalElasticity || newElasticity>=maximalElasticity) {
      
      //cerr << " ############# RESAMPLE-WARNING ################ " << endl;
      if (ruVerbosity>0)
	cout << " ElasticityResampling WARNING: kel=" << elasticity 
	     << ", kel_new=" << newElasticity 
	     << ", min=" << minimalElasticity
	     << ", max=" << maximalElasticity
	  //<< ". skipping"
	     << endl;
      //cerr << " ############# RESAMPLE-WARNING ################ " << endl;
      // return;
      //}
      if (newElasticity<minimalElasticity) newElasticity = minimalElasticity;
      if (newElasticity>maximalElasticity) newElasticity = maximalElasticity;
    
      //int iterativeLeadingParticleIndex = leadingParticleIndex; // first leading particle
      //bool done = false; // flag for iterative procedure
      //while(!done) {
      //ParticleBlockEntry& leadingParticle = particlesList.GetEntry(iterativeLeadingParticleIndex);
    
      const double EleadingTotNew = EleadingKin * (newElasticity/elasticity) + leading.GetMass();
      if (ruVerbosity>0) {
	cout << " Change elasticity : from " << elasticity << " to " << newElasticity 
	     << endl;
	cout << " Change leadingEtot: from " << EleadingTot << " to " << EleadingTotNew 
	     << endl; 
      }
    
      if (factor>1.0) { // increasing elasticity is straightforward 
      
	const double deltaE = EleadingTotNew - EleadingTot;
	const double sumEnergyKinScale = sumEnergyKin - EleadingKin;
      
	if (sumEnergyKinScale<=0.01) {
	  if (ruVerbosity>0) 
	    cout << " ElasticityResampling> WARNING: Not enough secondaries for resampling!"
		 << " skipping."
		 << endl;
	  continue;
	}
      
	const double secondaryModify = (sumEnergyKinScale - deltaE) / sumEnergyKinScale;
	// modify[id] = (sum[id].totEnergy - sumMass) / sumEkin;
      
	if (ruVerbosity>0)
	  cout << " Secondary kinetic energy scaling-factor: " << secondaryModify << endl;
	
	for (int iPart=0; iPart<nPart; ++iPart) { // loop all particles and modify energies
	
	  ParticleBlockEntry& secondary = particlesList.GetEntry(iPart);
	
	  if (!secondary.IsGood() ||
	      (particlesList.IsCenterOfMassSystem() &&
	       ((direction==ParticleBlock::eForward && !secondary.IsForward()) ||
		(direction==ParticleBlock::eBackward && secondary.IsForward())))) {
	    continue;
	  }
	
	  //const int id = abs((double)secondary.GetId());
	  const double mass  = secondary.GetMass();
	  const double energy = secondary.GetEnergy();
	  const double px = secondary.GetPx();
	  const double py = secondary.GetPy();
	  const double pz = secondary.GetPz();
	  const double Pi = sqrt(px*px + py*py + pz*pz);
	
	  const double modify = ( (particlesList.IsLeading(iPart, direction)) ?  
				  newElasticity/elasticity :                     // leading particle
				  secondaryModify );                             // rest of particles
	
	  const double Ekin = energy - mass;
	  const double energyNew = mass + Ekin*modify;
	  const double Pf = sqrt((energyNew+mass) * (energyNew-mass));
	
	  double newPx = 0;
	  double newPy = 0;
	  double newPz = Pf;
	  if (Pi<=0) {
	    if (ruVerbosity>0)
	      cout << " WARNING: P_f=" << Pf 
		   << " P_i=" << Pi
		   << " Ekin_i=" << Ekin
		   << " factor=" << factor
		   << " E_f=" << energyNew
		   << " [-> new momentum in z-direction]"
		   << endl;
	  } else {
	    if (Pf<=0) {
	      if (ruVerbosity>0)
		cout << " WARNING: P_f=" << Pf 
		     << " P_i=" << Pi
		     << " Ekin_i=" << Ekin
		     << " factor=" << factor
		     << " E_f=" << energyNew
		     << " [-> new momentum forecd to ZERO]"
		     << endl;
	      newPz = 0;
	    } else {
	      const double momentumFactor = Pf/Pi;
	      newPx = px * momentumFactor;
	      newPy = py * momentumFactor;
	      newPz = pz * momentumFactor;
	    }
	  }

	  secondary.SetPx(newPx);
	  secondary.SetPy(newPy);
	  secondary.SetPz(newPz);
	  secondary.SetEnergy(energyNew);
	
	  /*
	  // ## DEBUG ############################################
	  cout << setprecision(2);
	  cout << " summary-modify>    id=" << setw(5) << id 
	  << " N="     << setw(4) << sum[id].number+deltaN
	  << " C="     << setw(4) << sum[id].charge+deltaCharge;
	  if (ruVerbosity>3) 
	  cout << " Etot="  << setw(14)<< sumEtot
	  << " EtotIn=" << setw(14) << sum[id].totEnergy
	  << " sumMass=" << setw(14) << sumMass;
	  cout << " Ekin="  << setw(14)<< sumEkin
	  << " dN="    << setw(4) << deltaN
	  << " dC="    << setw(3) << deltaCharge
	  << " dEtot=" << setw(14)<< deltaEtot
	  << " dEkin=" << setw(14)<< deltaEkin
	  << " factor="<< setw(10) << setprecision(8) << modify[id] 
	  << endl;
	  //}
	  // ## end-DEBUG ############################################
	  */
      
	} // loop nPart
    
      } else { // for decreasing elasticity we have to work iteratively !

	bool reiterate = false;
	int countIterations = 0;
	vector<bool> exclude(nPart, false);
	do { 
	
	  ++countIterations;
	
	  // first check how much energy needs to be re-distributed
	  // (not only the leading particle may be affected by the elastictiy-mod!)
	  double deltaE = 0;
	  double sumE_others = 0;
	  double sumM_others = 0;
	  int nOthers = 0;
	  for (int iPart=0; iPart<nPart; ++iPart) { // loop all particles and modify energies
	
	    if (exclude[iPart])
	      continue;
	  
	    ParticleBlockEntry& secondary = particlesList.GetEntry(iPart);
	  
	    if (!secondary.IsGood() ||
		(particlesList.IsCenterOfMassSystem() &&
		 ((direction==ParticleBlock::eForward && !secondary.IsForward()) ||
		  (direction==ParticleBlock::eBackward && secondary.IsForward())))) {
	      continue;
	    }
	  
	    const double mass  = secondary.GetMass();
	    const double energy = secondary.GetEnergy();
	  
	    // cout << " id=" << setw(4) << iPart << " E=" << setw(15) << energy;

	    if (energy>=EleadingTotNew) {
	      exclude[iPart] = true;
	      deltaE += energy - EleadingTotNew;
	      secondary.SetEnergy(EleadingTotNew);
	      // cout << " leading detla=" << deltaE << endl;
	    } else {
	      ++nOthers;
	      sumE_others += energy;
	      sumM_others += mass;
	      // cout << " others sumE_others=" << sumE_others << " sumM=" << sumM_others << endl;
	    }
	  } // loop secondaries
	
	  if (nOthers == 0) {
	    if (ruVerbosity>0) 
	      cout << " ElasticityResampling> WARNING: No secondaries left for further resampling: break here "
		   << endl;
	    reiterate = false;
	    break;
	  }
	  
	  const double sumEkin_others = sumE_others - sumM_others;
	  const bool noKinE_others = (sumEkin_others <= 0);
	  if (noKinE_others) {
	    if (ruVerbosity>0) 
	      cout << " ElasticityResampling> INFO: The total secondary kinetic energy is: " << sumEkin_others 
		   << endl;
	  }
	  
	  const double secondaryModify = (!noKinE_others ? (sumEkin_others + deltaE) / sumEkin_others : 1);
	  // cout << " sumEkin_others=" << sumEkin_others << " deltaE=" << deltaE << endl;
	  // cout << " Secondary kinetic energy scaling-factor: " << secondaryModify << endl;
	
	  reiterate = false;
	  for (int iPart=0; iPart<nPart; ++iPart) { // loop all particles and modify energies
	  
	    if (exclude[iPart])
	      continue;
	  
	    ParticleBlockEntry& secondary = particlesList.GetEntry(iPart);
	  
	    if (!secondary.IsGood() ||
		(particlesList.IsCenterOfMassSystem() &&
		 ((direction==ParticleBlock::eForward && !secondary.IsForward()) ||
		  (direction==ParticleBlock::eBackward && secondary.IsForward())))) {
	      continue;
	    }
	  
	    const double mass  = secondary.GetMass();
	    const double energy = secondary.GetEnergy();
	  
	    if (energy<=EleadingTotNew) {
	      double energyNew = 0;
	      if (noKinE_others) {
		energyNew += deltaE / nOthers; // emergency sollution: if there is no kin-energy to scale, 
		                               // just add up the energy on top ...
	      } else {
		const double Ekin = energy - mass;
		energyNew = mass + Ekin * secondaryModify;
	      }
	      secondary.SetEnergy(energyNew);
	      // check if we scaled above the leading particle energy, and thus have to re-iterate
	      if (energyNew>EleadingTotNew) {
		reiterate = true;
	      }
	    }
	  } // loop secondaries
	} while(reiterate);
      
	// final loop to recalculate momenta
	for (int iPart=0; iPart<nPart; ++iPart) { // loop all particles and modify energies
	
	  ParticleBlockEntry& secondary = particlesList.GetEntry(iPart);
	
	  if (!secondary.IsGood() ||
	      (particlesList.IsCenterOfMassSystem() &&
	       ((direction==ParticleBlock::eForward && !secondary.IsForward()) ||
		(direction==ParticleBlock::eBackward && secondary.IsForward())))) {
	    continue;
	  }
	
	  //const int id = abs((double)secondary.GetId());
	  const double mass  = secondary.GetMass();
	  const double energy = secondary.GetEnergy();
	  const double px = secondary.GetPx();
	  const double py = secondary.GetPy();
	  const double pz = secondary.GetPz();
	  const double Pi = sqrt(px*px + py*py + pz*pz);
	  const double Pf = sqrt((energy+mass) * (energy-mass));
	
	  double newPx = 0;
	  double newPy = 0;
	  double newPz = Pf;
	  if (Pi<=0) {
	    if (ruVerbosity>0)
	      cout << " WARNING: P_f=" << Pf 
		   << " P_i=" << Pi
		   << " factor=" << factor
		   << " [-> new momentum in z-direction]"
		   << endl;
	  } else {
	    if (Pf<=0) {
	      if (ruVerbosity>0)
		cout << " WARNING: P_f=" << Pf 
		     << " P_i=" << Pi
		     << " factor=" << factor
		     << " [-> new momentum forecd to ZERO]"
		     << endl;
	      newPz = 0;
	    } else {
	      const double momentumFactor = Pf/Pi;
	      newPx = px * momentumFactor;
	      newPy = py * momentumFactor;
	      newPz = pz * momentumFactor;
	    }
	  }
	  
	  secondary.SetPx(newPx);
	  secondary.SetPy(newPy);
	  secondary.SetPz(newPz);
	
	} // loop nPart
      
	if (ruVerbosity>0) {
	  cout << " ElasticityResampling> INFO needed " << countIterations << " iterations to lower elasticity " << endl;
	}

      } // end increase/decrease elasticity
    
      /*
      // CHECK total/kinetic energy in forward/backward direction
      double sumEnergyTot2 = 0;
      double sumEnergyKin2 = 0;
    
      for (int iSecondary = 0; 
      iSecondary<particlesList.Size();
      ++iSecondary) {
      
      const ParticleBlockEntry& secondary = particlesList.GetEntry(iSecondary);
      
      if (!secondary.IsGood() ||
      (direction==ParticleBlock::eForward && !secondary.IsForward()) ||
      (direction==ParticleBlock::eBackward && secondary.IsForward())) {
      continue;
      }

      const double energy = secondary.GetEnergy();
      
      // - energy -
      sumEnergyTot2 += energy;
      sumEnergyKin2 += energy - secondary.GetMass();
      }

      cout << " CHECK Energy flow in " << (iMode==0 ? "forward" : "backward" ) << "-direction"
      << " total=" << sumEnergyTot2 << " kinetic=" << sumEnergyKin2 << endl; 
      */
    
    } // loop forward/backward
  
    map<int, SummaryInfo> sumCheck = particlesList.MakeSummary(ParticleBlock::eTypes, ParticleBlock::eIncludeLeading);
    double sumEnergyTotCheck = 0;
    double sumEnergyKinCheck = 0;
    for (map<int, SummaryInfo>::iterator i = sumCheck.begin(); i != sumCheck.end(); ++i) {
      sumEnergyTotCheck += i->second.totEnergy;
      sumEnergyKinCheck += i->second.kinEnergy;
      if (ruVerbosity>1) {
	cout << setprecision(3);
	cout << " summary-check>     id=" << setw(5) << i->first 
	     << " Etot=" << setw(14) << i->second.totEnergy
	     << " Ekin=" << setw(14) << i->second.kinEnergy
	     << " EkinMin=" << setw(14) << i->second.minKinEnergy
	     << " EkinMax=" << setw(14) << i->second.maxKinEnergy
	     << endl;
      }
    }
    
    if (ruVerbosity>0) {
      cout << setprecision(3);
      cout << " summary-check>      sumTotEnergy=" << sumEnergyTotCheck 
	   << " (delta=" << scientific << sumEnergyTotCheck-sumEnergyTotIni << fixed << ")"
	   << " sumKinEnergy=" << sumEnergyKinCheck 
	   << " (delta=" << sumEnergyKinCheck-sumEnergyKinIni << ")"
	   << endl;     
      // particlesList.FindLeadingParticle();
    }

  } // end of ElasticityResampling


  /*
   * EMRatioResampling
   * *****************
   */
  void
  EMRatioResampling(const double factor, ParticleBlock& particlesList) 
  {  
    const int nPart = particlesList.Size();

    // get leading particles in forward/backward direction
    vector<int> leadingParticleIndices;
    particlesList.IdentifyLeadingParticle(ParticleBlock::eForward);
    leadingParticleIndices.push_back(particlesList.GetLeadingParticleIndex(ParticleBlock::eForward));
    if (particlesList.IsCenterOfMassSystem()) {
      particlesList.IdentifyLeadingParticle(ParticleBlock::eBackward);
      leadingParticleIndices.push_back(particlesList.GetLeadingParticleIndex(ParticleBlock::eBackward));
    }
  
    // get overview of the particle block
    double sumEnergyTot = 0;
    double sumEnergyKin = 0;
    double sumEnergyTotEM = 0;
    double sumEnergyKinEM = 0;
    map<int, SummaryInfo> sum = particlesList.MakeSummary(ParticleBlock::eTypes, ParticleBlock::eExcludeLeading);
    for (map<int, SummaryInfo>::iterator i = sum.begin(); i != sum.end(); ++i) {
      sumEnergyTot += i->second.totEnergy;
      sumEnergyKin += i->second.kinEnergy;
      switch(i->first) {
      case CommonBlockParticleCONEX::eGamma:
      case CommonBlockParticleCONEX::eElectron:
      case CommonBlockParticleCONEX::ePositron:
      case CommonBlockParticleCONEX::ePi0:
	sumEnergyTotEM += i->second.totEnergy;
	sumEnergyKinEM += i->second.kinEnergy;
	break;
      }
      if (ruVerbosity>=1) {
	cout << setprecision(3);
	cout << " summary-particles> id=" << setw(5) << i->first 
	     << " Etot=" << setw(14) << i->second.totEnergy
	     << " Ekin=" << setw(14) << i->second.kinEnergy
	     << " EkinMin=" << setw(14) << i->second.minKinEnergy
	     << " EkinMax=" << setw(14) << i->second.maxKinEnergy
	     << endl;
      }
    }

    double EtotLeading = 0;
    for (unsigned int iLeading=0; 
	 iLeading != leadingParticleIndices.size(); 
	 ++iLeading) {
      const ParticleBlockEntry& leading = particlesList.GetEntry(leadingParticleIndices[iLeading]);
      EtotLeading += leading.GetEnergy();
      if (ruVerbosity>=1) {
	cout << " leading particle (" << (iLeading==0 ? "FORWARD": "BACKWARD" )
	     << ") " << leading.GetName() << "(" << leading.GetId() << ") "
	     << " E_tot=" << leading.GetEnergy()
	     << endl;
      }
    }

    const double sumMassEM = sumEnergyTotEM-sumEnergyKinEM;
    const double sumMass = sumEnergyTot-sumEnergyKin;
    const double sumEnergyKinOther = sumEnergyKin-sumEnergyKinEM;
  
    if (ruVerbosity>=1) {
      cout << setprecision(3);
      cout << " summary-particles    (all) > "
	   << " sumTotEnergy=" << sumEnergyTot 
	   << " sumKinEnergy=" << sumEnergyKin
	   << " sumMass=" << sumMass << endl;
      cout << " summary-particles   (e.m.) > "
	   << " sumTotEnergyEM=" << sumEnergyTotEM 
	   << " sumKinEnergyEM=" << sumEnergyKinEM
	   << " sumMassEM=" << sumMassEM << endl;
      cout << " summary-particles (others) > "
	   << " sumEnergyKinOther=" << sumEnergyKinOther
	   << " totalE_leading=" << EtotLeading
	   << endl;
    }
  
    if (sumEnergyKinOther<=0.00001) {
      if (ruVerbosity>0)
	cout << " EMRatioResampling WARNING cannot change em-ratio with sumEnergyKinOther=" << sumEnergyKinOther << endl;
      return;
    }

    const double Etot = sumEnergyTot;

    const double minimalRatio = sumMassEM/Etot + 0.00001;
    const double maximalRatio = (sumEnergyTotEM+sumEnergyKinOther)/Etot - 0.00001;
  
    const double ratioEM = sumEnergyTotEM/sumEnergyTot;
    double newRatioEM = ratioEM*factor;
    if (newRatioEM<minimalRatio) {
      if (ruVerbosity>0)
	cout << " EMRatioResampling WARING rem_mod=" << newRatioEM << " smaller than min(rem)=" << minimalRatio
	     << endl;
      newRatioEM = minimalRatio;
    }
    if (newRatioEM>maximalRatio) {
      if (ruVerbosity>0)
	cout << " EMRatioResampling WARNING rem_mod=" << newRatioEM << " larger than max(rem)=" << maximalRatio
	     << endl;
      newRatioEM = maximalRatio;
    }
    if (ratioEM<=0 || newRatioEM<=0) {
      if (ruVerbosity>0)
	cout << " EMRatioResampling ERROR rem=" << ratioEM << ". skipping" << endl;
      return;
    }
    const double sumEnergyTotEMnew = sumEnergyTot * newRatioEM;
    const double sumEnergyKinEMnew = sumEnergyTotEMnew - sumMassEM;
    const double deltaE = sumEnergyTotEMnew - sumEnergyTotEM;
    const double deltaEkin = sumEnergyKinEMnew - sumEnergyKinEM;
    if (ruVerbosity>0) {
      cout << " ratioEM=" << ratioEM << " -> " << newRatioEM << endl; 
      cout << " sumEtot=" << Etot << " ETotEM= " << sumEnergyTotEM << " -> " << sumEnergyTotEMnew 
	   << " delta=" << deltaE
	   << endl; 
      cout << " sumEkin=" << sumEnergyKin << " ETotEM(kin)= " << sumEnergyKinEM << " -> " << sumEnergyKinEMnew 
	   << " delta(kin)=" << deltaEkin
	   << endl; 
    }
  
    const double modifyEM = sumEnergyKinEMnew/sumEnergyKinEM;
    const double modifyOther = (sumEnergyKinOther-deltaEkin)/sumEnergyKinOther;
  
    //const double sumEnergyKinScale = sumEnergyKin-sumEnergyKinEM;
    const double otherModify = modifyOther;
    // modify[id] = (sum[id].totEnergy - sumMass) / sumEkin;
  
    if (ruVerbosity>0) {
      cout << "            e.m. Ekin-factor=" << modifyEM << endl;
      cout << " Other than e.m. Ekin-factor=" << otherModify << endl;
    }
  
    for (int iPart=0; iPart<nPart; ++iPart) {
    
      ParticleBlockEntry& secondary = particlesList.GetEntry(iPart);
    
      if (!secondary.IsGood() ||             // do not take BAD particles
	  particlesList.IsLeading(iPart, ParticleBlock::eBoth)) {  // do not harm leading particle
	continue;
      }
    
      const int id = (int)abs((double)secondary.GetId());
      const double mass = secondary.GetMass();
      const double energy = secondary.GetEnergy();
      const double px = secondary.GetPx();
      const double py = secondary.GetPy();
      const double pz = secondary.GetPz();
      const double Pi = sqrt(px*px + py*py + pz*pz);
      
      const double modifyKin = ((id==CommonBlockParticleCONEX::eGamma ||
				 id==CommonBlockParticleCONEX::eElectron ||
				 id==CommonBlockParticleCONEX::ePi0) ?
				modifyEM :
				modifyOther);
      
      const double Ekin = energy - mass;
      const double energyNew = mass + Ekin*modifyKin;
      const double Pf = sqrt((energyNew+mass) * (energyNew-mass));
      
      
      double newPx = 0;
      double newPy = 0;
      double newPz = Pf;
      if (Pi<=0) {
	if (ruVerbosity>0)
	  cout << " WARNING: P_f=" << Pf 
	       << " P_i=" << Pi
	       << " Ekin_i=" << Ekin
	       << " factor=" << factor
	       << " E_f=" << energyNew
	       << " id=" << id
	       << " [-> new momentum in z-direction]"
	       << endl;
      } else {
	if (Pf<=0) {
	  if (ruVerbosity>0)
	    cout << " WARNING: P_f=" << Pf 
		 << " P_i=" << Pi
		 << " Ekin_i=" << Ekin
		 << " factor=" << factor
		 << " E_f=" << energyNew
		 << " id=" << id
		 << " [-> new momentum forecd to ZERO]"
		 << endl;
	  newPz = 0;
	} else {
	  const double momentumFactor = Pf/Pi;
	  newPx = px * momentumFactor;
	  newPy = py * momentumFactor;
	  newPz = pz * momentumFactor;
	}
      }
      
      secondary.SetPx(newPx);
      secondary.SetPy(newPy);
      secondary.SetPz(newPz);
      secondary.SetEnergy(energyNew);
    
      /*
      // ## DEBUG ############################################
      cout << setprecision(2);
      cout << " summary-modify>    id=" << setw(5) << id 
      << " N="     << setw(4) << sum[id].number+deltaN
      << " C="     << setw(4) << sum[id].charge+deltaCharge;
      if (ruVerbosity>3) 
      cout << " Etot="  << setw(14)<< sumEtot
      << " EtotIn=" << setw(14) << sum[id].totEnergy
      << " sumMass=" << setw(14) << sumMass;
      cout << " Ekin="  << setw(14)<< sumEkin
      << " dN="    << setw(4) << deltaN
      << " dC="    << setw(3) << deltaCharge
      << " dEtot=" << setw(14)<< deltaEtot
      << " dEkin=" << setw(14)<< deltaEkin
      << " factor="<< setw(10) << setprecision(8) << modify[id] 
      << endl;
      //}
      // ## end-DEBUG ############################################
      */
    }
  
  
    map<int, SummaryInfo> sumCheck = particlesList.MakeSummary(ParticleBlock::eTypes, ParticleBlock::eExcludeLeading);
    double sumEnergyTotCheck = 0;
    double sumEnergyKinCheck = 0;
    double sumEnergyTotEMCheck = 0;
    double sumEnergyKinEMCheck = 0;
    for (map<int, SummaryInfo>::iterator i = sumCheck.begin(); 
	 i != sumCheck.end(); ++i) {
      sumEnergyTotCheck += i->second.totEnergy;
      sumEnergyKinCheck += i->second.kinEnergy;
      switch(i->first) {
      case CommonBlockParticleCONEX::eGamma: 
      case CommonBlockParticleCONEX::eElectron:
      case CommonBlockParticleCONEX::ePositron:
      case CommonBlockParticleCONEX::ePi0:
	sumEnergyTotEMCheck += i->second.totEnergy;
	sumEnergyKinEMCheck += i->second.kinEnergy;
      }
      if (ruVerbosity>1) {
	cout << setprecision(3);
	cout << " summary-check>     id=" << setw(5) << i->first 
	     << " Etot=" << setw(14) << i->second.totEnergy
	     << " Ekin=" << setw(14) << i->second.kinEnergy
	     << " EkinMin=" << setw(14) << i->second.minKinEnergy
	     << " EkinMax=" << setw(14) << i->second.maxKinEnergy
	     << endl;
      }
    }
    
    if (ruVerbosity>0) {
      cout << setprecision(3);
      cout << " summary-check>      sumTotEnergy=" << sumEnergyTotCheck << " (delta=" << scientific << sumEnergyTotCheck-sumEnergyTot << fixed << ")"
	   << " sumKinEnergy=" << sumEnergyKinCheck << " (delta=" << sumEnergyKinCheck-sumEnergyKin << ")"
	   << endl;
      cout << " summary-check>      sumTotEnergyEM=" << sumEnergyTotEMCheck << " (delta=" << sumEnergyTotEMCheck-sumEnergyTotEM << ")"
	   << " sumKinEnergyEM=" << sumEnergyKinEMCheck << " (delta=" << sumEnergyKinEMCheck-sumEnergyKinEM << ")"
	   << endl;
      // particlesList.FindLeadingParticle();
    }
  
  } // end of EMRatioResampling


  
  /* ChargeRatioResampling
   * *********************
   */
  void 
  ChargeRatioResampling(const double factor, const int projId, ParticleBlock& particlesList) 
  {    
    if (factor==1)
      return;
    
    const bool projectilePositive = (CommonBlockParticleCONEX::GetCharge(projId) > 0);
    const bool projectileNegative = (CommonBlockParticleCONEX::GetCharge(projId) < 0);
    
    const int nPart = particlesList.Size();
    
    // get leading particles in forward/backward direction
    vector<int> leadingParticleIndices;
    particlesList.IdentifyLeadingParticle(ParticleBlock::eForward);
    leadingParticleIndices.push_back(particlesList.GetLeadingParticleIndex(ParticleBlock::eForward));
    if (particlesList.IsCenterOfMassSystem()) {
      particlesList.IdentifyLeadingParticle(ParticleBlock::eBackward);
      leadingParticleIndices.push_back(particlesList.GetLeadingParticleIndex(ParticleBlock::eBackward));
    }
    
    // get overview of the particle block
    double sumEnergyTot = 0;
    int sumNumberPi0 = 0;
    int sumNumberPiM = 0;
    int sumNumberPiP = 0;
    map<int, SummaryInfo> sum = particlesList.MakeSummary(ParticleBlock::eExact, ParticleBlock::eExcludeLeading);
    for (map<int, SummaryInfo>::iterator i = sum.begin(); i != sum.end(); ++i) {
      sumEnergyTot += i->second.totEnergy;
      switch(i->first) {
      case CommonBlockParticleCONEX::ePiP:
	sumNumberPiP += i->second.number;
	break;
      case CommonBlockParticleCONEX::ePiM:
	sumNumberPiM += i->second.number;
	break;
      case CommonBlockParticleCONEX::ePi0:
	sumNumberPi0 += i->second.number;
	break;
      }
      if (ruVerbosity>=1) {
	cout << setprecision(3);
	cout << " summary-particles> id=" << setw(5) << i->first 
	     << " Etot=" << setw(14) << i->second.totEnergy
	     << " Ekin=" << setw(14) << i->second.kinEnergy
	     << " EkinMin=" << setw(14) << i->second.minKinEnergy
	     << " EkinMax=" << setw(14) << i->second.maxKinEnergy
	     << endl;
      }
    }
    
    double EtotLeading = 0;
    for (unsigned int iLeading=0; iLeading != leadingParticleIndices.size(); ++iLeading) {
      const ParticleBlockEntry& leading = particlesList.GetEntry(leadingParticleIndices[iLeading]);
      EtotLeading += leading.GetEnergy();
      if (ruVerbosity>=1) {
	cout << " leading particle (" << (iLeading==0 ? "FORWARD": "BACKWARD" )
	     << ") " << leading.GetName() << "(" << leading.GetId() << ") "
	     << " E_tot=" << leading.GetEnergy()
	     << endl;
      }
    }
    
    int maxPionTransfer = sumNumberPi0;  // factor<1  :  pi0  ->  pi+-
    if (factor>1) {                      // factor>1  :  pi+- ->  pi0
      maxPionTransfer = sumNumberPiM + sumNumberPiP;               // neutral projectiles
      if (projectileNegative) maxPionTransfer = sumNumberPiP;      // negative projectiles
      else if (projectilePositive) maxPionTransfer = sumNumberPiM; // positive projectiles
    }
    
    const int totalNumberOfPions = sumNumberPi0 + sumNumberPiM + sumNumberPiP;
    
    if (totalNumberOfPions == 0) {
      if (ruVerbosity>20)
	cout << " ChargeRatioResampling> No pions in list of secondaries after hadronic interaction! Skipping" << endl;
      return;      
    }
    
    const double chargeRatio = double(sumNumberPi0) / totalNumberOfPions;
    double chargeRatioNew = chargeRatio * factor;
    if (factor>1) {
      const double maximalChargeRatio = double(sumNumberPi0 + maxPionTransfer) / totalNumberOfPions;
      if (chargeRatioNew>maximalChargeRatio) {
	if (ruVerbosity>20)
	  cout << " ChargeRatioResampling WARNING r_ch_mod=" << chargeRatioNew 
	       << " larger than max(r_ch)=" << maximalChargeRatio
	       << endl;
	chargeRatioNew = maximalChargeRatio;
      }
    }
    
    if (maxPionTransfer==0) {
      if (ruVerbosity>50)
	cout << " ChargeRatioResampling> No pions to transfer for charge-ratio resampling! Skipping" << endl;
      return;
    }
    
    const double probConvert = (factor < 1 ? 
				1 - factor :
				(factor-1) * sumNumberPi0 / maxPionTransfer);
    
    if (ruVerbosity>0) {
      cout << " ChargeRatioResampling> "
	   << " projectile pos " << projectilePositive << " neg " << projectileNegative << " id " << projId      
	   << " npi0=" << sumNumberPi0
	   << " npi-=" << sumNumberPiM << " npi+=" << sumNumberPiP
	   << " Npitot=" << totalNumberOfPions
	   << endl;
      cout << "                        " 
	   << " c=" << chargeRatio << " f=" << factor << " maxNpi=" << maxPionTransfer << " " 
	   << " p=" << probConvert 
	   << " p*maxNpi=" << maxPionTransfer*probConvert
	   << endl;
    }
    
    for (int iPart=0; iPart<nPart; ++iPart) {
      
      ParticleBlockEntry& secondary = particlesList.GetEntry(iPart);      
      const int id = secondary.GetId();
      
      if (!secondary.IsGood() ||             // do not take BAD particles
	  particlesList.IsLeading(iPart, ParticleBlock::eBoth))
	continue;    
      
      if (factor<1 && id != CommonBlockParticleCONEX::ePi0) // for factor<1 only pi0 are relevant
	continue;

      if (factor>1 && 
	  ((projectilePositive && id != CommonBlockParticleCONEX::ePiM) || // for factor>1 and pos. proj. only pi-
	   (projectileNegative && id != CommonBlockParticleCONEX::ePiP) || // for factor>1 and neg. proj. only pi+
	   (id != CommonBlockParticleCONEX::ePiP && id != CommonBlockParticleCONEX::ePiM)))           // for factor>1 and neut. proj. only pi-/pi+
	continue;
      
      if (gRandom->Rndm() > probConvert)
	continue;
      
      //const double mass = secondary.GetMass();
      const double energy = secondary.GetEnergy();
      const double px = secondary.GetPx();
      const double py = secondary.GetPy();
      const double pz = secondary.GetPz();
      const double Pi = sqrt(px*px + py*py + pz*pz);
      
      // for factor>1
      int newId = CommonBlockParticleCONEX::ePi0;     // pi+ or pi-    ->     pi0
      double massNew = 0.1349766; // GeV
      if (factor<1) {
	massNew = 0.13957018;     // GeV
	if (projectileNegative) 
	  newId = CommonBlockParticleCONEX::ePiP;     // pi0   ->   pi+
	else if (projectilePositive)
	  newId = CommonBlockParticleCONEX::ePiM;     // pi0   ->   pi-
	else	                                      // pi0   ->   pi+   or   pi-
	  newId = (gRandom->Rndm()>0.5 ? CommonBlockParticleCONEX::ePiM : CommonBlockParticleCONEX::ePiP);
      }
      
      // const double deltaMass = massNew - mass;
      const double EkinNew = energy - massNew;
      if (EkinNew<=0){
	if (ruVerbosity>20) {
	  cout << " ChargeRatioResampling> cannot switch particle id, since EkinNew=" << EkinNew << endl;
	}
	continue;
      }
      const double Pf = sqrt((energy+massNew) * (energy-massNew));
      
      double newPx = 0;
      double newPy = 0;
      double newPz = Pf;
      if (Pi<=0) {
	if (ruVerbosity>0)
	  cout << " WARNING: P_f=" << Pf 
	       << " P_i=" << Pi
	       << " factor=" << factor
	       << " id=" << id
	       << " [-> new momentum in z-direction]"
	       << endl;
      } else {
	if (Pf<=0) {
	  if (ruVerbosity>0)
	    cout << " WARNING: P_f=" << Pf 
		 << " P_i=" << Pi
		 << " factor=" << factor
		 << " id=" << id
		 << " [-> new momentum forecd to ZERO]"
		 << endl;
	  newPz = 0;
	} else {
	  const double momentumFactor = Pf/Pi;
	  newPx = px * momentumFactor;
	  newPy = py * momentumFactor;
	  newPz = pz * momentumFactor;
	}
      }

      secondary.SetId(newId);
      secondary.SetMass(massNew);
      secondary.SetPx(newPx);
      secondary.SetPy(newPy);
      secondary.SetPz(newPz);
      
      if (ruVerbosity>50) 
	cout << " Switch pion: index=" << iPart << " id=" << id << " -> " << newId << endl;
    }

    // *********************************************
    double sumEnergyTotCheck = 0;
    int sumNumberPi0Check = 0;
    int sumNumberPiMCheck = 0;
    int sumNumberPiPCheck = 0;
    map<int, SummaryInfo> sumCheck = particlesList.MakeSummary(ParticleBlock::eExact, ParticleBlock::eExcludeLeading);
    for (map<int, SummaryInfo>::iterator i = sumCheck.begin(); i != sumCheck.end(); ++i) {
      sumEnergyTotCheck += i->second.totEnergy;
      switch(i->first) {
      case CommonBlockParticleCONEX::ePiP:
	sumNumberPiPCheck += i->second.number;
	break;
      case CommonBlockParticleCONEX::ePiM:
	sumNumberPiMCheck += i->second.number;
	break;
      case CommonBlockParticleCONEX::ePi0:
	sumNumberPi0Check += i->second.number;
	break;
      }
      if (ruVerbosity>=1) {
	cout << setprecision(3);
	cout << " summary-check> id=" << setw(5) << i->first 
	     << " Etot=" << setw(14) << i->second.totEnergy
	     << " Ekin=" << setw(14) << i->second.kinEnergy
	     << " EkinMin=" << setw(14) << i->second.minKinEnergy
	     << " EkinMax=" << setw(14) << i->second.maxKinEnergy
	     << endl;
      }
    }
    
    if (ruVerbosity>0) {
      cout << setprecision(3);
      cout << " summary-check>         #pi0=" << sumNumberPi0Check << " (delta=" << sumNumberPi0Check-sumNumberPi0 << ")"
	   << " #pi-=" << sumNumberPiMCheck << " (delta=" << sumNumberPiMCheck-sumNumberPiM << ")"
	   << " #pi+=" << sumNumberPiPCheck << " (delta=" << sumNumberPiPCheck-sumNumberPiP << ")"
	   << " total=" << sumNumberPiPCheck+sumNumberPiMCheck+sumNumberPi0Check 
	   << " (delta=" << sumNumberPiPCheck+sumNumberPiMCheck+sumNumberPi0Check-totalNumberOfPions << ")"
	   << endl;    
      /*
	cout << " summary-check>      sumTotEnergy=" << sumEnergyTotCheck << " (delta=" << scientific << sumEnergyTotCheck-sumEnergyTot << fixed << ")" 
	<< endl;
      */
      // particlesList.FindLeadingParticle();
    }

  } // end of ChargeRatioResampling

  
  void 
  Pi0SpectrumResampling(const double factor, resample::ParticleBlock& particlesList)
  {    
    if (factor == 1)
      return;
    
    const double slopeFactor = factor - 1;
    
    const int nPart = particlesList.Size();
    
    // get leading particles in forward/backward direction
    vector<int> leadingParticleIndices;
    particlesList.IdentifyLeadingParticle(ParticleBlock::eForward);
    leadingParticleIndices.push_back(particlesList.GetLeadingParticleIndex(ParticleBlock::eForward));
    if (particlesList.IsCenterOfMassSystem()) {
      particlesList.IdentifyLeadingParticle(ParticleBlock::eBackward);
      leadingParticleIndices.push_back(particlesList.GetLeadingParticleIndex(ParticleBlock::eBackward));
    }
    
    // get overview of the particle block
    double sumEnergyTot = 0;
    double sumEnergyPi0 = 0;
    int sumNPi0 = 0;
    map<int, SummaryInfo> sum = particlesList.MakeSummary(ParticleBlock::eExact, ParticleBlock::eExcludeLeading);
    for (map<int, SummaryInfo>::iterator i = sum.begin(); i != sum.end(); ++i) {
      sumEnergyTot += i->second.totEnergy;
      switch(i->first) {
      case CommonBlockParticleCONEX::ePi0:
	sumEnergyPi0 += i->second.kinEnergy;
	sumNPi0 += i->second.number;
	break;
      }
      if (ruVerbosity>1) {
	cout << setprecision(3);
	cout << " summary-particles> id=" << setw(5) << i->first 
	     << " Etot=" << setw(14) << i->second.totEnergy
	     << " Ekin=" << setw(14) << i->second.kinEnergy
	     << " EkinMin=" << setw(14) << i->second.minKinEnergy
	     << " EkinMax=" << setw(14) << i->second.maxKinEnergy
	     << endl;
      }
    }
    
    
    if (ruVerbosity>0) {
      cout << " Pi0SpectrumResampling> "
	   << " sumEnergyPi0=" << sumEnergyPi0
	   << " sumNPi0=" << sumNPi0
	   << " avgE=";
      if (sumNPi0>0)
	cout << sumEnergyPi0/sumNPi0;
      else
	cout << "------";
      cout << endl;
    }
    
    if (sumNPi0<1)
      return;
    
    const double avgEnergyPi0 = sumEnergyPi0 / sumNPi0; // use as reference energy
    
    // -- for debug --
    int dbgN = 0;
    double dbgEi = 0;
    double dbgEf = 0;
    double dbgE2i = 0;
    double dbgE2f = 0;
    // --
    
    double sumPx = 0;
    double sumPy = 0;
    double sumPz = 0;
    double sumEnergyPi0Mod = 0;
    for (int iPart=0; iPart<nPart; ++iPart) {
      
      ParticleBlockEntry& secondary = particlesList.GetEntry(iPart);      
      const int id = secondary.GetId();
      
      if (!secondary.IsGood() ||             // do not take BAD particles
	  particlesList.IsLeading(iPart, ParticleBlock::eBoth))
	continue;    
      
      if (id != CommonBlockParticleCONEX::ePi0) 
	continue;
      
      sumPx += secondary.GetPx();
      sumPy += secondary.GetPy();
      sumPz += secondary.GetPz();
      
      const double mass = secondary.GetMass();
      const double energy = secondary.GetEnergy();
      
      const double Ekin = energy - mass;
      const double slopeWeight = pow(Ekin, -slopeFactor);
      const double EkinNew = Ekin * slopeWeight;
      const double energyNew = mass + EkinNew;
      
      secondary.SetEnergy(energyNew);
      
      sumEnergyPi0Mod += EkinNew;
      
      // -- for dbg --
      ++dbgN;
      dbgEi += Ekin;
      dbgEf += EkinNew;
      dbgE2i += Ekin*Ekin;
      dbgE2f += EkinNew*EkinNew;

    } //  end loop particles
    
    const double sumP = sqrt(sumPx*sumPx + sumPy*sumPy + sumPz*sumPz);
    if (sumP<=0) {
      cout << " Pi0SpectrumResampling: no momentum to work with" << endl;
      return;
    }
    sumPx /= sumP;
    sumPy /= sumP;
    sumPz /= sumP;
    
    
    
    const double scaleEkin = sumEnergyPi0 / sumEnergyPi0Mod;
    /*
    dbgEi /= dbgN;
    dbgEf /= dbgN;
    dbgE2i = sqrt(dbgE2i/dbgN - dbgEi*dbgEi);
    dbgE2f = sqrt(dbgE2f/dbgN - dbgEf*dbgEf);

    cout << " sumEnergyPi0=" << sumEnergyPi0 
	 << " sumEnergyPi0Mod=" << sumEnergyPi0Mod
	 << " scaleEkin=" << scaleEkin << endl;
    cout << " dbgE2i=" << dbgE2i
	 << " dbgE2f=" << dbgE2f
	 << endl;
    */

    // loop again, to re-scale pi0 energies
    double sumPxnew = 0;
    double sumPynew = 0;
    double sumPznew = 0;
    for (int iPart=0; iPart<nPart; ++iPart) {
      
      ParticleBlockEntry& secondary = particlesList.GetEntry(iPart);      
      const int id = secondary.GetId();
      
      if (!secondary.IsGood() ||             // do not take BAD particles
	  particlesList.IsLeading(iPart, ParticleBlock::eBoth))
	continue;    
      
      if (id != CommonBlockParticleCONEX::ePi0) 
	continue;
      
      const double mass = secondary.GetMass();
      const double energyNew = secondary.GetEnergy();

      const double px = secondary.GetPx();
      const double py = secondary.GetPy();
      const double pz = secondary.GetPz();
      const double Pi = sqrt(px*px + py*py + pz*pz);
      
      const double EkinScaled = (energyNew - mass) * scaleEkin;
      const double energyScaled = mass + EkinScaled;
      const double Pf = sqrt((energyScaled+mass) * (energyScaled-mass));
      
      double newPx = 0;
      double newPy = 0;
      double newPz = Pf;
      if (Pi<=0) {
	if (ruVerbosity>0)
	  cout << " WARNING: P_f=" << Pf 
	       << " P_i=" << Pi
	       << " factor=" << factor
	       << " id=" << id
	       << " [-> new momentum in z-direction]"
	       << endl;
      } else {
	if (Pf<=0) {
	  if (ruVerbosity>0)
	    cout << " WARNING: P_f=" << Pf 
		 << " P_i=" << Pi
		 << " factor=" << factor
		 << " id=" << id
		 << " [-> new momentum forecd to ZERO]"
		 << endl;
	  newPz = 0;
	} else {
	  const double momentumFactor = Pf/Pi;
	  newPx = px * momentumFactor;
	  newPy = py * momentumFactor;
	  newPz = pz * momentumFactor;
	}
      }
      secondary.SetEnergy(energyScaled);
      secondary.SetPx(newPx);
      secondary.SetPy(newPy);
      secondary.SetPz(newPz);
      sumPxnew += newPx;
      sumPynew += newPy;
      sumPznew += newPz;
    } //  end loop particles
    
    const double sumPnew = sqrt(sumPxnew*sumPxnew + sumPynew*sumPynew + sumPznew*sumPznew);
    if (sumPnew<=0) {
      cout << " Pi0SpectrumResampling: no new momentum to work with" << endl;
      return;
    }
    sumPxnew /= sumPnew;
    sumPynew /= sumPnew;
    sumPznew /= sumPnew;
    
    
    //    const double alpha = sumPx / sumPnew; 
    
    // final loop to re-align momenta
    for (int iPart=0; iPart<nPart; ++iPart) {
      
      ParticleBlockEntry& secondary = particlesList.GetEntry(iPart);      
      const int id = secondary.GetId();
      
      if (!secondary.IsGood() ||             // do not take BAD particles
	  particlesList.IsLeading(iPart, ParticleBlock::eBoth))
	continue;    
      
      if (id != CommonBlockParticleCONEX::ePi0) 
	continue;
    }
 
    
    
    
    // *********************************************
    double sumEnergyTotCheck = 0;
    double sumEnergyPi0Check = 0;
    int sumNPi0Check = 0;
    map<int, SummaryInfo> sumCheck = particlesList.MakeSummary(ParticleBlock::eExact, 
							       ParticleBlock::eExcludeLeading);
    for (map<int, SummaryInfo>::iterator i = sumCheck.begin(); i != sumCheck.end(); ++i) {
      sumEnergyTotCheck += i->second.totEnergy;
      switch(i->first) {
      case CommonBlockParticleCONEX::ePi0:
	sumEnergyPi0Check += i->second.kinEnergy;
	sumNPi0Check += i->second.number;
	break;
      }
      if (ruVerbosity>=1) {
	cout << setprecision(3);
	cout << " summary-check> id=" << setw(5) << i->first 
	     << " Etot=" << setw(14) << i->second.totEnergy
	     << " Ekin=" << setw(14) << i->second.kinEnergy
	     << " EkinMin=" << setw(14) << i->second.minKinEnergy
	     << " EkinMax=" << setw(14) << i->second.maxKinEnergy
	     << endl;
      }
    }
    
    const double avgEnergyPi0Check = sumEnergyPi0Check/sumNPi0Check;

    if (ruVerbosity>0) {
      cout << setprecision(3);
      cout << " summary-check>       Ekin(pi0)=" 
	   << sumEnergyPi0Check
	   << " (delta=" << sumEnergyPi0Check-sumEnergyPi0 << ")"
	   << " deltaN=" << sumNPi0Check-sumNPi0 << " (must be zero!)"
	   << " avgEPi0=" << avgEnergyPi0Check << " (delta="<< avgEnergyPi0Check-avgEnergyPi0 <<")" 
	   << endl;    
    }

  } // end of Pi0SpectrumResampling



  /* Rho0Resampling
   * *********************
   */
  void 
  Rho0Resampling(const double factor, const int projId, ParticleBlock& particlesList) 
  {    
    if (factor==1) {
      return;
    }
    
    if (ruVerbosity>50) {
      cout << " Rho0Resampling> projId=" << projId << endl;;
    }
    
    // only for charged pion primaries
    if (projId != CommonBlockParticleCONEX::ePiP &&
	projId != CommonBlockParticleCONEX::ePiM) {
      return;
    }
    
    // get leading particles in forward/backward direction
    vector<int> leadingParticleIndices;
    particlesList.IdentifyLeadingParticle(ParticleBlock::eForward);
    leadingParticleIndices.push_back(particlesList.GetLeadingParticleIndex(ParticleBlock::eForward));
    if (particlesList.IsCenterOfMassSystem()) {
      particlesList.IdentifyLeadingParticle(ParticleBlock::eBackward);
      leadingParticleIndices.push_back(particlesList.GetLeadingParticleIndex(ParticleBlock::eBackward));
    }
    
    // get overview of the particle block
    double sumEnergyTot = 0;
    map<int, SummaryInfo> sum = particlesList.MakeSummary(ParticleBlock::eExact, ParticleBlock::eExcludeLeading);
    for (map<int, SummaryInfo>::iterator i = sum.begin(); i != sum.end(); ++i) {
      sumEnergyTot += i->second.totEnergy;
      if (ruVerbosity>=1) {
	cout << setprecision(3);
	cout << " summary-particles> id=" << setw(5) << i->first 
	     << " Etot=" << setw(14) << i->second.totEnergy
	     << " Ekin=" << setw(14) << i->second.kinEnergy
	     << " EkinMin=" << setw(14) << i->second.minKinEnergy
	     << " EkinMax=" << setw(14) << i->second.maxKinEnergy
	     << endl;
      }
    }
    
    
    double probConvert = 0;
    if (factor>1) {                      // factor>1  :  pi0 ->  rho0      
      probConvert = factor-1.0;
    }
    else {                      // factor<1  :  pi0 <-  rho0
      probConvert = 1.0-factor;
    }

    
    if (ruVerbosity>20) {
      cout << " Rho0Resampling> probConvert=" << probConvert << endl;
    }
    

    for (unsigned int iLeading=0; iLeading != leadingParticleIndices.size(); ++iLeading) {

      ParticleBlockEntry& leading = particlesList.GetEntry(leadingParticleIndices[iLeading]);
      const int id = (int)abs((double)leading.GetId());
      
      if (ruVerbosity>20) {
	cout << " Rho0Resampling> leading particle: ";
	leading.Dump();
      }
        
      double massNew = 0;
      int newId = 0;
    
      if (id == CommonBlockParticleCONEX::ePi0 && factor>1) {

	if (gRandom->Rndm() > probConvert)
	  continue; 
	

	if (ruVerbosity>20) {
	  cout << " Rho0Resampling> Pi0 -> Rho0" << endl;
	}

	// switch pi0 to rho0

	newId = CommonBlockParticleCONEX::eRho0;
	massNew = 0.77549; // GeV

      } // finished with pi0 -> rho0
      
      else if (id == CommonBlockParticleCONEX::eRho0 && factor<1) {
	
	if (gRandom->Rndm() > probConvert)
	  continue;

	// switch rho0 to pi0

	if (ruVerbosity>20) {
	  cout << " Rho0Resampling> Rho0 -> Pi0" << endl;
	}

	newId = CommonBlockParticleCONEX::ePi0;
	massNew = 0.1349766; // GeV

      } else {
	continue;
      }
	
      const double energy = leading.GetEnergy();
      const double px = leading.GetPx();
      const double py = leading.GetPy();
      const double pz = leading.GetPz();
      const double Pi = sqrt(px*px + py*py + pz*pz);
      
      
      const double EkinNew = energy - massNew;
      if (EkinNew<=0){
	if (ruVerbosity>20) {
	  cout << " ChargeRatioResampling> cannot switch particle id, since EkinNew=" << EkinNew << endl;
	}
	continue;
      }
      const double Pf = sqrt((energy+massNew) * (energy-massNew));
      
      double newPx = 0;
      double newPy = 0;
      double newPz = Pf;
      if (Pi<=0) {
	if (ruVerbosity>0)
	  cout << " WARNING: P_f=" << Pf 
	       << " P_i=" << Pi
	       << " factor=" << factor
	       << " id=" << id
	       << " [-> new momentum in z-direction]"
	       << endl;
      } else {
	if (Pf<=0) {
	  if (ruVerbosity>0)
	    cout << " WARNING: P_f=" << Pf 
		 << " P_i=" << Pi
		 << " factor=" << factor
		 << " id=" << id
		 << " [-> new momentum forecd to ZERO]"
		 << endl;
	  newPz = 0;
	} else {
	  const double momentumFactor = Pf/Pi;
	  newPx = px * momentumFactor;
	  newPy = py * momentumFactor;
	  newPz = pz * momentumFactor;
	}
      }
      
      leading.SetId(newId);
      leading.SetMass(massNew);
      leading.SetPx(newPx);
      leading.SetPy(newPy);
      leading.SetPz(newPz);

      if (ruVerbosity>50) 
	cout << " Switch pi0/rho0: index=" << iLeading << " id=" << id << " -> " << newId << endl;
      
    }
      
    // *********************************************
    double sumEnergyTotCheck = 0;
    map<int, SummaryInfo> sumCheck = particlesList.MakeSummary(ParticleBlock::eExact, ParticleBlock::eExcludeLeading);
    for (map<int, SummaryInfo>::iterator i = sumCheck.begin(); i != sumCheck.end(); ++i) {
      sumEnergyTotCheck += i->second.totEnergy;
      if (ruVerbosity>=1) {
	cout << setprecision(3);
	cout << " summary-check> id=" << setw(5) << i->first 
	     << " Etot=" << setw(14) << i->second.totEnergy
	     << " Ekin=" << setw(14) << i->second.kinEnergy
	     << " EkinMin=" << setw(14) << i->second.minKinEnergy
	     << " EkinMax=" << setw(14) << i->second.maxKinEnergy
	     << endl;
      }
    }
    
    if (ruVerbosity>0) {
      //      cout << setprecision(3);
      /*cout << " summary-check>         #pi0=" << sumNumberPi0Check << " (delta=" << sumNumberPi0Check-sumNumberPi0 << ")"
	   << " #pi-=" << sumNumberPiMCheck << " (delta=" << sumNumberPiMCheck-sumNumberPiM << ")"
	   << " #pi+=" << sumNumberPiPCheck << " (delta=" << sumNumberPiPCheck-sumNumberPiP << ")"
	   << " total=" << sumNumberPiPCheck+sumNumberPiMCheck+sumNumberPi0Check 
	   << " (delta=" << sumNumberPiPCheck+sumNumberPiMCheck+sumNumberPi0Check-totalNumberOfPions << ")"
	   << endl;    
      */
    }

  } // end of Rho0Resampling


  
  /* BaryonProductionResampling
   * *********************
   */
  void 
  BaryonProductionResampling(const double factor, const int projId, ParticleBlock& particlesList) 
  {    
    if (factor==1) {
      return;
    }
        
    // get leading particles in forward/backward direction
    vector<int> leadingParticleIndices;
    particlesList.IdentifyLeadingParticle(ParticleBlock::eForward);
    leadingParticleIndices.push_back(particlesList.GetLeadingParticleIndex(ParticleBlock::eForward));
    if (particlesList.IsCenterOfMassSystem()) {
      particlesList.IdentifyLeadingParticle(ParticleBlock::eBackward);
      leadingParticleIndices.push_back(particlesList.GetLeadingParticleIndex(ParticleBlock::eBackward));
    }
        
    // get overview of the particle block
    double sumEnergyTot = 0;
    map<int, SummaryInfo> sum = particlesList.MakeSummary(ParticleBlock::eExact, ParticleBlock::eExcludeLeading);
    for (map<int, SummaryInfo>::iterator i = sum.begin(); i != sum.end(); ++i) {
      sumEnergyTot += i->second.totEnergy;
      if (ruVerbosity>=1) {
	cout << setprecision(3);
	cout << " summary-particles> id=" << setw(5) << i->first 
	     << " Etot=" << setw(14) << i->second.totEnergy
	     << " Ekin=" << setw(14) << i->second.kinEnergy
	     << " EkinMin=" << setw(14) << i->second.minKinEnergy
	     << " EkinMax=" << setw(14) << i->second.maxKinEnergy
	     << endl;
      }
    }

    const int nPart = particlesList.Size();
    
    // count pi+/pi- pairs, sort by energy
    vector<int> allPiM;
    vector<int> allPiP;
    vector<int> allBaryons;
    vector<int> allAntiBaryons;
    for (int iPart=0; iPart<nPart; ++iPart) {
      
      ParticleBlockEntry& secondary = particlesList.GetEntry(iPart);      
      const int id = secondary.GetId();
      
      if (!secondary.IsGood() ||             // do not take BAD particles
	  particlesList.IsLeading(iPart, ParticleBlock::eBoth))
	continue;    
      
      switch (id) {
      case CommonBlockParticleCONEX::ePiP:
	allPiP.push_back(iPart);
	break;
	
      case CommonBlockParticleCONEX::ePiM:
	allPiM.push_back(iPart);
	break;
	
      case CommonBlockParticleCONEX::eProton:
      case CommonBlockParticleCONEX::eNeutron:
      case CommonBlockParticleCONEX::eLambda:
	allBaryons.push_back(iPart);
	break;

      case CommonBlockParticleCONEX::eProtonBar:
      case CommonBlockParticleCONEX::eNeutronBar:
      case CommonBlockParticleCONEX::eLambdaBar:
	allAntiBaryons.push_back(iPart);
	break;
      }
    }

    if (ruVerbosity>=10) {
      cout << " BaryonProductionResampling> npi+=" << allPiP.size() << " npi-=" << allPiM.size()
	   << " nb=" << allBaryons.size() << " nbbar=" << allAntiBaryons.size() 
	   << endl;
    }
    
    
    // calculate pi+pi- pair energy
    //vector<double> pairEnergy;
    int maxLambda=-1, maxNeutron=-1, maxProton=-1;
    for (unsigned int iPair=0; iPair<min(allPiP.size(), allPiM.size()); ++iPair) {
      const double e1 = particlesList.GetEntry(allPiM[iPair]).GetEnergy();
      const double e2 = particlesList.GetEntry(allPiP[iPair]).GetEnergy();
      //pairEnergy.push_back(e1+e2);
      if (e1>CommonBlockParticleCONEX::fgMassProton &&
	  e2>CommonBlockParticleCONEX::fgMassProton) 
	maxProton = iPair;
      if (e1>CommonBlockParticleCONEX::fgMassNeutron &&
	  e2>CommonBlockParticleCONEX::fgMassNeutron)
	maxNeutron = iPair;
      // no lambda, sibyll seems to not like it....
      /*if (e1>CommonBlockParticleCONEX::fgMassLambda &&
	  e2>CommonBlockParticleCONEX::fgMassLambda)
	  maxLambda = iPair;*/
    }
    
    if (ruVerbosity>=10) {
      cout << " BaryonProductionResampling> maxProton=" << maxProton << " maxNeutron=" << maxNeutron
	   << " maxLambda=" << maxLambda
	   << endl;
    }

    if (maxProton<0)
      return;
    
    const int nBaryons = allBaryons.size() + allAntiBaryons.size(); 
    const int nMesons = (maxProton + 1) * 2;
    
    const double ratio = double(nBaryons) / double(nMesons);
    const double ratioNew = ratio * factor;
    
    const double nDiff = (ratioNew * double(nMesons) - double(nBaryons)) / (1.0 - ratioNew);
    
    const int nDiffPair = int(nDiff/2 + 0.5);

    const int nDiffPairRnd = gRandom->Poisson(abs(nDiff)/2);
    
    if (ruVerbosity>=10) {
      cout << " BaryonProductionResampling> ratio=" << ratio << " ratioNew=" << ratioNew
	   << " nDiff=" << nDiff << " nDiffPair=" << nDiffPair
	   << " nDiffPairRnd=" << nDiffPairRnd
	   << endl;
	}

    if (nDiffPair>0) { // meson -> baryon
      
      set<int> usedMesons; // PiP, usedPiM, usedBaryons, usedAntiBaryons; 
      
      for (int iDiffPair=0; iDiffPair<nDiffPairRnd; ++iDiffPair) {
	
	CommonBlockParticleCONEX::EParticleId targetId = CommonBlockParticleCONEX::eProton;
	
	const int nTypes = maxProton + max(maxNeutron, 0) + max(maxLambda, 0);
	int index = 0; // index of meson pair
	int iSel = 0; // selected meson pair
	bool mesonPairGood = false;
	int notFoundCounter = 0;
	
	do { // select pair of mesons
	  
	  const int rndmType = gRandom->Rndm() * nTypes; // select target baryon	
	  
	  if (maxLambda>=0 && rndmType<maxLambda) {
	    
	    // -> lambda
	    
	    targetId = CommonBlockParticleCONEX::eLambda;
	    index = rndmType;
	    
	  } else if (maxNeutron>=0 && rndmType<maxNeutron+maxLambda) { // includes Lambdas
	    
	    // -> neutron
	    targetId = CommonBlockParticleCONEX::eNeutron;
	    index = rndmType - max(maxLambda, 0);
	    
	  } else { // includes everything
	    
	    // -> proton
	    targetId = CommonBlockParticleCONEX::eProton;
	    index = rndmType - max(maxLambda, 0) - max(maxNeutron, 0);
	    
	  } // 
	  
	  const double massTarget = CommonBlockParticleCONEX::GetMass(targetId);
	  
	  int search = 0;
	  bool found = false;
	  for (; iSel<maxProton; ++iSel) {
	    
	    const double e1 = particlesList.GetEntry(allPiM[iSel]).GetEnergy();
	    const double e2 = particlesList.GetEntry(allPiP[iSel]).GetEnergy();
	    if (e1>massTarget &&
		e2>massTarget) {
	    
	      if (usedMesons.find(iSel)==usedMesons.end()) {
		if (search >= index) {
		  usedMesons.insert(iSel);
		  found = true;
		  break;
		}
	      }
	      
	      ++search;
	    } // end mass if
	    
	  } // end for loop 
	  
	  if (found) {
	    mesonPairGood = true;
	    break; // ok, break this loop
	  }
	  // not found...
	  if (++notFoundCounter>maxProton) {
	    break; // stop here
	  }
 	  
	} while (true); // found meson-pair to change 
	
	if (!mesonPairGood)
	  continue;
	

	// CHANGE PIDs HERE
	
	ParticleBlockEntry& p1 = particlesList.GetEntry(allPiM[iSel]);
	ParticleBlockEntry& p2 = particlesList.GetEntry(allPiP[iSel]);
	
	if (ruVerbosity) {
	  cout << " BaryonProductionResampling> " << endl;
	  cout << " meson1: index=" << allPiM[iSel] << " "; p1.Dump();
	  cout << " meson2: index=" << allPiP[iSel] << " "; p2.Dump();
	}
	
	if (!resample::SetNewPId(targetId, p1))
	  continue;
	
	if (!resample::SetNewPId(CommonBlockParticleCONEX::GetIdOfAntiParticle(targetId), p2))
	  continue;

	if (ruVerbosity) {
	  cout << " BaryonProductionResampling> " << endl;
	  cout << " baryon1: "; p1.Dump();
	  cout << " baryon2: "; p2.Dump();
	}

      } // end loop change meson pairs

      

    } else { // baryon -> meson


      if (ruVerbosity>50) {	
	cout << " BaryonProductionResampling> baryons: " << endl;
	for (unsigned int i=0; i<allBaryons.size(); ++i) {
	  particlesList.GetEntry(allBaryons[i]).Dump();
	}
	cout << " BaryonProductionResampling> anti-baryons: " << endl;
	for (unsigned int i=0; i<allAntiBaryons.size(); ++i) {
	  particlesList.GetEntry(allAntiBaryons[i]).Dump();
	}
      }


      set<int> usedBaryons, usedAntiBaryons;
      
      for (int iDiffPair=0; iDiffPair<nDiffPairRnd; ++iDiffPair) {
	
	// pick one randomly and search for counterpart
	
	int pick1 = 0;
	int typeOne = 0;
	int countTry = 0;
	do {
	  pick1 = gRandom->Rndm() * allBaryons.size();
	  typeOne = particlesList.GetEntry(allBaryons[pick1]).GetId();
	  
	  if (usedBaryons.find(pick1) == usedBaryons.end()) {
	    break;
	  }
	  
	} while (countTry++<nBaryons*2);
	
	if (countTry >= nBaryons*2-1)
	  continue; // failed....
	
	int pick2 = 0;
	countTry = 0;
	do {
	    
	  pick2 = gRandom->Rndm() * allAntiBaryons.size();
	  const int typeTwo = particlesList.GetEntry(allAntiBaryons[pick2]).GetId();

	  if (usedAntiBaryons.find(pick2) == usedAntiBaryons.end()) {
	    if (typeTwo == -typeOne) {
	      break;
	    }
	  }

	} while(countTry++<nBaryons*2);
	
	if (countTry >= nBaryons*2-1)
	  continue; // failed....

	usedBaryons.insert(pick1);
	usedAntiBaryons.insert(pick2);

	// CHANGE PIDs HERE

	ParticleBlockEntry& p1 = particlesList.GetEntry(allBaryons[pick1]);
	ParticleBlockEntry& p2 = particlesList.GetEntry(allAntiBaryons[pick2]);

	if (ruVerbosity) {
	  cout << " BaryonProductionResampling> " << endl;
	  cout << " baryon1: index=" << allBaryons[pick1] << " "; p1.Dump();
	  cout << " baryon2: index=" << allAntiBaryons[pick2] << " "; p2.Dump();
	}
	
	if (!SetNewPId(CommonBlockParticleCONEX::ePiP, p1))
	  continue;
	
	if (!SetNewPId(CommonBlockParticleCONEX::ePiM, p2))
	  continue;

	if (ruVerbosity) {
	  cout << " BaryonProductionResampling> " << endl;
	  cout << " meson1: "; p1.Dump();
	  cout << " meson2: "; p2.Dump();
	}
	
      } // end loop change baryon pairs
      
    }

      
    // *********************************************
    double sumEnergyTotCheck = 0;
    map<int, SummaryInfo> sumCheck = particlesList.MakeSummary(ParticleBlock::eExact, ParticleBlock::eExcludeLeading);
    for (map<int, SummaryInfo>::iterator i = sumCheck.begin(); i != sumCheck.end(); ++i) {
      sumEnergyTotCheck += i->second.totEnergy;
      if (ruVerbosity>=1) {
	cout << setprecision(3);
	cout << " summary-check> id=" << setw(5) << i->first 
	     << " Etot=" << setw(14) << i->second.totEnergy
	     << " Ekin=" << setw(14) << i->second.kinEnergy
	     << " EkinMin=" << setw(14) << i->second.minKinEnergy
	     << " EkinMax=" << setw(14) << i->second.maxKinEnergy
	     << endl;
      }
    }
    
    if (ruVerbosity>0) {
      //      cout << setprecision(3);
      /*cout << " summary-check>         #pi0=" << sumNumberPi0Check << " (delta=" << sumNumberPi0Check-sumNumberPi0 << ")"
	   << " #pi-=" << sumNumberPiMCheck << " (delta=" << sumNumberPiMCheck-sumNumberPiM << ")"
	   << " #pi+=" << sumNumberPiPCheck << " (delta=" << sumNumberPiPCheck-sumNumberPiP << ")"
	   << " total=" << sumNumberPiPCheck+sumNumberPiMCheck+sumNumberPi0Check 
	   << " (delta=" << sumNumberPiPCheck+sumNumberPiMCheck+sumNumberPi0Check-totalNumberOfPions << ")"
	   << endl;    
      */
    }

  } // end of BaryonProductionResampling




  /* SecondaryResampling
   * *********************
   */
  void
  SecondaryResampling(const EResampleMode mode, const double factor, const int projId, 
		      ParticleBlock& particlesList) 
  {  
    switch(mode) {

    case eMultiplicity:
      MultiplicityResampling(factor, particlesList);
      break;

    case eEMRatio:
      EMRatioResampling(factor, particlesList);
      break;
        
    case eElasticity:
      ElasticityResampling(factor, particlesList);
      break;
	
    case eChargeRatio:
      ChargeRatioResampling(factor, projId, particlesList);
      break;

    case ePi0Spectrum:
      Pi0SpectrumResampling(factor, particlesList);
      break;

    case eLeadingRho0:
      Rho0Resampling(factor, projId, particlesList);
      break;

    case eBaryonProd:
      BaryonProductionResampling(factor, projId, particlesList);
      break;
      
    case eNone:
      cout << " ************ No resampling selected: continuing with CONEX ************* " << endl;
      break;
      
    default:
      cerr << "\n\n\n\n\n"
	   << " ************************************************************************** \n"
	   << " ************************************************************************** \n"
	   << " **                                                                      ** \n"
	   << " **                                                                      ** \n"
	   << " **                                                                      ** \n"
	   << " **     unknown resampling mode : " << mode << " (ignore)                             ** \n"
	   << " **                                                                      ** \n"
	   << " **                                                                      ** \n"
	   << " **                                                                      ** \n"
	   << " ************************************************************************** \n"
	   << " ************************************************************************** \n"
	   << endl;
    }
  } // end of SecondaryResampling



  bool
  SetNewPId(const CommonBlockParticleCONEX::EParticleId pid, ParticleBlockEntry& prt) 
  {
    prt.SetId(pid);
    
    const double massNew = CommonBlockParticleCONEX::GetMass(pid);
    
    const double energy = prt.GetEnergy();
    const double px = prt.GetPx();
    const double py = prt.GetPy();
    const double pz = prt.GetPz();
    const double Pi = sqrt(px*px + py*py + pz*pz);
    
    const double EkinNew = energy - massNew;
    if (EkinNew<=0){
      if (ruVerbosity>20) {
	cout << " SetNewPId> cannot switch particle id, since EkinNew=" << EkinNew << endl;
      }
      return false;
    }
    const double Pf = sqrt((energy+massNew) * (energy-massNew));
	
    double newPx = 0;
    double newPy = 0;
    double newPz = Pf;
    if (Pi<=0) {
      if (ruVerbosity>0)
	cout << " WARNING: P_f=" << Pf 
	     << " P_i=" << Pi
	  //	     << " factor=" << factor
	     << " id=" << prt.GetId()
	     << " [-> new momentum in z-direction]"
	     << endl;
    } else {
      if (Pf<=0) {
	if (ruVerbosity>0)
	  cout << " WARNING: P_f=" << Pf 
	       << " P_i=" << Pi
	    //<< " factor=" << factor
	       << " id=" << prt.GetId()
	       << " [-> new momentum forecd to ZERO]"
	       << endl;
	newPz = 0;
      } else {
	const double momentumFactor = Pf/Pi;
	newPx = px * momentumFactor;
	newPy = py * momentumFactor;
	newPz = pz * momentumFactor;
      }
    }
    
    prt.SetMass(massNew);
    prt.SetPx(newPx);
    prt.SetPy(newPy);
    prt.SetPz(newPz);
    
    return true;
  } // end SetNewPId
  
} // end namesapce resample

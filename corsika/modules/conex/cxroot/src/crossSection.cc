#include <crossSection.h>
#include <ConexDynamicInterface.h>
#include <conexExtensions.h>

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TSystem.h>
#include <TMath.h>

#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <string>

using namespace std;


const std::string cxroot::fgOutputPrefix = "   [crossSection]   ";


/*************************************************************
*
*  main()
*
**************************************************************/

int 
main(int argc, char **argv) 
{
  if (argc<2) {
    cout << "\n\nplease specify interaction model code\n\n"
	 << endl;
    return 1;
  }

  EHEModel HEModel = (EHEModel) atoi(argv[1]);

  if ((int)HEModel != 2 && 
      (int)HEModel != 4 && 
      (int)HEModel != 5 && 
      (int)HEModel != 6 && 
      (int)HEModel != 9) {
    cout << "\n\nUnknown interaction model: " << (int) HEModel << endl << endl;
    cout << "eQGSJet01 = 2\n"
      "eEposLHCR = 4\n"
      "eSibyll23 = 5\n"
      "eQGSJetII = 6\n"
      "eDPMJetIII = 9\n"
      "eQGSJetIII = 10\n"
	 << endl;

    exit(1);
  }

  const string conexRoot = (gSystem->Getenv("CONEX_ROOT") ?
			    gSystem->Getenv("CONEX_ROOT") :
			    string(gSystem->Getenv("PWD"))+string("/cfg"));
  cout << "CONEX: Using input steering files located at: \'" 
       << conexRoot << "\'" << endl;
  
  // set the RESAMPLING random generator seed
  // resample::InitResampling(gSeed + 100);

  const int crLength = conexRoot.length();
  
  int gSeed = 1;
  int nShower = 1;
  int maxDetail = -1;
#ifdef RU_CONEX_EXTENSIONS
  // turn off all modifications
  gParticleListMode = 0;
  gFactorCrossSection = 1;
  gFactorResampling = 1;
  gResamplingThreshold = 17 - 9;
  gResamplingMode = resample::eNone;
#endif

  int gSeedcx[3] = {0,0,0};
  gSeedcx[0] = gSeed;
  ConexDynamicInterface lib(HEModel);
  
  int targetCode = 0; // air

  lib.InitConex(nShower, gSeedcx, HEModel, maxDetail, 
#ifdef RU_CONEX_EXTENSIONS
		gParticleListMode,
#endif
#ifdef RU_CONEX_EXTENSIONS_HE_HIST
			     targetCode,
#endif
		conexRoot.c_str(), crLength);
  
/*
c---------------------------------------------------------------------
      double precision function rlam(np,ek,am)
c---------------------------------------------------------------------
c     inelastic interaction path for different models
c     projectile np = 1 : proton
c                     2 : charged pions
c                     3 : charged kaons
c                     4 : K-long
c                     5 : K-short
c                     6 : neutral pion
c                     7 : neutron
c                     9 : muon
c                   >10 : nuclei
c     ek=kinetic energy
c     am=mass
c     below 1 GeV, cross section is the one at 1 GeV ...
c---------------------------------------------------------------------
*/
  

  //double Emin = 15.-9;//log10(1.e2);      // GeV

  /*
  const double mp = 1.67262158e-27 * 1.e3; // [g]
  const double fgMeanMass = 14.5*mp;       // [g]
  const double A = fgMeanMass / (1.e-28*1.e4*1.e-3);
  */
  const double meanMass = 24158.575325530426; // CONEX-Number

  vector<int> primaries;
  primaries.push_back(1);  // proton
  primaries.push_back(2);  // charged pion
  primaries.push_back(3);  // charged kaon
  primaries.push_back(40);   // He
  primaries.push_back(160);  // O
  primaries.push_back(280);  // Si
  primaries.push_back(560);  // Fe
  
  
  // primary loop
  for (unsigned int iPrimary=0; iPrimary<primaries.size(); ++iPrimary) {
    
    int particleCode = primaries[iPrimary];
    int particleId = 0;
    switch (particleCode) {
    case 1: 
      particleId = 1120; // proton
      break;
    case 2:
      particleId = 120;
      break;
    case 3:
      particleId = 130;
      break;
    default:
      particleId = particleCode*10;
    }
    
    // const int nucleusA = max(1, particleCode / 10);
    
    // Hi-Lo threshold in conex is 80GeV
    //const double lgEmin = log10(80*nucleusA); //log10(1.e2);      // GeV
    const double lgEmin = log10(80);            //log10(1.e2);      // GeV
    const double lgEmax = 21.5 - 9;             //log10(1.e11);     // GeV
    const int n = 100;
    const double dE = (lgEmax-lgEmin) / (n-1);
       
    
    double mass = 0;
    lib.ParticleMass(particleId, mass);
    
    cout << "\n-------------------------------------------\n "
	 << " particleCode=" << particleCode << " particleId=" << particleId
	 << " mass=" << mass 
	 << " lgEmin=" << lgEmin << " lgEmax=" << lgEmax
	 << endl;
   
    TGraph gLambda(n);
    TGraph gSigma(n);
    
    // energy loop
    for (int i=0; i<n; i++) {
      
      double energy = pow(10., lgEmin + dE * i); // GeV
      
      const double lambda = lib.InteractionLength(particleCode, energy, mass);
      const double sigma = meanMass / lambda;
      
      if (TMath::IsNaN(lambda) || !TMath::Finite(lambda) ||
	  TMath::IsNaN(sigma) || !TMath::Finite(sigma))
	continue;
      
      const double eBinning = log10(energy) + 9;
      gLambda.SetPoint(i, eBinning, lambda);
      gSigma.SetPoint(i, eBinning, sigma);
      
      cout << " primary=" << particleId << " Ebin=" << i << " lgE="
	   << eBinning << " m=" << mass 
           << " lambda=" << lambda << " sigma=" << sigma << endl;
      
    } // end energy loop
    
    ostringstream name;
    
    switch (particleId) {      
       case 1120: name << "proton"; break;
       case 120: name << "pion"; break;
       case 130: name << "kaon"; break;
       case 400: name << "helium"; break;
       case 1200: name << "carbon"; break;
       case 1400: name << "nitrogen"; break;
       case 1600: name << "oxygen"; break;
       case 2800: name << "silicon"; break;
       case 5600: name << "iron"; break;
       default: name << "unknown"; break;
    }
    
    name << "_";
    
    switch(HEModel) {
      
        case 1: name << "NEXUS"; break;
        case 2: name << "QGSJET"; break;
        case 4: name << "EPOS"; break;
        case 5: name << "SIBYLL"; break;
        case 6: name << "QGSJETII"; break;
        case 9: name << "DPMJETIII"; break;
        case 9: name << "QGSJETIII"; break;
        default: name << "unknown"; break;
    }
    
    cout << " writing: " << name.str() << endl;
    
    ostringstream nameLambda, nameSigma;
    nameLambda << name.str() << "_Lambda";
    nameSigma << name.str() << "_Sigma";
    
    TFile outFile ("Sigma.root", "UPDATE");
    gLambda.Write (nameLambda.str().c_str(), TFile::kOverwrite);
    gSigma.Write (nameSigma.str().c_str(), TFile::kOverwrite);
    outFile.Close();
    
  } // end particle loop
  
  cout << endl << endl;

  return 0;
}

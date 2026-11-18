gROOT->Reset();
//
//  cx2cors()
//
//  ROOT macro to print list of particles produced by the first interaction
//  (to be used with STACKIN option in CORSIKA)
//
//    usage: .x cx2cors.C("conex.root",3) 
//
//       to print the list of particles from the first interaction of 
//       event 3 in file conex.root
//  
#include <iostream>
#include <iomanip>
using namespace std;

void cx2cors(char * file="", int ievt=-1) {

  if ( file=="" || ievt<0 ) {
    cout << "\n usage: .x plotProfile.C(\"conex.root\",3) \n"
	 << "         to print the list of particles from the first interaction of event 3 in file conex.root\n" << endl;
    return;
  }



  // open file

  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(file);
  if (!f) {
    f = new TFile(file);
  }
  if (f->IsZombie()) return;


  // set and read tree

  TTree * FirstInteraction = (TTree *) f->Get("FirstInteraction");
  int nevent=FirstInteraction->GetEntries();
  if ( ievt >= nevent ) {
    cout << " Error - requested event number exceeds number of events in file!"
	 << "\n         ievt=" << ievt << " nevent=" << nevent 
	 << "\n exit ..." << endl;
    return;
  }
  int gnInt1;           // number of interactions
  const int maxNInt = 100;
  double gkinel1[maxNInt];      // inelasticity
  int gpId1[maxNInt];           // parent id
  double gpEnergy1[maxNInt];    // parent energy
  int gMult1[maxNInt];          // multiplicity
  int gMatg1[maxNInt];          // target nucleus mass
  double gDepth1[maxNInt];      // depth
  double gHeight1[maxNInt];     // height
  
  int gNParticles1;
  const int maxNPart = 100000;
  double gEnergy1[maxNPart];    // size = 75000 from conex.incnex
  double gpx1[maxNPart];
  double gpy1[maxNPart];
  double gpz1[maxNPart];
  int gType1[maxNPart];
  int gIdInt1[maxNPart];
  
  FirstInteraction->SetBranchAddress("nInt",&gnInt1);                // number of interatctions
  FirstInteraction->SetBranchAddress("kinel",gkinel1);       // inelasticity
  FirstInteraction->SetBranchAddress("pId", gpId1);            // parent id
  FirstInteraction->SetBranchAddress("pEnergy", gpEnergy1);// parent energy
  FirstInteraction->SetBranchAddress("mult", gMult1);         // multiplicity
  FirstInteraction->SetBranchAddress("matg", gMatg1);         // mass of target nucleus
  FirstInteraction->SetBranchAddress("depth", gDepth1);      // slant depth
  FirstInteraction->SetBranchAddress("height", gHeight1);   // height
  
  FirstInteraction->SetBranchAddress("nPart",&gNParticles1);    // number of particles
  FirstInteraction->SetBranchAddress("Energy",gEnergy1);// energy
  FirstInteraction->SetBranchAddress("px",gpx1);// Momentum px
  FirstInteraction->SetBranchAddress("py",gpy1);// Momentum py
  FirstInteraction->SetBranchAddress("pz",gpz1);// Momentum pz
  FirstInteraction->SetBranchAddress("Type",gType1);      // type
  FirstInteraction->SetBranchAddress("idInt",gIdInt1);   // interaction number 
  
  FirstInteraction->GetEntry(ievt);
  
  cout << " " << gMult1[0] << "   " <<gpEnergy1[0] << endl ;
  /*  cout << " " << gMult1[1] << "   " <<gpEnergy1[1] << endl ;
  cout << " " << gMult1[2] << "   " <<gpEnergy1[2] << endl ;
  cout << " " << gMult1[3] << "   " <<gpEnergy1[3] << endl ; */
  int oldFormat = cout.flags(ios::scientific);
  int oldPrecision             = cout.precision(7);
  for (int i=0; i< gNParticles1; i++) {
    // Particles from the first interaction
    if ( gIdInt1[i] == 1 )          
    cout << setw(5) << i+1  << setw(5) << gType1[i]  << setw(16) << gEnergy1[i] << setw(16) << gpz1[i]<< setw(16) << gpx1[i]<< setw(16) << gpy1[i] << endl ;
    }
  cout << " Height (cm) :  " << gHeight1[0]*100 << endl ;
  cout << " Target :  " << gMatg1[0] << endl ;

  FirstInteraction->ResetBranchAddresses();

  return;

}


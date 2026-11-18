#include <conexExtensions.h>
#include <conexConfig.h>
//#include <leadingInteractionsData.h>

#ifdef CONEX_EXTENSIONS
#include <Verbosity.h>
#include <CommonBlockWrapperCONEX.h>
#include <CommonBlockParticleCONEX.h>
#include <CommonBlockWrapperSIBYLL.h>
#include <CommonBlockParticleSIBYLL.h>
#include <CommonBlockWrapperSIBYLL_LAB.h>
#include <resample.h>
#endif

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <map>
#include <sstream>

#include <TDirectory.h>
#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>

using namespace std;

#ifdef CONEX_EXTENSIONS

using namespace resample;

const int fVerbosity = 0;


CXParticle gCXParticle;
NXParticle gNXParticle;
CXParticle gCXProjectile;
COSSINS gCXGeometry;
CXInteraction gCXInteraction;
int gInteractionCounter;
int gCXParticleCounter;
int gCXInteractionCounter;

int gCXSeed1, gCXSeed2, gCXSeed3;

TFile *gListTreeFile = 0;
TTree *gListTree = 0;
TTree *gProjectileTree = 0;
TTree *gInteractionTree = 0;
TTree *gSeedTree = 0;

int gParticleListReadInteraction;




// ---------------------------------------------------------------------
void
CXParticle::Dump() const 
{
  std::cout << " cxParticle: (id,weight,generation) (" << (int) id 
	    << "," << weight << "," << generation << ")"
	    << " P(x,y,z): (" << Px << "," << Py << "," << Pz << ")"
	    << " E,m: " << Energy << "," << mass
	    << " pos(x,y,z): (" << x << "," << y << "," << height << ")"
	      << " slant (traversed/ToImpact): (" << slantTraversed << "," << slantToImpact << ")"
	    << std::endl;
}


// ---------------------------------------------------------------------
void 
NXParticle::Dump() const 
{
  std::cout << Px << " " << Py << " " << Pz << " " << Energy << " " << mass 
	    << " " << x << " " << y << " " << z 
	    << " " << time << " " << id << " " << t_formation << " " << t_destruction
	    << " " << id_origin << " " << id_father << " " << id_mother << " " << status
	    << std::endl;
}



// ---------------------------------------------------------------------
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
// RU Mon Oct 23 09:44:08 CEST 2006
void 
initListTree(const int i, 
	     const int particleListMode) 
{
  int mode = 0;
  if (particleListMode < 0) {
    gParticleListReadInteraction = -particleListMode;
    mode = 1;
  } else if (particleListMode > 0)
    mode = particleListMode + 2;  
  
  if (mode<1) 
    return;
  
  TDirectory *save = gDirectory;
  
  // open file, initialize seed-tree
  if (!gListTreeFile) { // do only once
    
    if (mode==1) { // READ file preparations
      
      gListTreeFile = new TFile(gParticleListModeFile.c_str(), "READ");
      
      if (!gListTreeFile || gListTreeFile->IsZombie()) {
	cerr << " initListTree: connot open file: " << gParticleListModeFile << endl;
	exit(1);
      }

    } else { // WRITE file preparations

      gListTreeFile = new TFile(gParticleListModeFile.c_str(), "RECREATE");

      if (!gListTreeFile || gListTreeFile->IsZombie()) {
	cerr << " initListTree: connot open file: " << gParticleListModeFile << endl;
	exit(1);
      }
    }
  }
  
  gListTreeFile->cd();

  ostringstream treeName, treeName2, treeName3, treeName4;
  treeName << "ParticleList_Event" << i;
  treeName2 << "InteractionList_Event" << i;
  treeName3 << "Projectile_Event" << i;
  treeName4 << "Seeds_Event" << i;
  gCXParticleCounter = 0;
  gCXInteractionCounter = 0;

  // initialize particle data tree
  if (mode==1) { // READing

    if (!gSeedTree) {
    
      gSeedTree = (TTree*) gListTreeFile->Get(treeName4.str().c_str());
      if (!gSeedTree) {
	cerr << " initListTree: connot read TTree seeds!" << endl;
	exit(1);
      }
      gSeedTree->SetBranchAddress("seed1a", &gCXSeed1); // seeda is before each new shower
      gSeedTree->SetBranchAddress("seed2a", &gCXSeed2);
      gSeedTree->SetBranchAddress("seed3a", &gCXSeed3);
      /*gSeedTree->SetBranchAddress("seed1b", &gCXSeed1b); // seedb is after the first interaction
	gSeedTree->SetBranchAddress("seed2b", &gCXSeed2b);
	gSeedTree->SetBranchAddress("seed3b", &gCXSeed3b);*/
    }
    
    
    if (!gListTree) {
      gListTree = (TTree*) gListTreeFile->Get(treeName.str().c_str());
      if (!gListTree) {
	cerr << " initListTree: connot read TTree: " << treeName.str() << endl;
	exit(1);
      }
      
      /*
	gListTree->SetBranchAddress("Px", &gCXParticle.Px);
	gListTree->SetBranchAddress("Py", &gCXParticle.Py);
	gListTree->SetBranchAddress("Pz", &gCXParticle.Pz);
	gListTree->SetBranchAddress("Energy", &gCXParticle.Energy);
	gListTree->SetBranchAddress("mass", &gCXParticle.mass);
	gListTree->SetBranchAddress("x", &gCXParticle.x);
	gListTree->SetBranchAddress("y", &gCXParticle.y);
	gListTree->SetBranchAddress("height", &gCXParticle.height);
	gListTree->SetBranchAddress("time", &gCXParticle.time);
	gListTree->SetBranchAddress("id", &gCXParticle.id);
	gListTree->SetBranchAddress("weight", &gCXParticle.weight);
	gListTree->SetBranchAddress("generation", &gCXParticle.generation);
	gListTree->SetBranchAddress("slantTraversed", &gCXParticle.slantTraversed);
	gListTree->SetBranchAddress("xShower", &gCXParticle.xShower);
	gListTree->SetBranchAddress("yShower", &gCXParticle.yShower);
	gListTree->SetBranchAddress("slantToImpact", &gCXParticle.slantToImpact);
	gListTree->SetBranchAddress("interactionCounter", &gInteractionCounter);
      */
      
      gListTree->SetBranchAddress("Px", &gNXParticle.Px);
      gListTree->SetBranchAddress("Py", &gNXParticle.Py);
      gListTree->SetBranchAddress("Pz", &gNXParticle.Pz);
      gListTree->SetBranchAddress("Energy", &gNXParticle.Energy);
      gListTree->SetBranchAddress("mass", &gNXParticle.mass);
      gListTree->SetBranchAddress("x", &gNXParticle.x);
      gListTree->SetBranchAddress("y", &gNXParticle.y);
      gListTree->SetBranchAddress("z", &gNXParticle.z);
      gListTree->SetBranchAddress("time", &gNXParticle.time);
      gListTree->SetBranchAddress("id", &gNXParticle.id);
      gListTree->SetBranchAddress("t_formation", &gNXParticle.t_formation);
      gListTree->SetBranchAddress("t_destruction", &gNXParticle.t_destruction);
      gListTree->SetBranchAddress("id_origin", &gNXParticle.id_origin);
      gListTree->SetBranchAddress("id_father", &gNXParticle.id_father);
      gListTree->SetBranchAddress("id_mother", &gNXParticle.id_mother);
      gListTree->SetBranchAddress("status", &gNXParticle.status);
      gListTree->SetBranchAddress("interactionCounter", &gInteractionCounter);
    }


    // read projectile
    if (!gProjectileTree) {
      gProjectileTree = (TTree*) gListTreeFile->Get(treeName3.str().c_str());
      if (!gProjectileTree) {
	cerr << " initListTree: connot read TTree: " << treeName3.str() << endl;
	exit(1);
      }
      
      gProjectileTree->SetBranchAddress("Px", &gCXProjectile.Px);
      gProjectileTree->SetBranchAddress("Py", &gCXProjectile.Py);
      gProjectileTree->SetBranchAddress("Pz", &gCXProjectile.Pz);
      gProjectileTree->SetBranchAddress("Energy", &gCXProjectile.Energy);
      gProjectileTree->SetBranchAddress("mass", &gCXProjectile.mass);
      gProjectileTree->SetBranchAddress("x", &gCXProjectile.x);
      gProjectileTree->SetBranchAddress("y", &gCXProjectile.y);
      gProjectileTree->SetBranchAddress("height", &gCXProjectile.height);
      gProjectileTree->SetBranchAddress("time", &gCXProjectile.time);
      gProjectileTree->SetBranchAddress("id", &gCXProjectile.id);
      gProjectileTree->SetBranchAddress("weight", &gCXProjectile.weight);
      gProjectileTree->SetBranchAddress("generation", &gCXProjectile.generation);
      gProjectileTree->SetBranchAddress("slantTraversed", &gCXProjectile.slantTraversed);
      gProjectileTree->SetBranchAddress("xShower", &gCXProjectile.xShower);
      gProjectileTree->SetBranchAddress("yShower", &gCXProjectile.yShower);
      gProjectileTree->SetBranchAddress("slantToImpact", &gCXProjectile.slantToImpact);
      gProjectileTree->SetBranchAddress("interactionCounter", &gInteractionCounter);
      // the geometry for rotations
      gProjectileTree->SetBranchAddress("s0xs", &gCXGeometry.s0xs);
      gProjectileTree->SetBranchAddress("c0xs", &gCXGeometry.c0xs);
      gProjectileTree->SetBranchAddress("s0s", &gCXGeometry.s0s);
      gProjectileTree->SetBranchAddress("c0s", &gCXGeometry.c0s);


      // ----------------------------------------------
      analyseInteractions();
    }

  } else { // WRITEing

    
    gSeedTree = new TTree(treeName4.str().c_str(), "seeds");
    gSeedTree->Branch("seed1a", &gCXSeed1, "seed1/I");
    gSeedTree->Branch("seed2a", &gCXSeed2, "seed2/I");
    gSeedTree->Branch("seed3a", &gCXSeed3, "seed3/I");
    /*gSeedTree->Branch("seed1b", &gCXSeed1b, "seed1b/I");
    gSeedTree->Branch("seed2b", &gCXSeed2b, "seed2b/I");
    gSeedTree->Branch("seed3b", &gCXSeed3b, "seed3b/I");*/


    gListTree = new TTree(treeName.str().c_str(), "particle list tree");

    /*
    gListTree->Branch("Px", &gCXParticle.Px, "Px/D");
    gListTree->Branch("Py", &gCXParticle.Py, "Py/D");
    gListTree->Branch("Pz", &gCXParticle.Pz, "Pz/D");
    gListTree->Branch("Energy", &gCXParticle.Energy, "Energy/D");
    gListTree->Branch("mass", &gCXParticle.mass, "mass/D");
    gListTree->Branch("x", &gCXParticle.x, "x/D");
    gListTree->Branch("y", &gCXParticle.y, "y/D");
    gListTree->Branch("height", &gCXParticle.height, "height/D");
    gListTree->Branch("time", &gCXParticle.time, "time/D");
    gListTree->Branch("id", &gCXParticle.id, "id/D");
    gListTree->Branch("weight", &gCXParticle.weight, "weight/D");
    gListTree->Branch("generation", &gCXParticle.generation, "generation/D");
    gListTree->Branch("slantTraversed", &gCXParticle.slantTraversed, "slantTraversed/D");
    gListTree->Branch("xShower", &gCXParticle.xShower, "xShower/D");
    gListTree->Branch("yShower", &gCXParticle.yShower, "yShower/D");
    gListTree->Branch("slantToImpact", &gCXParticle.slantToImpact, "slantToImpact/D");
    gListTree->Branch("interactionCounter", &gInteractionCounter, "interactionCounter/I");
    */

    gListTree->Branch("Px", &gNXParticle.Px, "Px/D");
    gListTree->Branch("Py", &gNXParticle.Py, "Py/D");
    gListTree->Branch("Pz", &gNXParticle.Pz, "Pz/D");
    gListTree->Branch("Energy", &gNXParticle.Energy, "Energy/D");
    gListTree->Branch("mass", &gNXParticle.mass, "mass/D");
    gListTree->Branch("x", &gNXParticle.x, "x/D");
    gListTree->Branch("y", &gNXParticle.y, "y/D");
    gListTree->Branch("z", &gNXParticle.z, "z/D");
    gListTree->Branch("time", &gNXParticle.time, "time/D");
    gListTree->Branch("id", &gNXParticle.id, "id/I");
    gListTree->Branch("t_formation", &gNXParticle.t_formation, "t_formation/D");
    gListTree->Branch("t_destruction", &gNXParticle.t_destruction, "t_destruction/D");
    gListTree->Branch("id_origin", &gNXParticle.id_origin, "id_origin/I");
    gListTree->Branch("id_father", &gNXParticle.id_father, "id_father/I");
    gListTree->Branch("id_mother", &gNXParticle.id_mother, "id_mother/I");
    gListTree->Branch("status", &gNXParticle.status, "status/I");
    gListTree->Branch("interactionCounter", &gInteractionCounter, "interactionCounter/I");


    gProjectileTree = new TTree(treeName3.str().c_str(), "projectile tree");

    gProjectileTree->Branch("Px", &gCXProjectile.Px, "Px/D");
    gProjectileTree->Branch("Py", &gCXProjectile.Py, "Py/D");
    gProjectileTree->Branch("Pz", &gCXProjectile.Pz, "Pz/D");
    gProjectileTree->Branch("Energy", &gCXProjectile.Energy, "Energy/D");
    gProjectileTree->Branch("mass", &gCXProjectile.mass, "mass/D");
    gProjectileTree->Branch("x", &gCXProjectile.x, "x/D");
    gProjectileTree->Branch("y", &gCXProjectile.y, "y/D");
    gProjectileTree->Branch("height", &gCXProjectile.height, "height/D");
    gProjectileTree->Branch("time", &gCXProjectile.time, "time/D");
    gProjectileTree->Branch("id", &gCXProjectile.id, "id/D");
    gProjectileTree->Branch("weight", &gCXProjectile.weight, "weight/D");
    gProjectileTree->Branch("generation", &gCXProjectile.generation, "generation/D");
    gProjectileTree->Branch("slantTraversed", &gCXProjectile.slantTraversed, "slantTraversed/D");
    gProjectileTree->Branch("xShower", &gCXProjectile.xShower, "xShower/D");
    gProjectileTree->Branch("yShower", &gCXProjectile.yShower, "yShower/D");
    gProjectileTree->Branch("slantToImpact", &gCXProjectile.slantToImpact, "slantToImpact/D");
    gProjectileTree->Branch("interactionCounter", &gInteractionCounter, "interactionCounter/I");
    // the geometry for rotations
    gProjectileTree->Branch("s0xs", &gCXGeometry.s0xs, "s0xs/D");
    gProjectileTree->Branch("c0xs", &gCXGeometry.c0xs, "c0xs/D");
    gProjectileTree->Branch("s0s", &gCXGeometry.s0s, "s0s/D");
    gProjectileTree->Branch("c0s", &gCXGeometry.c0s, "c0s/D");

    
    gInteractionTree = new TTree(treeName2.str().c_str(), "interaction list tree");
    /*
    gInteractionTree->Branch("pId", &gInteractionDataPartial.pId, "pId/I");
    gInteractionTree->Branch("pEnergy", &gInteractionDataPartial.pEnergy, "pEnergy/D");
    gInteractionTree->Branch("energyCM", &gInteractionDataPartial.energyCM, "energyCM/D");
    gInteractionTree->Branch("kinel", &gInteractionDataPartial.kinel, "kinel/D");
    gInteractionTree->Branch("mult", &gInteractionDataPartial.mult, "mult/I");
    gInteractionTree->Branch("matg", &gInteractionDataPartial.matg, "matg/I");
    gInteractionTree->Branch("Depth", &gInteractionDataPartial.Depth, "Depth/D");
    gInteractionTree->Branch("Height", &gInteractionDataPartial.Height, "Height/D");
    */
    
    gInteractionTree->Branch("idProj", &gCXInteraction.idProj, "idProj/I");
    gInteractionTree->Branch("idTarg", &gCXInteraction.idTarg, "idTarg/I");
    gInteractionTree->Branch("eProj", &gCXInteraction.eProj, "eProj/D");
    gInteractionTree->Branch("eCMS", &gCXInteraction.eCMS, "eCMS/D");
    gInteractionTree->Branch("mult", &gCXInteraction.mult, "mult/I");
    gInteractionTree->Branch("eProd", &gCXInteraction.eProd, "eProd/D");
    gInteractionTree->Branch("interactionCounter", &gCXInteraction.interactionCounter, "interactionCounter/I");
  }

  save->cd();
}


// ---------------------------------------------------------------------
// sort interactions for energy
void
analyseInteractions()
{
  const bool silent = true;
  
  map<double,int> interactionList;
  const int nP = gListTree->GetEntries();
  if (!silent)
    cout << " analyseInteractions " << nP << endl;
  int curInt = -1;
  int mult = 0;
  double eSum = 0;
  for (int iP=0; iP<nP; ++iP) {
    gListTree->GetEntry(iP);
    //cout << " curInt " << curInt << " " << gInteractionTree << " " << gNXParticle.Energy << endl;
    if (curInt==gInteractionCounter) {
      eSum += gNXParticle.Energy;
      ++mult;
    } else {
      if (curInt>=0) {
	interactionList[eSum] = curInt;
	if (!silent)
	  cout << "Interaction: " << curInt << ", eSum: " << eSum << ", mult: " << mult << endl;
      }
      curInt = gInteractionCounter;
      mult = 0;
      eSum = 0;
    }
  }

  if (!silent) {    
    cout << "Interaction: " << curInt << ", eSum: " << eSum << ", mult: " << mult << endl;
    int cI = 1;
    for (map<double,int>::const_iterator iI=interactionList.begin(); 
	 iI!=interactionList.end();
	 ++iI) {
      cout << " sorted: " << cI << " " << iI->second << " " << iI->first << endl;
      cI++;
    }
  }
  
  if (gParticleListReadInteraction>interactionList.size()) {
    cout << "Asked for interaction: " << gParticleListReadInteraction << " of a total of " 
	 << interactionList.size() << " found in input file " << endl;
    exit(5);
  }

  map<double,int>::const_iterator iI = interactionList.end();
  for(int i=0; i<gParticleListReadInteraction; ++i) {
    --iI;
  }
  cout << "AnalyseInteractions: search for " << gParticleListReadInteraction 
       << " highest energy interaction" << endl;
  gParticleListReadInteraction = iI->second; // select interaction
  cout << "Identified # " << iI->second << " with E=" << iI->first << " GeV, "
       << "lgE=" << log10(iI->first)+9 
       << endl;
  //  gSeedTree->GetEntry(gParticleListReadInteraction-1);
  
  if (!silent)
    cout << " THIS ONE: " << iI->second << " " << iI->first << endl;
}


// ---------------------------------------------------------------------
void 
finishListTree(const int mode) 
{
  if (mode==0)
    return;
  TDirectory *save = gDirectory;
  gListTreeFile->cd();
  if (mode<1) 
    return;
  if (mode>=2) {
    gListTree->Write();
    gProjectileTree->Write();
    gInteractionTree->Write();
    gSeedTree->Write();
  }
  delete gListTree;
  delete gInteractionTree;
  delete gProjectileTree;
  save->cd();
}

void 
presetSeed(int seed[3])
{
  gSeedTree->GetEntry(gParticleListReadInteraction-1);
  seed[0] = gCXSeed1;
  seed[1] = gCXSeed2;
  seed[2] = gCXSeed3;
}





#endif

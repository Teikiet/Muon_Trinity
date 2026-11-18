#include <conexExtensionsF.h>
#include <conexExtensions.h>
//#include <leadingInteractionsData.h>

#include <Verbosity.h>
#include <CommonBlockWrapperCONEX.h>
#include <CommonBlockParticleCONEX.h>
#include <CommonBlockWrapperSIBYLL.h>
#include <CommonBlockParticleSIBYLL.h>
#include <CommonBlockWrapperSIBYLL_LAB.h>
#include <resample.h>

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
using namespace resample;

#ifdef CONEX_EXTENSIONS


void 
get_particle_from_list_(CXParticle &p, 
			int &goon)
{
  if (fVerbosity>2) {
    cout << "io(" << gCXParticleCounter << ") " << flush;
  }
  gListTree->GetEntry(gCXParticleCounter++);
  p = gCXParticle;
  if (fVerbosity>2) {
    p.Dump();
  }
  goon = (gCXParticleCounter<gListTree->GetEntries());
} 


void 
get_particle_from_list_cx_(NXBlock& data, 
			   int& indexF,
			   int& goon)
{
  do {
    
    if (gCXParticleCounter>=gListTree->GetEntries()) {
      goon = 0;
      return;
    }
    
    gListTree->GetEntry(gCXParticleCounter);
    ++gCXParticleCounter;
    
  } while(gInteractionCounter<gParticleListReadInteraction); // forward to the right interaction
  
  if (gInteractionCounter!=gParticleListReadInteraction) { // next interaction
    goon = 0;
    return;
  }
  
  goon = 1;
  const int index = indexF - 1;
  
  data.xsptl[index][0] = gNXParticle.Px;         //c     xsptl(1,i) ..... x-component of particle momentum 
  data.xsptl[index][1] = gNXParticle.Py;         //c     xsptl(2,i) ..... y-component of particle momentum 
  data.xsptl[index][2] = gNXParticle.Pz;         //c     xsptl(3,i) ..... z-component of particle momentum 
  data.xsptl[index][3] = gNXParticle.Energy;     //c     xsptl(4,i) ..... particle energy 
  data.xsptl[index][4] = gNXParticle.mass;      // c     xsptl(5,i) ..... particle mass 
  data.xsorptl[index][0] = gNXParticle.x;        //  c     xsorptl(1,i) ... x-component of formation point
  data.xsorptl[index][1] = gNXParticle.y;        // c      xsorptl(2,i) ... y-component of formation point
  data.xsorptl[index][2] = gNXParticle.z;        // c      xsorptl(3,i) ... z-component of formation point
  data.xsorptl[index][3] = gNXParticle.time;     // c      xsorptl(4,i) ... formation time
  data.idptlxs[index] = gNXParticle.id;           //c      idptlxs(i) ...... particle id
  data.xstivptl[index][0] = gNXParticle.t_formation; // c   xstivptl(1,i) ... formation time (always in the pp-cms!)
  data.xstivptl[index][1] = gNXParticle.t_destruction;// c  xstivptl(2,i) ... destruction time (always in the pp-cms!)
  data.ityptlxs[index] = gNXParticle.id_origin;    // c     ityptlxs(i)  .... type of particles origin:
  data.iorptlxs[index] = gNXParticle.id_father;     // c    iorptlxs(i) ..... particle number of father (if .le. 0 : no father)
  data.jorptlxs[index] = gNXParticle.id_mother;    // c     jorptlxs(i) ..... particle number of mother (if .le. 0 : no mother)
  data.istptlxs[index] = gNXParticle.status;        //c     istptlxs(i) ..... status: 40 and 41 : Remnant
  
  /*
    cout << "p5=" << data.xsptl[index][0] << " " << data.xsptl[index][1] << " " << data.xsptl[index][2]
       << " " << data.xsptl[index][3]<< " " << data.xsptl[index][4]
       << endl;
  */
  
  //goon = (gCXParticleCounter<gListTree->GetEntries());
  //goon = gCXParticleCounter<100;
} 

void 
get_projectile_(CXParticle &p, 
		COSSINS& c)
{
  //cout << "Get_Projectile: interaction=" << gParticleListReadInteraction << endl;
  //gProjectileTree->GetEntry(gCXInteractionCounter++);
  gProjectileTree->GetEntry(gParticleListReadInteraction-1);
  p = gCXProjectile;
  c = gCXGeometry;
}


void 
write_particle_(int& interactionCounter, 
		CXParticle &p) 
{
  if (fVerbosity>2) {
    cout << "io(" << gCXParticleCounter++ << ") " << flush;
    p.Dump();
  }
  gCXParticle = p;
  gInteractionCounter = interactionCounter;
  gListTree->Fill();
}

void 
write_projectile_(int& interactionCounter, 
		  CXParticle &p,
		  COSSINS& c) 
{
  gCXProjectile = p;
  gCXGeometry = c;
  gInteractionCounter = interactionCounter;
  gProjectileTree->Fill();
}


void
write_interaction_(int& interactionCounter,
		   double& eProj,
		   double& eCMS,
		   double& eProd,
		   int& idProj,
		   int& idTarg,
		   int& mult) 
{
  gCXInteraction.eProj = eProj;
  gCXInteraction.eCMS = eCMS;
  gCXInteraction.eProd = eProd;
  gCXInteraction.idProj = idProj;
  gCXInteraction.idTarg = idTarg;
  gCXInteraction.mult = mult;
  gCXInteraction.interactionCounter = interactionCounter;
  gInteractionTree->Fill();
}


void
write_particle_cx_(int& interactionCounter, 
		   NXBlock& data, 
		   int& indexF)
{
  const int index = indexF - 1;
  gNXParticle.Px = data.xsptl[index][0];         //c     xsptl(1,i) ..... x-component of particle momentum 
  gNXParticle.Py = data.xsptl[index][1];         //c     xsptl(2,i) ..... y-component of particle momentum 
  gNXParticle.Pz = data.xsptl[index][2];         //c     xsptl(3,i) ..... z-component of particle momentum 
  gNXParticle.Energy = data.xsptl[index][3];     //c     xsptl(4,i) ..... particle energy 
  gNXParticle.mass = data.xsptl[index][4];      // c     xsptl(5,i) ..... particle mass 
  gNXParticle.x = data.xsorptl[index][0];        //  c     xsorptl(1,i) ... x-component of formation point
  gNXParticle.y = data.xsorptl[index][1];        // c      xsorptl(2,i) ... y-component of formation point
  gNXParticle.z = data.xsorptl[index][2];        // c      xsorptl(3,i) ... z-component of formation point
  gNXParticle.time = data.xsorptl[index][3];     // c      xsorptl(4,i) ... formation time
  gNXParticle.id = data.idptlxs[index];           //c      idptlxs(i) ...... particle id
  gNXParticle.t_formation = data.xstivptl[index][0]; // c   xstivptl(1,i) ... formation time (always in the pp-cms!)
  gNXParticle.t_destruction = data.xstivptl[index][1];// c  xstivptl(2,i) ... destruction time (always in the pp-cms!)
  gNXParticle.id_origin = data.ityptlxs[index];    // c     ityptlxs(i)  .... type of particles origin:
  gNXParticle.id_father = data.iorptlxs[index];     // c    iorptlxs(i) ..... particle number of father (if .le. 0 : no father)
  gNXParticle.id_mother = data.jorptlxs[index];    // c     jorptlxs(i) ..... particle number of mother (if .le. 0 : no mother)
  gNXParticle.status = data.istptlxs[index];        //c     istptlxs(i) ..... status: 40 and 41 : Remnant
  gInteractionCounter = interactionCounter;

  if (fVerbosity>2) {
    cout << "io(" << gCXParticleCounter++ << ") " << flush;
    gNXParticle.Dump();
  }

  gListTree->Fill();
}


void
write_seed_cx_(int seed[3]) 
{
  if (fVerbosity>2) {
    cout << " write_seed_cx_ " << seed[0] << " " << seed[1] << " " << seed[2] << endl;
  }
  gCXSeed1 = seed[0];
  gCXSeed2 = seed[1];
  gCXSeed3 = seed[2];
  gSeedTree->Fill();
}


void
set_seed_cx_(int seed[3])
{
  gSeedTree->GetEntry(gParticleListReadInteraction-1);
  seed[0] = gCXSeed1;
  seed[1] = gCXSeed2;
  seed[2] = gCXSeed3;
  if (fVerbosity>2) {
    cout << " set_seed_cx_ " << seed[0] << " " << seed[1] << " " << seed[2]
	 << endl;
  }
}

// ---------------------------------------------------------------------------------



#endif

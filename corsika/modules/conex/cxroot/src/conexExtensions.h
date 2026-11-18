#ifndef _include_conex_extensions_h_
#define _include_conex_extensions_h_

#include <conexConfig.h>
#include <conexExtensionsF.h>

#include <string>

class TFile;
class TTree;

#ifdef CONEX_EXTENSIONS

extern const int fVerbosity;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
// RU Mon Oct 23 10:52:35 CEST 2006
extern int gParticleListMode;
extern int gParticleListReadInteraction;
extern std::string gParticleListModeFile;



extern "C" {
  struct NXBlock {
#warning CHECK THIS HARDCODED NUMBER FOR POTENTIAL CHANGES:
    const static int mxptlxs = 100000;
    //parameter (nflavxs=6)           !max nr of flavors
    //parameter (mxptlxs=100000)   !max nr of particles in conex MC particle list
    //common/xscptl/
    double xsptl[mxptlxs][5];
    double xstivptl[mxptlxs][2];
    double xsorptl[mxptlxs][4];
    int ibptlxs[mxptlxs][4];
    int ifrptlxs[mxptlxs][2];
    int istptlxs[mxptlxs];
    int jorptlxs[mxptlxs];
    int idptlxs[mxptlxs];
    int iorptlxs[mxptlxs];
    int ityptlxs[mxptlxs];
    int nptlxs;
  };

  struct COSSINS {
    double s0xs;
    double c0xs;
    double s0s;
    double c0s;
  };
}


// ------------------------------------------------
//  function definitions
//


// *******************************
//  CONEX binary particle definition 
//
struct CXParticle {
  double Px;       // 5 momentum
  double Py;
  double Pz;
  double Energy;   // (total energy)
  double mass;
  double x;        // observer frame
  double y;
  double height;
  double time;     // [m]
  double id;       
  double weight;
  double generation;
  double slantTraversed;
  double xShower;
  double yShower;
  double slantToImpact;
  void Dump() const;
};

// int interactionCounter;

// *******************************
//  NEXUS binary particle definition 
//
struct NXParticle {
  double Px;         //c     xsptl(1,i) ..... x-component of particle momentum 
  double Py;         //c     xsptl(2,i) ..... y-component of particle momentum 
  double Pz;         //c     xsptl(3,i) ..... z-component of particle momentum 
  double Energy;     //c     xsptl(4,i) ..... particle energy 
  double mass;      // c     xsptl(5,i) ..... particle mass 
  double x;        //  c     xsorptl(1,i) ... x-component of formation point
  double y;        // c      xsorptl(2,i) ... y-component of formation point
  double z;        // c      xsorptl(3,i) ... z-component of formation point
  double time;     // c      xsorptl(4,i) ... formation time
  int id;           //c      idptlxs(i) ...... particle id
  double t_formation; // c   xstivptl(1,i) ... formation time (always in the pp-cms!)
  double t_destruction;// c  xstivptl(2,i) ... destruction time (always in the pp-cms!)
  int id_origin;    // c     ityptlxs(i)  .... type of particles origin:
  int id_father;     // c    iorptlxs(i) ..... particle number of father (if .le. 0 : no father)
  int id_mother;    // c     jorptlxs(i) ..... particle number of mother (if .le. 0 : no mother)
  int status;        //c     istptlxs(i) ..... status: 40 and 41 : Remnant
  void Dump() const;
};

struct CXInteraction {
  double eProj;
  double eCMS;
  double eProd;
  int idProj;
  int idTarg;
  int mult;
  int interactionCounter;
};


extern "C" {

  void presetSeed(int seed[3]);
  void finishListTree(const int mode);
  void initListTree(const int i, 
		    const int mode); 
  void analyseInteractions();
}

extern CXParticle gCXParticle;
extern NXParticle gNXParticle;
extern CXParticle gCXProjectile;
extern COSSINS gCXGeometry;
extern CXInteraction gCXInteraction;
extern int gInteractionCounter;
extern int gCXParticleCounter;
extern int gCXInteractionCounter;

extern int gCXSeed1, gCXSeed2, gCXSeed3;

extern TFile *gListTreeFile;
extern TTree *gListTree;
extern TTree *gProjectileTree;
extern TTree *gInteractionTree;
extern TTree *gSeedTree;

extern int gParticleListReadInteraction;



#endif

#endif // include protection

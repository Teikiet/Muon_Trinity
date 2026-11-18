#ifndef __include_CxRoot_h__
#define __include_CxRoot_h__

#include <conexConfig.h>
#include <conexHEModels.h>

#include <string>

#include <TStopwatch.h>


class TFile;
class TTree;


class ConexDynamicInterface;


class CxRoot {

public:
  CxRoot();
  ~CxRoot();
  bool init(int, char**);
  void run();

private:
  bool GetOptions(int, char**);
  void initFunctionPointers(EHEModel model);
  double GetEnergy(double) const;
  double GetTheta(double) const;
  double GetPhi(double) const;
  double GetImpactParameter(double);
  void rootOut(int, double e=0.);

  void interpolateProfile(const int nX, const float* X,
			  const float* N, float& Xmx, float& Nmx);
  bool invert3by3(double A[3][3]) const;
  bool solve3by3(double y[3], double A[3][3], double x[3]) const;

private:
  int fAutoSave;
  int fNShower;
  int fHadron;
  int fSeed;
  int fSeedcx[3];
  int fParticleID;
  double fAlpha;
  EHEModel fHEModel;
  double flgEmin;
  double flgEmax;
  double fTheta1;
  double fTheta2;
  double fMinImpactParameter;
  double fMaxImpactParameter;
  double fCurrImpactParameter;
  double fPhi;
  std::string fPrefix;

  int fMaxDetail;   // maximum number of interaction to be saved

  TStopwatch fStopwatch;
  double fTimeTotal;
  float fTimeSec;

  ConexDynamicInterface* fConexInterface;


  // ---- rootOut variables ----

  float fOutputVersion;

  TFile* fFile;
  TTree* fShower;
#ifdef __MC3D__
  TTree* fShower3D;
#endif

  static const int maxX = 20000;

  float X[maxX], H[maxX], D[maxX], N[maxX], dEdX[maxX],
    Mu[maxX], dMu[maxX], Gamma[maxX], Electrons[maxX], Hadrons[maxX],
    EGround[3], fitpars[13], currlgE, Xmx, Nmx, XmxdEdX, dEdXmx ;
  int nX;

#ifdef __MC3D__
  float SD1000[5], MuMIA, MuTr, HaEm, HaEs ;
#endif

#ifdef LEADING_INTERACTIONS_TREE
  TTree* fLeadingInteractions; // RU ++ Wed Jul  5 16:08:46 CEST 2006

  int gnInt1;           // number of interactions
  static const int maxNInt = 100;
  double gkinel1[maxNInt];      // inelasticity
  int gpId1[maxNInt];           // parent id
  double gpEnergy1[maxNInt];    // parent energy
  int gMult1[maxNInt];          // multiplicity
  int gMatg1[maxNInt];          // target nucleus mass
  double gDepth1[maxNInt];      // depth
  double gHeight1[maxNInt];     // height

  int gNParticles1;
  static const int maxNPart = 100000;
  double gEnergy1[maxNPart];    // size = 75000 from conex.incnex
  double gpx1[maxNPart];
  double gpy1[maxNPart];
  double gpz1[maxNPart];
  int gType1[maxNPart];
  int gIdInt1[maxNPart];
#endif

};

#endif

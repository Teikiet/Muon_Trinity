#ifndef _include_ConexDynamicInterface_h_
#define _include_ConexDynamicInterface_h_

#include <conexConfig.h>
#include <conexHEModels.h>

#include <string>
#include <vector>

class LeadingInteractionsData;
class LeadingInteractionsParticleData;

class ConexDynamicInterface {

  ConexDynamicInterface();
  ConexDynamicInterface(const ConexDynamicInterface&);
  ConexDynamicInterface& operator=(const ConexDynamicInterface&);

public:
  ConexDynamicInterface(const EHEModel model);
  ~ConexDynamicInterface();

  EHEModel GetHEModel() const { return fModel; }
  std::string GetHEModelName() const;

  // generic library symbol fetcher
  void* GetSymbol(const std::string& name) const;

  void Init(const int nShower, const int seeds[3],
	    const int maxDetail, const int particleListMode, const std::string& configPath);

  // function pointers to CONEX subroutines:
  
protected:
  // it is not recommended to directly access InitConex. Use Init instead.
  void (*InitConex)(int&, int*, EHEModel&, int&, 
#ifdef CONEX_EXTENSIONS
		    int&,
#endif
		    const char*, int);

public:
  void (*GetSeed)(int*);

  void (*EndMod)(int);

  void (*RunConex)(int*, double&, double&, double&, double&, int&);

  void (*GetHeaderVariables)(int&, float&, float&, int&, float&, int&, float&,  
                             int&, float&, int& , float&, float&, float&,
                             float&, float&, int&, float&);
#ifdef __MC3D__
  void (*GetHeaderVariables3D)(float&, float&, float&, float&, float&,
                               float&, float&, int&, float&, float&);
  void (*Get3DOutput)(int&, float&, float&, int&, float&, float&, float&);
#endif

  void (*GetShowerData)(int&, int&, int&, float&, float&,
			float&, float&, float&);
  void (*GetdEdXProfile)(int&, int&, float&, float&);
  void (*GetMuonProfile)(int&, int&, float&, float&);
  void (*GetGammaProfile)(int&, int&, float&);
  void (*GetElectronProfile)(int &, int&, float&);
  void (*GetHadronProfile)(int &, int&, float&);
  int (*GetNumberOfDepthBins)();

  double (*ConexRandom)(double& dummy);

  double (*InteractionLength)(int&, double&, double&);
  void (*ParticleMass)(int&, double&);

#ifdef __SIBYLL21__
  // prepare call to initialize (modified) nucleus-air cross sections (VERY SLOW)
  void (*DoSignucIniIni)(const int& init);
#endif      

#ifdef CONEX_EXTENSIONS
  void (*initListTree)(const int i, const int mode);
  void (*finishListTree)(const int mode); 
  void (*presetSeed)(int seed[3]);
#endif

#ifdef __CORSIKA8__
  void (*ConexRun)(int& ipart, double& energy, double& theta, double& phi, double& dimpact,
		   int ioseed[3]);
  void (*HadronCascade)(int&, int&, int&, int&);
  void (*SolveMomentEquations)(int&);
  void (*Shower)(int& iqi, double& ei, double& xmi, double& ymi, double& zmi, double& dmi,
                 double& xi, double& yi, double& zi, double& tmi, double& ui, double& vi,
                 double& wi, int& iri, double& wti, int& latchi);

#endif
  
  const std::vector<LeadingInteractionsData>& GetLeadingInteractionsData() const 
  { return *fInteractionsDataPtr; }
  
  const std::vector<LeadingInteractionsParticleData>& GetLeadingInteractionsParticleData() const 
  { return *fInteractionsParticleDataPtr; }
  

private:
  void initFunctionPointers(EHEModel model);  

protected:
  EHEModel fModel;

private:
  void* fConexLib;
  std::vector<LeadingInteractionsData>* fInteractionsDataPtr;
  std::vector<LeadingInteractionsParticleData>* fInteractionsParticleDataPtr;

};

#endif

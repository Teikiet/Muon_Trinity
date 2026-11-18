#include <ConexDynamicInterface.h>
#include <CxOutputPrefix.h>
#include <conexConfig.h>

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <dlfcn.h>

using namespace std;

class TTree;

#ifdef CONEX_EXTENSIONS
// the full relative path is needed here to allow compilation
// of ConexShowerGenerator module of Auger Offline
#include <../resampler/ResamplingMode.h> 

// RU: to change the CrossSection by a given factor, etc.
double gFactorCrossSection = 1;
double gFactorExtraMeson = 1;
double gFactorResampling = 1;
double gResamplingThreshold = 15-9;

resample::EResampleMode gResamplingMode = resample::eNone;

int gParticleListMode = 0;
std::string gParticleListModeFile = "unknown.dat";

#endif // extensions



ConexDynamicInterface::ConexDynamicInterface(const EHEModel model)
  : fModel(model),
    fConexLib(0) {
  initFunctionPointers(model);
}


ConexDynamicInterface::~ConexDynamicInterface()
{
  if (fConexLib)
    dlclose(fConexLib);
}


std::string
ConexDynamicInterface::GetHEModelName()
  const
{
  switch(fModel) {
  case eDPMJetIII:
    return "DPMJetIII.2017-1";
    break;
  case eSibyll23:
    return "Sibyll2.3e";
    break;
  case eQGSJet01:
    return "QGSJet01";
    break;
  case eQGSJetII:
    return "QGSJetII-04";
    break;
  case eQGSJetIII:
    return "QGSJetIII-01";
    break;
  case eEposLHC:
    return "EposLHC-R";
    break;
  default:
    return "unknonw";
  }
}

void
ConexDynamicInterface::Init(const int vnShower, const int vseeds[3], 
			    const int vmaxDetail, const int vparticleListMode,
			    const std::string& vconfigPath)
{
  int nShower = vnShower;
  int randomSeeds[3] = {vseeds[0], vseeds[1], vseeds[2]};
  EHEModel model = fModel;
  int maxDetail = vmaxDetail;
  int particleListMode = vparticleListMode;
  std::string configPath = vconfigPath;

  if (configPath.empty())
    configPath = CONEX_CONFIG_PATH;
  
#ifdef CONEX_EXTENSIONS
  this->InitConex(nShower,
		  randomSeeds,
		  model,
		  maxDetail,
		  particleListMode,
		  configPath.c_str(),
		  configPath.size());
#else
  particleListMode = particleListMode;   //to avoid compiler problem with [-Werror=unused-parameter] 
  this->InitConex(nShower,
		  randomSeeds,
		  model,
		  maxDetail,
		  configPath.c_str(),
		  configPath.size());
#endif
}


/*************************************************************
 *
 *************************************************************/
void
ConexDynamicInterface::initFunctionPointers(EHEModel model)
{
#ifdef _CONEX_PREFIX
  const string conex_prefix = "" _CONEX_PREFIX;
#else
  const string conex_prefix("");
#endif
#ifdef _CONEX_SYSTEM
  const string conex_system = "" _CONEX_SYSTEM;
#else
  const string conex_system("");
#endif  
  string sharedLibExt = ".so";
  if ( conex_system == "Darwin" )
    sharedLibExt = ".dynlib";

  ostringstream conexLibNameS;
  conexLibNameS << conex_prefix << "/lib/" << conex_system << "/";
  if ( model == eQGSJet01)
    conexLibNameS << "libCONEXqgsjet" << sharedLibExt;
  if ( model == eSibyll23)
    conexLibNameS << "libCONEXsibyll" << sharedLibExt;
  else if ( model == eDPMJetIII )
    conexLibNameS << "libCONEXdpmjetIII" << sharedLibExt;
  else if ( model == eQGSJetII )
    conexLibNameS << "libCONEXqgsjetII" << sharedLibExt;
  else if ( model == eQGSJetIII )
    conexLibNameS << "libCONEXqgsjetIII" << sharedLibExt;
  else if ( model == eEposLHC )
    conexLibNameS << "libCONEXepos" << sharedLibExt;

  const string conexLibName = conexLibNameS.str();

  cout << fgOutputPrefix << "loading library"
       << "\n" << setw(fgOutputPrefix.size()) << "" << "\'"
       << conexLibName << "\'" << endl;

  fConexLib = dlopen(conexLibName.c_str(), RTLD_LAZY);

  if ( ! fConexLib ) {
    ostringstream errMsg;
    errMsg << "\n can not open CONEX shared library " << conexLibName
	   << " at path \'" << conex_prefix << "/lib/" << conex_system << "\'\n\n"
	   << " Dynamic-link error:\n \"" << dlerror() << "\"\n";

    cout << errMsg.str() << endl;
    exit(1);
  }

  bool failed = false;

  fInteractionsDataPtr = (vector<LeadingInteractionsData>*) dlsym(fConexLib, "gInteractionData");
  fInteractionsParticleDataPtr = (vector<LeadingInteractionsParticleData>*) 
    dlsym(fConexLib, "gInteractionParticleData");

  InitConex =  (void(*)(int&, int*, EHEModel&, int&, 
#ifdef CONEX_EXTENSIONS
			int&,
#endif
			const char*, int)) dlsym(fConexLib, "initconex_");
  if (!InitConex) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function initconex_"
	 << endl;
    failed = true;
  }

  GetSeed  = (void (*)(int*)) dlsym(fConexLib, "getseed_");
  if (!GetSeed) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function getseed_" << endl;
    failed = true;
  }

  EndMod  = (void (*)(int)) dlsym(fConexLib, "endmod_");
  if (!EndMod) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function endmod_" << endl;
    failed = true;
  }

  RunConex  = (void (*)(int*, double&, double&, double&, double&, int&)) dlsym(fConexLib, "runconex_");
  if (!RunConex) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function runconex_" << endl;
    failed = true;
  }

  GetHeaderVariables = (void (*)(int&, float&, float&, int&, float&, int&, float&,  
                                 int&, float&, int& , float&, float&, float&,
                                 float&, float&, int&, float&)) dlsym(fConexLib, "get_data_card1_");
  if (!GetHeaderVariables) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_data_card1_"
	 << endl;
    failed = true;
  }

#if __MC3D__
  GetHeaderVariables3D = (void (*)(float&, float&, float&, float&, float&,
                                   float&, float&, int&, float&, float&)) dlsym(fConexLib, "get_data_card2_");
  if (!GetHeaderVariables3D) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_data_card2_"
	 << endl;
    failed = true;
  }
  
  Get3DOutput = (void (*)(int&, float&, float&, int&, float&, float&, float&)) dlsym(fConexLib, "get_shower_lat_");
  if (!Get3DOutput) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_shower_lat_"
	 << endl;
    failed = true;
  }
  
#endif

  GetShowerData =  (void (*)(int&, int&, int&, float&, float&,
                             float&, float&, float&)) dlsym(fConexLib, "get_shower_data_");
  if (!GetShowerData) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_shower_data_"
	 << endl;
    failed = true;
  }

  GetdEdXProfile = (void (*)(int&, int&, float&, float&)) dlsym(fConexLib, "get_shower_edep_");
  if (!GetdEdXProfile) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_shower_edep_"
	 << endl;
    failed = true;
  }

  GetMuonProfile = (void (*)(int&, int&, float&, float&))dlsym(fConexLib, "get_shower_muon_");
  if (!GetMuonProfile) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_shower_muon_"
	 << endl;
    failed = true;
  }

  GetGammaProfile = (void (*)(int&, int&, float&)) dlsym(fConexLib, "get_shower_gamma_");
  if (!GetGammaProfile) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_shower_gamma_"
	 << endl;
    failed = true;
  }

  GetElectronProfile = (void (*)(int &, int&, float&)) dlsym(fConexLib, "get_shower_electron_");
  if (!GetElectronProfile) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_shower_electron_"
	 << endl;
    failed = true;
  }

  GetHadronProfile = (void (*)(int &, int&, float&)) dlsym(fConexLib, "get_shower_hadron_");
  if (!GetHadronProfile) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_shower_hadron_"
	 << endl;
    failed = true;
  }

  GetNumberOfDepthBins = (int (*)()) dlsym(fConexLib, "get_number_of_depth_bins_");
  if (!GetNumberOfDepthBins) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function get_number_of_depth_bins_"
	 << endl;
    failed = true;
  }

  ConexRandom = (double (*)(double&)) dlsym(fConexLib, "drangen_");
  if (!ConexRandom) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function drangen_" << endl;
    failed = true;
  }

  InteractionLength = (double(*)(int&,double&,double&)) dlsym(fConexLib, "rlam_");
  if (!InteractionLength) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function rlam_"
	 << endl;
    failed = true;
  }

  ParticleMass = (void(*)(int&,double&)) dlsym(fConexLib, "cxidmass_");
  if (!ParticleMass) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX function cxidmass_"
	 << endl;
    failed = true;
  }

#ifdef __SIBYLL21__
  DoSignucIniIni = (void(*)(const int& init)) dlsym(fConexLib, "do_signuc_ini_ini_");
  if (!DoSignucIniIni) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find SIBYLL function do_signuc_ini_ini_"
	 << endl;
    failed = true;
  }
#endif

  
#ifdef CONEX_EXTENSIONS
  initListTree = (void(*)(const int, const int)) dlsym(fConexLib, "initListTree");
  if (!initListTree) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find SIBYLL function initListTree"
	 << endl;
    failed = true;
  }
  
  finishListTree = (void(*)(const int)) dlsym(fConexLib, "finishListTree");
  if (!finishListTree) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find SIBYLL function finishListTree"
	 << endl;
    failed = true;
  }
  
  presetSeed = (void(*)(int seed[3])) dlsym(fConexLib, "presetSeed");
  if (!presetSeed) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find SIBYLL function presetSeed"
	 << endl;
    failed = true;
  }
#endif

#ifdef __CORSIKA8__
  ConexRun = (void(*)(int& ipart, double& energy, double& theta, double& phi, double& dimpact,
		      int ioseed[3])) dlsym(fConexLib, "conexrun_");
  HadronCascade = (void(*)(int&, int&, int&, int&)) dlsym(fConexLib, "hadroncascade_");
  SolveMomentEquations = (void(*)(int&)) dlsym(fConexLib, "solvemomentequations_");
  Shower = (void(*)(int& iqi, double& ei, double& xmi, double& ymi, double& zmi, double& dmi,
		    double& xi, double& yi, double& zi, double& tmi, double& ui, double& vi,
		    double& wi, int& iri, double& wti, int& latchi)) = dlsym(fConexLib, "show_");
#endif

  
  if (failed) {
    cout << fgOutputPrefix << " ConexDynamicInterface failed!" << endl;
    exit(1);
  }
}

void* 
ConexDynamicInterface::GetSymbol(const string& name) 
  const
{
  void* ptr = (void*) dlsym(fConexLib, name.c_str());
  if (!ptr) {
    cout << fgOutputPrefix << " ConexDynamicInterface can not find CONEX symbol: \'" 
	 << name << "\";" 
	 << endl;
  }
  return ptr;
}

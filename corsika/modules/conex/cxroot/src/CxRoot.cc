/*************************************************************
 *
 *            CxRoot.cc
 *
 *  C++ interface to conex program
 *
 *
 *
 *  showers are diced with power law dN/dE ~ E^(-alpha)
 *
 *
 *  - command line options:
 *
 *    conex -s [random seed] -a [alpha] -e [log10(emin/eV)]  -E [log10(emax/eV)]
 *          -z [zenith angle]                                    RU 24.11.04
 *          -i [azimuth angle]                                   TP 02.09.05
 *          -n [nShower] -p [0=gam,100=p,5600=Fe,...]
 *          -m [2=QGSJET01, 4=EPOS LHC-R, 5=Sibyll 2.1, 6=QGSJETII-04, 10=QGSJETIII-04]
 *
 *  - steering files:
 *
 *    by default, the steering files
 *         conex_qgsjetIII.param or conex_qgsjetII.param or conex_epos.param or
 *         conex_qgsjet.param  or conex_sibyll.param
 *    are supposed to reside in $PWD/cfg. Alternatively the
 *    enviromental variable CONEX_ROOT can be set
 *
 *  - ROOT output:
 *
 *    by default, the ROOT output will be written to $PWD.
 *    Alternatively the enviromental variable ROOT_OUT can be set
 *
 *
 *  v1.0 08/07/03 Michael.Unger@ik.fzk.de
 *
 *  18/08/04 T.P.:
 *       - added dEdX Branch
 *
 *  08/09/04 M.U.:
 *
 *       - removed hard coded paths
 *         (use $CONEX_ROOT and $ROOT_OUT instead)
 *       - added interaction model switch
 *       - added some Header variables
 *
 *  21/01/05 F.S.:
 *        - added Emax as input option
 *
 *  24/01/05 M.U.:
 *        - added energy at ground level (EGround)
 *        - CPU timing
 *
 *  10/02/05 R.U.:
 *        - added zenith angle as input option
 *  10/02/05 T.P.:
 *        - update header information
 *        - change model identifier to match with Conex definition :
 *          1 - neXus
 *          2 - QGSJet
 *          3 - Gheisha
 *          4 - Sibyll
 *        - added muon profile
 *        - added seed information
 *  29/03/05 T.P.:
 *        - change number of decimal for energy printing
 *  31/03/05 T.P.:
 *        - change primary particle code to A*100 or PDG code
 *  01/04/05 T.P.:
 *        In Header :
 *        - add corresponding bin altitude above see level (X->H(m))
 *        - add corresponding bin distance from impact point (or parameter)
 *          (X->D(m))
 *          (OutputVersion -> 1.1)
 *  04/04/05 T.P.:
 *        - add an optional prefix to the file name
 *  22/06/05 T.P.:
 *        - add the muon production rate profile in the root file
 *          (OutputVersion -> 1.1)
 *
 *  29/06/05 RU:
 *        - added the electron/gamma profiles to root file
 *          (OutputVersion -> 1.1)
 *
 *  02/09/05 TP:
 *        - added QGSJETII
 *        - added angle phi as input and in the root file
 *          (OutputVersion -> 1.2)
 *
 *  06/09/05 RU:
 *        - added Inelasticity of 1st Interaction
 *          (OutputVersion -> 1.3)
 *
 *
 *  19/01/06 MU:
 *        - added zenith angle generation
 *          (OutputVersion -> 1.4)
 *
 *  05/05/06 TP:
 *        - update QGSJETII-2 to QGSJETII-3
 *        - new Makefile
 *
 *  09/06/06 TP:
 *        - change model identifier to match with NEXUS and CONEX
 *          definition :
 *          1 - neXus
 *          2 - QGSJet
 *          3 - Gheisha
 *          4 - EPOS
 *          5 - Sibyll
 *          6 - QGSJETII
 *
 *  05/07/06 RU:
 *        - added additional TTree for output,
 *          to store all secondaries of the first interaction
 *          USE PRE-COMPILER SWITCH IN MAKEFILE TO ACTIVATE !!!!!
 *
 *  25/07/06 TP:
 *        - Add EPOS
 *
 *  14/08/06 TP:
 *        - Add _EPOS preprocessing control
 *
 *  22/10/06 MU:
 *        - sort FirstInteraction Tree
 *        - added interpolated Xmax and Nmax
 *          (OutputVersion -> 1.5)
 *
 *  13/04/07 TP:
 *        - add momentum and height in First Interaction Tree
 *
 *  7/5/07 RU:
 *        - made momentum in First-Interaction-Tree optional
 *          you have to enable 'FIRST_INTERACTION_TREE_EXT'
 *        - changed help output from char[] to ostringstream
 *
 *  9/5/07 RU:
 *        - removed "trigger surface"
 *        ->  (OutputVersion -> 1.6)
 *
 *  7/8/08 SM+TP:
 *        - move Height and distance table in Shower Tree (because zenith
 *          angle dependence)
 *        ->  (OutputVersion -> 2.0)
 *
 *  7/8/08 SM:
 *        - Move X/slant depth table from Header to Shower Tree because
 *          the vector can differ in length in merged ROOT files!
 *        - Remove delX and nX from Header.
 *        -> Still OutputVersion 2.0!
 *
 *  19/9/08 RU:
 *       - Use of TTree auto-save concept in order to recover crashed run
 *         command line option -S is used to specify the save intervall
 *
 *  20/03/09 RU;
 *       - Added svn revision information to conex ROOT output
 *
 *  16/04/09 RU:
 *       - Added Xmax_dEdX and dEdXmax
 *       -> OutputVersion 2.1
 *
 *  12/04/10 TP:
 *       - Update input file handling with user defined starting slant depth
 *       -> OutputVersion 2.2 (change only if input file is used)
 *
 *  23/04/10 RU:
 *       - Added time per event to Shower()
 *       -> OutputVersion 2.3
 *
 *  04/05/10 MU:
 *       - refactor
 *       - move interaction TTree #defines to conexConfig.h
 *
 *  7/5/10 RU:
 *       - added random seed initialization from urandom in C++
 *
 *  8/5/10 RU/MU:
 *       - Every single HE model gets its own conex library
 *       - now a single cxroot binary can dynamically load
 *         conex libraries needed for the simulation
 *
 *  9/5/10 RU:
 *       - moved main part of code into new CxRoot class
 *         -> get rid of most global + static variables
 *
 * 17/6/10 RU:
 *       - renamed first-interaction tree into leading-interactions tree,
 *         to prevent confusion about its data content
 *
 * 16/12/10 Colin Baus
 *       -> OutputVersion 2.4 (since first interaction tree was renamed ...)
 *
 * 05/12/11 MU:
 *       -> add interaction length tables to Header (--> OutputVersion 2.41)
 *
 * 20/01/12 MU:
 *       -> add impact parameter to Shower and command line options
 *           (--> OutputVersion 2.42)
 *
 * 11/10/12 MU:
 *       -> bug fix, move XmxdEdX by half an X-bin
 *           (--> OutputVersion 2.43)
 *
 *  24/01/12 TP:
 *        - Remove NEXUS from model list and add UrQMD
 *          definition :
 *          2 - QGSJet
 *          3 - Gheisha
 *          4 - EPOS LHC
 *          5 - Sibyll 2.1
 *          6 - QGSJETII-04
 *          8 - UrQMD 1.3.1
 *        - Update the code to conex 3D but used in 1D
 *
 *  31/01/12 TP :
 *       -> Save Header at the beginning to have valid file even if run crashed
 *
 *  14/02/13 TP:
 *       - ADD S1000 (VEM) value for Auger SD signal estimator
 *         (total, and photon, electron, muon and hadron component)
 *       - add Truncated number of muons
 *       - add hadron profile and sum of hadron energy
 *         and maximum energy of the hadrons
 *       - uses emcut1 and hacut1 for simu (as before)
 *         but  emcut2 for gamma profile
 *              emcut3 for electron profile
 *              hacut2 for muon profiles
 *              hacut3 for hadron profile
 *         emcut2 and emcut3 do not have to be ordered but >= emcut1
 *         hacut2 and hacut3 do not have to be ordered but >= hacut1
 *       - add muon density of MIA experiment (mu threshold=0.85/cos(the))
 *           (--> OutputVersion 2.50 in 1D and 3.0 in 3D)
 *
 *  21/02/15 RU:
 *       - Added resampling information to Header output
 *       -> OutputVersion 2.51
 *
 *  19/12/15 TP and FR:
 *       - add Sibyll 2.3
 *  28/04/16 TP and FR :
 *       - update Sibyll 2.3 as released in CORSIKA
 *       - prepare release 5.3
 *  13/06/17 TP and FR :
 *       - new Sibyll 2.3c
 *  07/07/17 TP and AF :
 *       - add DPMJETIII
 *          9 - DPMJET III.2017-1
 *  08/10/19 TP :
 *       - add  option -H to run without the electromagnetic component
 *       - correct zenith angle distribution for upward going showers
 *        (from Massimo Mastrodicasa (Aquila, INFN) )
 *  20/05/20 TP and FR :
 *       - new Sibyll 2.3d
 *  09/03/20 RU;
 *       - removed svn revision information to conex ROOT output
 *  02/02/25 TP:
 *       - introduce EPOS LHC-R and QGSJETIII-01
 *  12/02/25 TP:
 *       - update Sibyll to version 2.3e (phi asymmetry corrected)
 **************************************************************/

/*************************************************************
 * conex CVS-branch: 'first_int_tree'
 *
 * Mon May  7 09:15:16 CEST 2007 RU:
 *       - introduced particle list model 'L' :
 *         reads/writes from extra ROOT file specified by 'F'
 *         L=0 OFF
 *         L=1 reads from 'F' and start shower with particle list
 *         L>1 writes to 'F' (all particles of generation L-2)
 *
 *       - new option 'X': modifyer of cross-section
 *         (see ... for details)
 *
 *       - option 'M': modifyer of multiplicity
 *         (see ... for details)
 *
 * 2009: renamed branch from 'first_int_tree' into 'resampling' on svn server
 *
 *************************************************************/

#include <CxRoot.h>
#include <CxOutputPrefix.h>

#include <conexConfig.h>
#include <conexFunctions.h>
#include <ConexDynamicInterface.h>
#ifdef LEADING_INTERACTIONS_TREE
#include <leadingInteractionsData.h>
#endif

#ifdef CONEX_EXTENSIONS
#include <conexExtensions.h>
#include <conexResampling.h>
#include <ResamplingMode.h>
#endif

#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>
#include <TDirectory.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <unistd.h>
using namespace std;

CxRoot::CxRoot()
  : fAutoSave(10),
    fNShower(1),
    fSeed(0),
    fParticleID(100),
    fAlpha(3.),
    fHEModel(eEposLHC),
    flgEmin(16.5),
    flgEmax(21.),
    fTheta1(60.),
    fTheta2(60.),
    fMinImpactParameter(0),
    fMaxImpactParameter(0),
    fPhi(0.),
    fHadron(0),
    fPrefix("conex"),
    fMaxDetail(0),   // maximum number of interaction to be saved
    fConexInterface(0),
#ifdef __MC3D__
    fOutputVersion(3.00)
#else
    fOutputVersion(2.51) // the version of the root-output
#endif
{
  fSeedcx[0] = 0;
  fSeedcx[1] = 0;
  fSeedcx[2] = 0;
}

CxRoot::~CxRoot()
{
  delete fConexInterface;
}

bool
CxRoot::init(int argc, char** argv)
{
  // read options from comand line
  if (!GetOptions(argc, argv))
    return false;

#ifdef CONEX_EXTENSIONS
  cout << fgOutputPrefix << "Resampling: f(CX)=" << gFactorCrossSection;
  if (gFactorExtraMeson!=1)
    cout << ", f(meson)=" << gFactorExtraMeson;
  cout << ", f(resample)=" << gFactorResampling << ", mode=";
  switch(gResamplingMode) {
  case resample::eMultiplicity:
    cout << "MULT";
    break;
  case resample::eElasticity:
    cout << "ELA";
    break;
  case resample::eEMRatio:
    cout << "EM-Ratio";
    break;
  case resample::eChargeRatio:
    cout << "Charge-Ratio";
    break;
  case resample::ePi0Spectrum:
    cout << "Pi0-Spectrum";
    break;
  case resample::eLeadingRho0:
    cout << "leading Rho0/Pi0";
    break;
  case resample::eBaryonProd:
    cout << "Baryon production";
    break;
  case resample::eNone:
  default:
    cout << "NONE";
    break;
  }
  cout << ", lgE_threshold/eV=" << gResamplingThreshold+9 << " " << endl;
  
  if (gParticleListMode!=0) {
    cout << fgOutputPrefix
	 << "Particle-List-Mode: \'";
    if (gParticleListMode<0) {
      cout << "read";
    } else {
      cout << "write generation " << gParticleListMode;
    }
    cout << "\', file: \'" << gParticleListModeFile << "\'"
	 << endl;
  }
#endif // end extensions


  const string conexRoot = (gSystem->Getenv("CONEX_ROOT") ?
			    gSystem->Getenv("CONEX_ROOT") :
			    string(gSystem->Getenv("PWD"))+string("/cfg"));
  //  const string dataPath = (gSystem->Getenv("CORSIKA_DATA") ?
  // gSystem->Getenv("CORSIKA_DATA") : "");
  cout << fgOutputPrefix << "Using input steering files located at"
       << "\n" << setw(fgOutputPrefix.size()) << "" << "\'"
       << conexRoot << "\'";
  cout << endl;

  const int crLength = conexRoot.length();
  fSeedcx[0] = fSeed;
  
  fConexInterface = new ConexDynamicInterface(fHEModel);
  
  int targetCode = 0; // air
  
#ifdef CONEX_EXTENSIONS
  // ru-fudge, to get seed of shower0:
  if (gParticleListMode<0) { // 1: read, >1 write
    fConexInterface->initListTree(0, gParticleListMode);
    fConexInterface->presetSeed(fSeedcx);
  }
#else
  const int gParticleListMode = 0;
#endif
  int NShow=fNShower*(1-2*fHadron);
  fConexInterface->Init(NShow, fSeedcx, fMaxDetail,
		        gParticleListMode,
                        conexRoot);

  // initialize ROOT output file
  rootOut(-1);

  return true;

}

void
CxRoot::run()
{
  fStopwatch.Stop();
  fTimeTotal = fStopwatch.CpuTime();
  fStopwatch.Start(false);
  fTimeSec = 0;

  // loop over showers
  for (int iShow=0; iShow<fNShower; ++iShow) {

    fConexInterface->GetSeed(fSeedcx);

    double energy = GetEnergy(iShow);
    double theta = GetTheta(iShow);
    double azimuth = GetPhi(iShow);
    double impact = GetImpactParameter(iShow);

    cout << "\n" << fgOutputPrefix << "generating shower number "
         << iShow+1 << "\n"
         << setw(fgOutputPrefix.size()) << ""
         << "lg(E/eV) = "
         << log10(energy)+9.
         << ", zenith = " << theta
	 << ", azimuth = " << azimuth << "\n"
         << setw(fgOutputPrefix.size()) << ""
	 << "impact parameter = " << impact/1000 << " km\n"
         << setw(fgOutputPrefix.size()) << ""
         << "random seeds: (" << fSeedcx[0] << ", "
         << fSeedcx[1] << ", " << fSeedcx[2] << ")"
         << endl;

#ifdef CONEX_EXTENSIONS
    fConexInterface->initListTree(iShow, gParticleListMode); // 1: read, >1 write
#endif
    fConexInterface->RunConex(fSeedcx, energy, theta, azimuth, impact, fParticleID);
#ifdef CONEX_EXTENSIONS
    fConexInterface->finishListTree(gParticleListMode);
#endif

    rootOut(0, energy);

  }

  rootOut(1);

  //  if (fHEModel==4)fConexInterface->EndMod(4);
  
  cout << "\n" << fgOutputPrefix << "successfully processed "
       << fNShower << " shower"
       << (fNShower > 1 ? "s \n" : " \n");

  fStopwatch.Stop();
  const double realTime = fStopwatch.RealTime();
  const double cpuTime = fStopwatch.CpuTime();

  cout << setw(fgOutputPrefix.size()) << ""
       << "in " << realTime << " sec, "
       << ((double)realTime/fNShower) << " sec/showers, with "
       << (cpuTime/realTime*100) << "% cpu usage. \n" << endl;
}


/*************************************************************
 *
 *  GetEnergy()  - returns energy in GeV distributed
 *                 according to power law
 *
 **************************************************************/
double
CxRoot::GetEnergy(double dummy)
  const
{
  double energy;

  double emin = pow(10., flgEmin);
  double emax = pow(10., flgEmax);

  // dice an energy

  if (emax != emin) {
    energy = 1.1 * emax;
    while (energy > emax) {
      const double rdm = fConexInterface->ConexRandom(dummy);
      if (fAlpha >= 0 && fAlpha <= 1)                   // option for flat in lgE
        energy=pow(10., flgEmin+rdm*(flgEmax-flgEmin));
      else if (fAlpha<0.) {                        // aequidistant in lgE
        const double dlgE = TMath::Abs(fAlpha);
        const double addlgE = double(int(rdm*(1.+(flgEmax-flgEmin)/dlgE)))*dlgE;
        energy = pow(10., flgEmin+addlgE);
      }
      else
        energy = pow(1.-rdm, 1./(1.-fAlpha) ) * emin;
    }
  }
  else
    energy = emin;

  return energy/1.e9; // convert to GeV
}


double
CxRoot::GetImpactParameter(double i)
{
  fCurrImpactParameter = fMinImpactParameter == fMaxImpactParameter ?
    fMinImpactParameter :
    fMinImpactParameter +
    fConexInterface->ConexRandom(i)* (fMaxImpactParameter-fMinImpactParameter);
  return fCurrImpactParameter;
}

/*************************************************************
 *
 *   GetTheta(): dice a zenith angle from isotropic flux
 *               on flat surface:
 *
 *         dN/dcos(theta) ~ cos(theta)
 *
 *  returns theta [deg.]
 *
 **************************************************************/
double
CxRoot::GetTheta(double i)
  const
{
  if (fTheta1 == fTheta2)
    return fTheta1;
  else {
    const double degrad = 180./TMath::Pi();
    const double c1 = cos(fTheta1/degrad);
    const double c2 = cos(fTheta2/degrad);
    const double denom = (c1*c1+c2*c2);
    const double u =
      c2*c2/denom + fConexInterface->ConexRandom(i)*(c1*c1-c2*c2)/denom;
    //    const double cTheta=sqrt(u*(c1*c1+c2*c2));
    double cTheta;

    if (fTheta1<90.) {
      cTheta=sqrt(u*(c1*c1+c2*c2));
    }else if (fTheta1>=90.){
      cTheta=-sqrt(u*(c1*c1+c2*c2));
    }
    if (fTheta1<90. && fTheta2>90.) {
      cout << " Error - wrong zenith angle range:  theta1<90° and theta2>90°!\n"
          << " exit ...\n" << endl;
      exit(1);
    }
    return acos(cTheta)*degrad;
  }
}

/*************************************************************
 *
 *   GetPhi(): dice an azimuth angle between 0 and 360 deg
 *
 *  returns phi [deg.]
 *
 **************************************************************/
double
CxRoot::GetPhi(double i)
  const
{
  if (fPhi >= 0.)
    return fPhi;
  else
    return fConexInterface->ConexRandom(i)*360.;
}

/*************************************************************
 *
 *   GetOptions()  - reads command line options
 *
 **************************************************************/
bool
CxRoot::GetOptions(int argc, char** argv)
{
  int c;
  ostringstream help;
  help << " conex2r -s [random seed] -S [autosave range] -a [alpha]\n"
       << "         -e [log10(emin/eV)] -E [log10(emax/eV)]\n"
       << "         -z [min zenith angle/deg] -Z [max zenith angle/deg]\n"
       << "         -i [azimuth angle/deg] -n [nShower]\n"
       << "         -p [0=gamma,100=p,5600=Fe,...] -x [prefix]\n"
       << "         -o [min. impact [m]] -O [max. impact [m]]\n"
       << "         -H [1 = no electromagnetic shower component]\n"
       << "         -m [" << eQGSJet01 <<"=QGSJET01, "
       << eEposLHC << "=EPOS LHC-R, "
       << eSibyll23 << "=SIBYLL 2.3e, "
       << eDPMJetIII << "=DPMJET III.2017-1, "
       << eQGSJetII << "=QGSJETII-04, "
       << eQGSJetIII << "=QGSJETIII-01]\n";
#ifdef LEADING_INTERACTIONS_TREE
  help << "         -K [maxDetailLevel]";
#endif
#ifdef CONEX_EXTENSIONS
  help << " -X [crossSectionFactor] -P [mesonExtraFactor] -M [resamplingFactor] -L [prtcleListMode] -F [prtcleListFile]"
       << " -T [modifications-threshold]"
       << " -R [resampling mode]\n"
       << "\n"
       << "     prtcle-list-mode: i (i>0) - write particles of generation i, 0 - OFF, <0 - read from file \n"
       << "      resampling-mode:  1 - multiplicity, 2 - elasticity, 3 - EMratio, 4 - ChargeRatio, 5 - Pi0spec, 6 - LeadingRho0, 7 - BaryonProduction \n";
  gParticleListMode = 0;
  gFactorCrossSection = 1;
  gFactorExtraMeson = 1;
  gFactorResampling = 1;
  gResamplingThreshold = 15 - 9; // GeV
  gResamplingMode = resample::eNone;
#endif


  string options = "S:s:a:e:E:o:O:z:Z:i:n:p:m:x:H:h::";

#ifdef LEADING_INTERACTIONS_TREE
  options = "K:" + options;
#endif
#ifdef CONEX_EXTENSIONS
  options = "X:P:M:T:L:F:R:" + options;
#endif
  while ((c = getopt (argc, argv, options.c_str())) != -1) {
    switch (c) {
#ifdef LEADING_INTERACTIONS_TREE
    case 'K':  // option to set generations in first interaction tree
      fMaxDetail = atoi(optarg);
      break;
#endif
#ifdef CONEX_EXTENSIONS
    case 'X':
      gFactorCrossSection = atof(optarg);
      break;
    case 'P':
      gFactorExtraMeson = atof(optarg);
      break;
    case 'M':
      gFactorResampling = atof(optarg);
      break;
    case 'T':
      gResamplingThreshold = atof(optarg) - 9; // GeV
      break;
    case 'L':
      {
	gParticleListMode = atoi(optarg);
	/*
	*/
	break;
      }
    case 'F':
      gParticleListModeFile = string(optarg);
      break;
    case 'R':
      {
	const int flag = atoi(optarg);
	switch (flag) {
	case 0: gResamplingMode = resample::eNone; break;
	case 1: gResamplingMode = resample::eMultiplicity; break;
	case 2: gResamplingMode = resample::eElasticity; break;
	case 3: gResamplingMode = resample::eEMRatio; break;
	case 4: gResamplingMode = resample::eChargeRatio; break;
	case 5: gResamplingMode = resample::ePi0Spectrum; break;
	case 6: gResamplingMode = resample::eLeadingRho0; break;
	case 7: gResamplingMode = resample::eBaryonProd; break;
	default: cerr << " unkown resampling mode : " << flag << endl; exit(1); break;
	}
	break;
      }
#endif // CONEX_EXTENSIONS
    case 'S':
      fAutoSave=atoi(optarg);
      break;
    case 's':
      fSeed=atoi(optarg);
      break;
    case 'a':
      fAlpha=atof(optarg);
      break;
    case 'e':
      flgEmin=atof(optarg);
      break;
    case 'E':
      flgEmax=atof(optarg);
      break;
    case 'o':
      fMinImpactParameter = atof(optarg);
      break;
    case 'O':
      fMaxImpactParameter = atof(optarg);
      break;
    case 'z':
      fTheta1=atof(optarg);
      break;
    case 'Z':
      fTheta2=atof(optarg);
      break;
    case 'i':
      fPhi=atof(optarg);
      break;
    case 'n':
      fNShower=atoi(optarg);
      break;
    case 'p':
      fParticleID=atoi(optarg);
      break;
    case 'm':
      fHEModel = (EHEModel)atoi(optarg);
      if (fHEModel<1 || fHEModel==3 || fHEModel==7 || fHEModel==8 || fHEModel>10) {
        cout << " unknown HE interaction model " << fHEModel << endl;
        cout << help.str() << endl;
        return false;
      }
      break;
    case 'x':
      fPrefix=optarg;
      break;
    case 'H':
      fHadron=min(1,abs(atoi(optarg)));
      break;
    case 'h':
      cout << help.str() << endl;
      return false;
    default:
      cout << help.str() << endl;
      return false;
    }
  }

  // check if random seed was provided, otherwise generate one
  const bool seedProvided = fSeed;
  if (fSeed == 0) {
    ifstream urandom("/dev/urandom", ios::in|ios::binary);
    urandom.read((char*)&fSeed, sizeof(fSeed)/sizeof(char));
    urandom.close();
    fSeed = abs(fSeed);
  }

  cout << "\n          >> conex2r <<\n\n"
       << "  prefix           : " << fPrefix << "\n"
       << "  random seed      : " << fSeed;
  if (seedProvided)
    cout << " (provided by user) \n";
  else
    cout << " (read from /dev/urandom) \n";
  cout << "  alpha            : " << fAlpha << "\n"
       << "  lgE              : [" << flgEmin << ", "
       << flgEmax << "]\n"
       << "  theta            : [" << fTheta1 << ", "
       << fTheta2 << "] deg \n";
  cout << "  impact parameter : [" << fMinImpactParameter << ", "
       << fMaxImpactParameter << "] m\n";
  if (fPhi<0)
    cout  << "  phi drawn from   : [0,360] deg. \n";
  else
    cout  << "  phi fixed at     : " << fPhi << " deg \n";
  cout << "  number of showers: " << fNShower << "\n"
       << "  HE model         : " << fHEModel;
  if (fHEModel == eQGSJet01)
    cout << " (QGSJET01) ";
  else if (fHEModel == eQGSJetII)
    cout << " (QGSJETII-04) ";
  else if (fHEModel == eEposLHC)
    cout << " (EPOS LHC-R) ";
  else if (fHEModel == eSibyll23)
    cout << " (Sibyll 2.3e) ";
  else if (fHEModel == eDPMJetIII)
    cout << " (DPMJET III.2017-1) ";
  else if (fHEModel == eQGSJetIII)
    cout << " (QGSJETIII-01) ";
  else
    cout << " (unknown) ";
  cout << "\n  particle type    : " << fParticleID << "\n"
       << endl;

#ifdef LEADING_INTERACTIONS_TREE
  if (fMaxDetail > 0)
    cout << "  Save leading interactions up to generation : " << fMaxDetail << "\n"
	 << endl;
#endif

  if (fTheta1 > fTheta2) {
    cout << " Error - theta1>theta2!\n"
	 << " exit ...\n" << endl;
    exit(1);
  }

#ifdef CONEX_EXTENSIONS
  if (gParticleListModeFile=="" && gParticleListMode!=0) {
    cerr << " Error - specify particle list file for input/output" << endl
	 << " exit ...\n" << endl << endl;
    exit(1);
  }
#endif

  cout.setf(ios::showpoint);
  cout.setf(ios::fixed);
  cout.precision(2);

  return true;
}


/*************************************************************
 *
 *   rootOut(int iopt)
 *
 *   opens ROOT file (iopt=-1), stores the data (iopt=0)
 *   and closes ROOT file (iopt=1)
 *
 **************************************************************/
void
CxRoot::rootOut(int iopt, double e)
{
  int icut = 1;
  int icutg = 2;
  int icute = 3;
  int icutm = 2;
  int icuth = 3;
  int iSec = 0;

  if (iopt == -1) { // open file and define tree

    // check if ROOT_OUT is set (if not assume PWD)

    const char* rootOutDir = gSystem->Getenv("ROOT_OUT");
    if (rootOutDir == NULL)
      rootOutDir = gSystem->Getenv("PWD");

    // open ROOT output file

    ostringstream rootFileName;
    rootFileName << rootOutDir << "/" << fPrefix << "_";
    switch (fHEModel) {
    case eQGSJetII: rootFileName << "qgsjetII04"; break;
    case eQGSJet01: rootFileName << "qgsjet01";   break;
    case eEposLHC: rootFileName << "eposlhcr";     break;
    case eSibyll23: rootFileName << "sibyll23e";   break;
    case eDPMJetIII: rootFileName << "dpmjetIII";   break;
    case eQGSJetIII: rootFileName << "qgsjetIII";   break;
    default:
      cout << " rootOut: error - unknown model " << fHEModel << endl;
      cout << "          exit ..." << endl;
      exit(1);
      break;
    }
    rootFileName << "_" << setw(9) << setfill('0') << fSeed
                 << "_" << fParticleID
                 << ".root";

#ifdef DEBUG_RESAMPLING
    const string name = rootFileName.str().substr(0, rootFileName.str().rfind(".root"))
      + ".resampling.root";
    fResamplingDebug = TFile::Open(name.c_str(), "RECREATE");
    tResamplingDebug = new TTree("interaction", "interaction");
    tResamplingDebug->Branch("primary", &resampling_primary, "primary/I");
    tResamplingDebug->Branch("Elab", &resampling_Elab, "Elab/D");
    tResamplingDebug->Branch("factor", &resampling_factor, "factor/D");
    tResamplingDebug->Branch("yboost", &resampling_yboost, "yboost/D");
    tResamplingDebug->Branch("inucleon", &resampling_inucleon, "inucleon/I");
    tResamplingDebug->Branch("nnucleon", &resampling_nnucleon, "nnucleon/I");
    tResamplingDebug->Branch("isCM", &resampling_cms, "isCM/I");
    resampling_before = new TClonesArray("particle", 100);
    resampling_after = new TClonesArray("particle", 100);
    tResamplingDebug->Branch("before", "TClonesArray", &resampling_before);
    tResamplingDebug->Branch("after", "TClonesArray", &resampling_after);
    resampling_before_n = 0;
    resampling_after_n = 0;
#endif // DEBUG_RESAMPLING


    cout << fgOutputPrefix << "opening output file"
         << "\n" << setw(fgOutputPrefix.size()) << "" << "\'"
         << rootFileName.str() << "\'" << endl;
    fFile = new TFile(rootFileName.str().c_str(),"RECREATE");

    // file header
    int LEModel ;
    float HiLowEgy, hadCut, emCut, muCut, gaCut, elCut, haCut, hadThr, emThr, muThr, Version;
    fConexInterface->GetHeaderVariables(icut, hadCut, emCut,
                                        icutm, muCut, icutg, gaCut, icute, elCut,
                                        icuth, haCut, hadThr,
                                        muThr, emThr, HiLowEgy, LEModel, Version);

    TTree* Header = new TTree("Header", "run header");

    Header->Branch("Seed1", &fSeedcx[0], "Seed1/I");    // random seed1
    Header->Branch("Particle", &fParticleID, "Particle/I");   // particle ID
    Header->Branch("Alpha", &fAlpha, "Alpha/D");         // spectral index
    Header->Branch("lgEmin", &flgEmin, "lgEmin/D");      // minimum energy
    Header->Branch("lgEmax", &flgEmax, "lgEmax/D");      // maximum energy
    Header->Branch("zMin", &fTheta1, "zMin/D");          // min. zenith angle
    Header->Branch("zMax", &fTheta2, "zMax/D");          // min. zenith angle
    Header->Branch("Version", &Version, "Version/F");   // conex version
    Header->Branch("OutputVersion", &fOutputVersion,
                   "OutputVersion/F");         // CxRoot output version
    Header->Branch("HEModel", &fHEModel, "HEModel/I");   // HE model flag
    Header->Branch("LEModel", &LEModel, "LEModel/I");   // LE model flag
    Header->Branch("HiLowEgy", &HiLowEgy, "HiLowEgy/F");   // HE-LE model energy transition
    Header->Branch("hadCut", &hadCut, "hadCut/F");      // hadron/muon cutoff
    Header->Branch("emCut", &emCut, "emCut/F");         // emag cutoff
    Header->Branch("hadThr", &hadThr, "hadThr/F");      // hadron relative threshold
    Header->Branch("muThr", &muThr, "muThr/F");         // muon relative threshold
    Header->Branch("emThr", &emThr, "emThr/F");         // emag relative threshold
    Header->Branch("haCut",&haCut,"haCut/F");     // hadron cutoff for profile
    Header->Branch("muCut",&muCut,"muCut/F");     // muon cutoff for profile
    Header->Branch("elCut",&elCut,"elCut/F");     // electron cutoff for profile
    Header->Branch("gaCut",&gaCut,"gaCut/F");     // gamma cutoff for profile

#ifdef CONEX_EXTENSIONS
    Header->Branch("resamplingMode", &gResamplingMode, "resamplingMode/I");
    Header->Branch("modThreshold", &gResamplingThreshold, "modThreshold/D");
    Header->Branch("f19_cx", &gFactorCrossSection, "f19_cx/D");
    Header->Branch("f19_meson", &gFactorExtraMeson, "f19_meson/D");
    Header->Branch("f19", &gFactorResampling, "f19/D");
#endif

    // fill some hadron-air interaction lengths used in this calculation
    const double lgELamMin = 15.;
    const double lgELamMax = 21.;
    const double dlgELam = 0.2;
    const int nLam = int((lgELamMax-lgELamMin)/dlgELam+1);
    double lambdaLgE[nLam];
    double lambdaProton[nLam];
    double lambdaHelium[nLam];
    double lambdaNitrogen[nLam];
    double lambdaIron[nLam];
    double lambdaPion[nLam];

    ostringstream branchSuffix;
    branchSuffix << "[" << nLam << "]/D";
    Header->Branch("lambdaLgE", lambdaLgE,
                   (string("lambdaLgE")+branchSuffix.str()).c_str());
    Header->Branch("lambdaProton", lambdaProton,
                   (string("lambdaProton")+branchSuffix.str()).c_str());
    Header->Branch("lambdaPion", lambdaPion,
                   (string("lambdaPion")+branchSuffix.str()).c_str());
    Header->Branch("lambdaHelium", lambdaHelium,
                   (string("lambdaHelium")+branchSuffix.str()).c_str());
    Header->Branch("lambdaNitrogen", lambdaNitrogen,
                   (string("lambdaNitrogen")+branchSuffix.str()).c_str());
    Header->Branch("lambdaIron", lambdaIron,
                   (string("lambdaIron")+branchSuffix.str()).c_str());

    double lgE = lgELamMin;
    for (int i=0; i<nLam; ++i) {
      int particleId = 1;
      double energyPerNucleonInGeV = pow(10.,lgE)/1e9;
      double mass = 0.938;
      lambdaProton[i] =
        fConexInterface->InteractionLength(particleId, energyPerNucleonInGeV, mass);

      mass = 0.13957;
      particleId = 2;
      lambdaPion[i] =
        fConexInterface->InteractionLength(particleId, energyPerNucleonInGeV, mass);

      mass = 0.938;
      energyPerNucleonInGeV = pow(10.,lgE)/1e9/4;
      particleId = 40;
      lambdaHelium[i] =
        fConexInterface->InteractionLength(particleId, energyPerNucleonInGeV, mass);

      energyPerNucleonInGeV = pow(10.,lgE)/1e9/14;
      particleId = 140;
      lambdaNitrogen[i] =
        fConexInterface->InteractionLength(particleId, energyPerNucleonInGeV, mass);

      energyPerNucleonInGeV = pow(10.,lgE)/1e9/56;
      particleId = 560;
      lambdaIron[i] =
        fConexInterface->InteractionLength(particleId, energyPerNucleonInGeV, mass);

      lambdaLgE[i] = lgE;
      lgE += dlgELam;
    }

#ifdef __MC3D__

    float WHmax,WMmax,WTmax,hadLow,muLow,emLow,zLow,Thin,HGround;
    int iThin;

    fConexInterface->GetHeaderVariables3D(WHmax,WMmax,WTmax,hadLow,muLow,emLow,zLow,iThin,Thin,HGround);

    Header->Branch("WHmax",&WHmax,"WHmax/F");// hadronic maximum relative weight
    Header->Branch("WMmax",&WMmax,"WMmax/F");// muon maximum relative weight
    Header->Branch("WTmax",&WTmax,"WTmax/F");// EM maximum relative weight
    Header->Branch("hadLow",&hadLow,"hadLow/F");// hadron minimum energy for CE
    Header->Branch("muLow",&muLow,"muLow/F");   // muon minimum energy for CE
    Header->Branch("emLow",&emLow,"emLow/F");   // EM minimum energy for CE
    Header->Branch("zLow",&zLow,"zLow/F");   // minimum depth for low energy MC
    Header->Branch("iThin",&iThin,"iThin/I");   // thinning flag for EM MC
    Header->Branch("Thin",&Thin,"Thin/F");  // thinning relative level for EM MC
    Header->Branch("HGround",&HGround,"HGround/F");  // ground height (m) at impact point

#endif

    Header->Fill();
    Header->AutoSave();


    // shower data

    fShower = new TTree("Shower", "shower info");

    fShower->Branch("lgE", &currlgE, "lgE/F");
    fShower->Branch("zenith", &fitpars[1], "zenith/F");
    fShower->Branch("azimuth", &fitpars[10], "azimuth/F");
    fShower->Branch("Seed2", &fSeedcx[1], "Seed2/I");    // random seed2
    fShower->Branch("Seed3", &fSeedcx[2], "Seed3/I");    // random seed3
    fShower->Branch("Xfirst", &fitpars[2], "Xfirst/F");
    fShower->Branch("Hfirst", &fitpars[12], "Hfirst/F");
    fShower->Branch("XfirstIn", &fitpars[11], "XfirstIn/F");
    fShower->Branch("altitude", &fCurrImpactParameter, "altitude/D");
    fShower->Branch("X0", &fitpars[4], "X0/F");
    fShower->Branch("Xmax", &fitpars[9], "Xmax/F");
    fShower->Branch("Nmax", &fitpars[3], "Nmax/F");
    fShower->Branch("p1", &fitpars[5], "p1/F");
    fShower->Branch("p2", &fitpars[6], "p2/F");
    fShower->Branch("p3", &fitpars[7], "p3/F");
    fShower->Branch("chi2", &fitpars[8], "chi2/F");
    fShower->Branch("Xmx", &Xmx, "Xmx/F");
    fShower->Branch("Nmx", &Nmx, "Nmx/F");
    fShower->Branch("XmxdEdX",&XmxdEdX, "XmxdEdX/F");
    fShower->Branch("dEdXmx", &dEdXmx, "dEdXmx/F");
    fShower->Branch("cpuTime", &fTimeSec, "cpuTime/F");   // cpu time (sec)

    // profile

    fShower->Branch("nX", &nX, "nX/I");    // number of points in slant depth
    fShower->Branch("X", X, "X[nX]/F");    // slant depth  (g/cm^2)
    fShower->Branch("N", N, "N[nX]/F");    // number of charged particles(X)
    fShower->Branch("H", H, "H[nX]/F");    // height     (m)
    fShower->Branch("D", D, "D[nX]/F");    // distance     (m)
    fShower->Branch("dEdX", dEdX, "dEdX[nX]/F"); // Energy deposit(X) in GeV/(g/cm^2)
    fShower->Branch("Mu", Mu, "Mu[nX]/F");    //  number of muons(X)
    fShower->Branch("Gamma", Gamma, "Gamma[nX]/F");    //  number of gammas
    fShower->Branch("Electrons", Electrons, "Electrons[nX]/F");// # of e+e-
    fShower->Branch("Hadrons", Hadrons, "Hadrons[nX]/F");// # Hadrons
    fShower->Branch("dMu", dMu, "dMu[nX]/F");    //  muon production rate (X)
    fShower->Branch("EGround", EGround, "EGround[3]/F");    // Energy of particles at Xmax
                                                            // 0=e+gamma; 1=hadrons;2=muons

#ifdef __MC3D__

    // 3D

    fShower->Branch("MuMIA",&MuMIA,"MuMIA/F");  // MIA density at 600 m
    fShower->Branch("SD1000",&SD1000[0],"SD1000/F");  // AUGER density at 1000 m
    fShower->Branch("SDPh",&SD1000[1],"SDPh/F");  // Photon AUGER density at 1000 m
    fShower->Branch("SDEl",&SD1000[2],"SDEl/F");  // Electron AUGER density at 1000 m
    fShower->Branch("SDMu",&SD1000[3],"SDMu/F");  // Mu AUGER density at 1000 m
    fShower->Branch("SDHa",&SD1000[4],"SDHa/F");  // Hadron AUGER density at 1000 m
    fShower->Branch("MuTr",&MuTr,"MuTr/F");  // KASCADE Truncated muon number at ground
    fShower->Branch("EHaMax",&HaEm,"EHaMax/F");  // Maximum energy of the hadrons at ground
    fShower->Branch("EHaSum",&HaEs,"EHaSum/F");  // Energy sum of the hadrons at ground
#endif

#ifdef LEADING_INTERACTIONS_TREE
    if (fMaxDetail>0) {
      fLeadingInteractions = new TTree("LeadingInteractions", "first interaction");

      fLeadingInteractions->Branch("nInt", &gnInt1, "nInt/I");                // number of interatctions
      fLeadingInteractions->Branch("kinel", gkinel1, "kinel[nInt]/D");       // inelasticity
      fLeadingInteractions->Branch("pId", gpId1, "pId[nInt]/I");            // parent id
      fLeadingInteractions->Branch("pEnergy", gpEnergy1, "pEnergy[nInt]/D");// parent energy
      fLeadingInteractions->Branch("mult", gMult1, "Mult[nInt]/I");         // multiplicity
      fLeadingInteractions->Branch("matg", gMatg1, "Matg[nInt]/I");         // mass of target nucleus
      fLeadingInteractions->Branch("depth", gDepth1, "depth[nInt]/D");      // slant depth
      fLeadingInteractions->Branch("height", gHeight1, "height[nInt]/D");   // height

      fLeadingInteractions->Branch("nPart", &gNParticles1, "nPart/I");    // number of particles
      fLeadingInteractions->Branch("Energy", gEnergy1, "Energy[nPart]/D");// energy
#ifdef LEADING_INTERACTIONS_TREE_EXT
      fLeadingInteractions->Branch("px", gpx1, "px[nPart]/D");// Momentum px
      fLeadingInteractions->Branch("py", gpy1, "py[nPart]/D");// Momentum py
      fLeadingInteractions->Branch("pz", gpz1, "pz[nPart]/D");// Momentum pz
#endif // LEADING_INTERACTIONS_TREE_EXT
      fLeadingInteractions->Branch("Type", gType1, "Type[nPart]/I");      // type
      fLeadingInteractions->Branch("idInt", gIdInt1, "idInt[nPart]/I");   // interaction number
    }
#endif // LEADING_INTERACTIONS_TREE


  }
  else if  (iopt == 0) { // fill shower

    // --------- shower timing ----------
    fStopwatch.Stop();
    const double cpuTime = fStopwatch.CpuTime();
    fTimeSec = (cpuTime - fTimeTotal);
    fTimeTotal = cpuTime;
    fStopwatch.Start(false);

#ifdef LEADING_INTERACTIONS_TREE
    if (fMaxDetail>0) {
      gnInt1 = gInteractionData.size();
      if (gnInt1>maxNInt) {
        cout << " ++++++++++++ ERROR ++++++++++" << endl;
        cout << " maximum number of interaction=" << maxNInt
             << " exceeded! (" << gnInt1 << ")" << endl;
        cout << " all interaction above " << maxNInt << " are lost! " << endl;
        cout << " ++++++++++++ ERROR ++++++++++" << endl;
        gnInt1 = maxNInt;
      }
      for (int i=0; i<gnInt1; ++i) {
	cout << setw(fgOutputPrefix.size()) << ""
	     << "leading interaction " << setw(2) << i << " : "
	     << " projectile id/E=" << gInteractionData[i].pId
	     << "/" << gInteractionData[i].pEnergy
	     << ", depth=" << gInteractionData[i].Depth <<" g/cm^2"
	     << ", mult=" << gInteractionData[i].mult
	     << ", kInel=" << gInteractionData[i].kinel
	     << endl;
        gkinel1[i] = gInteractionData[i].kinel;
        gpId1[i] = gInteractionData[i].pId;
        gpEnergy1[i] = gInteractionData[i].pEnergy;
        gMult1[i] = gInteractionData[i].mult;
        gMatg1[i] = gInteractionData[i].matg;
        gDepth1[i] = gInteractionData[i].Depth;
        gHeight1[i] = gInteractionData[i].Height;
      }


      // sort vector and copy into array
      sort(gInteractionParticleData.begin(),gInteractionParticleData.end());

      gNParticles1=gInteractionParticleData.size();
      if (gNParticles1>maxNPart) {
        cout << " ++++++++++++ ERROR ++++++++++\n"
             << " maximum number of particles=" << maxNPart
             << " exceeded! (" << gNParticles1 << ")\n"
             << " all particles above " << maxNPart << " are lost!\n"
             << " ++++++++++++ ERROR ++++++++++" << endl;
        gNParticles1 = maxNPart;
      }
      for (int i=0; i<gNParticles1; ++i) {
        gEnergy1[i]=gInteractionParticleData[gNParticles1-i-1].fEnergy;
#ifdef LEADING_INTERACTIONS_TREE_EXT
        gpx1[i]=gInteractionParticleData[gNParticles1-i-1].fpx;
        gpy1[i]=gInteractionParticleData[gNParticles1-i-1].fpy;
        gpz1[i]=gInteractionParticleData[gNParticles1-i-1].fpz;
#endif // LEADING_INTERACTIONS_TREE_EXT
        gType1[i]=gInteractionParticleData[gNParticles1-i-1].fID;
        gIdInt1[i]=gInteractionParticleData[gNParticles1-i-1].fnInt;
      }
      fLeadingInteractions->Fill();
      if (!(int(fLeadingInteractions->GetEntries()+1)%fAutoSave))
        fLeadingInteractions->AutoSave();

      gInteractionParticleData.clear();
      gInteractionData.clear();
    }
#endif // LEADING_INTERACTIONS_TREE

    nX = maxX;
    fConexInterface->GetShowerData(icut, iSec, nX, X[0], N[0], fitpars[0], H[0], D[0]);
    fConexInterface->GetdEdXProfile(icut, nX, dEdX[0], EGround[0]);
    fConexInterface->GetMuonProfile(icutm, nX, Mu[0], dMu[0]);
    fConexInterface->GetGammaProfile(icutg, nX, Gamma[0]);
    fConexInterface->GetElectronProfile(icute, nX, Electrons[0]);
    fConexInterface->GetHadronProfile(icuth, nX, Hadrons[0]);

    // maximum of N particles passing planes at X
    interpolateProfile(nX, X, N, Xmx, Nmx);
    // maximum if dE/dX within a volume between two X-planes
    const double dX = X[2] - X[1];
    interpolateProfile(nX, X, dEdX, XmxdEdX, dEdXmx);
    XmxdEdX += dX/2;

    currlgE = log10(e)+9.;

#ifdef __MC3D__
    fConexInterface->Get3DOutput(icuth, HaEs, HaEm, icut, MuTr, MuMIA, SD1000[0]);
#endif

    fShower->Fill();


    cout << setw(fgOutputPrefix.size()) << "" << "shower finished: "
	 << "Xmax(dEdX)=" << XmxdEdX << " g/cm^2 "
	 << "(cpuTime="  << fTimeSec << " sec)"
	 <<  endl;

    if (!(int(fShower->GetEntries()+1)%fAutoSave))
      fShower->AutoSave();
  }
  else{   // (iopt==1) close file

    fFile->cd();
    fFile->Write();
    fFile->Close();

#ifdef DEBUG_RESAMPLING
    fResamplingDebug->cd();
    fResamplingDebug->Write();
    fResamplingDebug->Close();
#endif
  }
}



/****************************************************
 *
 * inversion of 3-by-3 matrix A
 *
 * (returns false if matrix singular)
 ***************************************************/
bool
CxRoot::invert3by3(double A[3][3])
  const
{
  const double kSmall=1e-80;

  double determinant=
     A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
    -A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
    +A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);

  double absDet = fabs(determinant);

  if(absDet < kSmall) {
    cout << " invert3by3: Error-matrix singular (absDet="<<absDet<<")" << endl;
    return false;
  }

  double B[3][3];

  B[0][0] = A[1][1]*A[2][2]-A[1][2]*A[2][1];
  B[1][0] = A[1][2]*A[2][0]-A[2][2]*A[1][0];
  B[2][0] = A[1][0]*A[2][1]-A[1][1]*A[2][0];

  B[0][1] = A[0][2]*A[2][1]-A[2][2]*A[0][1];
  B[1][1] = A[0][0]*A[2][2]-A[2][0]*A[0][2];
  B[2][1] = A[0][1]*A[2][0]-A[0][0]*A[2][1];

  B[0][2] = A[0][1]*A[1][2]-A[1][1]*A[0][2];
  B[1][2] = A[0][2]*A[1][0]-A[1][2]*A[0][0];
  B[2][2] = A[0][0]*A[1][1]-A[0][1]*A[1][0];

  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      A[i][j]=B[i][j]/determinant;

  return true;
}

/****************************************************
 *
 * solves linear system
 *
 *   / y[0] \   / A[0][0] A[0][1] A[0][2] \   / x[0] \
 *   | y[1] | = | A[1][0] A[1][1] A[1][2] | * | x[1] |
 *   \ y[2] /   \ A[2][0] A[2][1] A[2][2] /   \ x[2] /
 *
 *
 * Input:  y[3] and A[3][3]
 * Output: returns true when succeded (i.e. A is not singular)
 *         A is overwritten with its inverse
 *
 * M.Unger 12/1/05
 *
 ****************************************************/
bool
CxRoot::solve3by3(double y[3], double A[3][3], double x[3])
  const
{
  if(invert3by3(A)) {

    for (int i=0;i<3;i++)
      x[i]= A[i][0]*y[0]+A[i][1]*y[1]+A[i][2]*y[2];

    return true;

  }
  else
    return false;

}


void
CxRoot::interpolateProfile(const int nX, const float* X,
			   const float* N, float& Xmx, float& Nmx)
{
  // no maximum found yet
  Xmx = 0.;
  Nmx = 0.;

  // search maximum point
  int iMax=-1;
  float nMax=-1.;

  for (int i=0; i< nX; ++i) {
    if(N[i] > nMax) {
      iMax = i;
      nMax = N[i];
    }
  }

  if (iMax < 2 && iMax > nX-2) // too close to border
    return;

  // quadratic "fit" around maximum to get dEdXmax and Xmax
  double x[3] = {X[iMax-1],X[iMax],X[iMax+1]};
  double y[3] = {N[iMax-1],N[iMax],N[iMax+1]};
  double A[3][3];
  A[0][0] = x[0]*x[0];   A[0][1]= x[0];   A[0][2]= 1.;
  A[1][0] = x[1]*x[1];   A[1][1]= x[1];   A[1][2]= 1.;
  A[2][0] = x[2]*x[2];   A[2][1]= x[2];   A[2][2]= 1.;

  double a[3];

  solve3by3(y,A,a);

  if(a[0]<0.) {
    Xmx = -a[1]/(2.*a[0]);
    Nmx = a[0]*Xmx*Xmx+a[1]*Xmx+a[2];
  }
}

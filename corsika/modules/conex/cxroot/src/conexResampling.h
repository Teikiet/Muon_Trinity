#ifndef _include_conex_resampling_h_
#define _include_conex_resampling_h_

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include <conexConfig.h>
#include <ResamplingMode.h>

class TFile;
class TTree;


#ifdef CONEX_EXTENSIONS

// RU: to change the CrossSection by a given factor
extern double gFactorCrossSection;
extern double gFactorExtraMeson;
extern double gFactorResampling;
extern double gResamplingThreshold;

extern resample::EResampleMode gResamplingMode;



// ------------------------------------------------
//  function definitions
//

class CommonBlockCONEX;
class CommonBlockSIBYLL;
class CommonBlockSIBYLL_LAB;

extern "C" {  
  
  void resampleconex_(CommonBlockCONEX& p, double& x, int& primaryId);
  void resamplesibyll_(CommonBlockSIBYLL& p, double& x, int& primaryId,
		       float& sqs, double& plab);
  void resamplesibylllab_(CommonBlockCONEX& p, CommonBlockSIBYLL_LAB& s,
			  double& x, int& primaryId);
  //				    float& sqs, double& plab);}
  
  void modifiercx_(double& factMod, const double& energy, const int& pid);
  void sibyllf19_(double& factModPP, const double& factModPair);

}



extern "C" {
  struct xschadron {
    double xsamproj;
    double xsamtarg;
    double xsypjtl;
    double xsyhaha;
    double xspnullx;
  };
}


// ---------------------------------------------
//  debugging
//
#ifdef DEBUG_RESAMPLING
class TClonesArray;
std::string gROOTfilename = "";
// extern "C" {void cxutlob5_(double& yboost, double& x1, double& x2, double& x3, double& x4,  double& x5);}
TFile* fResamplingDebug = 0;
TTree* tResamplingDebug = 0;
int resampling_primary = 0;
double resampling_Elab = 0;
double resampling_factor = 0;
double resampling_yboost = 0;
int resampling_inucleon = 0;
int resampling_nnucleon = 0;
int resampling_cms = false;
TClonesArray* resampling_before = 0;
TClonesArray* resampling_after = 0;
int resampling_before_n = 0;
int resampling_after_n = 0;
#endif // DEBUG_RESAMPLING
#endif



#endif

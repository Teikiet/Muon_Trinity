#ifndef _conexFunctions_h_
#define _conexFunctions_h_

#include <conexConfig.h>

extern "C" {
  
  void initconex_(int&, int*, int&, int&, const char*, int);
  void endmod_(int);
  void runconex_(int*, double&, double&, double&, double&, int&);
  
  void get_data_card1_(int&, float&, float&, int&, float&, int&, float&,  
                       int&, float&, int& , float&, float&, float&,
                       float&, float&, int&, float&);
#ifdef __MC3D__
  void get_data_card2_(float&, float&, float&, float&, float&,
                       float&, float&, int&, float&, float&);
  void get_shower_lat_(int&, float&, float&, int&, float&, float&, float&);
#endif
  void get_shower_data_(int&, int&, int&, float&, float&, 
			float&, float&, float&);
  void get_shower_edep_(int&, int&, float&, float&);
  void get_shower_muon_(int&, int&, float&, float&);
  void get_shower_gamma_(int&, int&, float&);
  void get_shower_electron_(int &, int&, float&);
  void get_shower_hadron_(int &, int&, float&);
  int get_number_of_depth_bins_();
  double drangen_(double dummy);
  
}

#endif

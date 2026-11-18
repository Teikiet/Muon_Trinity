#ifndef _resample_SummaryInfo_h_
#define _resample_SummaryInfo_h_

namespace resample {

  class SummaryInfo {
  public:
    double charge;
    double totEnergy;
    double kinEnergy;
    double minKinEnergy;
    double maxKinEnergy;
    int number;
    //  TH1D* hE;
  };

}

#endif


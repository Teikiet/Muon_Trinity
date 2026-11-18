#include <cmath>
using namespace std;

//
//  plotProfiles()
//
//  ROOT macro to display simulated profiles of several showers
//
//    usage: .x plotProfiles.C("<conex-file-name(s)>", startShower, endShower) 
//
//  

#include <TLine.h>

#include <string>
using namespace std;

void plotProfiles(const string file="", int iStart=0, int iStop=10000, int color=kRed, const string leg="",
		  const string fileB="", int iStartB=0, int iStopB=10000, int colorB=kBlue, const string legB="",
		  const string fileC="", int iStartC=0, int iStopC=10000, int colorC=kGreen+1, const string legC="") {
  
  if ( file=="" ) {
    cout << "\n usage: .x plotProfiles.C(\"<conex-file-name(s)>\", startShower, endShower) \n"
	 << endl;
    return;
  }
  
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleSize(0.06, "X");
  gStyle->SetTitleSize(0.06, "Y");
  gStyle->SetLabelSize(0.04, "X");
  gStyle->SetLabelSize(0.04, "Y");
  gStyle->SetTitleOffset(1.1,"X");
  gStyle->SetTitleOffset(1.2,"Y");


  

  TCanvas* cEl = gROOT->FindObject("cEl");
  if ( cEl != NULL ) delete cEl;
  cEl = new TCanvas("cEl", "profile for electrons", 800, 10, 700, 500);
  TLegend* legend = new TLegend(.6, .6, 0.99-gPad->GetRightMargin(), 0.99-gPad->GetTopMargin());
  legend->SetFillColor(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);


  TCanvas* cMu = gROOT->FindObject("cMu");
  if ( cMu != NULL ) delete cMu;
  cMu = new TCanvas("cMu", "profile for muons", 800, 10, 700, 500);

  bool firstGraph = true;
  if (file!="") plotProfilesFile(file, iStart, iStop, color, cEl, cMu, firstGraph, leg, legend);
  if (fileB!="") plotProfilesFile(fileB, iStartB, iStopB, colorB, cEl, cMu, firstGraph, legB, legend);
  if (fileC!="") plotProfilesFile(fileC, iStartC, iStopC, colorC, cEl, cMu, firstGraph, legC, legend);
  
  cEl->cd();
  legend->Draw();
}


void plotProfilesFile(const string& file, int iStart, int iStop, int color, 
		      TCanvas* cEl, TCanvas* cMu, bool& firstGraph,
		      const string& legStr, TLegend* legend) {

  int lineSize = 2;
  
  // set and read tree
  
  TChain Shower("Shower");
  TChain Header("Header");
  Shower.Add(file.c_str());
  Header.Add(file.c_str());
  
  const int maxX = 5000;
  float Xmax;
  float X[maxX], H[maxX], D[maxX], N[maxX],eN[maxX],dEdX[maxX], Mu[maxX], dMu[maxX];
  float fitpars[11],EGround[3],muThr,hadThr,emThr, outputVersion;
  int nX,iPart,iseed1,iseed2,iseed3,HEModel;
  
  //Shower.SetBranchAddress("lgE",&fitpars[0]);
  //Shower.SetBranchAddress("zenith",&fitpars[1]);
  //Shower.SetBranchAddress("azimuth",&fitpars[10]);
  //Shower.SetBranchAddress("Xfirst",&fitpars[2]);
  //Shower.SetBranchAddress("X0",&fitpars[4]);
  Shower.SetBranchAddress("Xmax",&Xmax);
  //Shower.SetBranchAddress("Nmax",&fitpars[3]);
  //Shower.SetBranchAddress("p1",&fitpars[5]);
  //Shower.SetBranchAddress("p2",&fitpars[6]);
  //Shower.SetBranchAddress("p3",&fitpars[7]);
  //Shower.SetBranchAddress("chi2",&fitpars[8]);
  Shower.SetBranchAddress("nX",&nX);    // number of points in slant depth     
  Shower.SetBranchAddress("N",N);       // particles(X)      
  //Shower.SetBranchAddress("dEdX",dEdX); // Energy deposit(X)      
  Shower.SetBranchAddress("Mu",Mu);       // muons(X)      
  Shower.SetBranchAddress("dMu",dMu);       // muon production(X)      
  //Shower.SetBranchAddress("EGround",EGround);
  //Header.SetBranchAddress("muThr",&muThr);  
  //Header.SetBranchAddress("emThr",&emThr);  
  //Header.SetBranchAddress("hadThr",&hadThr);  
  //Header.SetBranchAddress("HEModel",&HEModel);  // high energy model    
  //Header.SetBranchAddress("Particle",&iPart);  
  //Header.SetBranchAddress("Seed1",&iseed1);
  Header.SetBranchAddress("OutputVersion",&outputVersion);
  //Shower.SetBranchAddress("Seed2",&iseed2);
  //Shower.SetBranchAddress("Seed3",&iseed3);
  
  // output v2.0 moved H/D/X to Shower, removed nX from Header
  Header.GetEntry(0);
  if (outputVersion < 2.0) {
    Header.SetBranchAddress("H",H);       // height above see level  
    Header.SetBranchAddress("D",D);       // distance to obs. point
    Header.SetBranchAddress("X",X);       // slant depth     
  }
  else {
    Shower.SetBranchAddress("H",H);       // height above see level  
    Shower.SetBranchAddress("D",D);       // distance to obs. point
    Shower.SetBranchAddress("X",X);       // slant depth     
  }
  
  Header.GetEntry(0);

  bool firstGraphOfFile = true;
  //  for (int i=iStart; i<=(int)min((double)Shower.GetEntries(), (double)iStop); ++i) {
  for (int i=iStart; i<=iStop; ++i) {
    
    Shower.GetEntry(i);
    
    cEl->cd();
    TGraph* graphEl = new TGraph(nX, X, N);
    graphEl->SetMarkerStyle(20); 
    graphEl->SetMarkerSize(.5); 
    graphEl->SetMarkerColor(color); 
    graphEl->SetLineColor(color); 
    graphEl->SetLineWidth(lineSize); 
    graphEl->SetTitle("");
    graphEl->GetXaxis()->SetTitle("Depth   [g/cm^{2}]");
    graphEl->GetYaxis()->SetTitle("Number of electrons");
    //graphEl->GetYaxis()->SetTitleOffset(1.);
    //graphEl->Draw( (firstGraph ? "ALP" : "LP") ); 
    graphEl->Draw( (firstGraph ? "AL" : "L") ); 
    
    cMu->cd();
    TGraph* graphMu = new TGraph(nX, X, Mu);
    graphMu->SetMarkerStyle(20); 
    graphMu->SetMarkerSize(.5); 
    graphMu->SetMarkerColor(color); 
    graphMu->SetLineColor(color); 
    graphMu->SetLineWidth(lineSize); 
    graphMu->SetTitle("");
    graphMu->GetXaxis()->SetTitle("X [g/cm^{2}]");
    graphMu->GetYaxis()->SetTitle("number of muons");
    //graphMu->GetYaxis()->SetTitleOffset(1.);
    //graphMu->Draw( (firstGraph ? "ALP" : "LP") );
    graphMu->Draw( (firstGraph ? "AL" : "L") );

    if (firstGraph) { // to get the pad dimensions
      cEl->Modified();
      cEl->Update();
      //cMu->Modified();
      //cMu->Update();
    }

    if (firstGraphOfFile && legStr!="") { // legend
      legend->AddEntry(graphEl, legStr.c_str(), "l");
    }
    
    // draw Xmax-ticks
    cEl->cd();
    double yMin = gPad->GetUymin();
    double yMax = gPad->GetUymax();
    TLine* line = new TLine();
    line->SetLineColor(color);
    line->SetLineWidth(2);
    line->SetLineStyle(1);
    line->DrawLine(Xmax, yMax-(yMax-yMin)*0.025, Xmax, yMax);
    
    firstGraph = false;
    firstGraphOfFile = false;
  }  
}


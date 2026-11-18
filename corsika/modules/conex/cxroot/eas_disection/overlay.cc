#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TLegend.h>
#include <TTree.h>
#include <TH1.h>

#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <cmath>
using namespace std;

const int seed = 30;

// -----------------------------------------------------------
struct graphs {
  graphs() : lgE(0), dEdX(0), El(0), Mu(0), dEdXall(0), Elall(0), Muall(0) {}
  float lgE;
  TGraph* dEdX;
  TGraph* El;
  TGraph* Mu;
  TMultiGraph* dEdXall;
  TMultiGraph* Elall;
  TMultiGraph* Muall;
  void SetColor(const int c) {dEdX->SetLineColor(c); El->SetLineColor(c); Mu->SetLineColor(c);}
  void SetWidth(const int w) {dEdX->SetLineWidth(w); El->SetLineWidth(w); Mu->SetLineWidth(w);}
  void SetStyle(const int s) {dEdX->SetLineStyle(s); El->SetLineStyle(s); Mu->SetLineStyle(s);}
};


// -----------------------------------------------------------
graphs
getData(const string& prefix, 
	const int iFstart, const int iFstop, 
	const string& name)
{
  const int nMaxX = 10000;
  float Xsum[nMaxX], dEdXsum[nMaxX], Musum[nMaxX], Electronssum[nMaxX];
  for(int i=0; i<nMaxX; ++i) {
    Xsum[i] = 0;
    dEdXsum[i] = 0;
    Musum[i] = 0;
    Electronssum[i] = 0;
  }

  graphs gRet;
  gRet.dEdXall = new TMultiGraph();
  gRet.Elall = new TMultiGraph();
  gRet.Muall = new TMultiGraph();
  
  int nX_max = 0;
  for(int iF=iFstart; iF<iFstop; ++iF) {
    
    ostringstream fullname;
    fullname << prefix;
    if (iF>=0)
      fullname << iF;
    fullname << name;
    
    TFile* file = TFile::Open(fullname.str().c_str());
    cout << "Reading file: " << fullname.str() << endl;
    if (!file || file->IsZombie()) {
      cout << "Could not open file " << fullname.str() << endl;
      continue;
    }

    TTree* tree = (TTree*) file->Get("Shower");
    if (!tree) {
      cout << "Input file broken: Shower tree not found in " << fullname.str() << endl;
      continue;
    }
    int nX = 0;
    float X[nMaxX], dEdX[nMaxX], Mu[nMaxX], Electrons[nMaxX];    
    tree->SetBranchAddress("lgE", &gRet.lgE);
    tree->SetBranchAddress("nX", &nX);
    tree->SetBranchAddress("X", X);
    tree->SetBranchAddress("dEdX", dEdX);
    tree->SetBranchAddress("Mu", Mu);
    tree->SetBranchAddress("Electrons", Electrons);
    tree->GetEntry(0);

    
    nX_max = max(nX_max, nX);
    
    for(int i=0; i<nX; ++i) {
      Xsum[i] = X[i];
      dEdXsum[i] += dEdX[i];
      Musum[i] += Mu[i];
      Electronssum[i] += Electrons[i];
    }
    
    TGraph* gdEdX = new TGraph(nX-1, X, dEdX);
    gdEdX->SetLineStyle(2);
    gRet.dEdXall->Add(gdEdX);    
    
    TGraph* gEl = new TGraph(nX, X, Electrons);
    gEl->SetLineStyle(2);
    gRet.Elall->Add(gEl);    

    TGraph* gMu = new TGraph(nX, X, Mu);
    gMu->SetLineStyle(2);
    gRet.Muall->Add(gMu);    
    
    file->Close();
  }
  
  gRet.dEdX = new TGraph(nX_max-1, Xsum, dEdXsum);
  gRet.dEdX->GetHistogram()->SetXTitle("Depth [g/cm^{2}]");  
  gRet.dEdX->GetHistogram()->SetYTitle("Energy Deposit  [GeV cm^{2}/g]");
  gRet.dEdX->GetXaxis()->CenterTitle();
  gRet.dEdX->GetYaxis()->CenterTitle();
  gRet.El = new TGraph(nX_max, Xsum, Electronssum);
  gRet.El->GetHistogram()->SetXTitle("Depth [g/cm^{2}]");
  gRet.El->GetHistogram()->SetYTitle("Electrons");
  gRet.El->GetXaxis()->CenterTitle();
  gRet.El->GetYaxis()->CenterTitle();
  gRet.Mu = new TGraph(nX_max, Xsum, Musum);
  gRet.Mu->GetHistogram()->SetXTitle("Depth [g/cm^{2}]");
  gRet.Mu->GetHistogram()->SetYTitle("Muons");
  gRet.Mu->GetXaxis()->CenterTitle();
  gRet.Mu->GetYaxis()->CenterTitle();
  return gRet;
}


// -----------------------------------------------------------
void
overlay(const int seed, const int intMin, const int intMax)
{
  const int logY = 1;
  TCanvas* cShowerdEdX = new TCanvas("cShowerdEdX");
  cShowerdEdX->SetLogy(logY);
  TCanvas* cShowerMu = new TCanvas("cShowerMu");
  cShowerMu->SetLogy(logY);
  TCanvas* cShowerEl = new TCanvas("cShowerEl");
  cShowerEl->SetLogy(logY);

  TLegend* leg = new TLegend(0.55, 0.75, 0.7, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  
  ostringstream fullName;
  fullName << "sibyll_" << setw(9) << setfill('0') << seed << "_100.root";

  graphs gFull = getData("full_",-1, 0, fullName.str());
  gFull.SetColor(kBlue+1);
  gFull.SetWidth(3);
  ostringstream legFull;
  legFull << "Proton, 10^{" << gFull.lgE << "}eV";
  leg->AddEntry(gFull.dEdX, legFull.str().c_str(), "l");

  cShowerdEdX->cd();
  gFull.dEdX->Draw("al");

  cShowerEl->cd();
  gFull.El->Draw("al");

  cShowerMu->cd();
  gFull.Mu->Draw("al");

  ostringstream intName;
  intName << "_sibyll_" << setw(9) << setfill('0') << seed << "_100.root";
  
  graphs gInt = getData("int_", intMin, intMax, intName.str());
  gInt.SetColor(kRed);
  gInt.SetWidth(3);
  ostringstream legTxt;
  legTxt << intMax - intMin + 1 << " Highest Energy Interactions";
  leg->AddEntry(gInt.dEdX, legTxt.str().c_str(), "l");
  leg->AddEntry(gInt.dEdXall, "Individual Sub-Showers", "l");

  cShowerdEdX->cd();
  gInt.dEdX->Draw("l");
  gInt.dEdXall->Draw("l");
  leg->Draw();

  cShowerEl->cd();
  gInt.El->Draw("l");
  gInt.Elall->Draw("l");
  leg->Draw();

  cShowerMu->cd();
  gInt.Mu->Draw("l");
  gInt.Muall->Draw("l");
  leg->Draw();
}


// -----------------------------------------------------------
int
main(int argc, char** argv)
{
  gROOT->SetStyle("Plain");
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.055, "XY");
  gStyle->SetTitleFont(42, "XY");
  //  gStyle->SetHistTitleSize(0.06, "X");

  if (argc<4) {
    cout << " specify: [seed] [interaction min+max]" << endl;
    return 2;
  }
 
  const int seed = atoi(argv[1]);
  const int intMin = atoi(argv[2]);
  const int intMax = atoi(argv[3]);

  TApplication app("hold", 0, 0);
  overlay(seed, intMin, intMax);
  app.Run();
  return 0;
}

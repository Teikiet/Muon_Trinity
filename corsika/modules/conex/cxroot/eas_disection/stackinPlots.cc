// ROOT
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TFile.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TF1.h>
#include <THStack.h>

// c++
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <cmath>

// COAST
#include <crsRead/MCorsikaReader.h>
#include <crs/TSubBlock.h>
#include <crs/CorsikaConsts.h>
#include <crs/MRunHeader.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MParticleBlock.h>
#include <crs/MLongitudinalBlock.h>
#include <crs/MParticle.h>

using namespace std;


map<string, TH1*> 
readDAT(const string& filename, const int color, const int fill)
{
  const int nBins = 40;
  const double maxR = 5000;

  TH1D* hMu = new TH1D((filename+"_muons").c_str(), "Muon Density", 30, 0, maxR);
  hMu->SetXTitle("Distance / m");
  hMu->SetYTitle("Density / m^{2}");
  hMu->SetLineColor(color);
  hMu->SetFillColor(color);
  //hMu->SetFillStyle(0);
  hMu->SetLineWidth(3);
  hMu->SetFillStyle(fill);

  TH1D* hMuSpect = new TH1D((filename+"_muonsSpect").c_str(), "Muon Spectrum", 20, 0, 500); // GeV
  hMuSpect->SetXTitle("Energy / GeV");
  hMuSpect->SetYTitle("dN/dE / GeV^{-1}");
  hMuSpect->SetLineColor(color);
  hMuSpect->SetFillColor(color);
  //hMuSpect->SetFillStyle(0);
  hMuSpect->SetLineWidth(3);
  hMuSpect->SetFillStyle(fill);

  TH1D* hEl = new TH1D((filename+"_em").c_str(), "Electron Density", 30, 0, maxR);
  hEl->SetXTitle("Distance / m");
  hEl->SetYTitle("Density / m^{2}");
  hEl->SetLineColor(color);
  hEl->SetFillColor(color);
  hEl->SetLineWidth(3);
  hEl->SetFillStyle(fill);
  TH1D* hHad = new TH1D((filename+"_had").c_str(), "Hadron Density", 30, 0, maxR/100);
  hHad->SetXTitle("Distance / m");
  hHad->SetYTitle("Density / m^{2}");
  hHad->SetLineColor(color);
  hHad->SetFillColor(color);
  hHad->SetLineWidth(3);
  hHad->SetFillStyle(fill);
  TH1D* hGam = new TH1D((filename+"gam").c_str(), "Photon Density", 30, 0, maxR);
  hGam->SetXTitle("Distance / m");
  hGam->SetYTitle("Density / m^{2}");
  hGam->SetLineColor(color);
  hGam->SetFillColor(color);
  hGam->SetLineWidth(3);
  hGam->SetFillStyle(fill);

  map<string, TH1*> data;  
  data["em"] = hEl;
  data["mu"] = hMu;
  data["muspect"] = hMuSpect;
  data["had"] = hHad;
  data["gam"] = hGam;

  crsRead::MCorsikaReader cr(filename, 3);
  
  crs::MRunHeader Run;
  if (!cr.GetRun (Run)) {
    cout << "No valid RUNH found in DAT file" << endl;
    return data;
  }
  
  crs::MEventHeader Shower;
  if (!cr.GetShower(Shower)) {
    cout << "No valid EVTH found in DAT file" << endl;
    return data;
  }
  
  const double height = Shower.GetObservationHeight(0);
  cout << " init obs-level at h=" << height << endl;


  crs::TSubBlock Data;
  while (cr.GetData (Data)) {
        
    switch (Data.GetBlockType ()) {
          
    case crs::TSubBlock::ePARTDATA:
      {
	const crs::MParticleBlock& ParticleData = Data;
	crs::MParticleBlock::ParticleListConstIterator iEntry;
	for (iEntry = ParticleData.FirstParticle();
	     iEntry != ParticleData.LastParticle();
	     ++iEntry) {
	  
	  if (iEntry->IsParticle()) {
                  
	    crs::MParticle iPart(*iEntry);
                  
	    const int id    = iPart.GetParticleID();
	    const int level = iPart.GetObservationLevel();
	    const double w  = iPart.GetWeight();
	    const double e  = sqrt(pow(iPart.GetKinEnergy(),2) + pow(iPart.GetMass(),2));
	    const double x  = iPart.GetX() / crs::meter;
	    const double y  = iPart.GetY() / crs::meter;
	    const double r = sqrt(x*x + y*y);
	    
	    if (level!=1) 
	      continue;

	    switch (id) {
	      
	    case 2:
	    case 3:
	      hEl->Fill(r, w);
	      break;

	    case 5:
	    case 6:
	      hMu->Fill(r, w); // lateral
	      hMuSpect->Fill(e, w); // spectrum
	      break;
	      
	      //case 7: pi0
	    case 8:
	    case 9:
	    case 10:
	    case 11:
	    case 12:
	      // case 13: neutron 
	    case 14:
	    case 15:
	      hHad->Fill(r, w);
	      break;

	    case 1:
	    case 7:
	      hGam->Fill(r, w);
	      break;
	      
	    }
	    
	  }
                
	} // end particle loop
              
	break;
      }
      
    default:
      break;
    } // end data block
        
  } // loop data
      
  crs::MEventEnd ShowerSummary;
  cr.GetShowerSummary(ShowerSummary);
  const double Xmax = ShowerSummary.GetXmax();
  
  // bin-weights
  TF1* rScale = new TF1("rScale", "2.*3.141*x", 0, 100000);
  hMu->Divide(rScale, hMu->GetBinWidth(1));
  hHad->Divide(rScale, hHad->GetBinWidth(1));
  hEl->Divide(rScale, hEl->GetBinWidth(1));
  hGam->Divide(rScale, hGam->GetBinWidth(1));

  hMuSpect->Scale(1./(hMuSpect->GetBinWidth(1)));
  
  cout << "---------------------------------\n"
       << " Shower info:\n"
       << "  Xmax = " << Xmax << "\n";

  return data;
}



map<string, TH1*>
readLong(const string& filename, const int color, const int fill)
{  
  // parse long-file
  ifstream dataLong((filename+".long").c_str());
  
  char line[5000];
  dataLong.getline(line, 5000);

  string dummy;
  int nLong;
  istringstream liness(line);
  liness >> dummy >> dummy >> dummy >> nLong;
  
  dataLong.getline(line, 5000); 
  
  // DEPTH     GAMMAS   POSITRONS   ELECTRONS         MU+         MU-     HADRONS     CHARGED      NUCLEI   CHERENKOV
  double DEPTH, GAMMAS, POSITRONS, ELECTRONS, MUP, MUM, HADRONS, CHARGED, NUCLEI, CHERENKOV;
  vector<double> vecDepth, vecEm, vecMu, vecHad, vecCharged;
  vecDepth.push_back(0);
  for (int iLong=0; iLong<nLong; ++iLong) {
    dataLong >> DEPTH >> GAMMAS >> POSITRONS >> ELECTRONS >> MUP >> MUM >> HADRONS >> CHARGED >> NUCLEI >> CHERENKOV;    
    dataLong.ignore(9999, '\n');
    if (!dataLong.good()) {
      cout << "LONG FILE READ ERROR" << endl;
      break;
    }
    vecDepth.push_back(DEPTH);
    vecEm.push_back(/*GAMMAS+*/POSITRONS+ELECTRONS);
    vecMu.push_back(MUP+MUM);
    vecHad.push_back(HADRONS);
    vecCharged.push_back(CHARGED);
  } // end loop read file
  
  dataLong.close();
  
  if (vecDepth.size()<3) {
    vecDepth.push_back(0);
    vecEm.push_back(0);
    vecMu.push_back(0);
    vecHad.push_back(0);    
  }

  TH1* hEm = new TH1D((filename+"_longem").c_str(), "Electron Profile", vecDepth.size()-1, &vecDepth.front()); 
  for (int i=0; i<vecEm.size(); ++i) 
    hEm->SetBinContent(i+1, vecEm[i]);
  hEm->SetLineColor(color);
  hEm->SetFillColor(color);
  hEm->SetLineWidth(3);
  hEm->SetFillStyle(fill);

  TH1D* hMu = new TH1D((filename+"_longMu").c_str(), "Muon Profile", vecDepth.size()-1, &vecDepth.front());
  for (int i=0; i<vecMu.size(); ++i) 
    hMu->SetBinContent(i+1, vecMu[i]);
  hMu->SetLineColor(color);
  hMu->SetFillColor(color);
  hMu->SetLineWidth(3);
  hMu->SetFillStyle(fill);

  TH1D* hHad = new TH1D((filename+"_longHad").c_str(), "Hadron Profile", vecDepth.size()-1, &vecDepth.front());
  for (int i=0; i<vecHad.size(); ++i) 
    hHad->SetBinContent(i+1, vecHad[i]);
  hHad->SetLineColor(color);
  hHad->SetFillColor(color);
  hHad->SetLineWidth(3);
  hHad->SetFillStyle(fill);

  TH1D* hCharged = new TH1D((filename+"_longCh").c_str(), "Charged Profile", vecDepth.size()-1, &vecDepth.front());
  for (int i=0; i<vecCharged.size(); ++i) 
    hCharged->SetBinContent(i+1, vecCharged[i]);
  hCharged->SetLineColor(color);
  hCharged->SetFillColor(color);
  hCharged->SetLineWidth(3);
  hCharged->SetFillStyle(fill);
  
  map<string, TH1*> data;
  data["charged"] = hCharged;
  data["em"] = hEm;
  data["mu"] = hMu;
  data["had"] = hHad;
  return data;
}



void
plots()
{
  ifstream cfg("stackinPlots.cfg");

  if (!cfg.good()) {
    cout << "No config file found!" << endl;
    exit(1);
  }

  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05,"XY");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetPadBottomMargin(0.12);

  // long
  TCanvas* cLongEm = new TCanvas("longEm");
  TCanvas* cLongHad = new TCanvas("longHad");
  TCanvas* cLongMu = new TCanvas("longMu");
  //rad
  TCanvas* cRadEm = new TCanvas("radEm");
  TCanvas* cRadHad = new TCanvas("radHad");
  TCanvas* cRadMu = new TCanvas("radMu");
  TCanvas* cRadGam = new TCanvas("radGam");
  //spect
  TCanvas* cSpectMu = new TCanvas("spectMu");

  THStack* stackLongEm = new THStack("stackLongEm", "Electron Profile;Depth [g/cm^{2}];Number");
  THStack* stackLongHad = new THStack("stackLongHad", "Hadron Profile;Depth [g/cm^{2}];Number");
  THStack* stackLongMu = new THStack("stackLongMu", "Muon Profile;Depth [g/cm^{2}];Number");

  THStack* stackRadEm = new THStack("stackRadEm", "Electron Density;Distance [m];Density [1/m^{2}]");
  THStack* stackRadMu = new THStack("stackRadMu", "Muon Density;Distance [m];Density [1/m^{2}]");
  THStack* stackRadHad = new THStack("stackRadHad", "Hadron Density;Distance [m];Density [1/m^{2}]");
  THStack* stackRadGam = new THStack("stackRadGam", "Photon Density;Distance [m];Density [1/m^{2}]");

  THStack* stackSpectMu = new THStack("stackSpectMu", "Muon Spectrum;E [GeV];dN/dE [1/GeV]");
  
  cLongEm->SetLogy(1);
  cLongHad->SetLogy(1);
  cLongMu->SetLogy(1);
  cRadEm->SetLogy(1);
  cRadHad->SetLogy(1);
  cRadMu->SetLogy(1);
  cRadGam->SetLogy(1);

  cSpectMu->SetLogy(1);
  
  //TLegend* leg = new TLegend(0.13, 0.7, 0.3, 0.9);
  TLegend* leg = new TLegend(0.65, 0.65, 0.9, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  
  bool first = true;
  while (cfg.good()) {
    
    string filename;
    cfg >> filename;
    
    int color = 1;
    int fill = 0;
    cfg >> color;
    cfg >> fill;

    char comment[5000];
    cfg.getline(comment, 5000);
    
    if (!cfg.good())
      continue;
    
    cout << "ana: " << filename << " " << comment << endl;
    
    map<string, TH1*> plotsLong = 
      readLong(filename, color, fill); //, comment, "same", leg);
    
    if (plotsLong["mu"]->GetNbinsX()<10) // sanity check
      continue;

    map<string, TH1*> plotsDat = 
      readDAT(filename, color, fill); //, comment, "same", leg);
        
    leg->AddEntry(plotsLong["mu"], comment, "f");
    
    
    
    // long
    cLongMu->cd();
    stackLongMu->Add(plotsLong["mu"]);//->Draw( first?"alp":"lp" );

    cLongEm->cd();
    stackLongEm->Add(plotsLong["em"]);//->Draw( first?"alp":"lp" );

    cLongHad->cd();
    stackLongHad->Add(plotsLong["had"]);//->Draw( first?"alp":"lp" );

		     
    // radial
    cRadMu->cd();
    stackRadMu->Add(plotsDat["mu"]);//->Draw( first?"":"same" );

    cRadEm->cd();
    stackRadEm->Add(plotsDat["em"]); //->Draw( first?"":"same" );

    cRadHad->cd();
    stackRadHad->Add(plotsDat["had"]);//->Draw( first?"":"same" );

    cRadGam->cd();
    stackRadGam->Add(plotsDat["gam"]);//->Draw( first?"":"same" );

    // spectra
    cSpectMu->cd();
    stackSpectMu->Add(plotsDat["muspect"]);
    
    first = false;    
  }  
  
  cLongMu->cd();
  stackLongMu->Draw();
  //leg->Draw();
  cLongMu->SaveAs("LongMuons.eps");
  cLongMu->SaveAs("LongMuons.pdf");
  
  cLongHad->cd();
  stackLongHad->Draw();
  //leg->Draw();
  cLongHad->SaveAs("LongHadrons.eps");
  cLongHad->SaveAs("LongHadrons.pdf");
  
  cLongEm->cd();
  stackLongEm->Draw();
  //leg->Draw();
  cLongEm->SaveAs("LongElectrons.eps");
  cLongEm->SaveAs("LongElectrons.pdf");
  
  
  cRadMu->cd();
  stackRadMu->SetMaximum(stackRadMu->GetMaximum()*1.2);
  stackRadMu->Draw("c");
  leg->Draw();
  cRadMu->SaveAs("RadialMuons.eps");
  cRadMu->SaveAs("RadialMuons.pdf");
  
  cRadHad->cd();
  stackRadHad->Draw("c");
  leg->Draw();
  cRadHad->SaveAs("RadialHadrons.eps");
  cRadHad->SaveAs("RadialHadrons.pdf");
  
  cRadEm->cd();
  stackRadEm->Draw("c");
  leg->Draw();
  cRadEm->SaveAs("RadialElectrons.eps");
  cRadEm->SaveAs("RadialElectrons.pdf");
  
  cRadGam->cd();
  stackRadGam->Draw("c");
  leg->Draw();
  cRadGam->SaveAs("RadialGamma.eps");
  cRadGam->SaveAs("RadialGamma.pdf");

  cSpectMu->cd();
  stackSpectMu->Draw("c");
  leg->Draw();
  cSpectMu->SaveAs("SpectrumMuons.eps");
  cSpectMu->SaveAs("SpectrumMuons.pdf");
  

  // The S1000 plots

  const Int_t nx = 6;
  const char *label[nx] = {"Central","Endcap","Forward","CASTOR+T2","FSC","ZDC"};

  TCanvas* cS1000Mu = new TCanvas("cS1000Mu");  
  cS1000Mu->SetLeftMargin(0.175);

  TLegend* legS1000 = new TLegend(0.2, 0.75, 0.45, 0.9);
  legS1000->SetFillStyle(0);
  legS1000->SetBorderSize(0);
  legS1000->SetTextFont(42);
  legS1000->SetTextSize(0.045);
  
  

  TList* stackRadMuList = stackRadMu->GetHists();
  TH1D* hS1000Mu = new TH1D("hS1000Mu", "Density at 1000m; ; Density [1/m^{2}]", stackRadMuList->GetEntries(), 0, stackRadMuList->GetEntries());
  hS1000Mu->GetXaxis()->SetLabelSize(0.07);
  hS1000Mu->GetYaxis()->SetTitleOffset(1.7);
  hS1000Mu->SetLineWidth(4);
  hS1000Mu->SetLineColor(kRed);
  for (int i=0; i<stackRadMuList->GetEntries(); ++i) {
    TH1* hDet = (TH1*) stackRadMuList->At(i);
    const double n = hDet->Interpolate(1000.);
    hS1000Mu->SetBinContent(i+1, n);
    
    ostringstream hS1000MuName;
    hS1000MuName << "hS1000Mu_" << i;
    TH1D* hS1000MuDet = new TH1D(hS1000MuName.str().c_str(), "Density at 1000m; ; Density [1/m^{2}]", stackRadMuList->GetEntries(), 0, stackRadMuList->GetEntries());
    hS1000MuDet->GetXaxis()->SetLabelSize(0.07);
    hS1000MuDet->GetYaxis()->SetTitleOffset(1.7);
    hS1000MuDet->SetLineWidth(4);
    hS1000MuDet->SetLineColor(hDet->GetLineColor());
    hS1000MuDet->SetFillColor(hDet->GetFillColor());
    hS1000MuDet->SetFillStyle(hDet->GetFillStyle());
    //hS1000MuDet->
  }
  
  TList* stackRadElList = stackRadEm->GetHists();
  TH1D* hS1000El = new TH1D("hS1000El", "Density at 1000m; ; Density [1/m^{-2}]", stackRadElList->GetEntries(), 0, stackRadElList->GetEntries());
  hS1000El->GetXaxis()->SetLabelSize(0.07);
  hS1000El->GetYaxis()->SetTitleOffset(1.7);
  hS1000El->SetLineWidth(4);
  hS1000El->SetLineColor(kBlue);
  for (int i=0; i<stackRadElList->GetEntries(); ++i) {
    const double n = ((TH1*)stackRadElList->At(i))->Interpolate(1000.);
    hS1000El->SetBinContent(i+1, n);
  }
  
  for (int i=0; i<nx; ++i) {
    hS1000El->GetXaxis()->SetBinLabel(i+1, label[i]);
    hS1000Mu->GetXaxis()->SetBinLabel(i+1, label[i]);
  }
  
  TList* stackRadGamList = stackRadGam->GetHists();
  TH1D* hS1000Gam = new TH1D("hS1000Gam", "hS1000Gam", stackRadGamList->GetEntries(), 0, stackRadGamList->GetEntries());
  hS1000Gam->SetLineWidth(4);
  hS1000Gam->SetLineColor(kCyan+1);
  for (int i=0; i<stackRadGamList->GetEntries(); ++i) {
    const double n = ((TH1*)stackRadGamList->At(i))->Interpolate(1000.);
    hS1000Gam->SetBinContent(i+1, n);
  }
  hS1000El->Add(hS1000Gam);

  hS1000El->Draw();
  hS1000Mu->Draw("same");
  
  legS1000->AddEntry(hS1000Mu, "Muons", "l");
  legS1000->AddEntry(hS1000El, "Electrons+Photons", "l");
  legS1000->Draw();  

  cS1000Mu->SaveAs("S1000.eps");
  cS1000Mu->SaveAs("S1000.pdf");
}



int
main()
{
  TApplication app("holder", 0, 0);
  plots();
  app.Run();
  return 0;
}

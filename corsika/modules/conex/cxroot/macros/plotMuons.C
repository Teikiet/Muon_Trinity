//
//  plotProfiles()
//
//  ROOT macro to display simulated profiles of several showers
//
//    usage: .x plotProfiles.C("<conex-file-name(s)>", startShower, endShower) 
//
//  
void plotMuons(const string& file="", int iStart=0, int iStop=10000) {

  if ( file=="" ) {
    cout << "\n usage: .x plotMuons.C(\"<conex-file-name(s)>\", startShower, endShower) \n"
	 << endl;
    return;
  }
  
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  // set and read tree
  
  TChain* Shower = new TChain("Shower");
  TChain* Header = new TChain("Header");
  Shower->Add(file.c_str());
  Header->Add(file.c_str());
  
  const int maxX = 5000;
  float X[maxX], H[maxX], D[maxX], N[maxX],eN[maxX],dEdX[maxX], Mu[maxX], dMu[maxX];
  float fitpars[11],EGround[3],muThr,hadThr,emThr, outputVersion;
  int nX,iPart,iseed1,iseed2,iseed3,HEModel;
  
  Shower->SetBranchAddress("lgE",&fitpars[0]);
  Shower->SetBranchAddress("zenith",&fitpars[1]);
  Shower->SetBranchAddress("azimuth",&fitpars[10]);
  Shower->SetBranchAddress("Xfirst",&fitpars[2]);
  Shower->SetBranchAddress("X0",&fitpars[4]);
  Shower->SetBranchAddress("Xmax",&fitpars[9]);
  Shower->SetBranchAddress("Nmax",&fitpars[3]);
  Shower->SetBranchAddress("p1",&fitpars[5]);
  Shower->SetBranchAddress("p2",&fitpars[6]);
  Shower->SetBranchAddress("p3",&fitpars[7]);
  Shower->SetBranchAddress("chi2",&fitpars[8]);
  Shower->SetBranchAddress("nX",&nX);    // number of points in slant depth     
  Shower->SetBranchAddress("N",N);       // particles(X)      
  Shower->SetBranchAddress("dEdX",dEdX); // Energy deposit(X)      
  Shower->SetBranchAddress("Mu",Mu);       // muons(X)      
  Shower->SetBranchAddress("dMu",dMu);       // muon production(X)      
  Shower->SetBranchAddress("EGround",EGround);
  Header->SetBranchAddress("muThr",&muThr);  
  Header->SetBranchAddress("emThr",&emThr);  
  Header->SetBranchAddress("hadThr",&hadThr);  
  Header->SetBranchAddress("HEModel",&HEModel);  // high energy model    
  Header->SetBranchAddress("Particle",&iPart);  
  Header->SetBranchAddress("Seed1",&iseed1);
  Header->SetBranchAddress("OutputVersion",&outputVersion);
  Shower->SetBranchAddress("Seed2",&iseed2);
  Shower->SetBranchAddress("Seed3",&iseed3);
  
  tmp=gROOT->(TProfile *) gROOT->FindObject("allmu");
  if ( tmp ) delete tmp;
  TProfile * allmu  =  new TProfile("allmu","<m>",maxX,0,2010," ");

  // output v2.0 moved H/D/X to Shower, removed nX from Header
  Header->GetEntry(0);
  if (outputVersion < 2.0) {
    Header->SetBranchAddress("H",H);       // height above see level  
    Header->SetBranchAddress("D",D);       // distance to obs. point
    Header->SetBranchAddress("X",X);       // slant depth     
  }
  else {
    Shower->SetBranchAddress("H",H);       // height above see level  
    Shower->SetBranchAddress("D",D);       // distance to obs. point
    Shower->SetBranchAddress("X",X);       // slant depth     
  }
  
  Header->GetEntry(0);


  TCanvas* cMu = gROOT->FindObject("cMu");
  if ( cMu != NULL ) delete cMu;
  cMu = new TCanvas("cMu", "profile for muons", 800, 10, 700, 500);

  for (int i=iStart; i<=(int)TMath::Min((double)Shower->GetEntries(), (double)iStop); ++i) {
    
    Shower->GetEntry(i);
    
    if(fitpars[0]==19){

    for (int j=0;j<nX; j++) {
      allmu->Fill(X[j],Mu[j]);	
    }

    }

  }
  allmu->SetMarkerStyle(20); 
  allmu->SetMarkerSize(.5); 
  allmu->SetMarkerColor(kBlue); 
  allmu->SetLineColor(kBlue); 
  allmu->SetTitle("");
  allmu->GetXaxis()->SetTitle("X [g/cm^{2}]");
  allmu->GetYaxis()->SetTitle("number of muons");
  allmu->GetYaxis()->SetTitleOffset(1.);
  allmu->Draw( "LP" );
  
  
}


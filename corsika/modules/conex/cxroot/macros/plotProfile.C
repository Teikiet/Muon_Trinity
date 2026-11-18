//
//  plotProfile()
//
//  ROOT macro to display simulated profile
//
//    usage: .x plotProfile.C("conex.root",3)
//
//       to plot profile of event no. 3 in file conex.root
//
void plotProfile(char* file="", int ievt=-1) {

  if ( file=="" || ievt<0 ) {
    cout << "\n usage: .x plotProfile.C(\"conex.root\",3) \n"
	 << "         to plot profile of event no. 3 in file contex.root\n" << endl;
    return;
  }


  // open file

  TFile* f = (TFile*) gROOT->GetListOfFiles()->FindObject(file);
  if (!f)
    f = new TFile(file);
  if (!f || f->IsZombie())
    return;

   gStyle->SetOptStat(0);

   ostringstream postScriptFile;
   postScriptFile << "plotProfile_event" << ievt << ".ps";

   // set and read tree

   TTree* Shower = (TTree*) f->Get("Shower");
   TTree* Header = (TTree*) f->Get("Header");
   int nevent=Shower->GetEntries();
   if ( ievt >= nevent ) {
     cout << " Error - requested event number exceeds number of events in file!"
          << "\n         ievt=" << ievt << " nevent=" << nevent
          << "\n exit ..." << endl;
     return;
   }
   const int maxX=20000;
   float X[maxX], H[maxX], D[maxX],
     N[maxX],eN[maxX],hN[maxX],dEdX[maxX], Mu[maxX], dMu[maxX];
   float fitpars[11],EGround[3],muThr,hadThr,emThr, outputVersion;
   int nX,iPart,iseed1,iseed2,iseed3,HEModel;

   Shower->SetBranchAddress("lgE", &fitpars[0]);
   Shower->SetBranchAddress("zenith", &fitpars[1]);
   Shower->SetBranchAddress("azimuth", &fitpars[10]);
   Shower->SetBranchAddress("Xfirst", &fitpars[2]);
   Shower->SetBranchAddress("X0", &fitpars[4]);
   Shower->SetBranchAddress("Xmax", &fitpars[9]);
   Shower->SetBranchAddress("Nmax", &fitpars[3]);
   Shower->SetBranchAddress("p1", &fitpars[5]);
   Shower->SetBranchAddress("p2", &fitpars[6]);
   Shower->SetBranchAddress("p3", &fitpars[7]);
   Shower->SetBranchAddress("chi2", &fitpars[8]);
   Shower->SetBranchAddress("nX", &nX);    // number of points in slant depth
   Shower->SetBranchAddress("N",N);       // particles(X)
   Shower->SetBranchAddress("Electrons",eN);       // electrons(X)
   Shower->SetBranchAddress("Hadrons",hN);       // hadrons(X)
   Shower->SetBranchAddress("dEdX",dEdX); // Energy deposit(X)
   Shower->SetBranchAddress("Mu",Mu);       // muons(X)
   Shower->SetBranchAddress("dMu",dMu);       // muon production(X)
   Shower->SetBranchAddress("EGround",EGround);
   Header->SetBranchAddress("muThr", &muThr);
   Header->SetBranchAddress("emThr", &emThr);
   Header->SetBranchAddress("hadThr", &hadThr);
   Header->SetBranchAddress("HEModel", &HEModel);  // high energy model
   Header->SetBranchAddress("Particle", &iPart);
   Header->SetBranchAddress("Seed1", &iseed1);
   Header->SetBranchAddress("OutputVersion", &outputVersion);
   Shower->SetBranchAddress("Seed2", &iseed2);
   Shower->SetBranchAddress("Seed3", &iseed3);

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
   Shower->GetEntry(ievt);


   TCanvas* v1 = gROOT->FindObject("v1");
   if ( v1 != NULL )
     delete v1;
   v1 = new TCanvas("v1","profile for charged particles",800,10,700,500);

   const double Xstart = fitpars[4];
   const double Xstop = Xstart+2000.;

   TGraph* graph1 = new TGraph(nX,X,N);

   graph1->SetMarkerStyle(20);
   graph1->SetMarkerSize(.7);
   graph1->SetTitle("");
   graph1->GetXaxis()->SetTitle("X [g/cm^{2}]");
   graph1->GetYaxis()->SetTitle("number of charged particles");
   graph1->GetYaxis()->SetTitleOffset(1.);
   graph1->Draw("AP");

   // plot conex profile fit

   TF1* GH = gROOT->FindObject("GHf");
   if (GH != NULL)
     delete GH;
   GH = new TF1("GHf", GaisserHillas, Xstart, Xstop, 6);
   GH->SetParameters(fitpars[3],fitpars[9],fitpars[4],fitpars[5],fitpars[6],fitpars[7]);
   GH->SetLineColor(2);
   GH->SetLineWidth(0.1);
   GH->SetNpx(1000);
   GH->Draw("Same");

   // draw shower parameters

   TLatex* leg = new TLatex();
   leg->SetNDC();
   leg->SetTextAlign(12);
   leg->SetTextSize(0.03);


   float xfirst = 0.15;
   float yfirst = 0.85;
   float dy = 0.05;

   if(iPart==100 || iPart==2212)
     leg->DrawLatex(xfirst,yfirst,"    proton");
   else if (iPart==5600)
     leg->DrawLatex(xfirst,yfirst,"     iron");
   else if (iPart==0 || iPart==22)
     leg->DrawLatex(xfirst,yfirst,"    gamma");
   else {
     char primary[10];
     int idmod=iPart%100;
     if (idmod==0)
       sprintf(primary,"A = %d",iPart/100);
     else
       sprintf(primary,"PDG id = %d",iPart);
     leg->DrawLatex(xfirst,yfirst,primary);
   }

   yfirst-=dy;
   if(HEModel == 1)
     leg->DrawLatex(xfirst,yfirst,"   (NeXuS)");
   else if (HEModel == 2)
     leg->DrawLatex(xfirst,yfirst,"   (QGSJet)");
   else if (HEModel == 6)
     leg->DrawLatex(xfirst,yfirst,"   (QGSJet II)");
   else if (HEModel == 5)
     leg->DrawLatex(xfirst,yfirst,"   (SIBYLL)");
   else if (HEModel == 4)
     leg->DrawLatex(xfirst,yfirst,"   (EPOS)");

   float y = yfirst-0.07;
   char legend[100];
   sprintf(legend,"lg(E/eV)=%5.2f",fitpars[0]);
   leg->DrawLatex(xfirst, y, legend); y-=dy;
   sprintf(legend,"zenith=%5.2f",fitpars[1]);
   leg->DrawLatex(xfirst, y, legend); y-=dy;
   sprintf(legend,"azimuth=%5.2f",fitpars[10]);
   leg->DrawLatex(xfirst, y, legend); y-=dy;
   sprintf(legend,"N_{max}=  %7.1e ",fitpars[3]);
   leg->DrawLatex(xfirst, y, legend); y-=dy;
   sprintf(legend,"X_{max}=%5.0f g/cm^{2}",fitpars[9]);
   leg->DrawLatex(xfirst, y, legend); y-=dy;
   sprintf(legend,"X_{0}    =%5.0f g/cm^{2}",fitpars[4]);
   leg->DrawLatex(xfirst, y, legend); y-=dy;
   sprintf(legend,"X_{int }  =%5.0f g/cm^{2}",fitpars[2]);
   leg->DrawLatex(xfirst, y, legend); y-=dy;
   sprintf(legend,"#lambda_{0}     =%6.0f g/cm^{2}",fitpars[5]);
   leg->DrawLatex(xfirst, y, legend); y-=dy;
   sprintf(legend,"#chi^{2}/Ndf     =%5.1e ",fitpars[8]);
   leg->DrawLatex(xfirst, y, legend); y-=dy;
   double ratio = fitpars[3]/pow(10.,fitpars[0])*1.e9;
   sprintf(legend,"N_{max}/(E/GeV)=%5.2f",ratio);
   leg->DrawLatex(xfirst, y, legend); y-=dy;
   sprintf(legend,"Seed(1)=%9.0f",iseed1);
   leg->DrawLatex(xfirst, y, legend); y-=dy;
   sprintf(legend,"Seed(2)=%9.0f",iseed2);
   leg->DrawLatex(xfirst, y, legend); y-=dy;
   sprintf(legend,"Seed(3)=%9.0f",iseed3);
   leg->DrawLatex(xfirst, y, legend); y-=dy;

   v1->Print((postScriptFile.str()+"(").c_str());

   // plot simulated profile for muons

   TCanvas* v2 = gROOT->FindObject("v2");
   if ( v2 != NULL )
     delete v2;
   v2 = new TCanvas("v2","number of muons vs altitude ",800,600,700,500);

   TGraph* graph2 = new TGraph(nX,H,Mu);

   // log scale
   gPad->SetLogx(1);

   graph2->SetMarkerStyle(20);
   graph2->SetMarkerSize(.7);
   graph2->SetTitle("");
   graph2->GetXaxis()->SetTitle("H [m]");
   graph2->GetYaxis()->SetTitle("number of Muons");
   graph2->GetYaxis()->SetTitleOffset(1.);
   graph2->Draw("AP");

   v2->Print(postScriptFile.str().c_str());

   // plot simulated profile for muon production rate

   TCanvas* v5 = gROOT->FindObject("v5");
   if (v5 != NULL)
     delete v5;
   v5 = new TCanvas("v5","muon production rate vs core distance ",400,500,700,500);

   TGraph* graph5 = new TGraph(nX,H,dMu);

   // log scale
   gPad->SetLogx(1);

   graph5->SetMarkerStyle(20);
   graph5->SetMarkerSize(.7);
   graph5->SetTitle("");
   graph5->GetXaxis()->SetTitle("L [m]");
   graph5->GetYaxis()->SetTitle("Muon production");
   graph5->GetYaxis()->SetTitleOffset(1.);
   graph5->Draw("AP");
   v5->Print(postScriptFile.str().c_str());

   // plot simulated Edepo profile

   TCanvas* v3 = gROOT->FindObject("v3");
   if ( v3 != NULL )
     delete v3;
   v3 = new TCanvas("v3","profile for dE/dX",10,600,700,500);

   TH1D* h3 = gROOT->FindObject("h3");
   if (h3 != NULL)
     delete h3;


   // X[] denotes the depth of the planes that are crossed by N[],
   // dEdX[] is inside the volume at average depth (X[i]+X[i+1])/2
   float dEdXdepth[maxX];
   const float deltaX = X[2] - X[1];
   for (int i=0; i<nX; ++i)
     dEdXdepth[i] = X[i] + deltaX/2;

   TH1F* graph3 = new TH1F("h3", "h3" ,nX-1, dEdXdepth);

   const double PeV = 1e6;
   for (int i=0; i<nX-1; ++i) {
     graph3->SetBinContent(i,dEdX[i]/PeV);
     graph3->SetBinError(i, 0);
   }
   // log scale
   //  gPad->SetLogx(1);
   //  gPad->SetLogy(1);

   graph3->SetMarkerStyle(20);
   graph3->SetMarkerSize(.7);
   graph3->SetTitle("");
   graph3->GetXaxis()->SetTitle("X [g/cm^{2}]");
   graph3->GetYaxis()->SetTitle("Energy deposit [PeV/(g/cm^{2})]");
   graph3->GetYaxis()->SetTitleOffset(1.);
   graph3->Draw("P");

   // draw Energy parameters
   leg = new TLatex();
   leg->SetNDC();
   leg->SetTextAlign(12);
   leg->SetTextSize(0.03);

   xfirst = 0.15;
   yfirst = 0.85;
   dy = 0.05;
   float dx = 0.03;

   if(iPart == 100 || iPart == 2212)
     leg->DrawLatex(xfirst+dx,yfirst,"    proton");
   else if (iPart == 5600)
     leg->DrawLatex(xfirst+dx,yfirst,"     iron");
   else if (iPart == 0 || iPart == 22)
     leg->DrawLatex(xfirst+dx,yfirst,"    gamma");
   else {
     int idmod=iPart%100;
     if (idmod == 0)
       sprintf(primary,"A = %d",iPart/100);
     else
       sprintf(primary,"PDG id = %d",iPart);
     leg->DrawLatex(xfirst+dx,yfirst,primary);
   }

   yfirst-=dy;
   if(HEModel == 1)
     leg->DrawLatex(xfirst,yfirst,"   (NeXuS)");
   else if (HEModel == 2)
     leg->DrawLatex(xfirst,yfirst,"   (QGSJet)");
   else if (HEModel == 6)
     leg->DrawLatex(xfirst,yfirst,"   (QGSJet II)");
   else if (HEModel == 5)
     leg->DrawLatex(xfirst,yfirst,"   (SIBYLL)");
   else if (HEModel == 4)
     leg->DrawLatex(xfirst,yfirst,"   (EPOS)");

   dx = 0.01;
   y = yfirst-0.07;
   char legend[100];
   sprintf(legend,"lg(E/eV)=%5.2f",fitpars[0]);
   leg->DrawLatex(xfirst, y, legend); y-=dy;
   sprintf(legend,"zenith=%5.2f",fitpars[1]);
   leg->DrawLatex(xfirst, y, legend); y-=dy;
   sprintf(legend,"E (GeV) at %5.0f g/cm^{2} :",X[nX-1]);
   leg->DrawLatex(xfirst, y, legend); y-=dy;
   sprintf(legend,"> e/m      : %8.2e",EGround[0]);
   leg->DrawLatex(xfirst+dx, y, legend); y-=dy;
   sprintf(legend,"> hadrons  : %8.2e",EGround[1]);
   leg->DrawLatex(xfirst+dx, y, legend); y-=dy;
   sprintf(legend,"> muons    : %8.2e",EGround[2]);
   leg->DrawLatex(xfirst+dx, y, legend); y-=dy;
   sprintf(legend,"E_{max}(CE) / E_{max}(MC) :");
   leg->DrawLatex(xfirst, y, legend); y-=dy;
   sprintf(legend,"> had. :  %7.1e",hadThr);
   leg->DrawLatex(xfirst+dx, y, legend); y-=dy;
   sprintf(legend,"> e/m :  %7.1e",emThr);
   leg->DrawLatex(xfirst+dx, y, legend); y-=dy;
   sprintf(legend,"> mu  :  %7.1e",muThr);
   leg->DrawLatex(xfirst+dx, y, legend); y-=dy;

   v3->Print(postScriptFile.str().c_str());

   // plot geometry
   TCanvas* v4 = gROOT->FindObject("v4");
   if (v4 != NULL)
     delete v4;
   v4 = new TCanvas("v4","H vs D",1150,600,400,300);

   TGraph* graph4 = new TGraph(nX,D,H);

   gPad->SetLogy(1);

   graph4->SetMarkerStyle(20);
   graph4->SetMarkerSize(.4);
   graph4->SetTitle("");
   graph4->GetXaxis()->SetTitle("D [m]");
   graph4->GetYaxis()->SetTitle("H [m]");
   graph4->GetYaxis()->SetTitleOffset(1.);
   graph4->Draw("AP");

   v4->Print((postScriptFile.str()+")").c_str());

   TCanvas* v6 = gROOT->FindObject("v6");
   if (v6 != NULL)
     delete v6;
   v6 = new TCanvas("v6","e and h  vs X",150,300,700,600);

   TGraph* graph7 = new TGraph(nX,X,eN);

   gPad->SetLogy(1);

   graph7->SetMarkerStyle(20);
   graph7->SetMarkerSize(.7);
   graph7->SetTitle("");
   graph7->GetXaxis()->SetTitle("X [g/cm^{2}]");
   graph7->GetYaxis()->SetTitle("number of electrons");
   graph7->GetYaxis()->SetTitleOffset(1.);
   graph7->Draw("AP");

   TGraph* graph8 = new TGraph(nX,X,hN);

   graph8->SetMarkerStyle(24);
   graph8->SetLineColor(2);
   graph8->SetMarkerColor(2);
   graph8->SetMarkerSize(.7);
   graph8->SetTitle("");
   graph8->GetXaxis()->SetTitle("X [g/cm^{2}]");
   graph8->GetYaxis()->SetTitle("number of hadrons");
   graph8->GetYaxis()->SetTitleOffset(1.);
   graph8->Draw("PSame");


   return;

 }

 //
 //  double GaisserHillas(double *x, double *par)
 //
 //   returns number of electrons at depth x[0]
 //   with Gaisser-Hillas parameters
 //
 //        Nmax=par[0]
 //        Xmax=par[1]
 //        X0  =par[2]
 //        l0  =par[3]
 //        l1  =par[4]
 //        l2  =par[5]
 //
 //      (lambda=l0+l1*x+l2*x*x )
 //
 double GaisserHillas(double* x, double* par)
 {
   double t =x[0];
   double anmax   = par[0];
   double tmax   = par[1];
   double t0      = par[2];
   double lambda    = par[3]+par[4]*t+par[5]*t*t;

   double expo = (tmax-t0)/lambda;
   double f = anmax*pow((t-t0)/(tmax-t0),expo);
   expo = (tmax-t)/lambda;
   f = f*exp(expo);

   return f;
 }


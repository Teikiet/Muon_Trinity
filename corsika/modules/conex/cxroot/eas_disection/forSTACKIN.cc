#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>

#include <boost/multiprecision/mpfr.hpp>

#include <boost/format.hpp> 
#include <boost/tuple/tuple.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;
namespace ublas = boost::numeric::ublas;
namespace mp = boost::multiprecision;

extern "C" {
  int idtrafocx_(const char* in, const char* out, const int&);
}

int 
CorsikaPID(const int Id) 
{
  int IdCors = 0;
  if (Id%100 != 0) {
    IdCors = idtrafocx_("nxs", "cor", Id);
  } else if(Id<0) { // then         !strangelet
    IdCors = Id;
  } else { // nucleus
    IdCors = Id + int(double(Id/100)/2.15 + 0.7);
  }
  return IdCors;
}

ublas::vector<double>
cross(const ublas::vector<double>& a, const ublas::vector<double>& b) 
{
  if (a.size() != 3) {
    cout << "NON CARTESIAN VECTOR a" << endl;
    exit (1);
  }
  if (b.size() != 3) {
    cout << "NON CARTESIAN VECTOR b" << endl;
    exit (1);
  }
  ublas::vector<double> c(3);
  c(0) = a(1)*b(2) - a(2)*b(1);
  c(1) = a(2)*b(0) - a(0)*b(2);
  c(2) = a(0)*b(1) - a(1)*b(0);
  return c;
}


/* Matrix inversion routine.
   Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
bool 
InvertMatrix(const ublas::matrix<double>& input, 
	     ublas::matrix<double>& inverse) 
{
  using namespace ublas;
  typedef permutation_matrix<std::size_t> pmatrix;
  // create a working copy of the input
  matrix<double> A(input);
  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());
  
  // perform LU-factorization
  int res = lu_factorize(A,pm);
  if( res != 0 ) return false;
  
  // create identity matrix of "inverse"
  inverse.assign(ublas::identity_matrix<double>(A.size1()));
  
  // backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);
  
  return true;
}


ublas::matrix<double> 
findBase(const ublas::vector<double>& p) {

  if (p.size() != 3) {
    cout << "NON CARTESIAN VECTOR" << endl;
    exit (1);
  }
  
  ublas::vector<double> ez = p;

  ublas::vector<double> non_ez = ez;
  non_ez[0] += 0.1;
  non_ez[1] += 0.5;
  non_ez[2] -= 0.5;
  
  ublas::vector<double> ex = cross(ez, non_ez);
  ublas::vector<double> ey = cross(ez, ex);
  
  ex /= norm_2(ex);
  ey /= norm_2(ey);
  ez /= norm_2(ez);
  
  ublas::matrix<double> rot(3,3);
  row(rot, 0) = ex;
  row(rot, 1) = ey;
  row(rot, 2) = ez;
  
  return rot;
}



boost::tuple<double,double>
  LorentzBoost(const double yboost, const double E, const double m, 
	       const double px, const double py, const double pz) 
{
  // boost from cm to lab system
  double amt = m*m + px*px + py*py;
  if (amt>0) 
    amt = sqrt(amt);
  else 
    amt = std::max(1.e-15, sqrt(std::abs((E+pz)*(E-pz)))); 
  const double y = ( pz<0 ? -1.0 : 1.0) * log((E+std::abs(pz))/amt) - yboost;
  //const double px_lab = px;
  //const double py_lab = py;
  const double pz_lab = amt*sinh(y);
  const double E_lab = amt*cosh(y);
  return boost::tuple<double,double>(E_lab, pz_lab);
}



void
forSTACKIN(const string& inputName) 
{
  //TFile* in = TFile::Open("full-interactions.50.root");
  TFile* in = TFile::Open(inputName.c_str());//"full-interactions.90.root");
  
  int iEvent = 0;
  do { // loop events
    
    cout << "Processing Event " << iEvent << endl;
    
    ostringstream partListName, projListName, intListName;
    partListName << "ParticleList_Event" << iEvent;
    projListName << "Projectile_Event" << iEvent;
    intListName << "InteractionList_Event" << iEvent;
    
    double Px = 0, Py = 0, Pz = 0, Energy = 0, Mass = 0;
    int Id = 0, status = 0, interactionCounter = 0;
    TTree* partList = (TTree*) in->Get(partListName.str().c_str());
    if (!partList) {
      cout << "TTree " << partListName.str() << " not found" << endl;
      break;
    }       
    partList->SetBranchAddress("Px", &Px);
    partList->SetBranchAddress("Py", &Py);
    partList->SetBranchAddress("Pz", &Pz);
    partList->SetBranchAddress("Energy", &Energy);
    partList->SetBranchAddress("mass", &Mass);
    partList->SetBranchAddress("id", &Id);
    partList->SetBranchAddress("status", &status);
    partList->SetBranchAddress("interactionCounter", &interactionCounter);
    
    double firstHeight = 0, firstEnergy = 0, firstMass = 0;
    double firstPx = 0, firstPy = 0, firstPz = 0; // only used for reference frame transformation
    TTree* projList = (TTree*) in->Get(projListName.str().c_str());
    if (!projList) {
      cout << "TTree " << projListName.str() << " not found" << endl;
      break;
    }
    projList->SetBranchAddress("height", &firstHeight);
    projList->SetBranchAddress("Energy", &firstEnergy);  
    projList->SetBranchAddress("mass", &firstMass);  
    projList->SetBranchAddress("Px", &firstPx);  
    projList->SetBranchAddress("Py", &firstPy);  
    projList->SetBranchAddress("Pz", &firstPz);  
    
    int firstTarget = 0;
    int firstId = 0;
    double firstSqrtS = 0;
    TTree* intList = (TTree*) in->Get(intListName.str().c_str());
    if (!intList) {
      cout << "TTree " << intListName.str() << " not found" << endl;
      break;
    }
    intList->SetBranchAddress("idProj", &firstId);
    intList->SetBranchAddress("idTarg", &firstTarget);
    intList->SetBranchAddress("eCMS", &firstSqrtS);
    
    intList->GetEntry(0);
    projList->GetEntry(0);
    
    
    
    // Define the reference frames for transformations
    
    ublas::vector<double> firstDir(3);
    firstDir(0) = firstPx;
    firstDir(1) = firstPy;
    firstDir(2) = firstPz;
    
    const ublas::matrix<double> ref1 = findBase(firstDir);
    ublas::matrix<double> ref2(3,3);
    if (!InvertMatrix(ref1, ref2)) {
      cout << "Matrix inversion failed" << endl;
      exit (2);
    }
    
    firstDir /= norm_2(firstDir);
    const double thetaShw = acos(firstDir(2)) * 180./M_PI; 
    const double phiShw = atan2(firstDir(1), firstDir(0)) * 180./M_PI; 
    
    cout << " dir=" << firstDir << endl;
    cout << " thetaShw=" << thetaShw << " phiShw=" << phiShw << endl;
    cout  << "ref1=" << ref1 << endl;
    cout  << "ref2=" << ref2 << endl;
    cout << " ref1*ref2=" << prod(ref1,ref2) << endl;
    
    // const double gamma = firstEnergy / firstMass;
    // const double beta = sqrt((gamma-1.)*(gamma+1.))/gamma;
    
  
    
    // Calculate the Lorentz Boost
    
    // center-of-mass system
  
    
    
    typedef mp::number<mp::mpfr_float_backend<300> >  value_type;
    
    //  firstEnergy = 4331.93/2;
    
    const double mTarg = 0.93827998638153076;
    
    const value_type E(firstEnergy), M(firstMass);
    const value_type P = sqrt((E-M)*(E+M));
    const value_type RAP = 0.5 * log((E+P)/(E-P)) / 2;
    const double rap = RAP.convert_to<double>();
    
    // P    xspnll=sqrt(max(0.d0,xselab**2-xsamproj**2))
    // const value_type xsengy = sqrt(2*xselab*xsamtarg+xsamtarg**2+xsamproj**2 )
    const double xsengy = sqrt(2*firstEnergy*mTarg + mTarg*mTarg + firstMass*firstMass);
    //const value_type xsengy = sqrt(2*E*firstTarget);
    // xsecms=xsengy
    // xsekin=xselab-xsamproj
    
    // const value_type xsyhaha = log((sqrt(xspnll**2+s)+xspnll)/sqrt(s))
    const value_type xsyhaha = log((sqrt(P*P + xsengy*xsengy) + P) / xsengy);
    
    
    cout << " firstEnergy=" << firstEnergy << " firstMass=" << firstMass << " firstP=" << P << " E-P=" << (firstEnergy-P) << endl;
    //cout << " gamma=" << gamma << " beta=" << beta << endl;
    cout << " rap=" << rap << endl;
    
    cout << " firstSqrtS=" << firstSqrtS << endl;
    cout << " xsengy=" << xsengy << endl;
    cout << " xsyhaha=" << xsyhaha.convert_to<double>() << endl;
    
    
    
    // check how many secondaries there are 
    
    double esum = 0;
    int nSec = 0;
    for (int iSec=0; iSec<partList->GetEntries(); ++iSec) {
      partList->GetEntry(iSec);
      if (status!=0 || interactionCounter!=1 || Energy-Mass<0.1) { // first interaction ONLY
	continue;
      }
      esum += Energy;
      nSec++;
    }
    cout << "nPart=" << partList->GetEntries() << endl;
    cout << "nSec=" << nSec << endl;
    cout << "Esum=" << esum << endl;
    cout << "Esum/nSec=";
    if (nSec>0)
      cout << "esum/nSec";
    else 
      cout << "-";
    cout << endl;
    
    
    
    // loop and process secondaries
    
    
    ostringstream hEFlowName, hNFlowName, hENFlowName;
    hEFlowName << "hEFlowName_" << iEvent;
    hNFlowName << "hNFlowName_" << iEvent;
    hENFlowName << "hENFlowName_" << iEvent;

    const Int_t nx = 6;    
    TH1D* hEFlow = new TH1D(hEFlowName.str().c_str(), "Energy Flow", nx, 0, nx);
    TH1D* hNFlow = new TH1D(hNFlowName.str().c_str(), "Particle Flow", nx, 0, nx);
    TH1D* hENFlow = new TH1D(hENFlowName.str().c_str(), "Average Energy", nx, 0, nx);
    
    hEFlow->SetLineWidth(3);
    hEFlow->SetLineColor(kBlue);
    hEFlow->GetXaxis()->SetLabelSize(0.07);
    hEFlow->GetYaxis()->SetLabelSize(0.05);
    hEFlow->GetYaxis()->SetTitleOffset(0.8);
    hEFlow->GetYaxis()->SetTitleSize(0.06);
    hNFlow->SetLineWidth(3);
    hNFlow->SetLineColor(kBlue);
    hNFlow->GetXaxis()->SetLabelSize(0.07);
    hNFlow->GetYaxis()->SetLabelSize(0.05);
    hNFlow->GetYaxis()->SetTitleOffset(0.8);
    hNFlow->GetYaxis()->SetTitleSize(0.06);
    hENFlow->SetLineWidth(3);
    hENFlow->SetLineColor(kBlue);
    hENFlow->GetYaxis()->SetLabelSize(0.05);
    hENFlow->GetXaxis()->SetLabelSize(0.07);
    hENFlow->GetYaxis()->SetTitleOffset(0.8);
    hENFlow->GetYaxis()->SetTitleSize(0.06);
    
    const char *label[nx] = {"Central", "Endcap", "Forward", "CASTOR+T2", "FSC", "ZDC"};
    for (int i=0; i<nx; ++i) {
      hEFlow->GetXaxis()->SetBinLabel(i+1, label[i]);
      hNFlow->GetXaxis()->SetBinLabel(i+1, label[i]);
      hENFlow->GetXaxis()->SetBinLabel(i+1, label[i]);
    }
    
    
    
    map< int, vector<int> > mapping;
    map< int, double > mappingEnergy;
    
    
    int nSecCount = 0;
    for (int iSec=0; nSecCount<nSec; ++iSec) {
      
      partList->GetEntry(iSec);
      
      if (status!=0 || interactionCounter!=1 || Energy-Mass<0.1) { // first interaction ONLY
	continue;
      }
      
      nSecCount++;
      
      
      ublas::vector<double> secP(3);
      secP(0) = Px;
      secP(1) = Py;
      secP(2) = Pz;
      
      ublas::vector<double> secPcm(3);
      // secPcm = prod(ref1, secP);
      secPcm = secP;
      
      boost::tuple<double,double> ep = LorentzBoost(rap, Energy, Mass, secPcm(0), secPcm(1), secPcm(2));
      
      const double pTot = sqrt(secPcm(0)*secPcm(0) + secPcm(1)*secPcm(1) + ep.get<1>()*ep.get<1>());
      const double theta = acos(ep.get<1>()/pTot);
      const double prap = -log(tan(theta/2.));
      
      const double absprap = abs(prap);
      int idDet = -1;
      if (absprap<1)        // 0 - 1. 
	idDet = 0;
      else if (absprap<3.5) // 1. - 3.5
	idDet = 1;
      else if (absprap<5)   // 3.5 - 5
	idDet = 2;
      else if (absprap<6.6) // 5 - 6.6
	idDet = 3;
      else if (absprap<8)   // 6.6 - 8
	idDet = 4;
      else                  // 8 - infty
	idDet = 5;  
      
      
      // text OUTPUT
      
      const int IdCors = CorsikaPID(Id);

      if (idDet==5) {

	switch (IdCors) {
	      
	case 13:
	case 1:
	case 7:
	  idDet=6;
	  break;

	}
      }

      
      cout << boost::format("%|1$5i|%|2$5i| %|3$ 15.6E| %|4$ 15.6E| %|5$ 15.6E| %|6$ 15.6E|") 
	% nSecCount % IdCors % Energy % Pz % Px % Py 
	// << " "  << secPcm(2) << " " << ep.get<0>() << " " << ep.get<1>()
	   << ", pseudo=" << setw(10)<< prap
	   << ", det=" << setw(3) << idDet
	   << endl;
      // from CORSIKA docu: stackin format: (2I5, 4(1X,E15.7))        
      
      mapping[idDet].push_back(iSec);
      if (mappingEnergy.count(idDet)) {
	mappingEnergy[idDet] += Energy;
      } else {
	mappingEnergy[idDet] = Energy;
      }

      hEFlow->Fill(idDet+0.1, Energy);
      hNFlow->Fill(idDet+0.1);

    }
    
    hENFlow->Add(hEFlow);
    hENFlow->Divide(hNFlow);

    hEFlow->SetYTitle("Total Energy Flow  [GeV]");
    hNFlow->SetYTitle("Number of Particles");
    hENFlow->SetYTitle("Avg. Energy per Particle  [GeV]");
    
    // find scale
    const double maxY = hENFlow->GetMaximum();
    double scale = 1;
    if (maxY/1e3<1) { // GeV
    } else if (maxY/1e6<1) { // TeV
      scale = 1e3;
      hEFlow->SetYTitle("Total Energy Flow  [TeV]");
      hENFlow->SetYTitle("Avg. Energy per Particle  [TeV]");
    } else if (maxY/1e9<1) { // PeV 
      scale = 1e6;
      hEFlow->SetYTitle("Total Energy Flow  [PeV]");
      hENFlow->SetYTitle("Avg. Energy per Particle  [PeV]");
    } else if (maxY/1e12<1) { // EeV 
      hEFlow->SetYTitle("Total Energy Flow  [EeV]");
      hENFlow->SetYTitle("Avg. Energy per Particle  [EeV]");
      scale = 1e9;
    } else if (maxY/1e15<1) { // ZeV 
      hEFlow->SetYTitle("Total Energy Flow  [ZeV]");
      hENFlow->SetYTitle("Avg. Energy per Particle  [ZeV]");
      scale = 1e12;
    }
    
    hEFlow->Scale(1./scale);
    hENFlow->Scale(1./scale);
      

    
    // write output files
    
    for (map<int, vector<int> >::const_iterator iterDet = mapping.begin(); iterDet != mapping.end(); ++iterDet) {
      
      const int idDet = iterDet->first;
      const vector<int>& secList = iterDet->second;
      
      ostringstream steer_name, stack_name;
      steer_name << "stackin_evt" << iEvent << "_det" << idDet << ".inp";
      stack_name << "stackin_part_evt" << iEvent << "_det" << idDet << ".dat";
      
      ofstream outStack(stack_name.str().c_str());
      
      outStack << "   " << secList.size() << " " << firstId << " " << firstHeight << " " << mappingEnergy[idDet] << " " << firstTarget << " " << 0 << " " << 0 << " " << 0 << endl;    
      
      for (int iSec=0; iSec<secList.size(); ++iSec) {
	
	partList->GetEntry(secList[iSec]);
	
	const int IdCors = CorsikaPID(Id);
	
	outStack << boost::format("%|1$5i|%|2$5i| %|3$ 15.6E| %|4$ 15.6E| %|5$ 15.6E| %|6$ 15.6E|") 
	  % (iSec+1) % IdCors % Energy % Pz % Px % Py 
	  //<< " "  << secPcm(2) << " " << ep.get<0>() << " " << ep.get<1>()
	  //<< " " << prap
		 << endl;
	// from CORSIKA docu: stackin format: (2I5, 4(1X,E15.7))        
      }
      outStack.close();  
      
      
      ofstream outSteer(steer_name.str().c_str());
      
      outSteer << "RUNNR  " << idDet << endl;
      outSteer << "THETAP " << thetaShw << "     " << thetaShw << endl;
      outSteer << "PHIP   " << phiShw << "      " << phiShw << endl;
      outSteer << "INFILE " << stack_name.str().c_str() << endl;
      outSteer << "FIXHEI " << firstHeight*100 << "   0 " << endl; //  18765.43e2      0
      
      ifstream inSteer("stackin_template.inp");
      while (inSteer.good()) {
	char oneLine[10000];
	inSteer.getline(oneLine, 10000);
	outSteer << oneLine << endl;
      }
      inSteer.close();
      
      outSteer.close();
      
    } // end loop mapping
    
    ostringstream cEname, cNname, cNEname;
    cEname << "distE_evt" << iEvent;
    cNname << "distN_evt" << iEvent;
    cNEname << "distEN_evt" << iEvent;
    
    TCanvas* cE = new TCanvas(cEname.str().c_str());
    cE->SetBottomMargin(0.12);
    hEFlow->Draw();
    cE->SaveAs((cEname.str()+".eps").c_str());
    cE->SaveAs((cEname.str()+".pdf").c_str());
    
    TCanvas* cN = new TCanvas(cNname.str().c_str());
    cN->SetBottomMargin(0.12);
    hNFlow->Draw();
    cN->SaveAs((cNname.str()+".eps").c_str());
    cN->SaveAs((cNname.str()+".pdf").c_str());
    
    TCanvas* cNE = new TCanvas(cNEname.str().c_str());
    cNE->SetBottomMargin(0.12);
    hENFlow->Draw();
    cNE->SaveAs((cNEname.str()+".eps").c_str());
    cNE->SaveAs((cNEname.str()+".pdf").c_str());

    iEvent++;
    
    
  } while(true); // end loop events
  
}





int 
main(int argc, char** argv)
{
  gStyle->SetOptStat(0);

  if (argc<2) {
    cout << "Must specify input conex files, run with CONEX_EXTENSIONS and first interaction tree output" << endl;
    return 1;
  }

  const string inputName(argv[1]);

  forSTACKIN(inputName);
  
  return 0;
}


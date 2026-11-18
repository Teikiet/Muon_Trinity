#include <TChain.h>
#include <TFile.h>

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
using namespace std;

void help() {
  cout << "\n please specify input data files! \n "
       << " merge(\"input_data_file\", \"output=all.root\"), [firstInteractionTree], [miniCONEX/special-options]) \n"
       << "   [firstInteractionTree] = \"on\" or \"off\" (default=off)\n" 
       << "   [miniCONEX/special-options] = \"mini\" or \"none\" (default=none)\n" 
       << endl;
}

/*
  Function to merge several input CONEX files into one output file
*/
void merge(const string& input_files="", 
	   const string& output="merged.root",
	   const string& firstIntTree="off",
	   const string& miniCONEX="off") {

  if (input_files=="" ||
      (firstIntTree!="on" && firstIntTree!="ON" && firstIntTree!="On" &&
       firstIntTree!="off" && firstIntTree!="OFF" && firstIntTree!="Off")) {
    help();
    return;
  }
  
  const bool firstInteractionTree = (firstIntTree=="on" ||
				     firstIntTree=="On" ||
				     firstIntTree=="ON");
  
  const bool isMiniCONEX = (miniCONEX=="on" ||
			    miniCONEX=="On" ||
			    miniCONEX=="ON");
  
  cout << " > Merging of FirstInteractionTree is "
       << (firstInteractionTree ? "ON" : "OFF")
       << endl;

  cout << " > Merging to CONEX-format: " << (isMiniCONEX ? "miniCONEX (Profile data will not be stored)" : "std"  )<< endl; 
  
  TFile* file = new TFile(output.c_str(), "RECREATE");
  cout << " Creating merged file: " << output << endl;
  TChain shower("Shower");
  TChain header("Header");
  int nShower = shower.Add(input_files.c_str());
  int nHeader = header.Add(input_files.c_str());
  if (nShower != nHeader) {
    cout << " Check you input files, some are corrupt! "  << endl;
    cout << " nShower= " << nShower << " nHeader= " << nHeader << endl;
  } else {
    cout << " adding " << nShower << " data files " << endl;
  }
  
  if (isMiniCONEX) {
      if (header.GetBranch("nX"))
        header.SetBranchStatus("nX", 0);
      if (header.GetBranch("X"))
	header.SetBranchStatus("X", 0);
      if (header.GetBranch("H"))
	header.SetBranchStatus("H", 0);
      if (header.GetBranch("D"))
	header.SetBranchStatus("D", 0);
      shower.SetBranchStatus("nX", 0);
      if (shower.GetBranch("X"))
	shower.SetBranchStatus("X", 0);
      shower.SetBranchStatus("N", 0);
      if (shower.GetBranch("H"))
	shower.SetBranchStatus("H", 0);
      if (shower.GetBranch("D"))
	shower.SetBranchStatus("D", 0);
      shower.SetBranchStatus("dEdX", 0);
      shower.SetBranchStatus("Mu", 0);
      shower.SetBranchStatus("Gamma", 0);
      shower.SetBranchStatus("Electrons", 0);
      shower.SetBranchStatus("dMu", 0);
      shower.SetBranchStatus("EGround", 0);      
  }
  
  shower.Merge(file, 2000, "C");
  cout << " finished merging Shower branch. N=" << shower.GetEntries() << endl;
  
  file = new TFile(output.c_str(), "UPDATE");
  header.Merge(file, 2000, "C");
  cout << " finished merging Header branch. N=" << header.GetEntries() << endl;

  TChain ch2("FirstInteraction");
  if (firstInteractionTree) {
    int nCh2 = ch2.Add(input_files.c_str());
    file = new TFile(output.c_str(), "UPDATE");
    ch2.Merge(file, 2000, "C");
    cout << " finished merging FirstInteraction branch " << endl;
  }
  
  cout << " DONE:  N_header=" << header.GetEntries() << ", N_shower=" << shower.GetEntries();
  if (firstInteractionTree) cout << ", N_firstIntTree=" << ch2.GetEntries();
  cout << "\n -------------------------------------------------------------------\n"
       << endl;
  
}


/*
  entry point for compiled program
 */
int main(int argc, char** argv) {
  
  if (argc<3) {
    help();
    return 1;
  }
  
  const string input(argv[1]), output(argv[2]);  
  const string firstInteractionTree( (argc<4 ? "off" : argv[3] ) );
  const string specialOpts( (argc<5 ? "none" : argv[4] ) );
  
  if (firstInteractionTree!="on" && firstInteractionTree!="off") {
    help();
    return 1;
  }
  
  if (specialOpts=="none")
    merge(input, output, firstInteractionTree);
  else if (specialOpts=="mini")
    merge(input, output, firstInteractionTree, "on");
  else {
    cerr << "\n Special option: \"" << specialOpts << "\" unknown " << endl;
  }

  return 0;
}

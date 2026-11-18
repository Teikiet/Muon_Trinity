/*************************************************************
*
*           conex.cc 
*
*  Main program for conex/CxRoot. 
*
**************************************************************/

#include <CxRoot.h>

#include <iostream>
#include <string>
#include <iomanip>
#include <cstdlib>
using namespace std;



/*************************************************************
*
*  main()
*
*     run with command line option -h to get a list of options
*
**************************************************************/

int
main(int argc, char** argv) 
{
  CxRoot cxRoot;
  if (!cxRoot.init(argc, argv))
    return 1;
  cxRoot.run();  
  return 0;
}



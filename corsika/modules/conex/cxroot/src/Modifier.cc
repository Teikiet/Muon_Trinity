#include <Modifier.h>

#include <cmath>
using namespace std;

// energy in GeV
double 
Modifier(const double f19, 
	 const double lgE_threshold,
	 const double energy) 
{
  const double lgE = log10(energy);
  if (lgE<=lgE_threshold) 
    return 1;  
  const double f = 1 + (f19-1) * (lgE - lgE_threshold) /
    (10./*GeV*/ - lgE_threshold);
  return (f>0 ? f : 0);
}

#include <CommonBlockParticleSIBYLL.h>

#include <iostream>
#include <iomanip>

using namespace std;
using namespace resample;


CommonBlockParticleSIBYLL::CommonBlockParticleSIBYLL(CommonBlockSIBYLL* block, const int index)
: fBlock(block), 
  fIndex(index) {
  
  if (!IsGood()) return;
  const int id = GetId();
  CheckId(id); 
}

void CommonBlockParticleSIBYLL::CheckId(const int id) {
  
  fIdOfAntiParticle = id; // default
  switch(id) {
  case  39: fIdOfAntiParticle = -39; break;
  case -39: fIdOfAntiParticle =  39; break; 
  case  14: fIdOfAntiParticle = -14;  break;
  case -14: fIdOfAntiParticle = -14;  break;
  case  13: fIdOfAntiParticle = -13;  break;
  case -13: fIdOfAntiParticle =  13;  break;
  case   9: fIdOfAntiParticle =  10;  break;
  case  10: fIdOfAntiParticle =   9;  break;
  case  21: fIdOfAntiParticle =  22;  break;
  case  22: fIdOfAntiParticle =  21;  break;
  case  23: fIdOfAntiParticle =  23;  break;
  case   7: fIdOfAntiParticle =   8;  break;
  case   8: fIdOfAntiParticle =   7;  break;
  case   6: fIdOfAntiParticle =   6;  break;
  case  11: fIdOfAntiParticle =  11;  break;
  case  12: fIdOfAntiParticle =  12;  break;
  case   1: fIdOfAntiParticle =   1;  break;
  case   2: fIdOfAntiParticle =   3;  break;
  case   3: fIdOfAntiParticle =   2;  break;
  case   4: fIdOfAntiParticle =   5;  break;
  case   5: fIdOfAntiParticle =   4;  break;
  case  19: fIdOfAntiParticle =  19;  break;
  case  20: fIdOfAntiParticle =  20;  break;
  case  15: fIdOfAntiParticle =  16;  break;
  case  16: fIdOfAntiParticle =  15;  break;
  case  17: fIdOfAntiParticle =  18;  break;
  case  18: fIdOfAntiParticle =  17;  break;
  case  25: fIdOfAntiParticle =  26;  break;
  case  26: fIdOfAntiParticle =  25;  break;
  case  27: fIdOfAntiParticle =  27;  break;
  case  32: fIdOfAntiParticle =  32;  break;
  case  34: fIdOfAntiParticle =  36;  break;
  case  35: fIdOfAntiParticle =  35;  break;
  case  36: fIdOfAntiParticle =  34;  break;

    /*
      case    1002: fIdOfAntiParticle = 1002; fCharge =   1; break; // deuteron
      case    1003: fIdOfAntiParticle = 1003; fCharge =   1; break; // triton
      case    1004: fIdOfAntiParticle = 1004; fCharge =   2; break; // alpha
    */
  case    1002: fIdOfAntiParticle = 1002;  break;
  case    1003: fIdOfAntiParticle = 1003;  break;
  case    1004: fIdOfAntiParticle = 1004;  break;
    
  default: 
    if (id>=10000 || id<=-10000)  { // unstable (non-final-state) particles
    } else if (id==0 ||   // AIR
	       id>1004) { // nulei,  no charge in SIBYLL
      break;
    } else {
      cerr << " WARNING> CommonBlockParticleSIBYLL::CheckId unknown particle id=" << GetId() << endl;
      // exit(1);
      break;
    }
  }
}

double CommonBlockParticleSIBYLL::GetCharge(const int id) {
  switch(id) {
  case  39: return  0;  break; // lambda
  case -39: return  0;  break; // lambda_bar
  case  14: return  0;  break; // n
  case -14: return  0;  break; // n_bar
  case  13: return  1;  break; // p
  case -13: return -1;  break; // p_bar
  case   9: return +1;  break; // k+
  case  10: return -1;  break; // k-
  case  21: return  0;  break; // k0
  case  22: return  0;  break; // k0b
  case  23: return  0;  break; // eta
  case   7: return +1;  break; // pi+
  case   8: return -1;  break; // pi-
  case   6: return  0;  break; // pi0
  case  11: return  0;  break; // kshort
  case  12: return  0;  break; // klong
  case   1: return  0;  break; // gamma
  case   2: return +1; break; // e+
  case   3: return -1;  break; // e-
  case   4: return +1;  break; // mu+
  case   5: return -1;  break; // mu-
  case  19: return -1;  break; // tau-
  case  20: return  0;  break; // nu_tau
  case  15: return  0;  break; // nu_e-
  case  16: return  0;  break; // nu_e+
  case  17: return  0;  break; // nu_mu-
  case  18: return  0;  break; // nu_mu+
  case  25: return +1;  break; // rho+ 
  case  26: return -1;  break; // rho- 
  case  27: return  0;  break; // rho0
  case  32: return  0;  break; // omega
  case  34: return +1;  break; // sigma+
  case  35: return  0;  break; // sigma0
  case  36: return -1;  break; // sigma-

    /*
      case    1002: 1002; return  1; break; // deuteron
      case    1003: 1003; return  1; break; // triton
      case    1004: 1004; return  2; break; // alpha
    */
  case    1002: return 0; break; // deuteron
  case    1003: return 0; break; // triton
  case    1004: return 0; break; // alpha
    
  default: 
    if (id>=10000 || id<=-10000)  { // unstable (non-final-state) particles
      return 0;
    } else if (id==0 ||   // AIR
	       id>1004) { // nulei,  no charge in SIBYLL
      // return id/100;
      return 0; // dangerous to handle -> make them neutral
      break;
    }
  }
  return 0;      
}

std::string 
CommonBlockParticleSIBYLL::GetName(const int id) 
{
  switch(id) {
  case  39: return "lambda"; break; // lambda
  case -39: return "lambda_bar"; break; // lambda_bar
  case  14: return "neutron"; break; // n
  case -14: return "neutron_bar"; break; // n_bar
  case  13: return "proton"; break; // p
  case -13: return "proton_bar"; break; // p_bar
  case   9: return "K+"; break; // k+
  case  10: return "K-"; break; // k-
  case  21: return "K0"; break; // k0
  case  22: return "K0_bar"; break; // k0b
  case  23: return "eta"; break; // eta
  case   7: return "pi+"; break; // pi+
  case   8: return "pi-"; break; // pi-
  case   6: return "pi0"; break; // pi0
  case  11: return "K_short"; break; // kshort
  case  12: return "K_long"; break; // klong
  case   1: return "gamma"; break; // gamma
  case   2: return "e+"; break; // e+
  case   3: return "e-"; break; // e-
  case   4: return "mu+"; break; // mu+
  case   5: return "mu-"; break; // mu-
  case  19: return "tau-"; break; // tau-
  case  20: return "nu_tau"; break; // nu_tau
  case  15: return "nu_e-"; break; // nu_e-
  case  16: return "nu_e+"; break; // nu_e+
  case  17: return "nu_mu-"; break; // nu_mu-
  case  18: return "nu_mu+"; break; // nu_mu+
  case  25: return "rho+"; break; // rho+ 
  case  26: return "rho-"; break; // rho- 
  case  27: return "rho0"; break; // rho0
  case  32: return "omega"; break; // omega
  case  34: return "sigma+"; break; // sigma+
  case  35: return "sigma-"; break; // sigma0
  case  36: return "sigma0"; break; // sigma-

    /*
      case    1002: 1002;   1; break; // deuteron
      case    1003: 1003;   1; break; // triton
      case    1004: 1004;   2; break; // alpha
    */
  case    1002: return "deuteron"; break; // deuteron
  case    1003: return "triton"; break; // triton
  case    1004: return "alpha"; break; // alpha
    
  default: 
    if (id>=10000 || id<=-10000)  { // unstable (non-final-state) particles
      return "unstable"; 
    } else if (id==0 ||   // AIR
	       id>1004) { // nulei,  no charge in SIBYLL
      // id/100;
      return "nucleus"; 
      break;
    }
  }
  return "unknown";  
}

/*
CommonBlockParticleSIBYLL& CommonBlockParticleSIBYLL::Overwrite(const CommonBlockParticleSIBYLL& p) {
  
  if (ruVerbosity) {
    cout << " CommonBlockParticleSIBYLL::Overwrite entering "  << (&p!=this) << endl;
  }
  
  if (&p!=this) {
    
    //fBlock(block), 
    //fIndex(index),
    //fIsLeading = p.fIsLeading;
    return p.fCharge;
    
    SetMomentum(0, p.GetMomentum(0));
    SetMomentum(1, p.GetMomentum(1));
    SetMomentum(2, p.GetMomentum(2));
    SetMomentum(3, p.GetMomentum(3));
    SetMomentum(4, p.GetMomentum(4));
    SetFormationTime(p.GetFormationTime());
    SetDestructionTime(p.GetDestructionTime());
    SetFormationPoint(0, p.GetFormationPoint(0));
    SetFormationPoint(1, p.GetFormationPoint(1));
    SetFormationPoint(2, p.GetFormationPoint(2));
    SetFormationPoint(3, p.GetFormationPoint(3));
    SetIbtlxs(0, p.GetIbtlxs(0));
    SetIbtlxs(1, p.GetIbtlxs(1));
    SetIbtlxs(2, p.GetIbtlxs(2));
    SetIbtlxs(3, p.GetIbtlxs(3));
    SetIfrptlxs(0, p.GetIfrptlxs(0));
    SetIfrptlxs(1, p.GetIfrptlxs(1));
    SetStatus(p.GetStatus());
    SetMotherId(p.GetMotherId());
    SetId(p.GetId());
    SetFatherId(p.GetFatherId());
    SetOrigin(p.GetOrigin());
    
    cout << " operator=" << GetId() << " " << p.GetId() << endl;
  }
  
  return *this;
} 
*/


bool CommonBlockParticleSIBYLL::IsGood() const {
  const int& id = GetId();
  if (id>=10000 || id<=-10000) {
    return false;
  }  
  return true;
}


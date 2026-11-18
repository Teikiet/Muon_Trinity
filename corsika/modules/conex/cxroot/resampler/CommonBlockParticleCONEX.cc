#include <CommonBlockParticleCONEX.h>

#include <iostream>
#include <iomanip>

using namespace std;
using namespace resample;

// from CONEX subroutine cxidmass(idi,amass)
const double CommonBlockParticleCONEX::fgMassProton = 0.93828;
const double CommonBlockParticleCONEX::fgMassNeutron = 0.93957;
const double CommonBlockParticleCONEX::fgMassLambda = 1.1156;
const double CommonBlockParticleCONEX::fgMassPionPM = 0.13957018;
const double CommonBlockParticleCONEX::fgMassPion0 = 0.1349766;
const double CommonBlockParticleCONEX::fgMassEta0 = 0.547853;
const double CommonBlockParticleCONEX::fgMassKaon0 = 0.497614;
const double CommonBlockParticleCONEX::fgMassKaonPM = 0.493677;
const double CommonBlockParticleCONEX::fgMassElectron = 0.5109989e-3;
const double CommonBlockParticleCONEX::fgMassMuon = 0.105658;
const double CommonBlockParticleCONEX::fgMassTau = 1.777; // CHECK...


CommonBlockParticleCONEX::CommonBlockParticleCONEX(CommonBlockCONEX* block, int index)
: fBlock(block), 
  fIndex(index) {
  
  if (!IsGood()) 
    return;
  
  const int id = GetId();
  CheckId(id);
}

void 
CommonBlockParticleCONEX::CheckId(const int id) 
{  
  fIdOfAntiParticle = -id; // default
  switch(id) {
  case  1111: fIdOfAntiParticle = id; break; 
  case  1121: fIdOfAntiParticle = 2221; break; 
  case  1221: fIdOfAntiParticle = id; break; 
  case  2221: fIdOfAntiParticle = 1121; break; 
  case   330: fIdOfAntiParticle = id; break; 
  case   331: fIdOfAntiParticle = id; break; 
  case   220: fIdOfAntiParticle = id; break; // eta
  case   221: fIdOfAntiParticle = id; break; // omega
  case   111: fIdOfAntiParticle = id; break; // rho0
  case   110: fIdOfAntiParticle = id; break; // pi0
  case    10: fIdOfAntiParticle = id; break; // gamma
  case    17: fIdOfAntiParticle = id; break; // deuteron
  case    18: fIdOfAntiParticle = id; break; // triton
  case    19: fIdOfAntiParticle = id; break; // alpha
  default: 
    if (id>100 && id%100==0) { // nulei, no charge in CONEX
      break;
    }
  }
}


double 
CommonBlockParticleCONEX::GetCharge(const int id) 
{
  switch(id) {
  case  2130: return 0; break; // lambda
  case -2130: return  0; break; // lambda_bar
  case  1220: return  0; break; // n
  case -1220: return  0; break; // n_bar
  case  1120: return  1; break; // p
  case -1120: return -1; break; // p_bar
  case  1111: return 0/*+2*/; break; 
  case  1121: return +1; break; 
  case  1221: return  0; break; 
  case  2221: return -1; break; 
  case  1330: return  0; break; 
  case -1330: return  0; break; 
    /* not in conex
       case -1111: return -2; fName = "delta++_bar"; break; 
       case -1121: return -1; fName = "delta+_bar"; break; 
       case -1221: return  0; fName = "delta0_bar"; break; 
       case -2221: return +1; fName = "delta-_bar"; break; 
    */
  case  2330: return -1; break;
  case -2330: return -1; break;
  case  1230: return  0; break; //sigma0
  case -1230: return  0; break; //sigma0
  case  1130: return +1; break; //sigma+
  case -1130: return -1; break;
  case  2230: return -1; break; //sigma-
  case -2230: return +1; break; //sigma-_bar
  case   231: return  0; break; 
  case  -231: return  0; break; 
  case   330: return  0; break; 
  case   331: return  0; break; 
  case   131: return +1; break; 
  case  -131: return -1; break; 
  case   130: return +1; break; // k+
  case  -130: return -1; break; // k-
  case   230: return  0; break; // k0
  case  -230: return  0; break; // k0b
  case   220: return  0; break; // eta
  case   221: return 0/*-1*/; break; // omega
    // case  -221: return +1; fName = "Omega_bar"; break; // omega (NOT EXISTING IN CONEX)
  case   120: return +1; break; // pi+
  case  -120: return -1; break; // pi-
  case   111: return  0; break; // rho0
  case   121: return +1; break; // rho+ 
  case  -121: return -1; break; // rho- 
  case   110: return  0; break; // pi0
  case   -20: return  0; break; // kshort
  case    20: return  0; break; // klong
  case    10: return  0; break; // gamma
  case   -12: return +1; break; // e+
  case    12: return -1; break; // e-
  case   -14: return +1; break; // mu+
  case    14: return -1; break; // mu-
  case    16: return -1; break; // tau-
  case    15: return  0; break; // nu_tau
  case    11: 
  case   -11: return  0; break; // nu_e
  case    13: 
  case   -13: return  0; break; // nu_mu

    /*
      case    17: return  1; fName = ""; break; // deuteron
      case    18: return  1; fName = ""; break; // triton
      case    19: return  2; fName = ""; break; // alpha
    */
  case    17: return  0; break; // deuteron
  case    18: return  0; break; // triton
  case    19: return  0; break; // alpha

  default: 
    if (id>100 && id%100==0) { // nulei, no charge in CONEX
          // returnid/100;
      return 0; // dangerous to handle -> make them neutral
      break;
    }
  }
  return 0;
}

std::string 
CommonBlockParticleCONEX::GetName(const int id) 
{
  switch(id) {
  case  2130: return "lambda"; break; // lambda
  case -2130: return "lambda_bar"; break; // lambda_bar
  case  1220: return "neutron"; break; // n
  case -1220: return "neutron_bar"; break; // n_bar
  case  1120: return "proton"; break; // p
  case -1120: return "proton_bar"; break; // p_bar
  case  1111: return "delta++"; break; 
  case  1121: return "delta+"; break; 
  case  1221: return "delta0"; break; 
  case  2221: return "delta-"; break; 
  case  1330: return "xi0"; break; 
  case -1330: return "xi0_bar"; break; 
    /* not in conex
       case -1111: fCharge =  -2; return "delta++_bar"; break; 
       case -1121: fCharge =  -1; return "delta+_bar"; break; 
       case -1221: fCharge =   0; return "delta0_bar"; break; 
       case -2221: fCharge =  +1; return "delta-_bar"; break; 
    */
  case  2330: return "xi-"; break;
  case -2330: return "xi-_bar"; break;
  case  1230: return "sigma0"; break; //sigma0
  case -1230: return "sigma0_bar"; break; //sigma0
  case  1130: return "sigma+"; break; //sigma+
  case -1130: return "sigma+_bar"; break;
  case  2230: return "sigma-"; break; //sigma-
  case -2230: return "sigma-_bar"; break; //sigma-_bar
  case   231: return "K*0"; break; 
  case  -231: return "K*0_bar"; break; 
  case   330: return "etaprime"; break; 
  case   331: return "phi"; break; 
  case   131: return "K*+"; break; 
  case  -131: return "K*-"; break; 
  case   130: return "K+"; break; // k+
  case  -130: return "K-"; break; // k-
  case   230: return "K0"; break; // k0
  case  -230: return "K0_bar"; break; // k0b
  case   220: return "eta"; break; // eta
  case   221: return "Omega"; break; // omega
    // case  -221: fCharge =  +1; return "Omega_bar"; break; // omega (NOT EXISTING IN CONEX)
  case   120: return "pi+"; break; // pi+
  case  -120: return "pi-"; break; // pi-
  case   111: return "rho0"; break; // rho0
  case   121: return "rho+"; break; // rho+ 
  case  -121: return "rho-"; break; // rho- 
  case   110: return "pi0"; break; // pi0
  case   -20: return "K_short"; break; // kshort
  case    20: return "K_long"; break; // klong
  case    10: return "gamma"; break; // gamma
  case   -12: return "e+"; break; // e+
  case    12: return "e-"; break; // e-
  case   -14: return "mu+"; break; // mu+
  case    14: return "mu-"; break; // mu-
  case    16: return "tau-"; break; // tau-
  case    15: return "nu_tau"; break; // nu_tau
  case    11: 
  case   -11: return "nu_e"; break; // nu_e
  case    13: 
  case   -13: return "nu_mu"; break; // nu_mu

    /*
      case    17: fCharge =   1; return ""; break; // deuteron
      case    18: fCharge =   1; return ""; break; // triton
      case    19: fCharge =   2; return ""; break; // alpha
    */
  case    17: return "deuteron"; break; // deuteron
  case    18: return "triton"; break; // triton
  case    19: return "alpha"; break; // alpha

  default: 
    if (id>100 && id%100==0) { // nulei, no charge in CONEX
      return "nucleus"; 
    }
  }
  return "unknown";
}


double 
CommonBlockParticleCONEX::GetMass(const int id) 
{
  switch(id) {
  case  2130:  // lambda
  case -2130: return fgMassLambda; break; // lambda_bar
  case  1220: 
  case -1220: return fgMassNeutron; break; // n_bar
  case  1120: 
  case -1120: return fgMassProton; break; // p_bar
    // case  1111: return 0/*+2*/; break; 
    // case  1121: return +1; break; 
    //case  1221: return  0; break; 
    //case  2221: return -1; break; 
    //case  1330: return  0; break; 
    //case -1330: return  0; break; 
    /* not in conex
       case -1111: return -2; fName = "delta++_bar"; break; 
       case -1121: return -1; fName = "delta+_bar"; break; 
       case -1221: return  0; fName = "delta0_bar"; break; 
       case -2221: return +1; fName = "delta-_bar"; break; 
    */
    //case  2330: return -1; break;
    //case -2330: return -1; break;
    //case  1230: return  0; break; //sigma0
    //case -1230: return  0; break; //sigma0
    //case  1130: return +1; break; //sigma+
    //case -1130: return -1; break;
    //case  2230: return -1; break; //sigma-
    //case -2230: return +1; break; //sigma-_bar
    //case   231: return  0; break; 
    //case  -231: return  0; break; 
    //case   330: return  0; break; 
    //case   331: return  0; break; 
    //case   131: return +1; break; 
    //case  -131: return -1; break; 
  case   130: 
  case  -130: return fgMassKaonPM; break; // k-
  case   230: return fgMassKaon0; break; // k0
    //case  -230: return  0; break; // k0b
  case   220: return fgMassEta0; break; // eta
    //case   221: return 0/*-1*/; break; // omega
    //// case  -221: return +1; fName = "Omega_bar"; break; // omega (NOT EXISTING IN CONEX)
  case   120: 
  case  -120: return fgMassPionPM; break; // pi-
    //case   111: return  0; break; // rho0
    //case   121: return +1; break; // rho+ 
    //case  -121: return -1; break; // rho- 
  case   110: return fgMassPion0; break; // pi0
    //case   -20: return  0; break; // kshort
    //case    20: return  0; break; // klong
  case    10: return  0; break; // gamma
  case   -12:
  case    12: return fgMassElectron; break; // e-
  case   -14: 
  case    14: return fgMassMuon; break; // mu-
  case    16: 
  case    15: return fgMassTau; break; // nu_tau
  case    11: 
  case   -11: return  0; break; // nu_e
  case    13: 
  case   -13: return  0; break; // nu_mu

    /*
      case    17: return  1; fName = ""; break; // deuteron
      case    18: return  1; fName = ""; break; // triton
      case    19: return  2; fName = ""; break; // alpha
    */
  case    17: return  0; break; // deuteron
  case    18: return  0; break; // triton
  case    19: return  0; break; // alpha

  default: 
    if (id>100 && id%100==0) { // nulei, no charge in CONEX
          // returnid/100;
      return 0; // dangerous to handle -> make them neutral
      break;
    }
  }
  return 0;
}


resample::CommonBlockParticleCONEX::EParticleId
CommonBlockParticleCONEX::GetIdOfAntiParticle(const resample::CommonBlockParticleCONEX::EParticleId id) 
{
  int IdOfAntiParticle = -int(id); // default
  switch (int(id)) {
  case  1111: IdOfAntiParticle = id; break; 
  case  1121: IdOfAntiParticle = 2221; break; 
  case  1221: IdOfAntiParticle = id; break; 
  case  2221: IdOfAntiParticle = 1121; break; 
  case   330: IdOfAntiParticle = id; break; 
  case   331: IdOfAntiParticle = id; break; 
  case   220: IdOfAntiParticle = id; break; // eta
  case   221: IdOfAntiParticle = id; break; // omega
  case   111: IdOfAntiParticle = id; break; // rho0
  case   110: IdOfAntiParticle = id; break; // pi0
  case    10: IdOfAntiParticle = id; break; // gamma
  case    17: IdOfAntiParticle = id; break; // deuteron
  case    18: IdOfAntiParticle = id; break; // triton
  case    19: IdOfAntiParticle = id; break; // alpha
  default: 
    if (id>100 && id%100==0) { // nulei, no charge in CONEX
      break;
    }
  }
  return resample::CommonBlockParticleCONEX::EParticleId(IdOfAntiParticle);
}

/*
CommonBlockParticleCONEX& CommonBlockParticleCONEX::Overwrite(const CommonBlockParticleCONEX& p) {
  
  if (ruVerbosity) {
    cout << " CommonBlockParticleCONEX::Overwrite entering "  << (&p!=this) << endl;
  }
  
  if (&p!=this) {
    
    //fBlock(block), 
    //fIndex(index),
    //fIsLeading = p.fIsLeading;
    fCharge = p.fCharge;
    
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


bool 
CommonBlockParticleCONEX::IsGood() const 
{
  const int id = GetId();
  if (id==0 ||                                                     // crap
      GetStatus()!=0 ||                                            // no (valid) secondary of the current interaction
      // GetStatus()>1 ||                                             // no (valid) secondary of the current interaction
      ((id==1120 || id==1220) && (GetEnergy()-GetMass()<0.002))) { // target spectator     
    return false;
  }  
  return true;
}


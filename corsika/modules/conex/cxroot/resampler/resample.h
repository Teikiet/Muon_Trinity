/**

   Ralf Ulrich Wed May  9 15:59:11 CEST 2007
   
   started multiplicity resampling after HE-interaction
 
 */

#ifndef __INCLUDE_RESAMPLE_H__
#define __INCLUDE_RESAMPLE_H__

//#define DEBUG_RESAMPLING
#undef DEBUG_RESAMPLING

#include <ResamplingMode.h>

namespace resample {
  
  class ParticleBlock;

  void InitResampling(const int seed);

  // The resampling interface function
  void SecondaryResampling(const resample::EResampleMode mode, const double factor, const int projId,
			   resample::ParticleBlock& particlesList);
  
  // 
  // The internal function to really to specific resampling tasks
  void MultiplicityResampling(const double factor, resample::ParticleBlock& particlesList);
  void ElasticityResampling(const double factor, resample::ParticleBlock& particlesList);
  void EMRatioResampling(const double factor, resample::ParticleBlock& particlesList);
  void ChargeRatioResampling(const double factor, const int projId, resample::ParticleBlock& particlesList);
  void Pi0SpectrumResampling(const double factor, resample::ParticleBlock& particlesList);
  void Rho0Resampling(const double factor, const int projId, resample::ParticleBlock& particlesList);
  void BaryonProductionResampling(const double factor, const int projId, resample::ParticleBlock& particlesList);


  bool SetNewPId(const resample::CommonBlockParticleCONEX::EParticleId pid, resample::ParticleBlockEntry& prt);

}


#endif


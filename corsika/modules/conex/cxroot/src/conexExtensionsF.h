#ifndef _include_conex_extensionsf_h_
#define _include_conex_extensionsf_h_

#include <conexConfig.h>

#ifdef CONEX_EXTENSIONS

class NXParticle;
class CXParticle;
class NXBlock;
class COSSINS;

extern "C" {

  void write_seed_cx_(int seed[3]);
  void set_seed_cx_(int seed[3]);
  void get_particle_from_list_cx_(NXBlock& p, int& i, int& goon);
  void get_particle_from_list_(CXParticle& p, int& goon);
  void get_projectile_(CXParticle& p, COSSINS& c);
  void write_particle_(int& interactionCounter, CXParticle &p);
  void write_particle_cx_(int& interactionCounter, NXBlock& data, int& indexF);

  void write_projectile_(int& interactionCounter, CXParticle &p, COSSINS& c);
  void write_interaction_(int& interactionCounter,
			  double& eProj,
			  double& eCMS,
			  double& eProd,
			  int& idProj,
			  int& idTarg,
			  int& mult);
  
}

#endif

#endif // include protection

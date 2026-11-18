/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <string>

#include <rng_decl.hpp>
#include "epos-lhcr-types.hpp"

/**
 * \file epos-lhcr.hpp
 *
 * Interface file for the EPOS LHCR library.
 */

DECLARE_RNG(epos)

extern "C" {

// random number and seed functions
void ranfst_(int& seed);
void ranfgt_(int& seed);
void rmmard_(double rvec[], int const* lenv, int const* /*iseq*/);
void rmmaqd_(int seed[3], int&, char*, int);

// additional random number functions
void ranfini_(double&, int&, int&);
void ranfcv_(double&);

float rangen_();
double drangen_();

void aaset_(int&);
void atitle_();
void readidtable_();
void hnbcreate_();
void hnbdestroy_();
void ainit_();
void aepos_(int&);
void afinal_();
void alistf_(char* str, int str_length); // hidden str length
void hnbspd_(int&);
void hnbpajini_();
void conini_();
void psaini_();

void paramini_(int&);
void xsigma_();

//
//  cross section from tables
//
// c------------------------------------------------------------------------------
//   function eposcrse(ek,mapro,matar,id)
// c------------------------------------------------------------------------------
// c inelastic cross section of epos
// c (id=0 corresponds to air)
// c ek     - kinetic energy for the interaction in the lab
// c maproj - projec mass number     (1<maproj<64)
// c matarg - target mass number     (1<matarg<64)
// c------------------------------------------------------------------------------
float eposcrse_(float&, int&, int&, int&);
float eposelacrse_(float&, int&, int&, int&);

// calculate cross section
// c------------------------------------------------------------------------------
//       subroutine crseaaEpos(sigt,sigi,sigc,sige)
// c----------------------------------eposela--------------------------------------------
// c nucleus-nucleus (hadron) cross section of epos from simplified (realistic)
// c simulations
// c (id=0 corresponds to air)
// c  sigt = sig tot
// c  sigi = sig inelastic (cut + projectile diffraction)
// c  sigc = sig cut
// c  sige = sig elastic (includes target diffraction)
// c------------------------------------------------------------------------------
void crseaaepos_(float&, float&, float&, float&);

void emsaaa_(int&);
void gakfra_(int&, int&);
void utghost_(int&);
void bjinta_(int&);
void utresc_(int&);

void emsfrag_(int&);

// get particles hadron class: meson, baryon etc..???
void iclass_(int&, int&);
// get charge for id
void idchrg_(int&, int&);
// get isospin, spin, strangeness for id
void idspin_(int&, int&, int&, int&);
// get mass for id
void idmass_(int&, float&);
// convert id from one format to another
int idtrafo_(char[3], char[3], int&);

// common blocks as
// defined in epos.inc

extern EPOS_LHCR::CUNCERTMU cuncertmu_;
extern EPOS_LHCR::CHACAS chacas_;
extern EPOS_LHCR::LHCTUNE lhctune_;
extern EPOS_LHCR::CICNT cicnt_;
extern EPOS_LHCR::HADR6 hadr6_;
extern EPOS_LHCR::HADR16 hadr16_;
extern EPOS_LHCR::NUCL6 nucl6_;
extern EPOS_LHCR::CJINTI cjinti_;
extern EPOS_LHCR::NXSAIR nxsair_;
extern EPOS_LHCR::APPLI appli_;
extern EPOS_LHCR::XSAPPLI xsappli_;
extern EPOS_LHCR::EVENTS events_;
extern EPOS_LHCR::XSEVENT xsevent_;
extern EPOS_LHCR::METR1 metr1_;
extern EPOS_LHCR::OTHE2 othe2_;
extern EPOS_LHCR::METR7 metr7_;
extern EPOS_LHCR::HAD12 had12_;
extern EPOS_LHCR::NUCL1 nucl1_;
extern EPOS_LHCR::CHADRON chadron_;
extern EPOS_LHCR::HADR2 hadr2_;
extern EPOS_LHCR::HADR25 hadr25_;
extern EPOS_LHCR::LEPT1 lept1_;
extern EPOS_LHCR::ENRGY enrgy_;
extern EPOS_LHCR::HADR1 hadr1_;
extern EPOS_LHCR::CEVT cevt_;
extern EPOS_LHCR::CSEED cseed_;
extern EPOS_LHCR::OTHE1 othe1_;
extern EPOS_LHCR::FILES files_;
extern EPOS_LHCR::FNAME fname_;
extern EPOS_LHCR::NFNAME nfname_;
extern EPOS_LHCR::PRNT1 prnt1_;
extern EPOS_LHCR::PRNT3 prnt3_;
extern EPOS_LHCR::CPTL cptl_;
extern EPOS_LHCR::HADR5 hadr5_;
extern EPOS_LHCR::NODCY nodcy_;
extern EPOS_LHCR::HAD10 had10_;
}

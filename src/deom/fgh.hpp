#pragma once
/**
 * writtened by Zhijun Pan
 * 2018.3.26
 * -------------------
 * FGH-DVR for eigenstate 
 * of Morse potential
 */
#ifndef FGH_H_
#define FGH_H_
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include "armadillo"
#include "json11.hpp"
#include "deomConst.hpp"

static const double    cm2unit = 4.5554927e-6;
static const double    kt2unit = 3.1662e-6;
static const double    fs2unit = 41.34902;
static const double    am2unit = 1822.8885;
static const double    ar2unit = 1.88964;//Angstrom
static const int    mmax = 3;           //Maximum 3D potential

using namespace json11;
//Encode maybe faster
//and let's try
int encode (const uvec i_list,const uvec cum_prod,const int nlen); 
//Decode should be accelerate.
uvec decode (int i_rank,const fvec cum_divd,const uvec cum_prod,const int nlen); 

double get_pot_1d (const double m0,const double de,const double beta,const double r0,const double r); 
//double get_pot_1d (double m0,Json j,double q1,double q2,double de,double beta,double r0,double r){
 /*   double lamd=0;
    double r2=r*r;
    for (int idru=0; idru<(len(j)); ++idru)
        lamd +=j[idru]["gamd"].number_value()*j[idru]["lamd"].number_value()*cm2unit*cm2unit;
*/   

double get_pot_2d (const double gxy,const double  x,const double  y); 

mat get_ham_1d(const Json m);
sp_mat mk_T(const field<mat>& hams,const uvec& i_list,const uvec& nr,const uvec& cum,const int loc,const int nmds);
sp_mat fgh_dvr(const int ndim,const int nmds,const uvec& nr,const fvec& ri,const fvec& rf,const mat& vmn,const field<mat>& hams);
#endif

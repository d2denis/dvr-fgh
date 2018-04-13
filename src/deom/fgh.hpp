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

static const int    mmax = 3;           //Maximum 3D potential

using namespace json11;
//Encode: (i,j,k) -> (ijk)
int encode (const uvec i_list,const uvec cum_prod,const int nlen);
//Decode: (ijk) -> (i,j,k)
uvec decode (int i_rank,const fvec cum_divd,const uvec cum_prod,const int nlen);

double get_pot_1d (const double m0,const double de,const double beta,const double r0,const double r);
//double get_pot_1d (double m0,Json j,double q1,double q2,double de,double beta,double r0,double r){
 /*   double lamd=0;
    double r2=r*r;
    for (int idru=0; idru<(len(j)); ++idru)
        lamd +=j[idru]["gamd"].number_value()*j[idru]["lamd"].number_value()*cm2unit*cm2unit;
*/
double get_energy(const Json m,const uword ie=3);
double get_pot_2d (const double gxy,const double  x,const double  y);

mat get_ham_1d(const Json m);
sp_mat fgh_dvr(const int ndim,const int nmds,const uvec& nr,const fvec& ri,const fvec& rf,const mat& vmn,const field<mat>& hams);
void printime(const std::string & strout);
void wave_trun(mat& wavetrun,const mat& wavefun, const uword nmds, const uword ne,const fvec& Eie,const uvec& nr,const fvec& ri,const fvec& rf,const mat& vmn,const field<mat>& hams);
#endif

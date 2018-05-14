/**
 * writtened by Zhijun Pan
 * 2018.5.9
 * -------------------
 * My deom response function
 * no output to any screen/file
 */
#include <ctime>
#include <chrono>
#include <iomanip>

#include <iostream>
#include <sstream>
#include <string>
#include "deom.hpp"

using namespace std;
using namespace arma;
 static    const int npsd = 1;//
 static    const int pade = 2;//
 static    const double temp = 9.5e-4;//
 static    const int lmax = 5;//
 static    const int nmax = 500;//
 static    const double ferr = 2e-9;//
 static    const double w1_max = 0.01214;//
 static    const int    nt1 = 400;//
 static    const int num_ob =nt1/2;
 static    const double thedt = 10.0;//
 static    const double staticErr = 2e-5;//
 static    const int    nk = 32;//
 static      const string sch_hei = "s";

void resp1 (vec& fw_im, const double w_max, const int mynt, const double mydt,
              const double staticErr, const int nk,
              const mat& sdip, const cube& pdip, const vec& bdip,
              const syst& s, const bath& b, const hidx& h) {

    const double dt1 = 2.0*deom_pi/w_max;
    const int    mt  = floor(dt1/mydt);
    const double dt1_res = dt1-mydt*mt;

    deom d1(s,b,h);

    cx_cube rho_t0 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);

    const mat& exph= expmat(-real(d1.ham1)/d1.temperature);
    rho_t0.slice(0).set_real(exph/trace(exph));
    d1.equilibrium (rho_t0,mydt,staticErr,nk);

    cx_cube rho_t1 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);
    d1.oprAct(rho_t1,sdip,pdip,bdip,rho_t0,'c');

    cx_vec ft = zeros<cx_vec>(mynt);

        for (int it=0; it<mynt; ++it) {

            ft(it) = deom_ci*d1.Trace(sdip,pdip,bdip,rho_t1);

//            printf ("In sch-pic: it=%d, nddo=%d, lddo=%d\n", it, d1.nddo, d1.lddo);
            double t1 = it*dt1;
            for (int jt=0; jt<mt; ++jt) {
                d1.rk4 (rho_t1,t1,mydt);
                t1 += mydt;
            }
            d1.rk4 (rho_t1,t1,dt1_res);
        }
    ft(0) *= 0.0;//Infinity(nan) occur nddo>17, should be fixed
    const cx_vec& fw = ifft(ft)*mynt*dt1;
    fw_im = imag(fw.rows(0,mynt/2-1));
}

mat matprod(cube a, vec b){
    mat result = zeros<mat>(size(a.slice(0)));
    uword nmds = a.n_slices;
    for(uword i=0;i<nmds;++i){
        result += a.slice(i)*b(i);
    }
    return result;
}
cube cubprod(cube a, vec b){
    cube result = zeros<cube>(size(a));
    uword nmds = a.n_slices;

    for(uword i=0;i<nmds;++i){
        result.slice(i) = a.slice(i)*b(i);
    }
    return result;
}
cube cubprod2(cube a, vec b){
    cube result = zeros<cube>(size(a));
    uword nmds = a.n_slices;

    for(uword i=0;i<nmds;++i){
        result.slice(i) = a.slice(i)*a.slice(i)*b(i);
    }
    return result;
}

vec fit1d(const cx_mat& hm,const cube& dip,const cube& pol,vec d1,vec d2,vec q1,vec q2,vec lamd,vec gamd) {

    cx_cube qm =deom_c1*(cubprod(dip,q1)+cubprod(pol,q2));
//    cx_cube qm =deom_c1*(cubprod(dip,q1)+cubprod2(dip,q2));
    syst s(hm,qm);
    bath b(npsd,pade,temp,lamd,gamd);
    cx_vec index = b.expn_gam;
    hidx h(lmax,nmax,ferr,index);
	    

    mat  sdip = matprod(dip,d1)+matprod(pol,d2);;
    cube pdip = zeros<cube>(size(qm));
    vec  bdip = zeros<vec>(size(index));
    vec  fw_im(num_ob);
//    resp1st (w1_max, nt1, thedt, staticErr, nk, sdip, pdip, bdip, 's', s, b, h); 
    resp1 (fw_im,w1_max, nt1, thedt, staticErr, nk, sdip, pdip, bdip,  s, b, h);
    return fw_im;
}

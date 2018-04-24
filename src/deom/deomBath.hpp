/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#ifndef DEOMBATH_H_
#define DEOMBATH_H_

#include <map>
#include <string>
#include "armadillo"
#include "json11.hpp"
#include "deomConst.hpp"

using namespace std;
using namespace arma;
using namespace json11;

double jwdru (double omg,vec jdru);
double jwsdr (double omg,vec jdru);
double fBose(double x,vec& pole,vec& resi,double rn,double tn);
void  psd(vec& pole,vec& resi,double& rn,double& tn,int N,int BoseFermi,int pade);

class bath {

    public:

    double temperature;
    ivec   modLabel;
    cx_vec coef_lft;
    cx_vec coef_rht;
    vec    coef_abs;
    cx_vec expn_gam;
    vec    delt_res;

    bath (const Json& json) {

    uword npsd = json["npsd"].int_value();   
    uword pade = json["pade"].number_value();
    double temp = json["temp"].number_value();
    Json jomg = json["jomg"];
    uword nmod = json["nmod"].int_value();   

    uword nind =0;
//    nind =  0;
    uvec mode ;
    for (uword m=0;m< nmod; m++){
        uword ndru = 0;
        uword nsdr = 0;
        try{
	    for (uword nmode=0;!jomg[m]["jdru"][nmode].is_null();nmode++){
                ndru++;
	    }
        }catch (const std::exception& e){
            ndru = 0;
/*        }try{
            nsdr = jomg[m]["jsdr"].n_elem;
        }catch (const std::exception& e){
            nsdr = 0;
*/        }
        uword nper = ndru+2*nsdr+npsd;
        uvec minsert=m*ones<uvec>(nper);
        mode.insert_rows(nind,minsert);
        nind += nper;
    }
    double rn,tn;
    vec pole,resi;
    psd(pole, resi, rn, tn,npsd,1,pade);
    }

    bath (const bath& rhs): 
          temperature(rhs.temperature),
          modLabel(rhs.modLabel),
          coef_lft(rhs.coef_lft), 
          coef_rht(rhs.coef_rht), 
          coef_abs(rhs.coef_abs), 
          expn_gam(rhs.expn_gam), 
          delt_res(rhs.delt_res) {}

    bath (const double temp, const ivec& mlbl, const cx_vec& etal, const cx_vec& etar, const vec& etaa, 
          const cx_vec& expn, const vec& delr): 
          temperature(temp),
          modLabel(mlbl),
          coef_lft(etal), 
          coef_rht(etar), 
          coef_abs(etaa),
          expn_gam(expn), 
          delt_res(delr) {}

    bath& operator= (const bath& rhs) {
        if (this != &rhs) {
            temperature = rhs.temperature;
            modLabel = rhs.modLabel;
            coef_lft = rhs.coef_lft;
            coef_rht = rhs.coef_rht;
            coef_abs = rhs.coef_abs;
            expn_gam = rhs.expn_gam;
            delt_res = rhs.delt_res;
        }
        return *this;
    }

   ~bath () {};

};

#endif

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

cx_double jwdru (cx_double omg,cx_double lamda,cx_double gamma){
    return 2.0*lamda*gamma*omg/(omg*omg+gamma*gamma);
}
cx_double jwsdr (cx_double omg,vec jdru);
cx_double fBose(cx_double x,vec& pole,vec& resi,double rn,double tn);
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

        const string modeFile = json["modeFile"].string_value();
        const string etalFile = json["etalFile"].string_value();
        const string etarFile = json["etarFile"].string_value();
        const string etaaFile = json["etaaFile"].string_value();
        const string expnFile = json["expnFile"].string_value();
        const string delrFile = json["delrFile"].string_value();

        temperature = json["temp"].number_value();

        printf ("$InitBath\n");
        if (modLabel.load (modeFile, arma_ascii)) {
            modLabel.print("modLabel");
        } else {
            printf ("modLabel is not loaded!\n");
        }
        if (coef_lft.load (etalFile, arma_ascii)) {
            coef_lft.print("coef_lft");
        } else {
            printf ("coef_lft is not loaded!\n");
        }
        if (coef_rht.load (etarFile, arma_ascii)) {
            coef_rht.print("coef_rht");
        } else {
            printf ("coef_rht is not loaded!\n");
        }
        if (coef_abs.load (etaaFile, arma_ascii)) {
            coef_abs.print("coef_abs");
        } else {
            printf ("coef_abs is not loaded!\n");
        }
        if (expn_gam.load (expnFile, arma_ascii)) {
            expn_gam.print("expn_gam");
        } else {
            printf ("expn_gam is not loaded!\n");
        }
        if (delt_res.load (delrFile, arma_ascii)) {
            delt_res.print("delt_res");
        } else {
            printf ("delt_res is not loaded!\n");
        }
        printf ("$InitBath\n\n");
    }

    bath (const int npsd, const int pade, const double temp,const vec& lamd,const vec& gamd) {

	    int nind =0;
	    int nmod = gamd.n_rows;
	    for (uword m=0;m< nmod; m++){
		int ndru = 1;//special drue case
		int nsdr = 0;
		int nper = ndru+2*nsdr+npsd;
		ivec minsert=ones<ivec>(nper)*m;
		modLabel.insert_rows(nind,minsert);
		nind += nper;
	    }
	    double rn,tn;
	    vec pole,resi;
	    psd(pole, resi, rn, tn,npsd,1,pade);
            temperature = temp;
	    uword iind=0;
	    coef_lft = zeros<cx_vec>(nind);
	    coef_rht = zeros<cx_vec>(nind);
	    coef_abs = zeros<vec>(nind);
	    expn_gam = zeros<cx_vec>(nind);
	    delt_res = zeros<vec>(nmod);
	    vec delr = zeros<vec>(nmod);
	    for (uword m=0;m< nmod; m++){
	       for(uword idru=0;idru<1;idru++){
		expn_gam(iind) = gamd(m);
            	cx_double eta = -2*deom_ci*lamd(m)*gamd(m)*fBose(-1*deom_ci*gamd(m)/temp,pole,resi,rn,tn);
            	coef_lft(iind) = eta;
                coef_rht(iind) = conj(eta);
                coef_abs(iind) = abs(eta);
                delr(m) += 2*lamd(m)*gamd(m)/temp*rn;
		iind++;
               }
	       for(uword ipsd=0;ipsd<1;ipsd++){
		expn_gam(iind) = pole(ipsd)*temp;
		cx_double zomg = -1*deom_ci*expn_gam(iind);
		cx_double jsmd = 0;
	       	for(uword idru=0;idru<1;idru++){
			jsmd += jwdru(zomg,lamd(m),gamd(m));
                }
            	cx_double eta = -2*deom_ci*resi(ipsd)*temp*jsmd;
            	coef_lft(iind) = eta;
                coef_rht(iind) = conj(eta);
                coef_abs(iind) = abs(eta);
		iind++;
               }
            }
            delt_res = delr;
    }

    bath(const bath& rhs): 
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

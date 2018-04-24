/**
 * writtened by Zhijun Pan
 * 2018.3.26
 * -------------------
 * FGH-DVR for eigenstate 
 * of Morse potential
 */
#include <ctime>
#include <chrono>
#include <iomanip>

#include <iostream>
#include <sstream>
#include <string>
#include "fgh.hpp"
#include "deom.hpp"

using namespace std;
using namespace arma;
using namespace json11;
int main() {

    ifstream jsonFile("morse.json");
    stringstream strStream;
    strStream << jsonFile.rdbuf();
    string jsonStr = strStream.str();
    string err;

    const Json json = Json::parse(jsonStr,err);
    if (!err.empty()) {
        printf ("Error in parsing input file: %s\n", err.c_str());
        return 0;
    }

    const uword    ne = json["ne"].int_value();
    const uword    ie = json["ie"].int_value();
    const Json  mds = json["modes"];
    uword  nmds = 0;
    uword  ndim = 1;
    uvec nr (mmax);
    fvec ri (mmax);
    fvec rf (mmax);
    fvec Eie(mmax); 
    field<mat> hams(mmax);   // field for different nr in modes
    for (uword nmode=0;!mds[nmode].is_null();nmode++){
      nr[nmode] = mds[nmode]["nr"].int_value();
      ri[nmode] = mds[nmode]["ri"].number_value();
      rf[nmode] = mds[nmode]["rf"].number_value();
      ndim *= nr[nmode];
      hams(nmode) = get_ham_1d(mds[nmode]);
      Eie(nmode)  = get_energy(mds[nmode],ie);
      nmds++;
    }
    cout << "there are " << nmds << " modes" << '\n';
    mat vmn = zeros<mat>(nmds,nmds);// the couple coefficient, not consider now
    printime("Begin calculation wavefunction and energy:");
    mat wavefun;
    vec energy ;
    wave_energy( wavefun, energy, nmds, ne,Eie,nr,ri, rf,vmn,hams);
    printime("End loading wavefunction and energy:");

    fvec q1(nmds);
    fvec q2(nmds);
    vec lamd(nmds);
    vec gamd(nmds);
    for (uword nmode=0;nmode<nmds;nmode++){
      q1[nmode] = mds[nmode]["q1"].number_value();
      q2[nmode] = mds[nmode]["q2"].number_value();
      lamd[nmode] = mds[nmode][0]["lamd"].number_value();
      gamd[nmode] = mds[nmode][0]["gamd"].number_value();
    }
    ifstream jsonFile2("default.json");
    stringstream strStream2;
    strStream2 << jsonFile2.rdbuf();
    string jsonStr2 = strStream2.str();

    const Json deft = Json::parse(jsonStr2,err);
    if (!err.empty()) {
        printf ("Error in parsing input file: %s\n", err.c_str());
        return 0;
    }
//    const Json  mysyst = deft["syst"];
    const Json  mybath = deft["bath"];
    const Json  myhidx = deft["hidx"];
    int npsd = mybath["npsd"].int_value();
    int pade = mybath["pade"].int_value();
    double temp = mybath["npsd"].number_value();
    int lmax = myhidx["lmax"].int_value();
    int nmax = myhidx["nmax"].int_value();
    double ferr = myhidx["ferr"].number_value();
    cx_mat hm = deom_c1*diagmat(energy);
    cx_cube qm =deom_c1*qmod(q1,q2,nr,ri,rf,wavefun);
    syst s(hm,qm);

    bath b(npsd,pade,temp,lamd,gamd);
//        printime ("End parsing input file: ");
	    
    hidx h(lmax,nmax,ferr,b.expn_gam);
    return 0;
}

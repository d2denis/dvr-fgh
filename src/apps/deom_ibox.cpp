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

    fvec lamd,gamd,lams,gams;
    fvec q1(nmds);
    fvec q2(nmds);
    for (uword nmode=0;nmode<nmds;nmode++){
      q1[nmode] = mds[nmode]["q1"].number_value();
      q2[nmode] = mds[nmode]["q2"].number_value();
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
    const Json  mysyst = deft["syst"];
    const Json  mybath = deft["bath"];
    const Json  myhidx = deft["hidx"];
    Json::object inidic = mybath;
    inidic["jomg"]=Json::array{80,80};
    inidic["nmod"]=nmds;
    Json nimei = Json{inidic};
    bath(nimei);
	    cout << bath.modLabel.n_elem <<'\n';
    return 0;
}

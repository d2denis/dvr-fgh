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

    const int    ne = json["ne"].int_value();
    const int    ie = json["ie"].int_value();
    const Json  mds = json["modes"];
    int  nmds = 0;
    int  ndim = 1;
    uvec nr = ones<uvec>(mmax);
    fvec ri = ones<fvec>(mmax);
    fvec rf = ones<fvec>(mmax);
    field<mat> hams(mmax);   // field for different nr in modes
    for (int nmode=0;!mds[nmode].is_null();nmode++){
      nr[nmode] = mds[nmode]["nr"].int_value();
      ri[nmode] = mds[nmode]["ri"].number_value();
      rf[nmode] = mds[nmode]["rf"].number_value();
      ndim *= nr[nmode];
      hams(nmode) = get_ham_1d(mds[nmode]);
      nmds++;
    }
    cout << "there are " << nmds << " modes" << '\n';
    mat vmn = zeros<mat>(nmds,nmds);// the couple coefficient, not consider now
    /*
     * Fourier Grid Hamiltonian
     */
    printime("Begin Hamiltonian construction: ");
    sp_mat hamt = fgh_dvr(ndim,nmds,nr,ri,rf,vmn,hams);
    printime("Hamiltonian constructed: ");
    int sparse = hamt.n_nonzero;
    printf("Total nonzero: %12d\t Sparce: %6.4f\n",sparse,1.0*sparse/ndim/ndim);
    mat wavefun;
    vec energy ;
    eigs_sym(energy,wavefun,hamt,ne,"sm",.000000000000000000000001);
    printime("Diagonalized :") ;
    energy.save("energy.dat",raw_ascii);
    wavefun.save("wavefun.dat",raw_ascii);
    return 0;
}
    

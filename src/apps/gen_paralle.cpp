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
    time_t secb=time(0);
    tm *ti = localtime(&secb);
    printf("Begin Hamiltonian construction: %02d:%02d:%02d\n",ti->tm_hour,ti->tm_min,ti->tm_sec);
    sp_mat hamt(ndim,ndim);
    hamt = fgh_dvr(ndim,nmds,nr,ri,rf,vmn,hams);
    time_t sece=time(0);
    tm *tf = localtime(&sece);
    printf("Hamiltonian constructed: %02d:%02d:%02d\n",tf->tm_hour,tf->tm_min,tf->tm_sec);
    int sparse = nonzeros(hamt).n_rows;
    printf("Total nonzero: %12d\t Sparce: %6.4f\n",sparse,1.0*sparse/ndim/ndim);
    mat wavefun;
    vec energy ;
    eigs_sym(energy,wavefun,hamt,ne,"sm",.000000000000000000000001);
    time_t secu=time(0);
    tm *tu = localtime(&secu);
    printf("Diagonalized : %02d:%02d:%02d\n",tu->tm_hour,tu->tm_min,tu->tm_sec);
    energy.save("energy.dat",raw_ascii);
    wavefun.save("wavefun.dat",raw_ascii);
    /*
     * Read states and generate input 
     */
//    mat Psi(ie,ndim);
// sort the state according to <X^2>
    return 0;
}
    

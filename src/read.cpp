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

    const uword    ne = json["ne"].int_value();
    const uword    ie = json["ie"].int_value();
    const Json  mds = json["modes"];
    uword  nmds = 0;
    uword  ndim = 1;
    uvec nr = ones<uvec>(mmax);
    fvec ri = ones<fvec>(mmax);
    fvec rf = ones<fvec>(mmax);
    field<mat> hams(mmax);   // field for different nr in modes
    for (uword nmode=0;!mds[nmode].is_null();nmode++){
      nr[nmode] = mds[nmode]["nr"].int_value();
      ri[nmode] = mds[nmode]["ri"].number_value();
      rf[nmode] = mds[nmode]["rf"].number_value();
      ndim *= nr[nmode];
      nmds++;
    }
    cout << "there are " << nmds << " modes" << '\n';
/*    time_t sece=time(0);
    tm *tf = localtime(&sece);
    printf("Hamiltonian constructed: %02d:%02d:%02d\n",tf->tm_hour,tf->tm_min,tf->tm_sec);
    int sparse = nonzeros(hamt).n_rows;
    printf("Total nonzero: %12d\t Sparce: %6.4f\n",sparse,1.0*sparse/ndim/ndim);
*/
    mat wavefun;
    vec energy ;
    time_t secu=time(0);
    tm *tu = localtime(&secu);
    printf("Loading energy and wavefunction : %02d:%02d:%02d\n",tu->tm_hour,tu->tm_min,tu->tm_sec);
    energy.load("energy.dat",raw_ascii);
    wavefun.load("wavefun.dat",raw_ascii);
    for (uword m=0; m<ie; m++)
	    cout << energy(m) <<'\n';

    /*
     * Truncate states and generate input 
     */
    mat Psi(ne,ndim);
// Label the states according to single mode H_s
    ivec label = ones<ivec>(ne);
    for(uword i=0; i<ne; i++){
	for(uword m=0; m< nmds; m++){
	    for(uword j=0; j< ndim; j++){
		hams(m)
    return 0;
}
    

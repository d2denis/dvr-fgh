/**
 * writtened by Zhijun Pan
 * 2018.3.26
 * -------------------
 * FGH-DVR for eigenstate 
 * of Morse potential
 */
#include "fgh.hpp"

using namespace std;
using namespace arma;

double get_pot_1d (const double m0,const double de,const double beta,const double r0,const double r) {
    double morse = (1-exp(-beta*(r-r0)));
    return de*morse*morse;//+0.5*m0*lamd*(q1*r2*(1+0.5*q2*r2));
}
double get_pot_2d (const double gxy,const double  x,const double  y) {
    return gxy*x*y*(x+y);
}
mat get_ham_1d(Json m){
      double  m0 = m["m0"].number_value()*am2unit;
      double  ri = m["ri"].number_value()*ar2unit;
      double  rf = m["rf"].number_value()*ar2unit;
      double  de = m["de"].number_value()*cm2unit;
      double  be = m["be"].number_value()/ar2unit;
      double  re = m["re"].number_value()*ar2unit;
//      double  q1 = m["q1"].number_value();
//      double  q2 = m["q2"].number_value();
      int     nr = m["nr"].int_value();
      int     nk = (nr-1)/2;
      double  dr = (rf-ri)/nr;
      double  dk = 2.0*deom_pi/(rf-ri);
      double ttmp= dk*dk/(2.0*m0);
      mat ham    = zeros<mat>(nr,nr);
      for(int i=0;i<nr;++i){
	  for(int j=0;j<nr;++j){
              for(int k=0;k<nk;++k){
		  double tk = 2.0*deom_pi*(k+1)*(i-j)/nr;
		  ham(i,j) += cos(tk)*ttmp*(k+1)*(k+1)*2.0/nr;
	      }
	  }
//	  ham(i,i) += get_pot_1d(m0,m["jdru"],q1,q2,de,be,re,ri+i*dr);
	  ham(i,i) += get_pot_1d(m0,de,be,re,ri+i*dr);
      }
      return ham;
}
//sp_mat fgh_dvr(int ndim,int nmds,uvec nr,fvec ri,fvec rf,mat vmn,field<mat> hams){
sp_mat fgh_dvr(const int ndim,const int nmds,const uvec& nr,const fvec& ri,const fvec& rf,const mat& vmn,const field<mat>& hams){
    uvec cum_prod = zeros<uvec>(nmds);
    cum_prod(0) = nr[nmds-1];
    for (int m=1;m<nmds;m++)
	    cum_prod(m) = cum_prod(m-1)*nr(nmds-m-1);
    fvec cum_divd = ones<fvec>(nmds)/cum_prod;
// batch insertion to accelerate sparce matrix construct
//    sp_mat hamt(ndim,ndim);
    uword nelem = 0;
    uword nnz = cum_prod(nmds-1)*((nr(0)-1)*nmds+1);

    umat locations(2,nnz);
    vec values(nnz);
    for(int i=0; i<ndim; i++){
	uvec i_list = decode(i,cum_divd,cum_prod,nmds);
        for(int j=i; j<ndim; j++){
    	    uvec j_list = decode(j,cum_divd,cum_prod,nmds);
	    if(i != j){
		int num=0;
		int loc=-1;
		for( int m=0;m<nmds;m++){
		        if (i_list(m) != j_list(m)){
		    	    loc = m;
		    	    num +=1;
		    	    if (num>1){
		    		    num = 2;
		    		    loc =-1;
		    		    break;
		    	    }
		        }
		}
	        if(num == 1 && loc > -1){
		    mat ham0 =	hams(loc);
		locations(0,nelem)=locations(1,nelem+1)=i;
		locations(1,nelem)=locations(0,nelem+1)=j;
		    double tmp = ham0(i_list(loc),j_list(loc));
		values(nelem)=values(nelem+1)= tmp; 
		nelem += 2;
		}
	    }else{
		    double tmp = 0.0;
		    for( int m = 0; m < nmds; m++){
                        int rm = ri(m)+i_list(m)*(rf(m)-ri(m))/nr(m);
			mat ham0 = hams(m);
                        tmp += ham0(i_list(m),i_list(m));
		        for( int n = m; n < nmds; n++){
                            int rn = ri(n)+i_list(n)*(rf(n)-ri(n))/nr(n);
                            tmp += get_pot_2d(vmn(m,n),rm,rn);
			}
		    }
//		    hamt(i,j) = tmp;
		locations(0,nelem)=locations(1,nelem)=i;
		values(nelem)= tmp; 
		nelem += 1;
	    }
	}
    }
    sp_mat hamt(locations,values);
    return hamt;
}

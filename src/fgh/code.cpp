/**
 * writtened by Zhijun Pan
 * 2018.3.26
 * -------------------
 * FGH-DVR for eigenstate 
 * of Morse potential
 */
#include "fgh.hpp"
int encode (const uvec i_list,const uvec cum_prod,const int nlen) {
    int i_rank = 0;
    for (int i=0; i<(nlen-1); ++i)
        i_rank += i_list(i)*cum_prod(i);
    i_rank += i_list(nlen-1);
    return i_rank;
}
//Decode should be accelerate.
uvec decode (int i_rank,const fvec cum_divd,const uvec cum_prod,const int nlen) {
    uvec i_list(nlen);
    for (int i=0; i<(nlen-1); ++i){
        i_list(i)=floor(i_rank*cum_divd(nlen-i-2));
        i_rank = i_rank-i_list(i)*cum_prod(nlen-i-2);
    }
    i_list(nlen-1) = i_rank;
    return i_list;
}

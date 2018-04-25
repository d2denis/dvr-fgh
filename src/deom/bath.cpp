/**
 * git clone https://github.com/d2denis/dvr-fgh.git
 * ---
 * Written by Zhijun Pan
 */
#include <string>
#include "armadillo"
#include "deomConst.hpp"

using namespace std;
using namespace arma;

cx_double jwdru (cx_double omg,cx_double lamda,cx_double gamma){
    return 2.0*lamda*gamma*omg/(omg*omg+gamma*gamma);
}
cx_double jwsdr (cx_double omg,vec jdru){
    cx_double lamda = jdru[0];
    cx_double omega = jdru[1];
    cx_double gamma = jdru[2];
    cx_double omegamma = omg*gamma;
    cx_double omega2 = omega*omega;
    return 2.0*lamda*omegamma*omega2/((omg*omg-omega2)*(omg*omg-omega2)+omegamma*omegamma);
}

cx_double fBose(cx_double x,vec& pole,vec& resi,double rn,double tn){
        cx_double sumpole = 0.0;
        for (int i=0;i<pole.n_elem;i++)
            sumpole += 2.0*resi[i]*x/(x*x+pole[i]*pole[i]); 
	return 1.0/x+0.5+rn*x+tn*x*x*x+sumpole;
}
vec mytrunc(const vec& input, int N){
    uword size = input.n_elem;
    if(N>size){
        cout << "can't trunctate\n";
        return zeros<vec>(N);
    }
    vec output(N);
    for(int i=0;i<N;i++){
        output(i) = input(i);
    }
    return output;
}

vec tseig(const vec& D,const vec& E){
    mat H = diagmat(E,-1) + diagmat(D,0) + diagmat(E,1);
    vec eigval = eig_sym(H);
    return sort(eigval,"descend");
}
void MSD(vec& pole,vec& resi,double& rn,double& tn,int N,int BoseFermi){
    rn = tn = 0.0;
    resi = ones<vec>(N);
    if(BoseFermi == 1){
        pole = linspace(2,2,N)*deom_pi;
    }else if (BoseFermi == 2){
        pole = linspace(1,2,N)*deom_pi;
    }
}
vec sqrtoff(int N,double temp){
    vec fu(N);
    for(int i=0;i<N;i++)
        fu(i) = 1.0/sqrt((temp+2.0*i)*(temp+2.0*(i+1)));
    return fu;
}

void  psd(vec& pole,vec& resi,double& rn,double& tn,int N,int BoseFermi,int pade){
//Pan disable the pade=3 option
    if (N < 0 || BoseFermi<1 || BoseFermi>2 || pade<0 || pade>2){
        printf("N or BoseFermi or pade has wrong value!");
    }

    if (pade == 0){
        MSD(pole,resi,rn,tn,N,BoseFermi);
    }else if (pade == 1 || pade == 2){
        if(N > 0){
            int M = 2*N+pade/2;
            double temp = (BoseFermi==1)? 3.0: 1.0;
            vec diag = zeros<vec>(M);
            vec doff = sqrtoff((M-1),temp);
	    vec tseigen = tseig(diag,doff);
            pole = 2.0/mytrunc(tseigen,N);
            vec pol2 = square(pole);
            M -= 1;
            double tmp = (BoseFermi==1)? 5.0: 3.0;
            diag = zeros<vec>(M);
            doff = sqrtoff((M-1),tmp);
            M /= 2;
            vec eig2 = square(2.0/mytrunc(tseig(diag,doff),M));
            double   scaling = (pade == 1)? N*(2.0*N+temp) : 1.0/(4.0*(N+1.0)*(2.0*N+temp));
	    resi = zeros<vec>(N);
            for(int j=0;j<N;j++){
                if (pade == 2){
                    temp = 0.5*scaling*(eig2[j]-pol2[j]);
                }else if (pade == 1){
                    if (j == N-1){
                        temp = 0.5*scaling;
                    }else{
                        temp = 0.5*scaling*(eig2[j]-pol2[j])/(pol2[N-1]-pol2[j]);
		    }
		}
                for(int k=0;k<M;k++){
                    temp *= (k!=j) ? (eig2[k]-pol2[j])/(pol2[k]-pol2[j]): 1.0;
                }
                resi[j] = temp;
            }
        }
        rn = tn = 0.0;
        if(BoseFermi==1 && pade==2){
            rn = 1.0/(4.0*(N+1.0)*(2.0*N+3.0));
        }
    }
}
/*
 int main(){
  cout << "test psd\n" ;
    double rn=1.0;
    double tn;
    vec pole,resi;
    psd(pole,resi,rn,tn,4,1,2);
    printf( "psd of 4,1,2:%f\t%f\t%f\t%f\n", pole(3), resi(3), rn, tn);
}
*/

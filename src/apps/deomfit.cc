/*
 * 1D fit for one dimension  morse potential
 * by Zhijun Pan
 * 2018.5.3
 */
#include "ceres/ceres.h"
#include "glog/logging.h"

#include <ctime>
#include <chrono>
#include <iomanip>

#include <iostream>
#include <sstream>
#include <string>
#include "fgh.hpp"
#include "deom.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>


using ceres::AutoDiffCostFunction;
using ceres::NumericDiffCostFunction;
using ceres::CENTRAL;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using namespace std;
using namespace arma;
using namespace json11;

    const uword    nmds =1;
    const uword    ne = 30;//
    const uword    ie = 2;//
    const int npsd = 2;//
    const int pade = 2;//
    const double temp = 9.5e-4;//
    const int lmax = 5;//
    const int nmax = 50;//
    const double ferr = 2e-7;//my
    const double w1_max = 0.01;//
    const int    nt1 = 400;//
    const int num_ob =nt1/2;
    const double dt = 10.0;//
    const double staticErr = 2e-7;//
    const int    nk = 32;//
    const string sch_hei = "s";//
  const int  num_para = 4;

mat matprod(cube a, vec b){
    mat result = zeros<mat>(size(a.slice(0)));
    uword nmds = b.n_elem;
    for(uword i=0;i<nmds;++i){
        result += a.slice(i)*b(i);
    }
    return result;
}
cube cubprod(cube a, vec b){
    cube result = zeros<cube>(size(a));
    uword nmds = b.n_elem;

    for(uword i=0;i<nmds;++i){
        result.slice(i) = a.slice(i)*b(i);
    }
    return result;
}
cube cubprod2(cube a, vec b){
    cube result = zeros<cube>(size(a));
    uword nmds = b.n_elem;

    for(uword i=0;i<nmds;++i){
        result.slice(i) = a.slice(i)*a.slice(i)*b(i);
    }
    return result;
}

int json2para(const uword nmds,double* d1,double* d2,double* p1,double* p2,double* q1,double* q2,double* gamd){
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

    Json mds = json["modes"];

    for (uword nmode=0;nmode<nmds;nmode++){
      d1[nmode] = mds[nmode]["d1"].number_value();
      d2[nmode] = mds[nmode]["d2"].number_value();
      p1[nmode] = mds[nmode]["p1"].number_value();
      p2[nmode] = mds[nmode]["p2"].number_value();
      q1[nmode] = mds[nmode]["q1"].number_value();
      q2[nmode] = mds[nmode]["q2"].number_value();
//      lamd[nmode] = mds[nmode]["jdru"][0]["lamd"].number_value()*cm2unit;
      gamd[nmode] = mds[nmode]["jdru"][0]["gamd"].number_value()*cm2unit;
    }
    return 1;
}
int setjson(Json& mds,const uword nmds,vec& d1,vec& d2,vec& q1,vec& q2,vec& lamd,vec& gamd){
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

    mds = json["modes"];

    for (uword nmode=0;nmode<nmds;nmode++){
      d1[nmode] = mds[nmode]["d1"].number_value();
      d2[nmode] = mds[nmode]["d2"].number_value();
      q1[nmode] = mds[nmode]["q1"].number_value();
      q2[nmode] = mds[nmode]["q2"].number_value();
      lamd[nmode] = mds[nmode]["jdru"][0]["lamd"].number_value()*cm2unit;
      gamd[nmode] = mds[nmode]["jdru"][0]["gamd"].number_value()*cm2unit;
    }
    return 1;
}

void readinp(const Json mds,cx_mat& hm,cube& dip,cube& pol,const uword nmds=1){
    uword  ndim = 1;
    uvec nr (nmds);
    fvec ri (nmds);
    fvec rf (nmds);
    fvec Eie(nmds); 
    field<mat> hams(mmax);   // field for different nr in modes
    for (uword nmode=0;!mds[nmode].is_null();nmode++){
      nr[nmode] = mds[nmode]["nr"].int_value();
      ri[nmode] = mds[nmode]["ri"].number_value()*ar2unit;
      rf[nmode] = mds[nmode]["rf"].number_value()*ar2unit;
      ndim *= nr[nmode];
      hams(nmode) = get_ham_1d(mds[nmode]); //immutable Json, ar2unit don't affect ri,rf
      Eie(nmode)  = get_energy(mds[nmode],ie);
    }
    mat vmn = zeros<mat>(nmds,nmds);// the couple coefficient, not consider now
    mat wavefun;
    vec energy ;
    wave_energy( wavefun, energy, nmds, ne,Eie,nr,ri, rf,vmn,hams);
    hm = deom_c1*diagmat(energy);
    qmod(dip,pol,nr,ri,rf,wavefun,nmds);
}
void resp1 (vec& fw_im, const double w_max, const int mynt, const double dt,
              const double staticErr, const int nk,
              const mat& sdip, const cube& pdip, const vec& bdip,
              const syst& s, const bath& b, const hidx& h) {

    const double dw1 = w_max/mynt;
    const double dt1 = 2.0*deom_pi/w_max;
    const int    mt  = floor(dt1/dt);
    const double dt1_res = dt1-dt*mt;

    deom d1(s,b,h);

    cx_cube rho_t0 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);

    const mat& exph= expmat(-real(d1.ham1)/d1.temperature);
    rho_t0.slice(0).set_real(exph/trace(exph));
    d1.equilibrium (rho_t0,dt,staticErr,nk);

    cx_cube rho_t1 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nmax);
    d1.oprAct(rho_t1,sdip,pdip,bdip,rho_t0,'c');

    cx_vec ft = zeros<cx_vec>(mynt);

        for (int it=0; it<mynt; ++it) {

            ft(it) = deom_ci*d1.Trace(sdip,pdip,bdip,rho_t1);

//            printf ("In sch-pic: it=%d, nddo=%d, lddo=%d\n", it, d1.nddo, d1.lddo);
            double t1 = it*dt1;
            for (int jt=0; jt<mt; ++jt) {
                d1.rk4 (rho_t1,t1,dt);
                t1 += dt;
            }
            d1.rk4 (rho_t1,t1,dt1_res);
        }
    ft(0) *= 0.0;//Infinity(nan) occur nddo>17, should be fixed
    const cx_vec& fw = ifft(ft)*mynt*dt1;
    fw_im = imag(fw.rows(0,mynt/2-1));
}


vec fit1d(const cx_mat& hm,const cube& dip,const cube& pol,vec d1,vec d2,vec q1,vec q2,vec lamd,vec gamd) {

    cx_cube qm =deom_c1*(cubprod(dip,q1)+cubprod2(dip,q2));
    syst s(hm,qm);
    bath b(npsd,pade,temp,lamd,gamd);
    cx_vec index = b.expn_gam;
    hidx h(lmax,nmax,ferr,index);
	    

    mat  sdip = matprod(dip,d1)+matprod(pol,d2);;
    cube pdip = zeros<cube>(size(qm));
    vec  bdip = zeros<vec>(size(index));
    vec  fw_im(num_ob);
    resp1 (fw_im,w1_max, nt1, dt, staticErr, nk, sdip, pdip, bdip,  s, b, h);
    return fw_im;
}

// Read a Brownian_oscillator Absorption Large dataset.
class BALProblem {
 public:
  ~BALProblem() {
    delete[] observations_;
    delete[] parameters_;
  }

  int num_parameters()       const { return num_parameters_;               }
  double* observations() const { return observations_;                   }
  double* mutable_points()           { return parameters_;                     }


  bool LoadFile(const char* filename) {
    FILE* fptr = fopen(filename, "r");
    if (fptr == NULL) {
      return false;
    };


    observations_ = new double[num_ob];

/*    FscanfOrDie(fptr, "%d", &num_parameters_);
    parameters_ = new double[num_parameters_];
    for (int i = 0; i < num_parameters_; ++i) {
      FscanfOrDie(fptr, "%lf", parameters_ + i);
    }
*/
    num_parameters_ = nmds*num_para;
    parameters_ = new double[num_parameters_];
    double nimei[3];
    json2para(nmds,parameters_,nimei,nimei,nimei,parameters_+nmds,parameters_+2*nmds,parameters_+3*nmds);
    for (int i = 0; i < num_ob; ++i) {
        FscanfOrDie(fptr, "%lf", observations_ + i );
    }

    return true;
  }

 private:
  template<typename T>
  void FscanfOrDie(FILE *fptr, const char *format, T *value) {
    int num_scanned = fscanf(fptr, format, value);
    if (num_scanned != 1) {
      LOG(FATAL) << "Invalid DEOM data file.";
    }
  }

  int num_parameters_;

  double* observations_;
  double* parameters_;
};
//array to vector
vec a2v(const double* x){
    vec v(nmds);
    for(uword i=0;i<nmds;i++){
        v(i) = *(x+i);
    }
    return v;
}
struct Resp1Residual {
  Resp1Residual( double* y)
      :  y_(y) {}

//  template <typename T> bool operator()(const T* const m,
   bool operator()(const double* const d,const double* const q,const double* const qq,const double* const g,  double* residual) const {
  
  Json mds;
  vec  d1(nmds);
  vec  q1(nmds);
  vec  d2(nmds);
  vec  q2(nmds);
  vec  lamd(nmds);
  vec  gamd(nmds);
  if(setjson(mds,nmds,d1,d2,q1,q2,lamd,gamd)){
//parameter block
     d1=a2v(d);
     q1=a2v(q);
     q2=a2v(qq);
//     lamd(0)=m[1]*cm2unit; // fixed to avoid over fitted, the factor is q*lamd.
     gamd=a2v(g);
//system block
     cx_mat hm;
     cube dip, pol;
     readinp(mds,hm,dip,pol,nmds);

     vec fw_im(num_ob);
     vec nimei = fit1d(hm,dip,pol,d1,d2,q1,q2, lamd, gamd);
  nimei.save("fitted.dat",raw_ascii);
     for(int i=0;i<num_ob;i++){
         residual[i] = y_[i]-nimei(i);
     }
     
  }
    return true;
  }

 private:
  double* y_;
};

int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);
/*
  if (argc != 2) {
    std::cerr << "usage: mycf <bal_problem>\n";
    return 1;
  }
  */

  BALProblem bal_problem;
  if (!bal_problem.LoadFile("resp1st_im.w")) {
    std::cerr << "ERROR: unable to open file "<< argv[1]  << "\n";
    return 1;
  }
  double* datax = bal_problem.mutable_points();
  double* datay = bal_problem.observations();
  if(bal_problem.num_parameters()!=nmds*num_para){
    std::cerr << "ERROR: number of parameter not match"  << "\n";
    return 1;
  }
  std::cout <<bal_problem.num_parameters()<< " parameters, d: " << datax[0] << " q1: " << datax[1]<< " q2: " << datax[2]<< " gamd: " << datax[3]/cm2unit<< "\n";
  Problem problem;
    problem.AddResidualBlock(
//        new AutoDiffCostFunction<Resp1Residual, num_ob, 1>(
        new NumericDiffCostFunction<Resp1Residual,CENTRAL, num_ob, nmds,nmds,nmds,nmds>(
            new Resp1Residual(datay)),
        NULL,
        datax,datax+1,datax+2,datax+3);

  Solver::Options options;
  options.max_num_iterations = 500;
  options.linear_solver_type = ceres::DENSE_SCHUR;
//  options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
  options.minimizer_progress_to_stdout = true;

  Solver::Summary summary;
  Solve(options, &problem, &summary);
  std::cout << summary.BriefReport() << "\n";
  std::cout << "Final   d: " << datax[0] << " q1: " << datax[1]<< " q2: " << datax[2]<< " gamd: " << datax[3]/cm2unit<< "\n";
  return 0;
}

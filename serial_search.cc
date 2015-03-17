#include <stdio.h>
#include <math.h>
#include <cblas.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#define NCHANNELS 8

double sum_array(double* array,int n){
    int i;
    double result=0;
    for (i=0;i<n;i++){
        result += array[i];
    }
    return result;
}

double dot_prod(double* a, double* b, int n){
    int i;
    double result=0;
    for (i=0;i<n;i++){
        result += a[i] * b[i];
    }
    return result;
}

//recorded data
const double sr = 10; //samping rate in samples/ms
double *current_recorded;
double *time_recorded;

//Nernest potentials in mv
const double vNa=50; 
const double vK=-90; 
const double vL=-70; 
const double vH=-30;

// Maximal conductances
double gNa=450;
double gK=70;
double gL=2;
double gCa=19;
double gNap=1;
double gSK=8.5;
double gH=4;
double gA=0;
double gT=2.65;

//Capacitance
const double C=50;

//Time constants in ms
const double tauNbar=6;
const double tauHPbar=1000;
const double tauE=20;
const double tauH=1;
const double taurs=1500;
const double RTF=26.7;
const double Ca_ex=2.5;

//Function constants in mV
const double thetaM=-35;
const double thetaN=-30;
const double sigmaM=-5;
const double sigmaN=-5;
const double thetaS=-30;
const double sigmaS=-12;
const double thetaMP=-40;
const double sigmaMP=-6;
const double thetaHP=-48;
const double sigmaHP=6;
const double thetaA=-20;
const double sigmaA=-10;
const double thetaE=-60;
const double sigmaE=5;
const double thetaaT = -65;
const double sigmaaT = -7.8;
const double thetabT = 0.4;
const double sigmabT = -0.1;
const double thetarT = -67;
const double sigmarT = 2; 
const double taur0 = 40;
const double taur1 = 17.5;
const double thrT = 68;
const double sgmrt = 2.2;
const double phirT=0.2;
const double thetaRF=-105;
const double sigmaRF=5;
const double thetaRS=-105;
const double sigmaRS=25;

// Other constants (related to Ca dynamics)
const double kr=0.3;
const double prf=100;
double f=0.1;

const double eps=0.0015;
const double kca=0.3;
const double ks=0.7;

const double dt = .0005;
const double t_start = 0;
const double t_stop = 1000;

//const double iap = 70;
const double E[NCHANNELS] = {vNa, vK,vL,0,vK,vNa,vH,vK};
const double max_g[NCHANNELS] = {gNa, gK, gL, gCa, gSK, gNap, gH, gA};

// Initial conditions
double v = -70.14; 
double n = .0002;
double ca = 0.103;
double h = 0.01;
double rf = 0.03;
double rs = 0.03;
double e  = 0.05;
double rT = 0.01;
double hp = 0.1;

// other dynamic variables
 double minf;
double ninf;
double tauN;
double alphaH;
double betaH;
double hinf;
double s_inf;
double aTinf;
double bTinf;
double rTinf;
double taurT;
double kinf;
double mpinf;
double hpinf;
double tauhp;
double rinff;
double rinfs;
double taurf;
double ainf;
double einf;
double g[NCHANNELS];
double sum_g;
double v_inf;
double tau_v; 
double tau_ca;
double ca_inf;

double iNa;
double iK;
double iCa;
double iT;
double iSK;
double iNap;
double iH;
double iA;
double iL;

//interpolation 
gsl_interp *interpolation; 
gsl_interp_accel * accelerator;

int 
DiffEquations(double t, const double y[], double diff[], void *params)
{

    n=y[0];
    h=y[1];
    rf=y[2];
    rs=y[3];
    ca=y[4];
    v=y[5];
    e=y[6];
    rT=y[7];
    hp=y[8];

    // Na+ and K+ Equations and Currents
    minf = 1/(1+exp((v-thetaM)/sigmaM));
    ninf = 1/(1+exp((v-thetaN)/sigmaN));
    tauN = tauNbar/cosh((v-thetaN)/(2*sigmaN));

    alphaH = 0.128*exp(-(v+50)/18);
    betaH = 4/(1+exp(-(v+27)/5));
    hinf = alphaH/(alphaH+betaH);

    iNa = gNa*(pow(minf,3))*h*(v-vNa);
    iK = gK*(pow(n,4))*(v-vK);

    // L-Type Ca++ Equations and Current
    s_inf = 1/(1+exp((v-thetaS)/sigmaS));
    iCa = gCa*(pow(s_inf,2))*v*(Ca_ex/(1-exp((2*v)/RTF)));

    // T-Type Ca++ Equations and Current
    aTinf = 1/(1+exp((v-thetaaT)/sigmaaT));
    bTinf = 1/(1+exp((rT-thetabT)/sigmabT))-1/(1+exp(-thetabT/sigmabT));
    rTinf = 1/(1+exp((v-thetarT)/sigmarT));
    taurT = taur0+taur1/(1+exp((v-thrT)/sgmrt));
    iT=gT*pow(aTinf,3)*pow(bTinf,3)*v*(Ca_ex/(1-exp((2*v)/RTF)));

     // SK Equations and Current
    kinf = pow(ca,2)/(pow(ca,2)+pow(ks,2));
    iSK = gSK*kinf*(v-vK);

     // Na+ Persistant Current and Equations
    mpinf = 1/(1+exp((v-thetaMP)/sigmaMP));
    hpinf = 1/(1+exp((v-thetaHP)/sigmaHP));
    tauhp = tauHPbar/cosh((v-thetaHP)/(2*sigmaHP));
    iNap = gNap*mpinf*hp*(v-vNa);

    // Hyperpolarization activated inward current Ih
    rinff = 1/(1+exp((v-thetaRF)/sigmaRF));
    rinfs = 1/(1+exp((v-thetaRS)/sigmaRS));
    taurf = prf/(-7.4*(v+70)/(exp(-(v+70)/0.8)-1)+65*exp(-(v+56)/23));

    iH = gH*(kr*rf+(1-kr)*rs)*(v-vH);

    // A-type potassium current Ia
    ainf = 1/(1+exp((v-thetaA)/sigmaA));
    einf = 1/(1+exp((v-thetaE)/sigmaE));

    iA = gA*ainf*e*(v-vK);

    // Leak current
    iL = gL*(v-vL);

    // Applied current
    double iap;
    //int idx = t/.1;
    //    iap = current_recorded[idx];
    iap = gsl_interp_eval(interpolation, time_recorded, 
                          current_recorded, t, accelerator);
    // try {
    //     iap = gsl_interp_eval(interpolation, time_recorded, 
    //                           current_recorded, t, accelerator);
    // }
    // catch (int e) {
    //     // for (int idx; idx < sizeof(time_recorded)/sizeof(double); idx++){
    //     //     std::cout << time_recorded[idx] << "\n";
    //     // }
    //     std::cout << "time=" << t << "\n";
    //     // throw e;
    // }

    diff[0]=(ninf-n)/tauN;
    diff[1]=(hinf-h)/tauH;
    diff[2]=(rinff-rf)/taurf;
    diff[3]=(rinfs-rs)/taurs;
    diff[4]=(-f)*(eps*(iCa+iT)+ kca*(ca-0.1));
    diff[5]=(-iNa-iK-iCa-iT-iSK-iNap-iH-iA-iL+iap)/C;
    diff[6]=(einf-e)/tauE;
    diff[7]=phirT*(rTinf-rT)/taurT;
    diff[8]=(hpinf-hp)/tauhp;

    return GSL_SUCCESS;

}
// double iNa_v;
// double iK_v;
// double iCa_v;
// double iT_v;
// double iSK_v;
// double iNap_v;
// double iH_v;
// double iA_v;
// double iL_v;
// double ninf_v;
// double tauN;
// double alphaH_v;
// double beta

double dn_n;
double dn_h;
double dn_rf;
double dn_rs;
double dn_ca;
double dn_v;
double dn_e;
double dn_rT;
double dn_hp;
double dh_n;
double dh_h;
double dh_rf;
double dh_rs;
double dh_ca;
double dh_v;
double dh_e;
double dh_rT;
double dh_hp;
double drf_n;
double drf_h;
double drf_rf;
double drf_rs;
double drf_ca;
double drf_v;
double drf_e;
double drf_rT;
double drf_hp;
double drs_n;
double drs_h;
double drs_rf;
double drs_rs;
double drs_ca;
double drs_v;
double drs_e;
double drs_rT;
double drs_hp;
double dca_n;
double dca_h;
double dca_rf;
double dca_rs;
double dca_ca;
double dca_v;
double dca_e;
double dca_rT;
double dca_hp;
double dv_n;
double dv_h;
double dv_rf;
double dv_rs;
double dv_ca;
double dv_v;
double dv_e;
double dv_rT;
double dv_hp;
double de_n;
double de_h;
double de_rf;
double de_rs;
double de_ca;
double de_v;
double de_e;
double de_rT;
double de_hp;
double drT_n;
double drT_h;
double drT_rf;
double drT_rs;
double drT_ca;
double drT_v;
double drT_e;
double drT_rT;
double drT_hp;
double dhp_n;
double dhp_h;
double dhp_rf;
double dhp_rs;
double dhp_ca;
double dhp_v;
double dhp_e;
double dhp_rT;
double dhp_hp;

int 
jacobian(double t, const double y[], double *dfdy, 
         double dfdt[], void *params)
{

    double iNa_v = gNa*h*(-3*(exp((v-thetaM)/sigmaM))*(v-vNa)/(thetaM*pow(1+exp((v-thetaM)/sigmaM),4))
                          + 1/pow(1+exp((v-thetaM)/sigmaM),3));

    double iK_v = gK*pow(n,4);
                        
    double iCa_v = gCa*Ca_ex*(-2*v*exp((v-thetaS)/sigmaS)/(sigmaS*pow(1+exp((v-thetaS)/sigmaS),3)*(1-exp(2*v/RTF))) + (1-exp(2*v/RTF)+2*v*exp(2*v/RTF)/RTF)/pow((1+exp((v-thetaS)/sigmaS))*(1-exp(2*v/RTF)),2)); 
    bTinf = 1/(1+exp((rT-thetabT)/sigmabT))-1/(1+exp(-thetabT/sigmabT));
    double iT_v = gT*Ca_ex*pow(bTinf,3)*(-3*exp((v-thetaaT)/sigmaaT)*v/
                                         (sigmaaT*pow(1+exp((v-thetaaT)/sigmaaT),4)*(1-exp(2*v/RTF))) +
                                         (1-exp(2*v/RTF)+2*v*exp(2*v/RTF)/RTF)/
                                         (pow(1+exp((v-thetaaT)/sigmaaT),3)*pow(1+exp(2*v/RTF),2)));
    double iSK_v = gSK*pow(ca,2)/(pow(ca,2)+pow(ks,2));
    double iNap_v = gNap*hp*(-exp((v-thetaMP)/sigmaMP)*(v-vNa)/(sigmaMP*pow(1+exp((v-thetaMP)/sigmaMP),2)) +
                      1/(1+exp((v-thetaMP)/sigmaMP)));
    double iH_v = gH*(kr*rf+(1-kr)*rs);
    double iA_v = gA*e*(-exp((v-thetaA)/sigmaA)*(v-vK)/(sigmaA*pow(1+exp((v-thetaA)/sigmaA),2)) +
                 1/(1+exp((v-thetaA)/sigmaA)));
    double iL_v = gL;

    double dv_v = (-iNa_v-iK_v-iCa_v-iT_v-iSK_v-iNap_v-iH_v-iA_v-iL_v)/C;

    double dv_n = 4*gK*pow(n,3)*(v-vK)/C;
    double minf = 1/(1+exp((v-thetaM)/sigmaM))/C;
    double dv_h = gNa*pow(minf,3)*(v-vNa)/C;
    double dv_rf = gH*kr*(v-vH)/C;
    double dv_rs = gH*(1-kr)*(v-vH)/C;
    double dv_ca = gSK*(v-vK)*(2*ca*pow(ca,2+pow(ks,2))-2*pow(ca,3))/pow(pow(ca,2)+pow(ks,2),2)/C;
    double dv_e = gA*1/(1+exp((v-thetaA)/sigmaA))*(v-vK)/C;

    double aTinf = 1/(1+exp((v-thetaaT)/sigmaaT));
    double bTinf = 1/(1+exp((rT-thetabT)/sigmabT))-1/(1+exp(-thetabT/sigmabT));
    double bTinf_rT = -exp((rT-thetabT)/sigmabT)/(sigmabT*pow(1+exp((rT-thetabT)/sigmabT),2));
    double iT_rT = gT*pow(aTinf,3)*v*(Ca_ex/(1-exp((2*v)/RTF)))*3*pow(bTinf,2)*bTinf_rT;
    double dv_rT = -iT_rT/C;

    double mpinf = 1/(1+exp((v-thetaMP)/sigmaMP));
    double dv_hp = gNap*mpinf*(v-vNa)/C;

    //dn;
    double ninf = 1/(1+exp((v-thetaN)/sigmaN));
    double tauN = tauNbar/cosh((v-thetaN)/(2*sigmaN));
    double ninf_v = -exp((v-thetaN)/sigmaN)/(sigmaN*pow(1+exp((v-thetaN)/sigmaN),2));
    double tauN_v = -tauNbar*sinh((v-thetaN)/(2*sigmaN))/(2*sigmaN*pow(cosh((v-thetaN)/(2*sigmaN)),2));
    double dn_v = (tauN*ninf_v-(ninf-n)*tauN_v)/pow(tauN,2);
    double dn_n = -1/tauN;
    double dn_h = 0;
    double dn_rf = 0;
    double dn_rs = 0;
    double dn_ca = 0;
    double dn_e = 0;
    double dn_rT = 0;
    double dn_hp = 0;

    //dh;
    double alphaH = 0.128*exp(-(v+50)/18);
    double betaH = 4/(1+exp(-(v+27)/5));
    double hinf = alphaH/(alphaH+betaH);
    double alphaH_v = (-0.128/18)*exp(-(v+50)/18);
    double betaH_v = ((4./5)*exp(-(v+27)/5))/pow(1+exp(-(v+27)/5),2);
    double hinf_v = ((alphaH+betaH)*alphaH_v - alphaH*(alphaH_v+betaH_v))/pow(alphaH+betaH,2);
    double dh_v = hinf_v/tauH;
    double dh_n = 0;
    double dh_h = -1/tauH;
    double dh_rf = 0;
    double dh_rs = 0;
    double dh_ca = 0;
    double dh_e = 0;
    double dh_rT = 0;
    double dh_hp = 0;

    //drf;
    double rinff = 1/(1+exp((v-thetaRF)/sigmaRF));
    double taurf = prf/(-7.4*(v+70)/(exp(-(v+70)/0.8)-1)+65*exp(-(v+56)/23));
    double rinff_v = -exp((v-thetaRF)/sigmaRF)/(sigmaRF*pow(1+exp((v-thetaRF)/sigmaRF),2));
    double taurf_v = prf*(-1.25*(-7.4*v - 518.0)*exp(-1.25*v - 87.5)/pow(exp(-1.25*v - 87.5) - 1,2) + 65*exp(-v/23 - 56./23)/23 + 7.4/(exp(-1.25*v - 87.5) - 1))/pow((-7.4*v - 518.0)/(exp(-1.25*v - 87.5) - 1) + 65*exp(-v/23 - 56./23),2);
    double drf_v = (taurf*rinff_v-(rinff-rf)*taurf_v)/pow(taurf,2);
    double drf_n = 0;
    double drf_h = 0;
    double drf_rf = -1/taurf;
    double drf_rs = 0;
    double drf_ca = 0;
    double drf_e = 0;
    double drf_rT = 0;
    double drf_hp = 0;
    
    //drs;
    double rinfs = 1/(1+exp((v-thetaRS)/sigmaRS));
    double rinfs_v = -exp((v-thetaRS)/sigmaRS)/(sigmaRS*pow(1+exp((v-thetaRS)/sigmaRS),2));
    double drs_v = rinfs_v/taurs;
    double drs_n = 0;
    double drs_h = 0;
    double drs_rf = 0;
    double drs_rs = -1/taurs;
    double drs_ca = 0;
    double drs_e = 0;
    double drs_rT = 0;
    double drs_hp = 0;
    
    //dca;
    double dca_n = 0;
    double dca_h = 0;
    double dca_rf = 0;
    double dca_rs = 0;
    double dca_ca = -f*kca;
    double dca_v = -f*eps*(iCa_v+iT_v);
    double dca_e = 0;
    double dca_rT = -f*eps*iT_rT;
    double dca_hp = 0;
     //de;
    double einf = 1/(1+exp((v-thetaE)/sigmaE));
    double einf_v = -exp((v-thetaE)/sigmaE)/(sigmaE*pow(1+exp((v-thetaE)/sigmaE),2));
    double de_v = einf_v/tauE;
    double de_n = 0;
    double de_h = 0;
    double de_rf = 0;
    double de_rs = 0;
    double de_ca = 0;
    double de_e = -1/tauE;
    double de_rT = 0;
    double de_hp = 0;

    //drT;
    double rTinf = 1/(1+exp((v-thetarT)/sigmarT));
    double taurT = taur0+taur1/(1+exp((v-thrT)/sgmrt));
    double rTinf_v = -exp((v-thetarT)/sigmarT)/(sigmarT*pow(1+exp((v-thetarT)/sigmarT),2));
    double taurT_v = -taur1*exp((v-thrT)/sgmrt)/(sgmrt*pow(1+exp((v-thrT)/sgmrt),2));
    double drT_v = phirT*(taurT*rTinf_v-(rTinf-rf)*taurT_v)/pow(taurT,2);
    double drT_n = 0;
    double drT_h = 0;
    double drT_rf = 0;
    double drT_rs = 0;
    double drT_ca = 0;
    double drT_e = 0;
    double drT_rT = -phirT/taurT;
    double drT_hp = 0;

    //dhp;
    double hpinf = 1/(1+exp((v-thetaHP)/sigmaHP));
    double tauhp = tauHPbar/cosh((v-thetaHP)/(2*sigmaHP));
    double hpinf_v = -exp((v-thetaHP)/sigmaHP)/(sigmaHP*pow(1+exp((v-thetaHP)/sigmaHP),2));
    double tauhp_v = -tauHPbar*sinh((v-thetaHP)/(2*sigmaHP))/(2*sigmaN*cosh((v-thetaN)/(2*sigmaN)));
    double dhp_v = (tauhp*hpinf_v-(hpinf-hp)*tauhp_v)/pow(tauhp,2);
    double dhp_n = 0;
    double dhp_h = 0;
    double dhp_rf = 0;
    double dhp_rs = 0;
    double dhp_ca = 0;
    double dhp_e = 0;
    double dhp_rT = 0;
    double dhp_hp = -1/tauhp;

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 9, 9);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set(m,0,0,dn_n);
    gsl_matrix_set(m,0,1,dn_h);
    gsl_matrix_set(m,0,2,dn_rf);
    gsl_matrix_set(m,0,3,dn_rs);
    gsl_matrix_set(m,0,4,dn_ca);
    gsl_matrix_set(m,0,5,dn_v);
    gsl_matrix_set(m,0,6,dn_e);
    gsl_matrix_set(m,0,7,dn_rT);
    gsl_matrix_set(m,0,8,dn_hp);
    gsl_matrix_set(m,1,0,dh_n);
    gsl_matrix_set(m,1,1,dh_h);
    gsl_matrix_set(m,1,2,dh_rf);
    gsl_matrix_set(m,1,3,dh_rs);
    gsl_matrix_set(m,1,4,dh_ca);
    gsl_matrix_set(m,1,5,dh_v);
    gsl_matrix_set(m,1,6,dh_e);
    gsl_matrix_set(m,1,7,dh_rT);
    gsl_matrix_set(m,1,8,dh_hp);
    gsl_matrix_set(m,2,0,drf_n);
    gsl_matrix_set(m,2,1,drf_h);
    gsl_matrix_set(m,2,2,drf_rf);
    gsl_matrix_set(m,2,3,drf_rs);
    gsl_matrix_set(m,2,4,drf_ca);
    gsl_matrix_set(m,2,5,drf_v);
    gsl_matrix_set(m,2,6,drf_e);
    gsl_matrix_set(m,2,7,drf_rT);
    gsl_matrix_set(m,2,8,drf_hp);
    gsl_matrix_set(m,3,0,drs_n);
    gsl_matrix_set(m,3,1,drs_h);
    gsl_matrix_set(m,3,2,drs_rf);
    gsl_matrix_set(m,3,3,drs_rs);
    gsl_matrix_set(m,3,4,drs_ca);
    gsl_matrix_set(m,3,5,drs_v);
    gsl_matrix_set(m,3,6,drs_e);
    gsl_matrix_set(m,3,7,drs_rT);
    gsl_matrix_set(m,3,8,drs_hp);
    gsl_matrix_set(m,4,0,dca_n);
    gsl_matrix_set(m,4,1,dca_h);
    gsl_matrix_set(m,4,2,dca_rf);
    gsl_matrix_set(m,4,3,dca_rs);
    gsl_matrix_set(m,4,4,dca_ca);
    gsl_matrix_set(m,4,5,dca_v);
    gsl_matrix_set(m,4,6,dca_e);
    gsl_matrix_set(m,4,7,dca_rT);
    gsl_matrix_set(m,4,8,dca_hp);
    gsl_matrix_set(m,5,0,dv_n);
    gsl_matrix_set(m,5,1,dv_h);
    gsl_matrix_set(m,5,2,dv_rf);
    gsl_matrix_set(m,5,3,dv_rs);
    gsl_matrix_set(m,5,4,dv_ca);
    gsl_matrix_set(m,5,5,dv_v);
    gsl_matrix_set(m,5,6,dv_e);
    gsl_matrix_set(m,5,7,dv_rT);
    gsl_matrix_set(m,5,8,dv_hp);
    gsl_matrix_set(m,6,0,de_n);
    gsl_matrix_set(m,6,1,de_h);
    gsl_matrix_set(m,6,2,de_rf);
    gsl_matrix_set(m,6,3,de_rs);
    gsl_matrix_set(m,6,4,de_ca);
    gsl_matrix_set(m,6,5,de_v);
    gsl_matrix_set(m,6,6,de_e);
    gsl_matrix_set(m,6,7,de_rT);
    gsl_matrix_set(m,6,8,de_hp);
    gsl_matrix_set(m,7,0,drT_n);
    gsl_matrix_set(m,7,1,drT_h);
    gsl_matrix_set(m,7,2,drT_rf);
    gsl_matrix_set(m,7,3,drT_rs);
    gsl_matrix_set(m,7,4,drT_ca);
    gsl_matrix_set(m,7,5,drT_v);
    gsl_matrix_set(m,7,6,drT_e);
    gsl_matrix_set(m,7,7,drT_rT);
    gsl_matrix_set(m,7,8,drT_hp);
    gsl_matrix_set(m,8,0,dhp_n);
    gsl_matrix_set(m,8,1,dhp_h);
    gsl_matrix_set(m,8,2,dhp_rf);
    gsl_matrix_set(m,8,3,dhp_rs);
    gsl_matrix_set(m,8,4,dhp_ca);
    gsl_matrix_set(m,8,5,dhp_v);
    gsl_matrix_set(m,8,6,dhp_e);
    gsl_matrix_set(m,8,7,dhp_rT);
    gsl_matrix_set(m,8,8,dhp_hp);

    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    dfdt[2] = 0.0;
    dfdt[3] = 0.0;
    dfdt[4] = 0.0;
    dfdt[5] = 0.0;
    dfdt[6] = 0.0; 
    dfdt[7] = 0.0;
    dfdt[8] = 0.0;

    return GSL_SUCCESS;
   
}

int
main(void)
{
    //FILE *fp;
    //fp = fopen("/home/pmalonis/parameter_fitting/test.txt","w");
    
    //setting up interpollation for current

    gsl_odeiv2_system sys = {DiffEquations, jacobian, 9};
    gsl_odeiv2_driver * d = 
        gsl_odeiv2_driver_alloc_y_new (&sys,gsl_odeiv2_step_rk8pd,1e-6,1e-9,0.0);
    double t = 0.0, t1 = 999.9;
    double y[9] = {n,h,rf,rs,ca,v,e,rT,hp};
    double npoints = 9999;

    // reading injected current vector
    std::string line;
    std::ifstream myfile ("current.csv");
    int nlines = 0;
    while ( getline(myfile,line)) {
        nlines++;
    }
    myfile.clear();
    myfile.seekg(0, std::ios::beg);
    current_recorded = (double*) malloc((nlines)*sizeof(double));
    int i = 0;
    if (myfile.is_open()) {
        while ( getline (myfile,line)) {
            std::string::size_type sz; 
            current_recorded[i] = std::stof(line, &sz);
            i++;
            
        }
    }
    else std::cout << "Unable to open file"; 
    myfile.close();

    // getting recorded time
    time_recorded = (double*) malloc((nlines)*sizeof(double));
    for (i=0; i < nlines; i++) {
        time_recorded[i] = i/sr;
    }

    //setting up interpolation
    interpolation = gsl_interp_alloc (gsl_interp_linear,nlines);
    gsl_interp_init(interpolation, time_recorded, current_recorded, nlines);
    accelerator =  gsl_interp_accel_alloc();

    double param_min_array[] = {450, 70};
    double param_max_array[] = {500,85};
    double param_step_array[] = {10,5};
    int nparams = sizeof(param_min_array)/sizeof(double);
    printf("%d\n", nparams);
    std::vector<double> param_min (param_min_array, 
                                   param_min_array + nparams); 
    std::vector<double> param_max (param_max_array, 
                                   param_max_array + nparams);
    std::vector<double> param_step (param_step_array, 
                                    param_step_array + nparams);
    std::vector<int> ones (param_min.size(), 1);
    std::vector<double> nsteps;
    for (i = 0; i < nparams; i++) {
        printf("max: %f\n", param_max[i]);
        printf("min: %f\n", param_min[i]);
        printf("%f\n", 1 + (param_max[i] - param_min[i])/param_step[i]);
        nsteps.push_back(1 + (param_max[i] - param_min[i])/param_step[i]);
    }

    double n_sets = std::accumulate(nsteps.begin(), 
                                 nsteps.end(), 1, std::multiplies<int>());
    printf("%f\n", n_sets);
   
    // std::vector<int> nsteps = ones + (param_max - param_min)/param_step;
    for (int j = 0; j < 2; j++) {
        gNa += j*10;
        //printf("%f\n", gNa);
        gsl_odeiv2_system sys = {DiffEquations, jacobian, 9};
        gsl_odeiv2_driver * d = 
        gsl_odeiv2_driver_alloc_y_new (&sys,gsl_odeiv2_step_rk8pd,1e-6,1e-9,0.0);
        t=0.0;
        for (i = 0; i <= npoints; i++)
            {
                double ti = i * t1 / npoints;
                int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

                if (status != GSL_SUCCESS)
                    {
                        printf ("error, return value=%d\n", status);
                        break;
                    }

                //fprintf (fp,"%.5e\n", y[5]);
            }
        gsl_odeiv2_driver_free (d);
        
    }    

    //fclose(fp);
    return 0;
}


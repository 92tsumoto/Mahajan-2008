//#ifndef __SYSPARA_H_INCLUDE 
//#define __SYSPARA_H_INCLUDE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "mkl.h"
#include "./lib/xhplot.h"

#define NN 26
#define BUF 200
#define NUM 20

//#define R 8314.472
//#define F 96485.33771638995
//#define T 310
#define R 8.314472		//(J/M*K)
#define F 96.4853415	//(coulomb/mM)
#define T 308		//(K)

#define dvm 5
#define Emax 2000
#define Emin -2000
#define VNMAX (Emax-Emin)*dvm+1

struct varstruct {

	int datas;
	int line_wid[NUM];
	
	int n;
	double Istim,dIstim;
	double coef,dist;
	double s1,s2;
	double RTonF,RTon2F;

	// Drug concentration for Nifekalant
	double drug_conc,ic50,ec50,hill,ehill;
	double normal_block,fraction_facil;
	double block_rate;

	// Model & Cell type & Drug type
	// Cell tupe
	// 0: Endo, 1: Mid, 2: Epi
	int model_type,celltype,drug_type;

	// Cell Geometry
	double RGC,wca;

	// Ion Valences 
	double zna,zk,zca;

	// time ratio
	double ndis;

	// Ion concentrations
	double nao,ki,ko,cao;

	// Change rates for K channel conductances
	double ikrf,iksf,ikxf,itof;
	
	// reversal potential
	double Ena,Ek,Eca,Eks;

	// Slowly activating potassium current
	double prnak;

	// Total Ion currents 
	double INa_total;
	double IK_total;
	double ICa_total;
	
	// Difference total Ion current 
	double dCa_JSR,dICa_total;
	
	// Base Currnt Stimulus
	double Istim_base;

	// Sttimulus parameters
	double BCL;  // Base cycle length = stimulus period
	double dt;  // time step
	int beat; // Number of stimulus

	// debug variable
	double ca_pre,dca_now;

    int m;
    int l;

	double x0[NUM][NN];
    double tsign[NUM];
    double tend[NUM];

	int pflag;
	int write, graph;
	int write0;
	int half;
	int deb;
	int pswitch, sswitch;
	int out_data;

} var;

// Fast sodium current
struct inastruct {
	
	double *Tam,*Tbm,*Tah,*Tbh,*Taj,*Tbj;
	double ina,gna,am,bm,ah,bh,aj,bj;

} ina;

// L-type calcium current
struct ilcastruct {
	
	double *TexpCa;
	double *Tpoinf,*Tpoi,*TPr,*TPs,*Trecov;
	double poinf,poi,Pr,Ps,recov;
	double alpha,beta;
	double fca,po,rxa,cat;
	double s1,s2,s1t,s2t;
	double k1,k2,k3,k4,k5,k6;
	double k1t,k2t,k3t,k4t,k5t,k6t;
	double r1,r2;
	double tau3,taupo,tau_ca,tauca,tauba;
	double pca,gca,gacao,cpt,tca;
	double jca,ilca;
	double s6,vx,sx,vy,sy,vth,vyr,syr;

	double ilcatot;
	double gca_rate;

} ilca;

// inward-rectifired potassium current 
struct ik1struct {
	
	double *Taki,*Tbki;
	double aki,bki;
	double kin;
	double iki,gkix,gki;

} ik1;

// Delayed rectifier potassium current (rapid component)
struct ikrstruct {
	
	double *Txkrv1,*Txkrv2,*Txkrinf,*Trg;
	double xkrv1,xkrv2,xkrinf,rg;
	double gkr,gkr_max,taukr;
	double ikr;

} ikr;

// Delayed rectifier potassium current (slow component)
struct iksstruct {
	
	double *Txs1ss,*Ttauxs1;
	double xs1ss,tauxs1;
	double xs2ss,tauxs2;
	double gks_max,gksx,gks;
	double iks;

} iks;

// Transient outward potassium current
struct itostruct {
	
	double *Txtos_inf,*Tytos_inf,*Ttxs,*Ttys,*Ttxf,*Ttyf;
	double xtos_inf,ytos_inf,txs,tys,txf,tyf;
	double xtof_inf,ytof_inf;
	double rs_inf;
	double gtos,gtof;
	double ito,itos,itof;

} ito;

// Na-K pump current
struct inakstruct {
	
	double *Tfnak;
	double fnak,sigma,gnak,kmnai,kmko;
	double inak;
	
} inak;

// Sodium-Calcium Exchander current
struct incxstruct {
	
	double *Tzw4;
	double zw3,zw4,zw8;
	double aloss,yz1,yz2,yz3,yz4;
	double xkdna;
	double xmcao,xmnao,xmnai,xmcai;
	double csm,jnaca,gnaca,ca;

} incx;

// CICR Irel 
struct jrelstruct {
	
	double *TsparkV,*Texirp;
	double sparkV,exirp;
	double cstar;
	double gryr,gbarsr,gdyad;
	double ax,ay,av,bv;
	double taua,taur;
	double Qr0,Qr,spark_rate;
	double xirp,xicap,xiryr;

} jrel;

// SR uptake,leak
struct jcastruct {

	double vup,kmup,up;
	double gleak,kj,leak;

} jca;

// Ca concentration
struct concstruct {

	double bpxs,bcal,xkcal;
	double spxs,srmax,srkd;
	double mempxs,bmem,kmem;
	double sarpxs,bsar,ksar;
	double dcsib,dciib;
	double bpxi,spxi,mempxi,sarpxi;
	double jd,taud,taups;
	double xbi,xbs;
	double xkon,xkoff,btrop;

} conc;

// for main
void make_ExPTable();
void eular(int n, double h, double x[],double t);
void function(double x[],double f[],double t);
//void input_para(FILE *);
void static_paras(FILE *);
void mem();
void close_mem();

void eventloop(FILE *, int *mode, int *P, double m[]);
void orbit(int *mode, double m[], double x2);
void draw_p(int *mode, int P, double x[], double x2);
void mouse(int *mode, double x[], double x2);

// for dataout
void data_out(FILE *, double t, double u[]);
void current(FILE *,FILE *,FILE *, FILE *, double t,double x[]);

void out_ikr (FILE *, double time, double p[]);
void out_iks (FILE *, double time, double p[]);
void out_ical(FILE *, double time, double p[]);
void out_inaca (FILE *, double time, double p[]);
void out_inak (FILE *, double time, double p[]);
void out_cicr (FILE *, double time, double p[]);

// for calculation of Ion currents
void comp_rev(double x[]);
void comp_ina(double x[]);
void comp_ical(double x[]);
void comp_ikr(double x[]);
void comp_iki(double x[]);
void comp_iks(double x[]);
void comp_ito(double x[]);
void comp_inak(double x[]);
void comp_inaca(double x[]);
void comp_jrel(double x[]);
void comp_jca(double x[]);
void comp_conc (double x[]);

//void main(int argc, char **argv);


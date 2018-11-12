
#include "syspara.h"

inline void function(double x[], double f[], double t)
{
	
	comp_rev(x);
	comp_ina(x);
	comp_ical(x);
	comp_iki(x);
	comp_ikr(x);
	comp_iks(x);
	comp_ito(x);
	comp_inak(x);
	comp_inaca(x);
	comp_jrel(x);
	comp_jca(x);
	comp_conc(x);

	//var.INa_total = ina.ina + 3.0*inak.inak + 3.0*incx.ca;
	//var.IK_total = ikr.ikr + iks.iks + ik1.iki + ito.ito - 2.0*inak.inak + var.Istim;
	//var.ICa_total = ilca.ilca - 2.0*incx.ca;
	var.dCa_JSR = -x[17]+jca.up-jca.leak;

	f[0] = -(ina.ina+ikr.ikr+iks.iks+ik1.iki+ito.ito+ilca.ilca+inak.inak+incx.ca+var.Istim);

	f[1] = ina.am*(1.0 - x[1]) - ina.bm*x[1]; // m
	f[2] = ina.ah*(1.0 - x[2]) - ina.bh*x[2]; // h
	f[3] = ina.aj*(1.0 - x[3]) - ina.bj*x[3]; // j
	
	f[4] = ilca.alpha*x[5] + ilca.k2*x[6] + ilca.k2t*x[7] + ilca.r2*ilca.po - x[4]*(ilca.beta+ilca.r1+ilca.k1t+ilca.k1);	// c1
	f[5] = ilca.beta*x[4] + ilca.k5*x[8] + ilca.k5t*x[9] - x[5]*(ilca.k6+ilca.k6t+ilca.alpha);	// c2
	f[6] = ilca.k1*x[4] + ilca.k4*x[8] + ilca.s1*ilca.po - x[6]*(ilca.k3+ilca.k2+ilca.s2);	// xilca
	f[7] = ilca.k1t*x[4] + ilca.k4t*x[9] + ilca.s1t*ilca.po - x[7]*(ilca.k3t+ilca.k2t+ilca.s2t);	// xilba
	f[8] = ilca.k3*x[6] + ilca.k6*x[5] - x[8]*(ilca.k5+ilca.k4);	// xi2ca
	f[9] = ilca.k3t*x[7] + ilca.k6t*x[5] - x[9]*(ilca.k5t+ilca.k4t);	// xi2ba

	f[10] = (ikr.xkrinf - x[10])/ikr.taukr;	//xr

	f[11] = (iks.xs1ss - x[11])/iks.tauxs1;	//xs1ss
	f[12] = (iks.xs2ss - x[12])/iks.tauxs2;	//xs2ss

	f[13] = (ito.xtos_inf - x[13])/ito.txs; // xtos
	f[14] = (ito.ytos_inf - x[14])/ito.tys; // ytos
	f[15] = (ito.xtof_inf - x[15])/ito.txf; // xtof
	f[16] = (ito.ytof_inf - x[16])/ito.tyf; // ytof

	f[17] = jrel.spark_rate*jrel.Qr - (x[17]/jrel.taur)*(1.0-jrel.taur*var.dCa_JSR/x[22]);	// xir
	f[18] = conc.xbi;	// tropi
	f[19] = conc.xbs;	// trops
	f[20] = conc.dcsib*(50.0*(x[17] - conc.jd - ilca.jca + incx.jnaca) - conc.xbs);	//Ca_submem
	f[21] = jrel.xiryr - (x[21] - x[20])/conc.taups;	// Ca_dyad
	f[22] = -x[17] + jca.up - jca.leak;	// Ca_NSR
	f[23] = (x[22]-x[23])/jrel.taua;	// Ca_JSR
	f[24] = conc.dciib*(conc.jd-jca.up+jca.leak-conc.xbi);	// Ca_i
	f[25] = -(ina.ina+3.0*inak.inak+3.0*incx.ca)/(1000.0*var.wca);	// Na_i

}

void comp_rev(double x[])
{
	incx.csm = x[20]/1000.0;	// cellml org
	var.Ena = var.RTonF*log(var.nao/x[25]);
	//var.Ek = var.RTonF*log(var.ko/var.ki);
	var.Eks = var.RTonF*log((var.ko+var.prnak*var.nao)/(var.ki+var.prnak*x[25]));
	//printf("Ena=%lf,Ek=%lf,Eks=%lf\n",var.Ena,var.Ek,var.Eks);

}

// Fast Sodium Current (time dependant) */
inline void comp_ina(double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ina.am = ina.Tam[iV]*d2 + ina.Tam[iV+1]*d1;
	ina.bm = ina.Tbm[iV]*d2 + ina.Tbm[iV+1]*d1;
	ina.ah = ina.Tah[iV]*d2 + ina.Tah[iV+1]*d1;
	ina.bh = ina.Tbh[iV]*d2 + ina.Tbh[iV+1]*d1;
	ina.aj = ina.Taj[iV]*d2 + ina.Taj[iV+1]*d1;
	ina.bj = ina.Tbj[iV]*d2 + ina.Tbj[iV+1]*d1;

	ina.ina = ina.gna*x[1]*x[1]*x[1]*x[2]*x[3]*(x[0]-var.Ena);

}

// L-type calcium current
inline void comp_ical(double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	double expCa,expNa,expK;

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ilca.poinf = ilca.Tpoinf[iV]*d2 + ilca.Tpoinf[iV+1]*d1;
	ilca.poi = ilca.Tpoi[iV]*d2 + ilca.Tpoi[iV+1]*d1;
	ilca.Pr = ilca.TPr[iV]*d2 + ilca.TPr[iV+1]*d1;
	ilca.Ps = ilca.TPs[iV]*d2 + ilca.TPs[iV+1]*d1;
	ilca.recov = ilca.Trecov[iV]*d2 + ilca.Trecov[iV+1]*d1;
	expCa = ilca.TexpCa[iV]*d2 + ilca.TexpCa[iV+1]*d1;

	ilca.alpha = ilca.poinf/ilca.taupo;
	ilca.beta = (1.0-ilca.poinf)/ilca.taupo;

	ilca.fca = 1.0/(1.0+pow(ilca.cat/x[21],3.0));
	
	//ilca.s1 = 0.02*ilca.fca; // original paper setting
	//ilca.k1 = 0.03*ilca.fca; // original paper setting
	ilca.s1 = 0.0182688*ilca.fca; // cellML setting
	ilca.k1 = 0.024168*ilca.fca; // cellml setting

	ilca.s2 = ilca.s1*(ilca.r1/ilca.r2)*(ilca.k2/ilca.k1);
	ilca.s2t = ilca.s1t*(ilca.r1/ilca.r2)*(ilca.k2t/ilca.k1t);

	ilca.k3 = (1.0-ilca.poi)/ilca.tau3;
	ilca.k3t = ilca.k3;

	ilca.tau_ca = ilca.tca/(1.0+pow((x[21]/ilca.cpt),4.0))+0.1; // cellML setting
	//ilca.tau_ca = (ilca.tca+0.1*pow((1.0+x[21]/ilca.cpt),4))/(1.0+pow((x[21]/ilca.cpt),4.0));	// BJ paper setting
	ilca.tauca = (ilca.recov - ilca.tau_ca)*ilca.Pr + ilca.tau_ca;
	ilca.tauba = (ilca.recov - 450.0)*ilca.Pr + 450.0;

	ilca.k5 = (1.0-ilca.Ps)/ilca.tauca;
	ilca.k5t = (1.0-ilca.Ps)/ilca.tauba;
	ilca.k6 = ilca.fca*ilca.Ps/ilca.tauca;
	ilca.k6t = ilca.Ps/ilca.tauba;

	ilca.k4 = ilca.k3*(ilca.alpha/ilca.beta)*(ilca.k1/ilca.k2)*(ilca.k5/ilca.k6);
	ilca.k4t = ilca.k3t*(ilca.alpha/ilca.beta)*(ilca.k1t/ilca.k2t)*(ilca.k5t/ilca.k6t);

	ilca.po = 1.0 - x[4] - x[5] - x[6] - x[7] - x[8] - x[9];

	//csm = Ca_submem/1000.00;
	if(fabs(x[0]/var.RTon2F) < 0.001){
		ilca.rxa = 4.0*ilca.pca*F/var.RTonF*((incx.csm*expCa-0.341*var.cao)*var.RTon2F); // cellML setting
		//ilca.rxa = 4.0*ilca.pca*F/var.RTonF*((x[20]*expCa-0.341*var.cao)*var.RTon2F); // BJ paper setting
	} else {
		ilca.rxa = 4.0*ilca.pca*x[0]*F/var.RTonF*((incx.csm*expCa-0.341*var.cao)/(expCa-1.0)); // cellML setting
		//ilca.rxa = 4.0*ilca.pca*x[0]*F/var.RTonF*((x[20]*expCa-0.341*var.cao)/(expCa-1.0)); // BJ paper setting
	}

	ilca.jca = ilca.gca*ilca.po*ilca.rxa;

	ilca.ilca = 2.0*var.wca*ilca.jca;

}

// Potassium Current (time-independant) Ik1
inline void comp_iki (double x[])
{
	
	MKL_INT iV=0;
	double V1,V2,d1,d2;

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ik1.aki = ik1.Taki[iV]*d2 + ik1.Taki[iV+1]*d1;
	ik1.bki = ik1.Tbki[iV]*d2 + ik1.Tbki[iV+1]*d1;

	ik1.kin = ik1.aki/(ik1.aki+ik1.bki);

	ik1.iki = ik1.gki*ik1.kin*(x[0]-var.Ek);

}

// Rapidly Activating Potassium Current 
inline void comp_ikr (double x[])
{

	
	MKL_INT iV=0;   
	double V1,V2,d1,d2;
//	double Ikr,Ikr2;

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ikr.xkrv1 = ikr.Txkrv1[iV]*d2 + ikr.Txkrv1[iV+1]*d1;
	ikr.xkrv2 = ikr.Txkrv2[iV]*d2 + ikr.Txkrv2[iV+1]*d1;
	ikr.xkrinf = ikr.Txkrinf[iV]*d2 + ikr.Txkrinf[iV+1]*d1;
	ikr.rg = ikr.Trg[iV]*d2 + ikr.Trg[iV+1]*d1;

	ikr.taukr = 1.0/(ikr.xkrv1 + ikr.xkrv2);

	ikr.ikr = ikr.gkr*x[10]*ikr.rg*(x[0]-var.Ek);
	
	// Drug effects
/*	if(var.drug_type == 0){ // without facilitation effect

		var.ikr_without = var.normal_block*Ikr;
		var.ikr_with    = 0.0;

	} else if(var.drug_type == 1){ // with facilitation effect

		var.ikr_without = var.normal_block*var.fraction_facil*Ikr;
		var.ikr_with    = var.normal_block*(1.0-var.fraction_facil)*Ikr2;
	}
	
	var.ikr = var.ikr_with + var.ikr_without;
*/
}

// Slowly Activating Potassium Current 
inline void comp_iks (double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	iks.xs1ss = iks.Txs1ss[iV]*d2 + iks.Txs1ss[iV+1]*d1;
	iks.xs2ss = iks.xs1ss;
	iks.tauxs1 = iks.Ttauxs1[iV]*d2 + iks.Ttauxs1[iV+1]*d1;
	iks.tauxs2 = 4.0*iks.tauxs1;

	iks.gksx = 1.0+0.8/(1.0+pow((0.5/x[24]),3.0));
	iks.gks = iks.gks_max*iks.gksx;

	//iks.iks = var.gks_rate*iks.gks*x[11]*x[12]*(x[0]-var.Eks);
	iks.iks = iks.gks*x[11]*x[12]*(x[0]-var.Eks);

}

// Transient Outward Current
inline void comp_ito (double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ito.xtos_inf = ito.Txtos_inf[iV]*d2 + ito.Txtos_inf[iV+1]*d1;
	ito.ytos_inf = ito.Tytos_inf[iV]*d2 + ito.Tytos_inf[iV+1]*d1;
	ito.txs = ito.Ttxs[iV]*d2 + ito.Ttxs[iV+1]*d1;
	ito.tys = ito.Ttys[iV]*d2 + ito.Ttys[iV+1]*d1;
	ito.txf = ito.Ttxf[iV]*d2 + ito.Ttxf[iV+1]*d1;
	ito.tyf = ito.Ttyf[iV]*d2 + ito.Ttyf[iV+1]*d1;

	ito.xtof_inf = ito.xtos_inf;
	ito.ytof_inf = ito.ytos_inf;
	ito.rs_inf = ito.ytos_inf;

	ito.itos = ito.gtos*x[13]*(x[14]+0.5*ito.rs_inf)*(x[0]-var.Ek);
	ito.itof = ito.gtof*x[15]*x[16]*(x[0]-var.Ek);

	ito.ito = ito.itos+ito.itof;

}

// Sodium-Potassium Pump NaK ATPase
inline void comp_inak (double x[])
{

	MKL_INT iV=0,iNa=0;
	double V1,V2,d1,d2;
	double Na1,Na2,k1,k2;
	double tab_nak; 

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	inak.fnak=inak.Tfnak[iV]*d2 + inak.Tfnak[iV+1]*d1;

	tab_nak = x[25]/(x[25]+inak.kmnai);

	inak.inak = inak.gnak*inak.fnak*tab_nak*(var.ko/(var.ko+inak.kmko));

}

// Sodium-Calcium Exchanger V-S
inline void comp_inaca (double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	double exp2,exp3; 

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	//incx.csm = x[20];	// Mahajan org
	incx.zw4 = incx.Tzw4[iV]*d2 + incx.Tzw4[iV+1]*d1;
	incx.zw3 = x[25]*x[25]*x[25]*var.cao*exp(0.35*x[0]/var.RTonF) - var.nao*var.nao*var.nao*incx.csm*exp(x[0]*(0.35-1.0)/var.RTonF);
	
	incx.aloss = 1.0/(1.0 + pow(incx.xkdna/x[20],3.0));
	incx.yz1 = incx.xmcao*x[25]*x[25]*x[25] + incx.xmnao*incx.xmnao*incx.xmnao*incx.csm;
	incx.yz2 = incx.xmnai*incx.xmnai*incx.xmnai*var.cao*(1.0+incx.csm/incx.xmcai);
	incx.yz3 = incx.xmcai*var.nao*var.nao*var.nao*(1.0+((x[25]*x[25]*x[25])/(incx.xmnai*incx.xmnai*incx.xmnai)));
	incx.yz4 = x[25]*x[25]*x[25]*var.cao+var.nao*var.nao*var.nao*incx.csm;
	incx.zw8 = incx.yz1 + incx.yz2 + incx.yz3 + incx.yz4;

	incx.jnaca = (incx.gnaca*incx.aloss*incx.zw3)/(incx.zw4*incx.zw8);

    incx.ca = var.wca*incx.jnaca;

}

// Irel Current 
inline void comp_jrel (double x[])
{
	
	MKL_INT iV=0;
	double V1,V2,d1,d2;

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	jrel.sparkV = jrel.TsparkV[iV]*d2 + jrel.TsparkV[iV+1]*d1;
	jrel.exirp = jrel.Texirp[iV]*d2 + jrel.Texirp[iV+1]*d1;

	if(x[23]>50.0 && x[23] < jrel.cstar){
		jrel.Qr0 = x[23] - 50.0;
	} else if(x[23]>= jrel.cstar){
		jrel.Qr0 = jrel.av*x[23] + jrel.bv;
	} else {
		jrel.Qr0 = 0.0;
	}

	jrel.Qr = x[22]*jrel.Qr0/jrel.cstar;

	jrel.spark_rate = jrel.gryr*ilca.po*fabs(ilca.rxa)*jrel.sparkV;
	
	jrel.xirp = jrel.gbarsr*ilca.po*jrel.Qr*fabs(ilca.rxa)*jrel.exirp;	// cellML setting
	//jrel.xirp = jrel.gbarsr*ilca.po*jrel.Qr0*fabs(ilca.rxa)*jrel.exirp;	// BJ paper setting
	jrel.xicap = ilca.po*jrel.gdyad*fabs(ilca.rxa);

	jrel.xiryr = jrel.xirp + jrel.xicap;

}

// SR uptake and leakage
inline void comp_jca (double x[])
{

	jca.up = jca.vup*x[24]*x[24]/(x[24]*x[24]+jca.kmup*jca.kmup);

	jca.leak = jca.gleak*x[22]*x[22]/(x[22]*x[22]+jca.kj*jca.kj)*(x[22]*16.667 - x[24]);

}

// intracellular and subcellular Ca concentration
inline void comp_conc (double x[])
{

	conc.bpxs = conc.bcal*conc.xkcal/(conc.xkcal+x[20])/(conc.xkcal+x[20]);
	conc.spxs = conc.srmax*conc.srkd/(conc.srkd+x[20])/(conc.srkd+x[20]);
	conc.mempxs = conc.bmem*conc.kmem/(conc.kmem+x[20])/(conc.kmem+x[20]);
	conc.sarpxs = conc.bsar*conc.ksar/(conc.ksar+x[20])/(conc.ksar+x[20]);

	conc.dcsib = 1.0/(1.0+conc.bpxs+conc.spxs+conc.mempxs+conc.sarpxs);

	conc.bpxi = conc.bcal*conc.xkcal/(conc.xkcal+x[24])/(conc.xkcal+x[24]);
	conc.spxi = conc.srmax*conc.srkd/(conc.srkd+x[24])/(conc.srkd+x[24]);
	conc.mempxi = conc.bmem*conc.kmem/(conc.kmem+x[24])/(conc.kmem+x[24]);
	conc.sarpxi = conc.bsar*conc.ksar/(conc.ksar+x[24])/(conc.ksar+x[24]);

	conc.dciib = 1.0/(1.0+conc.bpxi+conc.spxi+conc.mempxi+conc.sarpxi);

	conc.jd = (x[20]-x[24])/conc.taud;
	
	conc.xbi = conc.xkon*x[24]*(conc.btrop - x[18])-conc.xkoff*x[18];	// tropi
	conc.xbs = conc.xkon*x[20]*(conc.btrop - x[19])-conc.xkoff*x[19];	// trops
}


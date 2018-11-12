#include "syspara.h"

void static_paras(FILE *fp1)
{
	int i,j,k;

	// general
		var.RTonF = R*T/F;
		var.RTon2F = R*T/(2.0*F);

	// Cell Geometry */
		var.wca = 8.0;	//component cell (mV/microM)

	// Ion Valences
		var.zna = 1;  // Na valence
		var.zk = 1;   // K valence
		var.zca = 2;  // Ca valence
	
	// ion concentrations
		var.ki = 140.0;      // K (mM)
		var.ko = 5.4;      // K (mM)
		var.cao = 1.8;     // Ca (mM)
		var.nao = 136.0;     // Na (mM) Correct value
		var.Ek = var.RTonF*log(var.ko/var.ki);
		printf("EK=%lf\n",var.Ek);

	if(var.celltype==0){ // Endo cell
		var.iksf = 1.0;
	} else if(var.celltype==1){	// Mid cell
		var.iksf = 1.0;
	} else if(var.celltype==2){ // Epi cell
		var.iksf = 1.0;
	}	

	// Sodium channel current
		ina.gna = 12.0;	//(mS/uF);

	// L-type calcium current
		ilca.gca = 182.0;
		ilca.pca = 0.00054;	// Permiability of membrane to Ca (cm/s)
		ilca.vth = 0.0;

		ilca.s6 = 8.0;
		ilca.vx = -40.0;
		ilca.sx = 3.0;
		ilca.vy = -40.0;
		ilca.sy = 4.0;
		ilca.vyr = -40.0;
		ilca.syr = 11.32;

		ilca.cat = 3.0;
		ilca.cpt = 6.09365;

		ilca.k2 = 1.03615E-4;
		ilca.k1t = 0.00413;
		ilca.k2t = 0.00224;
		ilca.r1 = 0.3;
		ilca.r2 = 3.0;
		ilca.s1t = 0.00195;

		ilca.tca = 78.0329;
		ilca.taupo = 1.0;
		ilca.tau3 = 3.0;
		ilca.gacao = 0.341;	// Activity coefficient of Ca

	// Inward rectifier potassium current: IK1
		ik1.gkix = 0.3;	// (microS/nF)
		ik1.gki = ik1.gkix*sqrt(var.ko/5.4);

	// Rapid Activated Potassium Current: IKr
		ikr.gkr_max = 0.0125;
		ikr.gkr = ikr.gkr_max*sqrt(var.ko/5.4); //(microS/nF)
	
	// Slow Activated Potassium Current: IKs
		var.prnak = 0.01833;     // Na/K Permiability Ratio
		iks.gks_max = 0.1386; // (microS/nF)

	// Ito Transient Outward Current
		ito.gtos = 0.04;		// (microS/nF);
		ito.gtof = 0.11;		// (microS/nF);
		
	// Sodium-Potassium Pump (NaK ATPase)
		inak.gnak = 1.5;   // Max. current through Na-K pump (nA/nF)
		inak.kmko = 1.5;    // Half-saturation concentration of NaK pump (mM)
		inak.kmnai = 12.0;    // Half-saturation concentration of NaK pump (mM)
		inak.sigma = (exp(var.nao/67.3)-1.0)/7.0;

	// Sodium-Calcium Exchanger V-S (NCX)
		incx.xkdna = 0.3;   // (microM)
		incx.xmcao = 1.3;   // (mM)
		incx.xmnao = 87.5;  // (mM)
		incx.xmnai = 12.3;	// (mM)
		incx.xmcai = 0.0036;	// (mM)
		incx.gnaca = 0.84;	// CONSTANTS[41] is gNaCa in component INaCa (uM_per_ms).
	
	// CICR, Irel
		jrel.cstar = 90.0; // CONSTANTS[48] is cstar in component Irel (uM).
		jrel.gryr = 2.58079;	// CONSTANTS[49] is gryr in component Irel (per_ms).
		jrel.gbarsr = 26841.8;	// CONSTANTS[50] is gbarsr in component Irel (dimensionless).
		jrel.gdyad = 9000.0;	// CONSTANTS[51] is gdyad in component Irel (mmole_per_coulomb_cm).
		jrel.ax = 0.3576;	// CONSTANTS[52] is ax in component Irel (per_mV).
		jrel.ay = 0.05;	// CONSTANTS[53] is ay in component Irel (per_mV).
		jrel.av = 11.3;	// CONSTANTS[54] is av in component Irel (per_ms).
		jrel.bv = (1.0 - jrel.av)*jrel.cstar - 50.0;	// CONSTANTS[78] is bv in component Irel (uM_per_ms).
		jrel.taua = 100.0;	// CONSTANTS[54] is taua in component Irel (ms).
		jrel.taur = 30.0;	// CONSTANTS[55] is taur in component Irel (ms).

	// NSR Ca Ion Concentration Changes 
		jca.vup = 0.4; //CONSTANTS[58] is vup in component Iup (uM_per_ms).
		jca.kmup = 0.5;	// CONSTANTS[56] is cup in component Iup (uM).
		jca.gleak = 2.069E-5; // CONSTANTS[59] is gleak in component Ileak (per_ms).
		jca.kj = 50.0;	//CONSTANTS[57] is kj in component Ileak (uM).

	// Intracellular and subcellular Ca Ion Concentration Changes 
		conc.bcal = 24.0;	// CONSTANTS[60] is bcal in component Ca (uM).
		conc.xkcal = 7.0;	// CONSTANTS[61] is xkcal in component Ca (uM).
		conc.srmax = 47.0;	// CONSTANTS[62] is srmax in component Ca (uM).
		conc.srkd = 0.6;	// CONSTANTS[63] is srkd in component Ca (uM).
		conc.bmem = 15.0;	// CONSTANTS[64] is bmem in component Ca (uM).
		conc.kmem = 0.3;	// CONSTANTS[65] is kmem in component Ca (uM).
		conc.bsar = 42.0;	// CONSTANTS[66] is bsar in component Ca (uM).
		conc.ksar = 13.0;	// CONSTANTS[67] is ksar in component Ca (uM).
		conc.xkon = 0.0327;	// CONSTANTS[68] is xkon in component Ca (per_uM_per_ms).
		conc.xkoff = 0.0196;	// CONSTANTS[69] is xkoff in component Ca (per_ms).
		conc.btrop = 70.0;	// CONSTANTS[70] is btrop in component Ca (uM).
		conc.taud = 4.0;	// CONSTANTS[71] is taud in component Ca (ms).
		conc.taups = 0.5;	// CONSTANTS[72] is taups in component Ca (ms).
}


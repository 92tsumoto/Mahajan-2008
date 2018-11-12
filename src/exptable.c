#include "syspara.h"

void make_ExpTable()
{

	int vindex,naiindex,kiindex,tindex;
	double v,nai,ki,ton;
    
	for(vindex=0;vindex<VNMAX;vindex++){

		v = (double)vindex/dvm-(double)Emax;

		
        /** for ina **/
		if(fabs(v+47.13) > 0.001){
        	ina.Tam[vindex] = 0.32*(v+47.13)/(1.0-exp(-0.1*(v+47.13)));
		} else {
        	ina.Tam[vindex] = 3.2;
		}
		ina.Tbm[vindex] = 0.08*exp(-v/11.0);
        
		if ( v >= -40.0 ) { // V>=40mV
            ina.Tah[vindex] = 0.0;
            ina.Tbh[vindex] = 1.0/(0.13*(1.0+exp((v+10.66)/-11.1)));
            ina.Taj[vindex] = 0.0;
            ina.Tbj[vindex] = (0.3*exp(-0.0000002535*v))/(1.0+exp(-0.1*(v+32.0)));
        } else { // V < -40mV
            ina.Tah[vindex] = 0.135 * exp(-(80.0+v)/6.8);
            ina.Tbh[vindex] = 3.56*exp(0.079*v)+3.1e5 * exp(0.35*v);
            ina.Taj[vindex] = (-127140.0*exp(0.2444*v)-0.00003474*exp(-0.04391*v))*(v+37.78)/(1.0+exp(0.311*(v+79.23)));
            ina.Tbj[vindex] = (0.1212*exp(-0.01052*v))/(1.0+exp(-0.1378*(v+40.14)));
        }

		// for ical
        ilca.Tpoinf[vindex] = 1.0/(1.0+exp(-(v-ilca.vth)/ilca.s6)); // 1/(1+exp(-(v-vth)/s6)); vth=0.0,s6=8.0
        ilca.Tpoi[vindex] = 1.0/(1.0+exp(-(v-ilca.vx)/ilca.sx)); // 1/(1+exp(-(v-vx)/sx)); vx=-40.0,sx=3.0
        ilca.TPr[vindex] = 1.0 - 1.0/(1.0+exp(-(v-ilca.vy)/ilca.sy)); // 1-(1/(1+exp(-(v-vy)/sy))); vy=-40.0,sy=4.0
        ilca.Trecov[vindex] = 10.0 + 4954.0*exp(v/15.6); 
        ilca.TPs[vindex] = 1.0/(1.0+exp(-(v-ilca.vyr)/ilca.syr)); // 1/(1+exp(-(v-vyr)/syr)); vx=-40.0,sx=11.32

        ilca.TexpCa[vindex] = exp(v/var.RTon2F);

		// for ik1
		ik1.Taki[vindex] = 1.02/(1.0+exp(0.2385*(v-var.Ek-59.215)));
		ik1.Tbki[vindex] = (0.49124*exp(0.08032*(v-var.Ek+5.476))+exp(0.06175*(v-var.Ek-594.31)))/(1.0+exp(-0.5143*(v-var.Ek+4.753)));

        // for ikr 
		if(fabs(v+7.0)>0.001){
        	ikr.Txkrv1[vindex] = 0.001381*(v+7.0)/(1.0-exp(-0.123*(v+7.0)));
		} else {
        	ikr.Txkrv1[vindex] = 0.00138/0.123;
		}
		if(fabs(v+10.0)>0.001){
        	ikr.Txkrv2[vindex] = 0.000611*(v+10.0)/(exp(0.145*(v+10.0))-1.0);
		} else {
        	ikr.Txkrv2[vindex] = 0.00061/0.145;
		}
		ikr.Txkrinf[vindex] = 1.0/(1.0+exp(-(v+50.0)/7.5));
        ikr.Trg[vindex] = 1.0/(1.0+exp((v+33.0)/22.4));
        
        // for iks 
        iks.Txs1ss[vindex] = 1.0/(1.0+exp(-(v-1.5)/16.7));

		if(fabs(v+30.0) < 0.001/0.0687) {    // NaN will be occured at v=-30
			iks.Ttauxs1[vindex] = 1.0 /( 0.0000719/0.148 + 0.000131/0.0687 );
	    } else {
			iks.Ttauxs1[vindex] = 1.0/(0.0000719*(v+30.0)/(1.0-exp(-0.148*(v+30.0)))+0.000131*(v+30.0)/(exp(0.0687*(v+30.0))-1.0));
		}

		// ito
		ito.Txtos_inf[vindex] = 1.0/(1.0+exp(-(v+3.0)/15.0));
		ito.Tytos_inf[vindex] = 1.0/(1.0+exp((v+33.5)/10.0));
		ito.Ttxs[vindex] = 0.5+9.0/(1.0+exp((v+3.0)/15.0));
		ito.Ttys[vindex] = 30.0+3000.0/(1.0+exp((v+60.0)/10.0));
		ito.Ttxf[vindex] = 1.5+3.5*exp(-(v/30.0)*(v/30.0));
		ito.Ttyf[vindex] = 20.0+20.0/(1.0+exp((v+33.5)/10.0));
		
		// inak 
        inak.Tfnak[vindex] = 1.0/(1.0+0.1245*exp(-0.1*v/var.RTonF)+0.0365*inak.sigma*exp(-v/var.RTonF));

		// inaca
        incx.Tzw4[vindex] = 1.0+0.2*exp((0.35-1.0)*v/var.RTonF);
		
		// irel
		jrel.TsparkV[vindex] = exp(-0.05*(v+30.0))/(1.0+exp(-0.05*(v+30.0)));
		jrel.Texirp[vindex] = exp(-0.3576*(v+30.0))/(1.0+exp(-0.3576*(v+30.0)));

	}
	
	// Ena and Nak ATPase
	//for(naiindex=0;naiindex<NAIMAX;naiindex++){
	//	nai = (double)naiindex/dvm;
	//	var.TEna[naiindex] = var.RTonF*log(var.nao/nai);
	//	var.Tnak[naiindex] = 1.0/(1.0+pow(var.kmnai/nai,2.0));
	//}       
	// Ek
	//for(kiindex=0;kiindex<KIMAX;kiindex++){
	//	ki = (double)kiindex/dvm;
	//	var.TEk[kiindex] = var.RTonF*log(var.ko/ki);
		// for Ik1      
	//		var.T5Ki[kiindex]=exp(-0.2385*var.RTonF*log(var.ko/ki));
	//		var.T6Ki[kiindex]=exp(-0.08032*var.RTonF*log(var.ko/ki));
	//		var.T7Ki[kiindex]=exp(-0.06175*var.RTonF*log(var.ko/ki));
	//		var.T8Ki[kiindex]=exp(0.5143*var.RTonF*log(var.ko/ki));
	//}       



}

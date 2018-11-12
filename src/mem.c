#include "syspara.h"

typedef double Number;
typedef long long Lint;

void mem()
{
	int i,k,m,l,z;

	// initialized tablization memorys for Exp functions
	ina.Tam=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Tbm=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Tah=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Tbh=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Taj=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Tbj=(Number *)calloc(VNMAX,sizeof(Number));
	if( ina.Tam==NULL || ina.Tah==NULL || ina.Taj==NULL || ina.Tbm==NULL || ina.Tbh==NULL || ina.Tbj==NULL ) exit(1);
	
	ilca.Tpoinf=(Number *)calloc(VNMAX,sizeof(Number));
	ilca.Tpoi=(Number *)calloc(VNMAX,sizeof(Number));
	ilca.TPr=(Number *)calloc(VNMAX,sizeof(Number));
	ilca.TPs=(Number *)calloc(VNMAX,sizeof(Number));
	ilca.Trecov=(Number *)calloc(VNMAX,sizeof(Number));
	ilca.TexpCa=(Number *)calloc(VNMAX,sizeof(Number));
	if( ilca.Tpoinf==NULL || ilca.Tpoi==NULL || ilca.TPr==NULL || ilca.TPs==NULL || ilca.Trecov==NULL || ilca.TexpCa==NULL ) exit(1);

	ik1.Taki=(Number *)calloc(VNMAX,sizeof(Number));
	ik1.Tbki=(Number *)calloc(VNMAX,sizeof(Number));
	if( ik1.Taki==NULL || ik1.Tbki==NULL) exit(1);

	ikr.Txkrv1=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.Txkrv2=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.Txkrinf=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.Trg=(Number *)calloc(VNMAX,sizeof(Number));
	if( ikr.Txkrv1==NULL || ikr.Txkrv2==NULL || ikr.Txkrinf==NULL || ikr.Trg==NULL) exit(1);

	iks.Txs1ss=(Number *)calloc(VNMAX,sizeof(Number));
	iks.Ttauxs1=(Number *)calloc(VNMAX,sizeof(Number));
	if( iks.Txs1ss==NULL || iks.Ttauxs1==NULL ) exit(1);

	ito.Txtos_inf=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Tytos_inf=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Ttxs=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Ttys=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Ttxf=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Ttyf=(Number *)calloc(VNMAX,sizeof(Number));
	if( ito.Txtos_inf==NULL || ito.Tytos_inf==NULL || ito.Ttxs==NULL || ito.Ttys==NULL || ito.Ttxf==NULL || ito.Ttyf==NULL ) exit(1);

	inak.Tfnak=(Number *)calloc(VNMAX,sizeof(Number));
	if( inak.Tfnak==NULL ) exit(1);

	incx.Tzw4=(Number *)calloc(VNMAX,sizeof(Number));
	if( incx.Tzw4==NULL ) exit(1);

	jrel.TsparkV=(Number *)calloc(VNMAX,sizeof(Number));
	jrel.Texirp=(Number *)calloc(VNMAX,sizeof(Number));
	if( jrel.TsparkV==NULL || jrel.Texirp==NULL ) exit(1);

}
		
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
void close_mem()
{

	free(ina.Tam);free(ina.Tah);free(ina.Taj);free(ina.Tbm);free(ina.Tbh);free(ina.Tbj);
	free(ilca.Tpoinf);free(ilca.Tpoi);free(ilca.TPr);free(ilca.TPs);free(ilca.Trecov); free(ilca.TexpCa);
	free(ik1.Taki); free(ik1.Tbki);
	free(ikr.Txkrv1);free(ikr.Txkrv2);free(ikr.Txkrinf);free(ikr.Trg);
	free(iks.Txs1ss);free(iks.Ttauxs1);
	free(ito.Txtos_inf);free(ito.Tytos_inf);free(ito.Ttxs);free(ito.Ttys);free(ito.Ttxf);free(ito.Ttyf);
	free(inak.Tfnak);
	free(incx.Tzw4);
	free(jrel.TsparkV);free(jrel.Texirp);

}

#include <string.h>
#include <stdlib.h>
#include "syspara.h"

FILE *fopen(), *fpin, *fp0, *fp1, *fp2, *fp3, *fp4, *fp5;
FILE *fp6, *fp7, *fp8, *fp9, *fp10, *fp11, *fp12, *fp13, *fp14, *fp15;
int mode = 1;
int P = 8;
int beats = 10;

typedef double Number;
typedef long long Lint;

main(int argc, char **argv)
{
	int i,w;
	int ii=0;
	double x[NN];
	double t = 0.0;
	double tt,ttt;
	double time=0.0;
	double ttime=0.0;
	double h;
	double v_old,dvdt,dvdt_new;
	double t_stok;
	double ikr_peak=-100.0;
	char *tmpname;
	char cmd[BUFSIZ];
	double tend;
	double apd,rmbp,vmax,toneapd,ttwoapd,dvdtmax;

	tmpname = "temp";

	sprintf(cmd, "/usr/bin/cpp -P %s > %s", argv[1],tmpname);
	if(system(cmd) == -1){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if((fpin=fopen(tmpname,"r"))==NULL){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if ((fp1 = fopen("status.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp2 = fopen("data.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp3 = fopen("ndata.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}

	// parameter inputs
	input_para(fpin);

	if (var.write){
		if ((fp0 = fopen(argv[2],"w"))==NULL){
			fprintf(stderr, "%s cannot open.\n",argv[2]);
			exit(-1);
		}
	}
	if (var.out_data){
		if ((fp4 = fopen("ikr.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp5 = fopen("iks.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp6 = fopen("ical.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
	//	if ((fp7 = fopen("incx.out","w")) == NULL){
	//		printf("Can't open File\n");
	//		exit(1);
	//	}
	//	if ((fp8 = fopen("inak.out","w")) == NULL){
	//		printf("Can't open File\n");
	//		exit(1);
	//	}
		if ((fp9 = fopen("jrel.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}

	}

	xhplot(WINDOW, 700.0, 700.0, WHITE);
	xhplot(DIRECT, 0.0, 0.0, WHITE);

	for (ii = 0; ii < var.datas; ii++){
		long j;
		time = 0.0;
		tend = var.tend[ii];

		for (i = 0; i < NN; i++){
			x[i] = var.x0[ii][i];
		}

		// Eular's step size
		h = 1.0 / var.m;
		h *= var.tsign[ii];
		xddp.line_wid = var.line_wid[ii];
		xhplot(LINEATT,0,0,WHITE);

		// input static parameters
		static_paras(fp1);
		printf("finished input of static parameters\n");

		if(var.model_type == 0){
			var.normal_block = 1.0/(1.0+pow((var.drug_conc/var.ic50),var.hill));
			var.block_rate = var.normal_block;
		} else if(var.model_type == 1){
			var.normal_block=var.block_rate;
			var.drug_conc = var.ic50*pow(1.0/var.block_rate-1.0,1/var.hill);
			fprintf(fp1,"concentration=%e\n",var.drug_conc);
		}
		if(var.drug_type==0){ // without facilitation
			var.fraction_facil = 0.0;
			printf("without facilitation model\n");
		} else if(var.drug_type==1){ // with facilitation
			var.fraction_facil = 1.0/(1.0+pow((var.drug_conc/var.ec50),var.ehill));
			printf("with facilitation model\n");
		}
		printf("fraction=%lf,noaffected=%lf,affected=%lf,total=%lf\n",
		var.fraction_facil,var.block_rate*var.fraction_facil,var.block_rate*(1.0-var.fraction_facil),var.block_rate);
		fprintf(fp1,"fraction=%lf,noaffected=%lf,affected=%lf,total=%lf\n",
		var.fraction_facil,var.block_rate*var.fraction_facil,var.block_rate*(1.0-var.fraction_facil),var.block_rate);


		// initialized memory
		mem();
		printf("exit memory initialization\n");

		printf("Istim=%lf\n",var.Istim_base);

		// Tablize exp functions.   
		printf("start tablization\n");
		make_ExpTable();
		printf("finished tablization\n");

		// Initialization time
		var.beat = 0;
		var.dt = h;

		if(var.deb==1){
			for(i=0;i<NN;i++){
				printf("%lf ", x[i]);
			} printf("\n");
			printf("time=%lf,Istim=%lf\n",time,var.Istim);
			printf("dvdtmax[%d]=%lf\n",var.beat,dvdt);
		//	printf("ENa=%lf, EK=%lf ECa=%lf\n", var.Ena,var.Ek,var.Eca);
		//	printf("vr1=%lf, vr2=%lf vr3=%lf\n", var.vr1,var.vr2,var.vr3);
		}

		tt = var.ndis*(double)var.m*var.BCL;
		ttt = (1.0-var.ndis)*(double)var.m*var.BCL;
		printf("tt=%lf,ttt=%lf\n",tt,ttt);

		while(1){
			eventloop(fp1,&mode,&P,x);
			time=var.beat*var.BCL;

			apd = 0.0; toneapd = 0.0; ttwoapd = 0.0; rmbp = x[0]; vmax = -90.0; dvdtmax = 0.0, v_old = x[0];

			for (j = 0; j< (int)tt; j++){
				t = h*(double)j;
				v_old = x[0];
				var.Istim = var.Istim_base;
				eular(NN,h,x,t);
				dvdt_new = (x[0]-v_old)/(h);
				if(x[0] > vmax){vmax = x[0];}
				if(dvdt_new > dvdtmax){
					dvdtmax = dvdt_new;
					toneapd = time;
				}
				if(dvdt_new < 0 && x[0] >= (vmax - 0.9*(vmax - rmbp))) ttwoapd = time;

				if(ikr_peak < ikr.ikr){
					ikr_peak = ikr.ikr;
				}

				if (var.pflag) orbit(&mode,x,dvdt_new);
				if (var.pswitch==1){
					if(j%10==0){
						data_out(fp2,ttime,x);
						if(var.out_data){
							current(fp4,fp5,fp6,fp9,ttime,x);
						}
					}
				}

				time += h;
				ttime=time;

			} // end 1st loop

			for (j = 0; j< (int)ttt; j++){
				t = h*(double)j;
				v_old = x[0];
				var.Istim = 0.0;
				eular(NN,h,x,t);
				dvdt_new = (x[0]-v_old)/(h); // LRd -> dvdtnew
				if(x[0] > vmax){vmax = x[0];}
				if(dvdt_new > dvdtmax){
					dvdtmax = dvdt_new;
					toneapd = time;
				}
				if(dvdt_new < 0 && x[0] >= (vmax - 0.9*(vmax - rmbp))) ttwoapd = time;

				if(ikr_peak < ikr.ikr){
					ikr_peak = ikr.ikr;
				}

				if (var.pflag) orbit(&mode,x,dvdt_new);
				if (var.pswitch==1){
					if(j%10==0){
						data_out(fp2,ttime,x);
						if(var.out_data){
							current(fp4,fp5,fp6,fp9,ttime,x);
						}
					}
				}

				time += h;
				ttime=time;
			} // end 2nd loop

			fprintf(fp3,"#beats=%d\n",var.beat);
			for(w=0;w<NN;w++){
				fprintf(fp3,"%16.15e\n",x[w]);
			}

			if(var.sswitch==1){
				printf("%d %lf ",var.beat,time);
				for(w=0;w<NN+1;w++){
					if(w!=NN){
						printf("%10.9lf ",x[w]);
					} else {
						printf("/ %10.9lf %10.9e %10.9lf %10.9e\n",ttwoapd - toneapd,var.drug_conc,var.block_rate,ikr_peak);
						fprintf(fp1,"APD90=%lf,Drug_conc=%lf,Block_rate=%lf,IKr_max=%e\n",ttwoapd-toneapd,var.drug_conc,var.block_rate,ikr_peak);
					}
				}
			}
			draw_p(&mode,P,x,dvdt);
			mouse(&mode,x,dvdt);
			if (fabs(time) > tend &&  tend != 0.0) break;
			var.beat++;
			//printf("boolien=%d\n",var.boolien);

			} // end while loop

		} // end for ii-loop


		//Termination and release of memory. */
		fclose(fp1); fclose(fp2); fclose(fp3);
		
		if(var.write0){
			fclose(fp4); fclose(fp5); fclose(fp6); fclose(fp7); fclose(fp9);
		}

		close_mem();


}


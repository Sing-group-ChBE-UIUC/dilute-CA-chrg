#include <stdio.h>
#include <stdlib.h>

int main(){
	int i,j,samp_max,printProps,ntraj,sample,it,n;
	double to,strain;
	unsigned long t,tmax,equil_start;
	double lambda_d,barrier,barrier2,extx,exty,extz;
	int N,tr,s,restart,over_write,flowType,HI_type,e;
	double epsilon,kappas,kappab,dt,Wi,tdtau_max,qmax,tau;
	char tplgy[100],dcmp_meth[100],EV[100],spring[100]; // input strings
	char *txt,*data,*outp;
	FILE *txtfile,*datafile,*outputfile;
	txt = malloc(100*sizeof(char));
	data = malloc(100*sizeof(char));
	outp = malloc(100*sizeof(char));
	N = 100;
	lambda_d = 10.0;
	barrier = 3;
	barrier2 = 3;
	ntraj = 10; //number of trajectories
	it = 1 ;
	int nWi = 10 ;
	int nc = 1 ;
	int Ncharges[1] = {1} ;
	double edot[10] = {1e-2,2e-2,5e-2,1e-1,0.25,0.5,1,2,5,10} ;
	for (e=0;e<nWi;e++){
		for (n=0;n<nc;n++){
			for (tr=1;tr<=ntraj;tr++){



				// sprintf(txt,"txt/P%d_%d_%.4lf_%.2lf_%.2lf_%.5lf_1_2.txt",N,Ncharges[n],lambda_d,barrier,barrier2,edot[e],tr,it);
				// txtfile = fopen(txt,"r");
				// if(!txtfile){
				// 	printf("Error - couldn't find txtfile %s, exiting\n",txt);
				// 	exit(1);
				// }
				// fscanf(txtfile,"tplgy = %s\n",&tplgy);
				// fscanf(txtfile,"HI_type = %d\n",&HI_type);
				// fscanf(txtfile,"dcmp_meth = %s\n",&dcmp_meth);
				// fscanf(txtfile,"EV = %s\n",&EV);
				// fscanf(txtfile,"spring = %s\n",&spring);
				// fscanf(txtfile,"epsilon = %lf\n",&epsilon);
				// fscanf(txtfile,"qmax = %lf\n",&qmax);
				// fscanf(txtfile,"kappas = %lf\n",&kappas);
				// fscanf(txtfile,"kappab = %lf\n",&kappab);
				// fscanf(txtfile,"srpy = %d\n",&s);
				// fscanf(txtfile,"N = %d\n",&N);
				// fscanf(txtfile,"tau = %lf\n",&tau);
				// fscanf(txtfile,"edot = %lf\n",&edot);
				// fscanf(txtfile,"flowType %d\n",&flowType);
				// fscanf(txtfile,"dt = %lf\n",&dt);
				// fscanf(txtfile,"equil_start = %lu\n",&equil_start);
				// fscanf(txtfile,"tmax = %lu\n",&tmax);
				// fscanf(txtfile,"tdtau_max = %lf\n",&tdtau_max);
				// fscanf(txtfile,"restart = %d\n",&restart);
				// fscanf(txtfile,"over_write = %d\n",&over_write);
				// fscanf(txtfile,"trace = %d\n",&tr);
				// fclose(txtfile);
				printProps = 1000;
				tmax = 199999000;
				samp_max = tmax/printProps + 1;
				double *ens_avg_x = calloc(samp_max,sizeof(double));
				long int *ens_count = calloc(samp_max,sizeof(long int));

				for(tr=1;tr<=ntraj;++tr){
					sprintf(data,"prop/E%d_%d_%.4f_%.2f_%.2f_%.5f_%d_%d.txt",N,Ncharges[n],lambda_d,barrier,barrier2,edot[e],tr,it);
					datafile = fopen(data,"r");
					if(!datafile){
						printf("Error - couldn't find datafile %s, exiting\n",data);
						exit(1);
					}
					sample = 0;
					while(!feof(datafile)){
						fscanf(datafile,"%lu %lf %lf %lf %lf\n",&t,&strain,&extx,&exty,&extz);
						ens_avg_x[sample] += extx;
						ens_count[sample]++;
						sample++;
					}
					fclose(datafile);
				}
				sprintf(outp,"ens/E%d_%d_%.4f_%.2f_%.2f_%.5f_%d_%d.txt",N,Ncharges[n],lambda_d,barrier,barrier2,edot[e],1,it);
				outputfile = fopen(outp,"w");
				fprintf(outputfile,"t t*dt <x/L> <x/L>^2\n");
				for(i=0;i<samp_max;++i){
					if(ens_count[i]>0){
						ens_avg_x[i] /= (double)ens_count[i];
						fprintf(outputfile,"%d %lf %lf %lf\n",i*printProps,i*printProps*dt,ens_avg_x[i],ens_avg_x[i]*ens_avg_x[i]);
					}
				}
				fclose(outputfile);
			}
		}
	}
	return 0;
}

#include "main.h" // function definitions and global variable delcarations are here

int main(int argc,const char *argv[]){
	readInput(argc,argv);
	// printf("barrier: %.2lf\n", barrier);
	// printf("barrier2: %.2lf\n", barrier2);
	allocate();
	wall_time = omp_get_wtime();
	for(it=it_start;it<it_max;++it){
		initHI();
		printOutput();
		initChains();
		if(it == it_start){
			cnumber=0;
			for(i=0;i<Ncharges;i++){
					cnumber=ran1(idum)*N;
					if(Charge[cnumber]==0){
						Charge[cnumber]=1;
						Charge_indices[i]=cnumber;
						Charge_track[i]=cnumber;
						Charge_old[i]=cnumber;
				}
				else
				{
					i=i-1;
				}
			}

			// for(i=0;i<Ncharges;i++){
			// 	printf("Charge_track: %d \n",Charge_track[i]);
			// }
		}
		verletList();
		//printf("tstart: %d tmax: %d\n",tstart,tmax);
		tcount=tstart%printProps;
		for(t=tstart;t<tmax;++t){
			strain = ((double)t - (double)equil_start)*dt*edot;
			memset(f,0.0,3*N*sizeof(double));
			bondForce();
			EVForce();
			coulombForce();
			if(t%s==0){
				blockNoise();
			}
			updatePos();
			if(chargehop_type==1){
				chargeHop() ;
			}
			if(t%sample_HI==0){
				sampleHI();
			}
			if(t%print_HI==0){
				printHI();
			}
			checkVerlet();
			calcStress();

			if(t%printProps==0){
				calcConf();
			}
			if(t%printProps==0){
				calcCH();
				for(i=0;i<Ncharges*printProps;i++){
					dxt[i] = 0;
					dyt[i] = 0;
					dzt[i] = 0;
				}
				tcount=0;
			}
			if(t%printxyz==0){
				printTraj();
			}
			if((t%1000==0) && (restart==0)){
				printRestart(); // print the final frame for restarting simulations eg from steady to relax
			}
			calcDis();

			// When developing on the desktop, it is helpful to print some information to the terminal
			// This shows the time, strain, execution time of various functions, and the position of one chain end
			// IMPORTANT: Comment this out when compiling on the cluster
			// Otherwise you will get large output files that will clog up the head node
			// if(t%printxyz==0){
			// 	printf("Ts - %lu strain - %lf 1e4 ts in %lf seconds DC time - %lf seconds Df time - %lf seconds rx[0] - %lf\n",t,strain,omp_get_wtime()-wall_time,decomp_time,update_time,rb[0]);
			// 	wall_time = omp_get_wtime();
			// 	decomp_time = 0.0;
			// 	update_time = 0.0;
			// }
		}
		resetAverage();
		if(calcmsd==1){
			calcMSD();
		}
		//printf("c1: %d c2: %d c3: %d c4: %d c5: %d c6: %d\n",counter1,counter2,counter3,counter4,counter5,counter6);
	}
	return 0;
}
void readInput(int argc,const char *argv[]){
	// Read inputs from the command line and "Input.txt"
	// The command line variables are changed often, whereas variables in the input file are infrequently changed
	// This is convenient for making bash scripts in which you simply change the command line arguments, but only one
	// input file is needed
	FILE *inputfile;
	sscanf(argv[1],"%d",&N);
	sscanf(argv[2],"%d",&Ncharges);
	sscanf(argv[3],"%lf",&Wi);
	sscanf(argv[4],"%lf",&dt);
	sscanf(argv[5],"%lf",&lambda_d);
	sscanf(argv[6],"%lf",&barrier);
	sscanf(argv[7],"%lf",&barrier2);
	sscanf(argv[8],"%d",&tr);
	sscanf(argv[9],"%d",&it_start);
	sscanf(argv[10],"%d",&it_max);

	// There are some print statements throughout the code for testing and development. Generally they can by ignored
	// and commented out
	inputfile = fopen("Input.txt", "r");
	fscanf(inputfile, "tplgy = %s\n", &tplgy);
	fscanf(inputfile, "flowType = %d\n", &flowType);
	fscanf(inputfile, "EV = %s\n", &EV);
	// printf("tplgy %s flowType %d dcmp_meth %s EV %s\n",tplgy,flowType,dcmp_meth,EV);
	if(strcmp(EV,"WCA")==0){
		rc2 = 5.03968; // short cutoff distance - advantage of using WCA over LJ
	}
	else if(strcmp(EV,"LJ")==0){
		rc2 = 25.0; // cutoff distance for LJ is longer because it is not truncated
	}
	else if(strcmp(EV,"NoEV")==0){
		epsilon = 0.0;
		rc2 = 0.0;
	}
	else{
		printf("Error - did not recognize EV type, exiting\n");
		exit(1);
	}
	if(Ncharges>N){
		printf("Error - Ncharges exceeds N, exiting\n");
		exit(1);
	}
	// For WCA use epsilon = 1.0 (athermal solvent). If theta or poor solvent is desired, use LJ
	// For LJ theta solvent is epsilon = 0.31, good solvent is epsilon = 0.1. epsilon > 0.31 is a poor solvent
	// I recommend against turning off EV. It can have weird effects
	fscanf(inputfile, "epsilon = %lf\n", &epsilon);
	fscanf(inputfile, "spring = %s\n", &spring);
	// For Hookean qmax is not needed. For Kremer-Grest qmax = 3.0.
	fscanf(inputfile, "qmax = %lf\n", &qmax);
	qmax2 = qmax*qmax;
	// For Kremer-Grest (FENE), typically kappas = 7.5
	// For Hookean, CES typically uses kappas = 200.0
	fscanf(inputfile, "kappas = %lf\n", &kappas);
	// Change kappab to a positive number to include a bending potential
	// Note: currently bending potential has a small issue for rings, will correct later
	fscanf(inputfile, "kappab = %lf\n", &kappab);
	fscanf(inputfile, "srpy = %d\n", &s);
	fscanf(inputfile, "over_write = %d\n", &over_write);
	fscanf(inputfile, "tmax = %lu\n", &tmax);
	fscanf(inputfile, "lambda_b = %lf\n", &lambda_b);
	fscanf(inputfile, "chargehop_type = %d\n", &chargehop_type);
	fscanf(inputfile, "calcmsd = %d\n", &calcmsd);
	fscanf(inputfile, "MStep = %d\n", &MStep);
	fscanf(inputfile, "tau = %lf\n", &tau);
	fscanf(inputfile, "restart = %d\n", &restart);
	fclose(inputfile);

	edot = Wi/tau;
	// printf("N %d edot %lf tr %d strain_cess %lf\n",N,edot,tr,strain_cess);
	// exit(1);
}
void allocate(){
	// Here I allocate the arrays and pointers for variables and files
	// There are also some random calculations for tmax, etc
	double tau_FD;
	printProps = 2000;//1/dt;
	txt = malloc(sizeof(char)*100);
	xyz = malloc(sizeof(char)*100);
	rg = malloc(sizeof(char)*100);
	ree = malloc(sizeof(char)*100);
	cm = malloc(sizeof(char)*100);
	ext = malloc(sizeof(char)*100);
	strs = malloc(sizeof(char)*100);
	//charge stuff
	str2 = malloc(sizeof(char)*300);
	outp = malloc(sizeof(char)*400);
	outp1 = malloc(sizeof(char)*400);
	rbi = calloc(3*Ncharges, sizeof(double));

	dxt = calloc(Ncharges*printProps, sizeof(double));
	dyt = calloc(Ncharges*printProps, sizeof(double));
	dzt = calloc(Ncharges*printProps, sizeof(double));

	/////
	rst = malloc(sizeof(char)*100);
	Davg = malloc(sizeof(char)*100);
	idum = malloc(sizeof(long));
	rb = calloc(3*N,sizeof(double));
	ro = calloc(3*N,sizeof(double));
	drb = calloc(3*N,sizeof(double));
	nlist = calloc(N,sizeof(int));
	list = calloc(N*N,sizeof(int));
	f = calloc(3*N,sizeof(double));
	Dt = calloc(9*N*N,sizeof(double));
	D = calloc(9*N*N,sizeof(double));
	chol = calloc(9*N*N,sizeof(double));
	Z = calloc(3*N*s,sizeof(double));
	p = sqrt(2.0*dt);
	pc = sqrt(2.0/dt);
	txx=0; //average stress tensor tauxx
	txy=0; //average stress tensor tauxy
	txz=0; //average stress tensor tauxz
	tyy=0; //average stress tensor tauyy
	tyz=0; //average stress tensor tauyz
	tzz=0; //average stress tensor tauzz
	// Heuristic printing interval. 1/dt prints every 1 bead diffusion time

	// printProps = 10000;
	printxyz = 1000;//1/dt;
	// printxyz = 10000;
	sample_HI = 5/dt;
	print_HI = 10*sample_HI;
	// Initialize the random seed
	*idum = initRan();
	// You can also fix the random seed for testing purposes, but be very careful
	//*idum = 1;
	if(*idum==1){
		printf("SEED IS FIXED\n");
	}
	ah = 1.0; // Hydrodynamic radius. For hard potentials, ah is essentially equal to the bead radius
	rc = sqrt(rc2);
	rv = 3.0*2.0; // Heuristic for Verlet list. See Frenkel & Smit p 545 for more details
	rnew = 0.5*(rv-rc);
	omp_set_num_threads(1); // Use only one thread, eg serial simulations
	// tau_FD = 2.0*0.5295*pow(N,1.7592)/8.0;
	tau_FD = 2.0*0.5295*pow(N,1.7592)/4.0; // Esimate of the freely draining (non-HI) relaxation time for a linear chain
	//tmax = (unsigned long)tdtau_max*tau/dt;
	// How long to run the simulation in terms of polymer relaxation times
	// To specify in terms on strain, use strain = edot*t = Wi*t/tau
	// Then tdtau_max = strain/Wi
	// printf("%lu\n",tmax);
	Charge  = calloc(N, sizeof(int));
	Charge_track = calloc(Ncharges, sizeof(int));
	Charge_old = calloc(Ncharges, sizeof(int));
	Charge_indices = calloc(Ncharges, sizeof(int));
	i_values = calloc(N, sizeof(int));
}
void printOutput(){
	// Write the parameters out to a file for documentation
	sprintf(txt,"txt/P%d_%d_%.4lf_%.2lf_%.2lf_%.5lf_%d_%d.txt",N,Ncharges,lambda_d,barrier,barrier2,Wi,tr,it);

	txtfile = fopen(txt,"w");
	fprintf(txtfile,"tplgy = %s\n",tplgy);
	fprintf(txtfile,"EV = %s\n",EV);
	fprintf(txtfile,"spring = %s\n",spring);
	fprintf(txtfile,"epsilon = %lf\n",epsilon);
	fprintf(txtfile,"qmax = %lf\n",qmax);
	fprintf(txtfile,"kappas = %lf\n",kappas);
	fprintf(txtfile,"kappab = %lf\n",kappab);
	fprintf(txtfile,"srpy = %d\n",s);
	fprintf(txtfile,"N = %d\n",N);
	///charge stuff
	fprintf(txtfile,"Ncharges = %d\n",Ncharges);
	fprintf(txtfile,"lambda_d = %lf\n",lambda_d);
	fprintf(txtfile,"lambda_b = %lf\n",lambda_b);
	fprintf(txtfile,"MStep = %d\n",MStep);
	fprintf(txtfile,"barrier = %lf\n",barrier);
	fprintf(txtfile,"barrier2 = %lf\n",barrier2);
	fprintf(txtfile,"chargehop_type = %d\n",chargehop_type);
	//
	fprintf(txtfile,"tau = %lf\n",tau);
	fprintf(txtfile,"edot = %lf\n",edot);
	fprintf(txtfile,"Wi = %lf\n",Wi);
	fprintf(txtfile,"flowType %d\n",flowType);
	fprintf(txtfile,"dt = %lf\n",dt);
	fprintf(txtfile,"equil_start = %lu\n",equil_start);
	fprintf(txtfile,"tmax = %lu\n",tmax);
	fprintf(txtfile,"restart = %d\n",restart);
	fprintf(txtfile,"over_write = %d\n",over_write);
	fprintf(txtfile,"trace = %d\n",tr);
	fprintf(txtfile,"SEED = %ld\n",*idum);
	fclose(txtfile);

	sprintf(xyz,"xyz/R%d_%d_%.4lf_%.2lf_%.2lf_%.5lf_%d_%d.xyz",N,Ncharges,lambda_d,barrier,barrier2,Wi,tr,it);
	sprintf(rg,"prop/RG%d_%d_%.4lf_%.2lf_%.2lf_%.5lf_%d_%d.txt",N,Ncharges,lambda_d,barrier,barrier2,Wi,tr,it);
	sprintf(ree,"prop/REE%d_%d_%.4lf_%.2lf_%.2lf_%.5lf_%d_%d.txt",N,Ncharges,lambda_d,barrier,barrier2,Wi,tr,it);
	sprintf(cm,"prop/CM%d_%d_%.4lf_%.2lf_%.2lf_%.5lf_%d_%d.txt",N,Ncharges,lambda_d,barrier,barrier2,Wi,tr,it);
	sprintf(ext,"prop/E%d_%d_%.4lf_%.2lf_%.2lf_%.5lf_%d_%d.txt",N,Ncharges,lambda_d,barrier,barrier2,Wi,tr,it);
	sprintf(strs,"prop/V%d_%d_%.4lf_%.2lf_%.2lf_%.5lf_%d_%d.txt",N,Ncharges,lambda_d,barrier,barrier2,Wi,tr,it);
	sprintf(str2,"prop/CH%d_%d_%.4lf_%.2lf_%.2lf_%.5lf_%d_%d.txt",N,Ncharges,lambda_d,barrier,barrier2,Wi,tr,it);
	sprintf(outp,"prop/MSD%d_%d_%.4lf_%.2lf_%.2lf_%.5lf_%d_%d.txt",N,Ncharges,lambda_d,barrier,barrier2,Wi,tr,it);
	//sprintf(outp1,"prop/MSD_new_%d_%d_%.4lf_%.2lf_%.2lf_%.5lf_%d_%d.csv",N,Ncharges,lambda_d,barrier,barrier2,Wi,tr,it);
	sprintf(rst,"xyz/RST%d_%d_%.4lf_%.2lf_%.2lf_%.5lf_%d.xyz",N,Ncharges,lambda_d,barrier,barrier2,Wi,it);
}
void initChains(){
	// Initialize chain conformations
	int i,ans,j,test,Ntemp,startcount;
	double dx,dy,dz,r,theta,phi,init_spacing,yrot,xrot,zrot,rxi,ryi,rzi;
	// unsigned long tstart;
	unsigned long told,ttemp;
	startcount = 0;
	init_spacing = 2.05;
	// Start from a random walk, previous trajectory, or a stretched conformation
	// if(flowType==0){
		// After the starting iteration generate a new initial conformation
		if(restart==0 || it>it_start){
			tstart = 0;
			// Initialize linear chains as a random walk
			if(strcmp(tplgy,"linear")==0){
				rb[0] = 0.0;
				rb[1] = 0.0;
				rb[2] = 0.0;
				for(i=1;i<N;++i){
					// printf("%d\n",i);
					test = 0;
					while(test==0){
						test = 1;
						theta = ran1(idum)*2.0*3.14159;
						phi = acos(2.0*ran1(idum)-1.0);
						rb[3*i] = rb[3*(i-1)] + init_spacing*sin(phi)*cos(theta);
						rb[3*i + 1] = rb[3*(i-1) + 1] + init_spacing*sin(phi)*sin(theta);
						rb[3*i + 2] = rb[3*(i-1) + 2] + init_spacing*cos(phi);
						// printf("i - %d test - %d\n",i,test);
						for(j=0;j<i;++j){
							dx = rb[3*i] - rb[3*j];
							dy = rb[3*i + 1] - rb[3*j + 1];
							dz = rb[3*i + 2] - rb[3*j + 2];
							r = sqrt(dx*dx + dy*dy + dz*dz);
							// printf("%d %d %lf %lf %lf %lf %lf %lf %lf\n",i,j,r,rb[3*i],rb[3*i+1],rb[3*i+2],rb[3*j],rb[3*j+1],rb[3*j+2]);
							if(r<init_spacing){
								test = 0;
							}
						}
					}
				}
			}
			// Initialize rings as circular ellipses
			else if(strcmp(tplgy,"ring")==0){
				theta=2*M_PI/(double)N;//how far apart beads should be spaced. This is not in terms of a
				r=(double)N/M_PI;//Circumference in terms of a divided by 2*pi (twos cancel)
				// printf("%lf\n",theta);
				for(i=0;i<N;++i){
					rb[3*i] = r*cos(theta*i);
					rb[3*i+1] = r*sin(theta*i);
					rb[3*i+2] = 0.0;
				}
				// Randomly rotate the ring
				yrot = M_PI*(ran1(idum)-0.5);
				for(i=0;i<N;++i){
					rxi = rb[3*i];
					rzi = rb[3*i+2];
					rb[3*i] = cos(yrot)*rxi + sin(yrot)*rzi;
					rb[3*i+2] = -sin(yrot)*rxi + cos(yrot)*rzi;
				}
				xrot = M_PI*(ran1(idum)-0.5);
				for(i=0;i<N;++i){
					ryi = rb[3*i+1];
					rzi = rb[3*i+2];
					rb[3*i+1] = cos(xrot)*ryi - sin(xrot)*rzi;
					rb[3*i+2] = sin(xrot)*ryi + cos(xrot)*rzi;
				}
				zrot = M_PI*(ran1(idum)-0.5);
				for(i=0;i<N;++i){
					rxi = rb[3*i];
					ryi = rb[3*i+1];
					rb[3*i] = cos(zrot)*rxi - sin(zrot)*ryi;
					rb[3*i+1] = sin(zrot)*rxi + cos(zrot)*ryi;
				}
				// for(i=0;i<N;++i){
				// 	printf("%lf %lf %lf\n",rb[3*i],rb[3*i+1],rb[3*i+2]);
				// }
			}
		}
		// If restarting from an existing simulation, read back in the trajectory
		else if((restart==1)&&(it==it_start)){
			rstfile = fopen(rst,"r");
			if(!rstfile){
				printf("Error - could not find restart file %s, exiting\n",rst);
				exit(1);
			}
			printf("restarting from file %s\n",rst);
			while(!feof(rstfile)){
				fscanf(rstfile,"%d\n%lu\n",&Ntemp,&tstart);
				for(i=0;i<Ntemp;++i){
					fscanf(rstfile,"A %lf %lf %lf\n",&rb[3*i],&rb[3*i+1],&rb[3*i+2]);
				}
			}
			fclose(rstfile);
			tstart++;
			// if(over_write==1){
			// 	printf("Warning - Setting over_write to 0");
			// 	over_write=0;
			// }

		}
		// Start from stretched trajectory. Useful for finding relaxation time
		// Not recommended for dilute CA code, which doesn't handle transience
		else if(restart==2){
			tstart = 0;
			// This start option is not currently available for rings
			if(strcmp(tplgy,"ring")==0){
				printf("Start from stretched conf. not available for rings, exiting\n");
				exit(1);
			}
			rb[0] = 0.0;
			rb[1] = 0.0;
			rb[2] = 0.0;
			for(i=1;i<N;++i){
				rb[3*i] = rb[3*(i-1)] + 0.7*qmax;
				rb[3*i+1] = 0.0;
				rb[3*i+2] = 0.0;
			}
		}
		else{
			printf("Restart code %d not recognized for flowType %d, exiting\n",restart,flowType);
			exit(1);
		}
	// }
	// Don't bother with equilibration for dilute CA simulations
	// Equilibration is most important for getting the correct transient response on startup flow
	// The current dilute CA code doesn't do transience, only steady state, so equilibration is unimportant
	// else{
	// 	// Restarting from previously equilibrated conformation or a previous flow conformation
	// 	// Starting from a random walk is not an option in flow, as this may not be the eqm conf
	// 	if(restart==0 || restart==1){
	// 		if(restart==0){
	// 			sprintf(rst,"xyz/RST%d_%.5f_%d.xyz",N,0.0,tr);
	// 		}
	// 		rstfile = fopen(rst,"r");
	// 		if(!rstfile){
	// 			printf("Error - could not find restart file %s, exiting\n",rst);
	// 			exit(1);
	// 		}
	// 		while(!feof(rstfile)){
	// 			fscanf(rstfile,"%d\n%lu\n",&Ntemp,&ttemp);
	// 			for(i=0;i<Ntemp;++i){
	// 				fscanf(rstfile,"A %lf %lf %lf\n",&rb[3*i],&rb[3*i+1],&rb[3*i+2]);
	// 			}
	// 		}
	// 		fclose(rstfile);
	// 		if(restart==0){
	// 			// start flow from time 0 for a new simulation
	// 			tstart = 0;
	// 			sprintf(rst,"xyz/RST%d_%.5f_%d.xyz",N,edot,tr);
	// 		}
	// 		else{
	// 			tstart = ttemp;
	// 		}
	// 	}
	// 	else{
	// 		printf("Restart code %d not recognized for flowType %d, exiting\n",restart,flowType);
	// 		exit(1);
	// 	}
	// }
	// Advance one time step if restarting as the stopped simulation would have done
	// This prevents double printing for the restart time step
	// If writing over an old trajectory, clean the files first
	if(over_write==1){

		xyzfile = fopen(xyz,"w");
		fprintf(xyzfile,"");
		fclose(xyzfile);
		extfile = fopen(ext,"w");
		fprintf(extfile,"");
		fclose(extfile);
		stressfile = fopen(strs,"w");
		fprintf(stressfile,"");
		fclose(stressfile);
		rgfile = fopen(rg,"w");
		fprintf(rgfile,"");
		fclose(rgfile);
		//charge stuff///
		chfile = fopen(str2,"w");
		fprintf(chfile,"");
		fclose(chfile);
		///////////
		if(strcmp(tplgy,"Linear")==0){
			reefile = fopen(ree,"w");
			fprintf(reefile,"");
			fclose(reefile);
		}
		if(flowType==0){
			cmfile = fopen(cm,"w");
			fprintf(cmfile,"");
			fclose(cmfile);
		}
	}
}
void initHI(){
	int i,j,k,ind1,ind2,info,index;
	unsigned long ttemp;
	if(it>0){
		sprintf(Davg,"Davg/D%d_%.5lf_%d_%d.txt",N,Wi,tr,it-1);
		Davgfile = fopen(Davg,"r");
		if(!Davgfile){
			printf("Error - couldn't find average diffusion tensor file %s to load in, exiting\n",Davg);
			exit(1);
		}
		printf("Reading in avg diffusion tensor from %s\n",Davg);
		// Read the running average diffusion tensor from the previous iteration into D
		// D is used to update positions in the current iteration
		fscanf(Davgfile,"%lu\n%d\n",&ttemp,&count_HI);
		for(i=0;i<9*N*N;++i){
			fscanf(Davgfile,"%lf ",&D[i]);
		}
		fclose(Davgfile);
		// Check the diffusion tensor by printing out
		// Warning - only do this for small chains, the matrix is large
		// printf("Avg diffusion tensor\n");
		// for(i=0;i<3*N;++i){
		// 	for(j=0;j<3*N;++j){
		// 		printf("%lf ",D[3*N*i + j]);
		// 	}
		// 	printf("\n");
		// }
		// printf("\n");
		// Decompose D to get the correlated noise matrix
		info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'L',3*N,D,3*N);
		// The correlated noise matrix "B" is located in the lower triangular
		// Copy over to a new matrix for use in updating positions
		// Set the lower triangular of D to 0. Not strictly necessary because we will use symv, but it is a good practice
		for(i=0;i<3*N;++i){
			// printf("i %d\n",i);
			for(j=0;j<i+1;++j){
				// index = 3*N*i + j;
				// printf("j %d index %d\n",j,index);
				chol[3*N*i + j] = D[3*N*i + j];
				D[3*N*i + j] = 0.0;
			}
		}
		// Check the diffusion tensor by printing out
		// Warning - only do this for small chains, the matrix is large
		// printf("Avg diffusion tensor\n");
		// for(i=0;i<3*N;++i){
		// 	for(j=0;j<3*N;++j){
		// 		printf("%lf ",D[3*N*i + j]);
		// 	}
		// 	printf("\n");
		// }
		// printf("\n");
		// Check the decomposted tensor by printing out
		// Warning - only do this for small chains, the matrix is large
		// printf("Avg decomposed tensor\n");
		// for(i=0;i<3*N;++i){
		// 	for(j=0;j<3*N;++j){
		// 		printf("%lf ",chol[3*N*i + j]);
		// 	}
		// 	printf("\n");
		// }
		// printf("\n");
		count_HI = 0;
		// Set the diagonal of D back to 1.0 for Stokes drag
		for(i=0;i<3*N;++i){
			D[3*N*i+i] = 1.0;
		}
	}
	if(restart==1){
		sprintf(Davg,"Davg/D%d_%.5lf_%d_%d.txt",N,Wi,tr,it);
		Davgfile = fopen(Davg,"r");
		if(!Davgfile){
			printf("Warning - couldn't find average diffusion tensor file %s to load in, continuing without\n",Davg);
			// continue;
			return;
		}
		// Read the running average diffusion tensor from the current iteration in and continue averaging
		fscanf(Davgfile,"%lu\n%d\n",&ttemp,&count_HI);
		for(i=0;i<9*N*N;++i){
			fscanf(Davgfile,"%lf ",&Dt[i]);
			Dt[i] *= count_HI;
		}
	}
	sprintf(Davg,"Davg/D%d_%.5lf_%d_%d.txt",N,Wi,tr,it);
}
void verletList(){
	// Neighbor list for accelerating EV calculation. See Frenkel & Smit p 545 for details
	int i,j;
	double dx,dy,dz,r;
	for(i=0;i<N;++i){
		nlist[i] = 0;
	}
	for(i=0;i<3*N;++i){
		ro[i] = rb[i];
	}
	for(i=0;i<N;++i){
		for(j=i+1;j<N;++j){
			dx = rb[3*i] - rb[3*j];
			dy = rb[3*i + 1] - rb[3*j + 1];
			dz = rb[3*i + 2] - rb[3*j + 2];
			r = sqrt(dx*dx + dy*dy + dz*dz);
			if(r<rv){
				list[i*N + nlist[i]] = j;
				nlist[i]++;
			}
		}
	}
	// for(i=0;i<N;++i){
	// 	printf("%d %d\n",i,list[i]);
	// }
	// exit(1);
}
void checkVerlet(){
	// If bead displacements exceed the skin radius, remake the neighbor list
	int i;
	double dx,dy,dz,r;
	for(i=0;i<N;++i){
		dx = rb[3*i] - ro[3*i];
		dy = rb[3*i + 1] - ro[3*i + 1];
		dz = rb[3*i + 2] - ro[3*i + 2];
		r = sqrt(dx*dx + dy*dy + dz*dz);
		if(r > rnew){
			verletList();
			break;
		}
	}
}
void bondForce(){
	// Spring and bending potentials
	double dx1,dy1,dz1,rr,r,Fs,dx2,dy2,dz2,amp1,amp2,var,theta,coeff,dxn,dyn,dzn,dx,dy,dz; // ,rrn,rn,dx1n,dy1n,dz1n,rr1n,r1n,dx2n,dy2n,dz2n,rr2n,r2n,dx3n,dy3n,dz3n,rr3n,r3n;
	double rr1,r1,rr2,r2,rr3,r3,theta1,theta2,theta3,cost1,cost2,cost3;
	double dxm2,dym2,dzm2,dxm1,dym1,dzm1,dxp1,dyp1,dzp1,rrm2,rm2,rrm1,rm1,rrp1,rp1,thetam1,thetap1,costm1,cost,costp1;
	double dxm2n,dym2n,dzm2n,dxm1n,dym1n,dzm1n,dxp1n,dyp1n,dzp1n,rrn,rn;
	int i,j,monind;
	for(i = 1; i < N; ++i){
		dx = rb[3*i] - rb[3*(i-1)]; dy = rb[3*i+1] - rb[3*(i-1)+1]; dz = rb[3*i+2] - rb[3*(i-1)+2];
		rr = dx*dx + dy*dy + dz*dz;
		r = sqrt(rr);
		dxn = dx/qmax; dyn = dy/qmax; dzn = dz/qmax;
		if(strcmp(spring,"FENE")==0){
			if(rr>qmax2){
				// printTraj();
				printf("FENE overstretch at time step %lu bead %d r %lf\n",t,i,r);
				exit(1);
			}
			Fs = -kappas/(1.0-rr/qmax2);
		}
		else if(strcmp(spring,"Hookean")==0){
			Fs = -kappas*(r-2.0)/r;
		}
		else{
			printf("Error - spring type %s not recognized, exiting\n",spring);
			exit(1);
		}
		f[3*i] += Fs*dx;
		f[3*i+1] += Fs*dy;
		f[3*i+2] += Fs*dz;
		f[3*(i-1)] -= Fs*dx;
		f[3*(i-1)+1] -= Fs*dy;
		f[3*(i-1)+2] -= Fs*dz;
	}
	// Spring between the first and last beads for rings
	if(strcmp(tplgy,"ring")==0){
		dx = rb[0] - rb[3*(N-1)]; dy = rb[1] - rb[3*(N-1)+1]; dz = rb[2] - rb[3*(N-1)+2];
		rr = dx*dx + dy*dy + dz*dz;
		r = sqrt(rr);
		if(strcmp(spring,"FENE")==0){
			if(rr>qmax2){
				// printTraj();
				printf("FENE overstretch ring end time step %lu r %lf\n",t,r);
				exit(1);
			}
			Fs = -kappas/(1.0-rr/qmax2);
		}
		else if(strcmp(spring,"Hookean")==0){
			Fs = -kappas*(r-qmax);
		}
		else{
			printf("Error - spring type %s not recognized, exiting\n",spring);
			exit(1);
		}
		f[0] += Fs*dx;
		f[1] += Fs*dy;
		f[2] += Fs*dz;
		f[3*(N-1)] -= Fs*dx;
		f[3*(N-1)+1] -= Fs*dy;
		f[3*(N-1)+2] -= Fs*dz;
	}
	if(kappab>0.0){
		for(i=1;i<N-1;++i){
			dx1 = rb[3*i] - rb[3*(i-1)]; dy1 = rb[3*i+1] - rb[3*(i-1)+1]; dz1 = rb[3*i+2] - rb[3*(i-1)+2];
			dx2 = rb[3*i] - rb[3*(i+1)]; dy2 = rb[3*i+1] - rb[3*(i+1)+1]; dz2 = rb[3*i+2] - rb[3*(i+1)+2];
			amp1 = sqrt(dx1*dx1+dy1*dy1+dz1*dz1);
			amp2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
			var = (dx1*dx2 + dy1*dy2 + dz1*dz2)/(amp1*amp2);
			theta = M_PI - acos(var); // Bond angle, Flory coordinates convention
			coeff = -kappab*theta/sqrt(1-var*var); // theta_0 = 0 is the relaxed state
			f[3*(i-1)] += coeff*(var*dx1/(amp1*amp1)-dx2/(amp1*amp2));
			f[3*(i-1)+1] += coeff*(var*dy1/(amp1*amp1)-dy2/(amp1*amp2));
			f[3*(i-1)+2] += coeff*(var*dz1/(amp1*amp1)-dz2/(amp1*amp2));
			f[3*i] += coeff*((dx1+dx2)/(amp1*amp2)-var*(dx1/(amp1*amp1)+dx2/(amp2*amp2)));
			f[3*i+1] += coeff*((dy1+dy2)/(amp1*amp2)-var*(dy1/(amp1*amp1)+dy2/(amp2*amp2)));
			f[3*i+2] += coeff*((dz1+dz2)/(amp1*amp2)-var*(dz1/(amp1*amp1)+dz2/(amp2*amp2)));
			f[3*(i+1)] += coeff*(var*dx2/(amp2*amp2)-dx1/(amp1*amp2));
			f[3*(i+1)+1] += coeff*(var*dy2/(amp2*amp2)-dy1/(amp1*amp2));
			f[3*(i+1)+2] += coeff*(var*dz2/(amp2*amp2)-dz1/(amp1*amp2));
		}
	}
}
void EVForce(){
	double dx,dy,dz,rr,r,coeff,r6,ratio;
	int i,j,k;
	if(strcmp(EV,"NoEV")==0){
		return;
	}
	for(i=0;i<N-1;++i){
		for(j=0;j<nlist[i];++j){
			k = list[i*N + j];
			dx = rb[3*i] - rb[3*k];
			dy = rb[3*i + 1] - rb[3*k + 1];
			dz = rb[3*i + 2] - rb[3*k + 2];
			rr = dx*dx + dy*dy + dz*dz;
			if(rr<rc2){
				ratio = 4.00/rr;
				r6 = ratio*ratio*ratio;
				if(strcmp(EV,"WCA")==0){
					coeff = (48.0*epsilon/rr)*(r6*r6-0.5*r6);
				}
				else if(strcmp(EV,"LJ")==0){
					coeff = (12.0*epsilon/rr)*(r6*r6-r6);
				}
				f[3*i] += coeff*dx;
				f[3*i + 1] += coeff*dy;
				f[3*i + 2] += coeff*dz;
				f[3*k] -= coeff*dx;
				f[3*k + 1] -= coeff*dy;
				f[3*k + 2] -= coeff*dz;
			}
		}
	}
}
void coulombForce(){
	double dx,dy,dz,rr,r,coeff,r6,ratio,Fc;
	int i,j,k,l,count;

	// for(i=0;i<N;i++){
	// 	printf("Charge[i]: %d \n", Charge[i]) ;
	// }
	// printf("lambda_d: %lf \n", lambda_d);
	kappa_d = 1.0/lambda_d ;

	for(i=0;i<N;++i){
		for(j=i+1;j<N;++j){
			dx = rb[3*i] - rb[3*j];
			dy = rb[3*i + 1] - rb[3*j + 1];
			dz = rb[3*i + 2] - rb[3*j + 2];
			rr = dx*dx + dy*dy + dz*dz;
			r = sqrt(rr) ;
			//printf("k: %d l: %d \n", k,l) ;
			count++ ;

			if (r<(5*lambda_d) && (Charge[i]!=0) && (Charge[j]!=0)) {// Only between charges aka when bead i and bead j are charged are the only interactions we care for{

				Fc = lambda_b*Charge[i]*Charge[j]*exp(-kappa_d*r)*(1+r*kappa_d)/rr;
				//printf("Fc: %lf \n",Fc);

				f[3*i] += Fc*dx/r; //fx
				f[3*i + 1] += Fc*dy/r; //fy
				f[3*i + 2] += Fc*dz/r; //fz
				f[3*j] -= Fc*dx/r; //fx
				f[3*j + 1] -= Fc*dy/r; //fy
				f[3*j + 2] -= Fc*dz/r; //fz
			}
		}
	}

	// for(i=0;i<Ncharges;++i){
	// 	for(j=i+1;j<Ncharges;++j){
	// 		k = Charge_track[i];
	// 		l = Charge_track[j];
	// 		dx = rb[3*k] - rb[3*l];
	// 		dy = rb[3*k + 1] - rb[3*l + 1];
	// 		dz = rb[3*k + 2] - rb[3*l + 2];
	// 		rr = dx*dx + dy*dy + dz*dz;
	// 		r = sqrt(rr) ;
	// 		//printf("k: %d l: %d \n", k,l) ;
	// 		count++ ;
	//
	// 		if (r<(5*lambda_d)) {// Only between charges aka when bead i and bead j are charged are the only interactions we care for{
	//
	// 			Fc = lambda_b*Charge[k]*Charge[l]*exp(-kappa_d*r)*(1+r*kappa_d)/rr;
	// 			//printf("Fc: %lf \n",Fc);
	//
	// 			f[3*k] += Fc*dx/r; //fx
	// 			f[3*k + 1] += Fc*dy/r; //fy
	// 			f[3*k + 2] += Fc*dz/r; //fz
	// 			f[3*l] -= Fc*dx/r; //fx
	// 			f[3*l + 1] -= Fc*dy/r; //fy
	// 			f[3*l + 2] -= Fc*dz/r; //fz
	// 		}
	// 	}
	// }
}
void chargeHop(){
	int i, j, k,c,l, test, test2, cnumber, cnumber2, count, NCharges,num_of_particles;
	double E_initial, E_final, dx, dy, dz, r, rr , mm, dice_roll;

	NCharges = Ncharges;
	num_of_particles = N ;

	if (t%MStep==0) // Do Monte Carlo Move every MSteps
	{
		for (j=0;j<NCharges;++j)
		{
			for (i=0; i<NCharges;++i)
			{
				if (Charge_track[j]==Charge_indices[i] && i==j)
				{i_values[i]=j;}
			}
		}

		for(i=0;i<NCharges;++i)
		{	Charge_indices[i]=Charge_track[i];}

	for (k=0; k<NCharges; ++k) // Do MC_Move for each charge
	{
		for(i=0;i<num_of_particles;++i)
		{
			Charge[i]=0;
		}

		for(i=0;i<num_of_particles;++i)
		{
			for(j=0;j<NCharges;++j)
			{
				if(i==Charge_track[j])
				Charge[i]=1;

			}
		}

		E_initial=0.0;
		for (i=0; i<num_of_particles; ++i)
		{
			for (j=i+1; j<num_of_particles; ++j)
					{
						if (Charge[i]!=0 && Charge[j]!=0)/// If Bead is charged and testing the next bead neighboring bead
						{
							dx = rb[3*i] - rb[3*j];
							dy = rb[3*i + 1] - rb[3*j + 1];
							dz = rb[3*i + 2] - rb[3*j + 2];
							rr = dx*dx + dy*dy + dz*dz;
							r = sqrt(rr) ;
							E_initial += (lambda_b*Charge[i]*Charge[j]/r)*exp(-kappa_d*r); // Calculating for Initial Energy
				    }
					}
		}
		// printf("this is E_initial=");
		// printf("%lf\n",E_initial);
		mm=ran1(idum);
			//printf("mm=%f\n",mm);

			cnumber=Charge_indices[k]; //find the next charge on the shuffled list
			test=0;// test for am I close?
			test2=0;

//////////// Perform swap if NOT close enough, only to a neighbor bead///////////////

			if (mm>.5) //Move to the right
			{

			cnumber2=cnumber+1; //right hand neighbor
			//printf("t=%d R cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);

				for(i=0;i<NCharges;++i)
				{

					if(cnumber==(num_of_particles-1) || cnumber%(N)==N-1 || cnumber==N-1 || Charge_indices[i]==cnumber2) // don't check if right hand neighbor is charged
					{
						test=3;
					}
				}
				if (test!=0)
				{
					//printf("R Skip\n");
						continue;
				}
				for (i=0;i<NCharges;++i)
				{
					if (Charge_indices[i]==cnumber )
						{
							Charge_indices[i]=cnumber2;
							//printf("R Swap\n");
							test2=0;
							test=0;
						}
				}
			}
			else //move to the left - if random # is less than .5
			{

				cnumber2=cnumber-1; //left hand neighbor
				//printf("t=%d L cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);
				for (i=0;i<NCharges;++i)
				{
					if(cnumber2<0 || cnumber==0 || cnumber%N==0 || Charge_indices[i]==cnumber2)// don't attempt to go left at bead 0
					{
						test=3;
					}
				}
				if (test!=0)
				{
					//printf("L Skip\n");
						continue;
				}

				for (i=0;i<NCharges;++i)
				{
					if (Charge_indices[i]==cnumber)
					{
						Charge_indices[i]=cnumber2;
						//printf("L Swap\n");
						test2=0;
						test=0;

					}
				}
			}


		// Update the charge positions to test new energy
		for(i=0;i<num_of_particles;++i)
		{
			Charge[i]=0;
		}

		for(i=0;i<num_of_particles;++i)
		{
			for(j=0;j<NCharges;++j)
			{
				if(i==Charge_indices[j])
				Charge[i]=1;

			}
		}

		if (test==0 || test==1) // Perform MC but under neighbor conditions - barrier
		{
			//Calculate New Energy
			E_final=0;
			for (i=0; i<num_of_particles; ++i)
			{
				for (j=i+1; j<num_of_particles; ++j)
						{
							if (Charge[i]!=0 && Charge[j]!=0)/// If Bead is charged and testing the next bead neighboring bead
							{
								dx = rb[3*i] - rb[3*j];
								dy = rb[3*i + 1] - rb[3*j + 1];
								dz = rb[3*i + 2] - rb[3*j + 2];
								rr = dx*dx + dy*dy + dz*dz;
								r = sqrt(rr) ;
								E_final += (lambda_b*Charge[i]*Charge[j]/r)*exp(-kappa_d*r); // Calculating for Initial Energy
						 }
						}
			}
			// printf("this is Efinal=");
			// printf("%lf\n",E_final);
			// printf("this delta E=(%lf)\n", (E_final-E_initial));

			//Actual Monte Carlo Test
				dice_roll=ran1(idum);
				if (dice_roll<=exp((-barrier-0.5*(E_final - E_initial)))) // if greater, accept swap
				{	count=0;
											for(i = 0; i<NCharges; ++i)
											{
													if(Charge_track[i]==cnumber && i==i_values[i] && count==0)
													{
														count+=1;
														Charge_old[i]=cnumber;
														Charge_indices[i]=cnumber2;
															Charge_track[i]=cnumber2;
										// 				Output2 = fopen(str2, "a");
										// 				fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
										// fclose(Output2);
										// printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);
													}
											}
											counter2++ ;
				}
				else
				{	count=0;
					for(i=0;i<NCharges; ++i)
					{
						if (Charge_indices[i]==cnumber2 && i==i_values[i] && count==0)
						{	count+=1;
							// printf("count=%d tt=%d Charge_indices=%d\n",count,t,cnumber2);
							Charge_indices[i]=cnumber;
							Charge_old[i]=cnumber;
							}
					}
					counter3++ ;
				}
			}

			///Update position

			for(i=0;i<num_of_particles;++i)
			{
				Charge[i]=0;
			}

			for(i=0;i<num_of_particles;++i)
			{
				for(j=0;j<NCharges;++j)
				{
					if(i==Charge_indices[j])
					Charge[i]=1;

				}
			}

		///////// If there's a possibility of a bead (not neighbor) close by to hop a charge
		cnumber=Charge_indices[k];
		for(j = 0; j<num_of_particles; ++j) //for loop to check distance between cnumber and every other particle
					{
						if (j!=cnumber && j!=cnumber+1 && j!=cnumber-1)
						{
							dx = rb[3*cnumber] - rb[3*j];
							dy = rb[3*cnumber + 1] - rb[3*j + 1];
							dz = rb[3*cnumber + 2] - rb[3*j + 2];
								if(dx*dx+dy*dy+dz*dz<4.41)
								{
									cnumber2=j;
					//printf("t=%d P cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);
					for(i=0;i<NCharges;++i)
					{
						if (Charge_indices[i]==cnumber2)
						{
							test2=1;
						}
					}
					if (test2!=0)
					{
						//printf("Pop Skip\n");
							continue;
					}
					for(i=0;i<NCharges;++i)
					{
						if (Charge_indices[i]==cnumber)
						{
							//printf("Pop\n");
							Charge_indices[i]=cnumber2;
							test2=0;
						}
					}
					test=2;
					break;
								}

						}
				}

		if(test==2) ///Perform MC under NOT neighbor conditions - barrier2
		{
			//Calculate New Energy
			E_final=0;
			for (i=0; i<num_of_particles; ++i)
			{
				for (j=i+1; j<num_of_particles; ++j)
						{
							if (Charge[i]!=0 && Charge[j]!=0)/// If Bead is charged and testing the next bead neighboring bead
							{
								dx = rb[3*i] - rb[3*j];
								dy = rb[3*i + 1] - rb[3*j + 1];
								dz = rb[3*i + 2] - rb[3*j + 2];
								rr = dx*dx + dy*dy + dz*dz;
								r = sqrt(rr) ;
								E_final += (lambda_b*Charge[i]*Charge[j]/r)*exp(-kappa_d*r); // Calculating for Initial Energy
						 }
						}
			}

			//Actual Monte Carlo Test
				dice_roll=ran1(idum);
				if (dice_roll<=exp((-barrier2-0.5*(E_final - E_initial)))) // if greater, accept swap
				{
					count=0;
											for(i = 0; i<NCharges; ++i)
											{
													if(Charge_track[i]==cnumber && i==i_values[i] && count==0)
													{count+=1;
														Charge_old[i]=cnumber;
															Charge_track[i]=cnumber2;
															Charge_indices[i]=cnumber2;
										// 					Output2 = fopen(str2, "a");
										// 					fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
										// fclose(Output2);
										// printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);
													}
											}
											counter5++ ;
				}
				else
				{counter6++ ; //reject and swap back

					count=0;
					for(i=0;i<NCharges; ++i)
					{
						if (Charge_indices[i]==cnumber2 && i==i_values[i] && count==0)
						{count+=1;
							// printf("t= %d Charge_indices=%d\n",t,Charge_indices[i]);
						Charge_indices[i]=cnumber;
						Charge_old[i]=cnumber;
						}
					}


			}
		}
		else{continue;}

		//Update position
		for(i=0;i<num_of_particles;++i)
		{
			Charge[i]=0;
		}

		for(i=0;i<num_of_particles;++i)
		{
			for(j=0;j<NCharges;++j)
			{
				if(i==Charge_track[j])
				Charge[i]=1;

			}
		}


	}//end of k loop
}//end of MC step

	for(i=0;i<num_of_particles;++i)
	{
		Charge[i]=0;
	}

	for(i=0;i<num_of_particles;++i)
	{
		for(j=0;j<NCharges;++j)
		{
			if(i==Charge_track[j])
			Charge[i]=1;
		}
	}

}//end of function
void blockNoise(){
	// Get the random noise vectors
	// This returns a 3*N*s noise matrix
	// There are s columns, one for each time step before updating the RPY tensor
	// Each column holds 3*N Guassian random variables, one for each bead
	// So the structure of the array is:
	// [Zx(0,0) Zx(0,1) ... Zx(0,s-1)]
	// [Zy(0,0) Zy(0,1) ... Zy(0,s-1)]
	// [Zz(0,0) Zz(0,1) ... Zz(0,s-1)]
	// [Zx(1,0) Zx(1,1) ... Zx(1,s-1)]
	// [...........................]
	// [Zz(N-1,0) Zx(N-1,1) ... Zx(0,s-1)]
	// where (x,y,z) are the directions, the first index is the bead #, and the second is the time step
	// Eg Zy(5,3) is the y-component of the random force for the 5th bead during the 3rd time step
	int i,j,k,ind1,ind2,info;
	for(i=0;i<s;i++){
		for(j=0;j<3*N;j++){
			Z[s*j + i] = p*gasdev(idum);
		}
	}
	if(it>0){
		// printf("ts %lu updating with trmm\n",t);
		double start_time = omp_get_wtime();
		// Triangular matrix-matrix multiplication. The decomposed matrix B is lower triangular
		cblas_dtrmm(CblasRowMajor,CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,3*N,s,1.0,chol,3*N,Z,s);
		// cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,3*N,s,3*N,1.0,chol,3*N,Z,s,0.0,Z,s);
		decomp_time += omp_get_wtime() - start_time;
	}
}
void updatePos(){
	//Save previous particle positions of all charges//
	for(i=0;i<Ncharges;i++){
		j = Charge_track[i];
		rbi[3*i] = rb[3*j];
		rbi[3*i+1] = rb[3*j+1];
		rbi[3*i+2] = rb[3*j+2];
	}
	// Update particle positions
	int i,j;
	double comx,comy,comz;
	double start_time = omp_get_wtime();
	// Applied solvent flow rate
	// flowType = 0 -> equilibrium
	if(flowType==0){}
	else if(flowType==1){ // PEF
		for(i=0;i<N;++i){
			rb[3*i] += dt*edot*rb[3*i];
			rb[3*i+1] -= dt*edot*rb[3*i+1];
		}
	}
	else if(flowType==2){
		for(i=0;i<N;++i){ // PSF
			rb[3*i] += dt*edot*rb[3*i+1];
		}
	}
	else if(flowType==3){
		for(i=0;i<N;++i){ // PRF
			rb[3*i] += dt*edot*rb[3*i+1];
			rb[3*i+1] -= dt*edot*rb[3*i];
		}
	}
	else{
		printf("Error - flowType %d not recognized, exiting\n",flowType);
		exit(1);
	}
	// Deterministic part, D.F in the Langevin equation
	if(it==0){
		for(i=0;i<3*N;i++){
			drb[i] = dt*f[i];
		}
	}
	else{
		// printf("ts %lu updating with dsymv\n",t);
		cblas_dsymv(CblasRowMajor,CblasUpper,3*N,dt,D,3*N,f,1,0.0,drb,1);
	}
	// Add the determinstic part and random part stored in Z to the current position
	// Recall the form of the random vector array to understand the indexing
	for(i=0;i<3*N;++i){
		rb[i] += drb[i] + Z[t%s + s*i];
	}


	// Subtract off the center off mass when in flow
	if(flowType!=0){
		comx = 0.0; comy = 0.0; comz = 0.0;
		for(i=0;i<N;++i){
			comx += rb[3*i];
			comy += rb[3*i+1];
			comz += rb[3*i+2];
		}
		comx /= (double)N;
		comy /= (double)N;
		comz /= (double)N;
		for(i=0;i<N;++i){
			rb[3*i] -= comx;
			rb[3*i+1] -= comy;
			rb[3*i+2] -= comz;
		}
	}
	update_time += omp_get_wtime() - start_time;
}
void sampleHI(){
	int i,j,ind1,ind2;
	double dx,dy,dz,r,rr,C1,C2;
	for(i=0;i<3*N;++i){
		Dt[3*N*i + i] += 1.0;
	}
	for(i=0;i<N;++i){
		ind1 = 9*N*i;
		for(j=i+1;j<N;++j){
			ind2 = 3*j;
			dx = rb[3*i] - rb[3*j];
			dy = rb[3*i + 1] - rb[3*j + 1];
			dz = rb[3*i + 2] - rb[3*j + 2];
			rr = dx*dx + dy*dy + dz*dz;
			r = sqrt(rr);
			if(r>=2.0*ah){
				// C1 = 0.75/r+0.5/(r*rr);
				// C2 = (0.75-1.5/rr)/(r*rr);
				C1 = 3.0*ah/(4.0*r) + ah*ah*ah/(2.0*r*rr);
				C2 = 3.0*ah/(4.0*r*rr)*(1-2.0*ah*ah/rr);
			}
			else{
				// C1 = 1.0-0.28125*r;
				// C2 = 0.09375/r;
				C1 = 1.0 - 9.0*r/(32.0*ah);
				C2 = 3.0/(32.0*ah*r);
			}
			Dt[ind1 + ind2] += C1 + C2*dx*dx;
			Dt[ind1 + ind2 + 1] += C2*dx*dy;
			Dt[ind1 + ind2 + 2] += C2*dx*dz;
			Dt[3*N + ind1 + ind2] += C2*dx*dy;
			Dt[3*N + ind1 + ind2 + 1] += C1 + C2*dy*dy;
			Dt[3*N + ind1 + ind2 + 2] += C2*dy*dz;
			Dt[6*N + ind1 + ind2] += C2*dx*dz;
			Dt[6*N + ind1 + ind2 + 1] += C2*dy*dz;
			Dt[6*N + ind1 + ind2 + 2] += C1 + C2*dz*dz;
		}
	}
	count_HI++;
	// Only one half of the diffusion tensor needs to be calculated
	// The other half can be filled in by symmetry
	for(i=0;i<3*N;++i){
		for(j=i;j<3*N;++j){
			Dt[3*N*j + i] = Dt[3*N*i + j];
		}
	}
}
void printHI(){
	int i,j;
	sprintf(Davg,"Davg/D%d_%.5lf_%d_%d.txt",N,Wi,tr,it);
	Davgfile = fopen(Davg,"w");
	fprintf(Davgfile,"%lu\n%d\n",t,count_HI);
	for(i=0;i<9*N*N;++i){
		fprintf(Davgfile,"%lf ",Dt[i]/count_HI);
	}
	fprintf(Davgfile,"\n");
	fclose(Davgfile);
}
void resetAverage(){
	int i,j;
	printHI();
	for(i=0;i<9*N*N;++i){
		Dt[i] = 0.0;
		chol[i] = 0.0;
	}
	count_HI = 0;
}
void printTraj(){
	int i;
	////charge stuff////
	xyzfile = fopen(xyz,"a");
	fprintf(xyzfile,"%d\n%lu\n",N,t);
	for(i=0;i<N;++i){
		if(Charge[i] == 0){
			fprintf(xyzfile,"A %lf %lf %lf\n",rb[3*i],rb[3*i+1],rb[3*i+2]);
		}
		else if(Charge[i] == 1){
			fprintf(xyzfile,"B %lf %lf %lf\n",rb[3*i],rb[3*i+1],rb[3*i+2]);
		}
	}
	fclose(xyzfile);
	////////////
	if(restart==0){
		rstfile = fopen(rst,"w");
		fprintf(rstfile,"%d\n%lu\n",N,t);
		for(i=0;i<N;++i){
			fprintf(rstfile,"A %lf %lf %lf\n",rb[3*i],rb[3*i+1],rb[3*i+2]);
		}
		fclose(rstfile);
	}

}
void calcStress(){
	// Do this every time step for better averaging. Stress is a very noisy quantity
	int i;
	// Stress tensor is symmetric, 6 unique components
	for(i=0;i<N;++i){
		stress[0] -= rb[3*i]*f[3*i]; // tau_xx
		stress[1] -= rb[3*i]*f[3*i+1]; // tau_xy
		stress[2] -= rb[3*i]*f[3*i+2]; // tau_xz
		stress[3] -= rb[3*i+1]*f[3*i+1]; // tau_yy
		stress[4] -= rb[3*i+1]*f[3*i+2]; // tau_yz
		stress[5] -= rb[3*i+2]*f[3*i+2]; // tau_zz
	}
	if(t%printProps==0){
		stressfile = fopen(strs,"a");
		fprintf(stressfile,"%lu %lf %lf %lf %lf %lf %lf %lf\n",t,strain,stress[0]/printProps,stress[1]/printProps,stress[2]/printProps,stress[3]/printProps,stress[4]/printProps,stress[5]/printProps);
		fclose(stressfile);
		for(i=0;i<6;++i){
			stress[i] = 0.0;
		}
	}
}
void calcConf(){
	// Calculate and print the conformational properties of the ring
	int i,j;
	double comx,comy,comz,rg2,rgt,reex,reey,reez,dx,dy,dz,maxx,minx,extx,extz,maxz,minz,exty,maxy,miny;
	double stress[6],rgp[6];
	reex = rb[3*0] - rb[3*(N-1)];
	reey = rb[3*0+1] - rb[3*(N-1)+1];
	reez = rb[3*0+2] - rb[3*(N-1)+2];
	maxx = rb[0]; minx = rb[0];
	for(i=0;i<N;++i){
		maxx = max(maxx,rb[3*i]);
		minx = min(minx,rb[3*i]);
	}
	extx = (maxx - minx);
	maxy = rb[1]; miny = rb[1];
	for(i=0;i<N;++i){
		maxy = max(maxy,rb[3*i+1]);
		miny = min(miny,rb[3*i+1]);
	}
	exty = (maxy - miny);
	maxz = rb[2]; minz = rb[2];
	for(i=0;i<N;++i){
		maxz = max(maxz,rb[3*i+2]);
		minz = min(minz,rb[3*i+2]);
	}
	extz = (maxz - minz);
	// Calculation of the gyration tensor components
	// Gyration tensor is symmetric, so there are 6 unique components
	for(i=0;i<6;++i){
		rgp[i] = 0.0;
	}

	//////
	comx = 0.0; comy = 0.0; comz = 0.0;
	for(i=0;i<N;++i){
		comx += rb[3*i];
		comy += rb[3*i+1];
		comz += rb[3*i+2];
	}
	comx /= (double)N;
	comy /= (double)N;
	comz /= (double)N;
	//////
	for(i=0;i<N;++i){
		rgp[0] += (rb[3*i]-comx)*(rb[3*i]-comx); // rgxx
		rgp[1] += (rb[3*i]-comx)*(rb[3*i+1]-comy); // rgxy
		rgp[2] += (rb[3*i]-comx)*(rb[3*i+2]-comz); // rgxz
		rgp[3] += (rb[3*i+1]-comy)*(rb[3*i+1] - comy); // rgyy
		rgp[4] += (rb[3*i+1]-comy)*(rb[3*i+2]-comz); // rgyz
		rgp[5] += (rb[3*i+2]-comz)*(rb[3*i+2]-comz); // rgzz
	}
	rgp[0] /= (double)N; rgp[1] /= (double)N; rgp[2] /= (double)N;
	rgp[3] /= (double)N; rgp[4] /= (double)N; rgp[5] /= (double)N;
	// End-to-end vector is zero for a ring
	if(strcmp(tplgy,"linear")==0){
		reefile = fopen(ree,"a");
		fprintf(reefile,"%lu %lf %lf %lf %lf\n",t,strain,reex,reey,reez);
		fclose(reefile);
	}
	extfile = fopen(ext,"a");
	if(strcmp(tplgy,"linear")==0){
		fprintf(extfile,"%lu %lf %lf %lf %lf\n",t,strain,extx/((N-1)*qmax),exty/((N-1)*qmax),extz/((N-1)*qmax));
	}
	else if(strcmp(tplgy,"ring")==0){ // Ring polymer contour length is half that of the linear chains
		fprintf(extfile,"%lu %lf %lf %lf %lf\n",t,strain,extx/(0.5*(N-1)*qmax),exty/(0.5*(N-1)*qmax),extz/(0.5*(N-1)*qmax));
	}
	fclose(extfile);
	rgfile = fopen(rg,"a");
	fprintf(rgfile,"%lu %lf %lf %lf %lf %lf %lf %lf\n",t,strain,rgp[0],rgp[1],rgp[2],rgp[3],rgp[4],rgp[5]);
	fclose(rgfile);
	if(flowType==0){
		comx = 0.0;
		comy = 0.0;
		comz = 0.0;
		for(i=0;i<N;++i){
			comx += rb[3*i];
			comy += rb[3*i+1];
			comz += rb[3*i+2];
		}
		comx /= (double)N;
		comy /= (double)N;
		comz /= (double)N;
		cmfile = fopen(cm,"a");
		fprintf(cmfile,"%lu %lf %lf %lf\n",t,comx,comy,comz);
		fclose(cmfile);
	}
}
void calcMSD()
{
  int nsamp,ttemp,ind1,MSDstart,ind2,count,m,l, wmax,Tmax,q,k,ctemp,n_it , n_nc,nc,ttemp1, i, j, jcount, nb1, n_nb1, nb2, n_nb2,nedot,n_nedot,n_nit,nit, dint,dint2;
  double dcomx,dcomy,dcomz,avg,dx,dy,dz,rxtemp,rytemp,rztemp,MSDtemp,avgmsd;

	outputfile = fopen(outp,"w");
	fprintf(outputfile,"");
	fclose(outputfile);
	
	tmax = tmax - tstart + 1;
  Tmax = tmax/printProps;
	wmax = Tmax+1;
	nsamp = Ncharges*wmax;


	// if(restart==0){
	// 	MSDstart=1;
	// }
	// else if (restart==1){
	// 	MSDstart = (int)(tstart-1)/printProps+1;
	// }
		// printf("MSDstart: %d\n",MSDstart);

  double *msd = calloc(nsamp, sizeof(double));
  double *allmsd = calloc(nsamp, sizeof(double));
  //double *avgmsd = calloc(nsamp, sizeof(double));
  double *msdcount = calloc(nsamp, sizeof(double));
	double *vx = calloc(nsamp, sizeof(double));
	double *vy = calloc(nsamp, sizeof(double));
	double *vz = calloc(nsamp, sizeof(double));
	double *Ree2 = calloc(Ncharges, sizeof(double));
	double *Rcount = calloc(Ncharges, sizeof(double));
	FILE *datafile;

	datafile = fopen(str2,"r");
	if(!datafile)
	{
			printf("Error - datafile not found, check data string %s\n",str2);
			exit(1);
	}
	//printf("before read");
	for(i=0;i<nsamp;i++)
	{
			fscanf(datafile, "%d %d %d %d %lf %lf %lf\n", &ttemp, &ctemp,&dint2, &dint, &rxtemp, &rytemp, &rztemp);

			vx[i] = rxtemp;

			vy[i] = rytemp;

			vz[i] = rztemp;

	}
	fclose(datafile);

		// for(i=0;i<Ncharges*tmax;i++){
		// 	if(i%1000==0){
		// 		printf("i: %d dxt: %lf\n",i,dxt[i]);
		// 	}
		//
		// }

		// for(i=0;i<Ncharges*Tmax;i++){
		// 	printf("vx_new: %lf\n",vx[i]);
		// }


  //Calculate MSD//


	outputfile=fopen(outp,"a");
	for(i = 1; i<wmax-1; ++i) //i is shift
	{
		for(l=0;l<Ncharges;l++){
			Ree2[l] = 0 ;
			Rcount[l] = 0;
		}
		avg = 0 ;
		for(q=0;q<Ncharges;q++){
			for(j = 1; j<wmax-i; ++j)
			{
				dx = 0.0; dy = 0.0; dz = 0.0;
				for(k = 0; k<i; ++k)
				{
					dx += vx[j+k+q*Tmax];
					dy += vy[j+k+q*Tmax];
					dz += vz[j+k+q*Tmax];
				}
				Ree2[q] += dx*dx+dy*dy+dz*dz;
				Rcount[q]++;
			}
			avg+=Ree2[q]/(double)Rcount[q];
		}
		//average MSD over charges//
		//print average to file//
		fprintf(outputfile,"%f,%lf\n",i*dt*printProps,avg/Ncharges);

	}
	fclose(outputfile);
	//printf("MSD End\n");
}
void calcDis()
{

	for(i=0;i<Ncharges;i++){
		j = Charge_track[i];
		dxt[i*printProps+tcount] = rb[3*j] - rbi[3*i]  ;
		dyt[i*printProps+tcount] = rb[3*j+1] - rbi[3*i+1]  ;
		dzt[i*printProps+tcount] = rb[3*j+2] - rbi[3*i+2]  ;
	}
	tcount++ ;
}
void calcCH()
{
	int i,j,ind2;
	double dx,dy,dz;
	///charge stuff///
	chfile = fopen(str2,"a");
	for(j=0;j<Ncharges;j++){
		dx=0;
		dy=0;
		dz=0;
		for(i=0;i<printProps;i++){
				ind2 = j*printProps+i ;
				dx += dxt[ind2] ;
				dy += dyt[ind2] ;
				dz += dzt[ind2] ;
			}
		fprintf(chfile,"%lu %d %d %d %lf %lf %lf\n",t,j,Charge_old[j],Charge_track[j],dx,dy,dz);
	}
	fclose(chfile);
	//////////////////
}
void printRestart(){
	rstfile = fopen(rst,"w");
	fprintf(rstfile,"%d\n%lu\n",N,t);
	for(i=0;i<N;++i){
		fprintf(rstfile,"A %lf %lf %lf\n",rb[3*i],rb[3*i+1],rb[3*i+2]);
	}
	fclose(rstfile);
}
long initRan(){
    //time_t seconds;
    //time(&seconds);
    //return -1*(unsigned long)(seconds/12345); This is bad.  :(

    //This will hopefully allow us to have a unique seed even if executed multiple times a second-Got from Mike
    //http://stackoverflow.com/questions/322938/recommended-way-to-initialize-srand
    unsigned long a = clock();
    unsigned long b = time(NULL);
    unsigned long c = getpid();
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c%1000000000; //careful here.  Another 0 might break the ran1 (long long instead of just long)
}
float ran1(long *idum){
	int j;
	long k;
	static long idum2 = 123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if(*idum <= 0)
	{
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2 = (*idum);
		for(j=NTAB+7;j>=0;--j)
		{
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if(*idum<0) *idum+=IM1;
			if(j<NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if(*idum<0) *idum += IM1;
	k=idum2/IQ2;
	if(*idum<0) idum2+= IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if(iy<1) iy += IMM1;
	if((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
float gasdev(long *idum){
	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if (*idum < 0) iset = 0;
	if (iset == 0)
	{
		do
		{
			v1 = 2.0*ran1(idum)-1.0;
			v2 = 2.0*ran1(idum)-1.0;
			rsq = v1*v1+v2*v2;
		}
		while(rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	}
	else
	{
		iset = 0;
		return gset;
	}
}

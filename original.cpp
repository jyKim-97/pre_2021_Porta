#include<cmath>
#include<cstdlib>
#include<iostream>
#include<cstring>
#include<cstdio>


using namespace std;


/*Sender = Population 2 = Alpha oscillations*/
/*Receiver = Population 1 = Gamma oscillations*/

/*The indexes from 0 to 399 stand for excitatory neurons and from 400-499 for inhibitory*/

#define dt 0.05 /*step of integration ms*/
#define gES 0.8 /*Sender internal excitatory synaptic weight (ns)*/
#define gIS 16.4 /*Sender internal inhibitory synaptic weight (ns)*/
// #define gES 0.0 /*Sender internal excitatory synaptic weight (ns)*/
// #define gIS 0. /*Sender internal inhibitory synaptic weight (ns)*/

#define gER 3.0 /*Receiver internal excitatory synaptic weight (ns)*/
#define gIR 16.0 /*Receiver internal inhibitory synaptic weight (ns)*/
// #define gER 3.0 /*Receiver internal excitatory synaptic weight (ns)*/
// #define gIR 16.0 /*Receiver internal inhibitory synaptic weight (ns)*/

#define gSR 4.0 /*Sender to Receiver excitatory synaptic weight (ns)*/
#define gRS 0.15 /*Receiver to Sender excitatory synaptic weight (ns)*/
// #define gSR 0.0 /*Sender to Receiver excitatory synaptic weight (ns)*/
// #define gRS 0.0 /*Receiver to Sender excitatory synaptic weight (ns)*/

// #define gEPoisson 0.0 /*Poisson excitatory synaptic weight (ns)*/
#define gEPoisson 0.6 /*Poisson excitatory synaptic weight (ns)*/
#define totaltime 200000 /*steps of simulation*/
#define transient 40000 /*transient steps*/


float ran2(long *idum);
long seed;


int main(int argc, char *argv[]){


  if(argc != 2){
    cout << "Enter with the following parameters:" << endl;
    cout << "Seed for random number;" << endl;
    cout << "Example: ./exec 318446" << endl;
    return 1;
  }

  seed= (long) atoi(argv[1]);

  FILE *outf;
  FILE *outf2;
  FILE *out_monitor = fopen("./data/original_monitor.dat", "wb");

  char datamempotential[500];
  char dataspikes[500];

  sprintf(datamempotential,"./data/FeedBackFeedForward_MeanMemPot_gES%.3lf_gIS%.3lf_gER%.3lf_gIR%.3lf_gSR%.3lf_gRS%.3lf_seed%ld.mp",gES,gIS,gER,gIR,gSR,gRS,seed);
  sprintf(dataspikes,"./data/SpikeTrain_FeedBackFeedForward_gES%.3lf_gIS%.3lf_gER%.3lf_gIR%.3lf_gSR%.3lf_gRS%.3lf_seed%ld.st",gES,gIS,gER,gIR,gSR,gRS,seed);

	
   /*Warm up*/
   for(int i=0;i<100000;i++) double auxseed=ran2(&seed);

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  /*&&                    CONNECTIVITY MATRIX                   &&*/
  /*&&*/int neigh_send_internal[500][50];
  /*&&*/int neigh_send_external[500][20];
  /*&&*/int neigh_receiv_internal[500][50];
  /*&&*/int neigh_receiv_external[500][20];
  /*&&*/
  /*&&*/for(int j=0;j<500;j++){ /*50 random connec within Sender*/
  /*&&*/  for(int jj=0;jj<50;jj++){
  /*&&*/    neigh_send_internal[j][jj]=(int)(499.0*ran2(&seed));
  /*&&*/    if(neigh_send_internal[j][jj]==j){
  /*&&*/      do neigh_send_internal[j][jj]=(int)(499.0*ran2(&seed));
  /*&&*/      while(neigh_send_internal[j][jj]==j);
  /*&&*/    }
  /*&&*/  }
  /*&&*/}
  /*&&*/for(int j=0;j<500;j++){ /*50 random connec within Receiver*/
  /*&&*/  for(int jj=0;jj<50;jj++){
  /*&&*/    neigh_receiv_internal[j][jj]=(int)(499.0*ran2(&seed));
  /*&&*/    if(neigh_receiv_internal[j][jj]==j){
  /*&&*/      do neigh_receiv_internal[j][jj]=(int)(499.0*ran2(&seed));
  /*&&*/      while(neigh_receiv_internal[j][jj]==j);
  /*&&*/    }
  /*&&*/  }
  /*&&*/}
  /*&&*/for(int j=0;j<500;j++){ /*20 random exc connec from Send to Receiv to each neuron in Receiv*/
  /*&&*/  for(int jj=0;jj<20;jj++){
  /*&&*/    neigh_receiv_external[j][jj]=(int)(399.*ran2(&seed));
  /*&&*/    neigh_send_external[j][jj]=(int)(399.0*ran2(&seed));
  /*&&*/  }
  /*&&*/}
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  /*&&                HETEROGENEITY PARAMETERS                  &&*/
  /*&&*/double a_sender[500];
  /*&&*/double b_sender[500];
  /*&&*/double c_sender[500];
  /*&&*/double d_sender[500];
  /*&&*/double a_receiv[500];
  /*&&*/double b_receiv[500];
  /*&&*/double c_receiv[500];
  /*&&*/double d_receiv[500];
  /*&&*/for(int j=0;j<500;j++){
  /*&&*/  double auxrand1=ran2(&seed);
  /*&&*/  double auxrand2=ran2(&seed);
  /*&&*/  if(j<400){
  /*&&*/    a_sender[j]=0.02;
  /*&&*/    b_sender[j]=0.2;
  /*&&*/    c_sender[j]=-65.0+(15.0*auxrand1*auxrand1);
  /*&&*/    d_sender[j]=8.0-(6.0*auxrand1*auxrand1);
  /*&&*/    a_receiv[j]=0.02;
  /*&&*/    b_receiv[j]=0.2;
  /*&&*/    c_receiv[j]=-65.0+(15.0*auxrand2*auxrand2);
  /*&&*/    d_receiv[j]=8.0-(6.0*auxrand2*auxrand2);
  /*&&*/  }
  /*&&*/  else{
  /*&&*/    a_sender[j]=0.02+(0.08*auxrand1);
  /*&&*/    b_sender[j]=0.25-(0.05*auxrand1);
  /*&&*/    c_sender[j]=-65.0;
  /*&&*/    d_sender[j]=2.0;
  /*&&*/    a_receiv[j]=0.02+(0.08*auxrand2);
  /*&&*/    b_receiv[j]=0.25-(0.05*auxrand2);
  /*&&*/    c_receiv[j]=-65.0;
  /*&&*/    d_receiv[j]=2.0;
  /*&&*/  } 
  /*&&*/} 
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  /*&&*/double v_sender[500];
  /*&&*/double v_receiver[500];
  /*&&*/double u_sender[500];
  /*&&*/double u_receiver[500];
  /*&&*/double rsynE_sender[500];
  /*&&*/double rsynI_sender[500];
  /*&&*/double rsynE_receiver[500];
  /*&&*/double rsynI_receiver[500];
  /*&&*/for(int j=0;j<500;j++){
  /*&&*/  v_sender[j]=v_receiver[j]=-60.0;
  /*&&*/  u_sender[j]=u_receiver[j]=0.2*(-60.0);
  /*&&*/  rsynE_sender[j]=rsynE_receiver[j]=0.001;
  /*&&*/  rsynI_sender[j]=rsynI_receiver[j]=0.001;
  /*&&*/}
  /*&&*/const double poisson_sender=(1.0-exp(-2.4*dt));
  /*&&*/const double poisson_receiver=(1.0-exp(-3.0*dt));
  /*&&*/int spikes_R[500];
  /*&&*/int spikes_S[500];
  /*&&*/
  /*&&*/double IsynExcSender[500];
  /*&&*/double IsynInhSender[500];
  /*&&*/double IsynExcReceiver[500];
  /*&&*/double IsynInhReceiver[500];
  /*&&*/double IsynExcReceiver_ext[500];
  /*&&*/double IsynExcSender_ext[500];
  double Iexc1, Iinh1, Iexc2, Iinh2, rexc1, rinh1, rexc2, rinh2;
  // /*&&*/double *Iexc1=new double[totaltime];
  // /*&&*/double *Iinh1=new double[totaltime];
  // /*&&*/double *Iexc2=new double[totaltime];
  // /*&&*/double *Iinh2=new double[totaltime];
  /*&&*/double *memV_Sender=new double[totaltime]; /*LFP SIGNAL*/
  /*&&*/double *memV_Receiver=new double[totaltime];/*LFP SIGNAL*/
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/




  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  int num_monitor_neurons=200;
  int monitor_neurons[num_monitor_neurons];
  for(int i=0;i<100;i++){
  	monitor_neurons[i]=(int)(399.0*ran2(&seed));
  	for(int j=0;j<i;j++){
  		do monitor_neurons[i]=(int)(399.0*ran2(&seed));
  		while(monitor_neurons[i]==monitor_neurons[j]);
  	}
  }
  for(int i=100;i<num_monitor_neurons;i++){
  	monitor_neurons[i]=(int)(399.0*ran2(&seed));
  	for(int j=0;j<i;j++){
  		do monitor_neurons[i]=(int)(399.0*ran2(&seed));
  		while(monitor_neurons[i]==monitor_neurons[j]);
  	}
  }

  int auxspikes[num_monitor_neurons];
  for(int i=0;i<num_monitor_neurons;i++) auxspikes[i]=0;
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
 
  /**/outf2 = fopen(dataspikes, "a");


  /*HERE STARTS THE ITERATIONS - LOOPING OVER TIME*/
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  for(int t=0;t<totaltime;t++){ /*Looping over time*/


    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  	if(t*dt*0.001>=9 and t*dt*0.001<20){ /*Save 11s of spike trains*/
	  	//*SPIKE TRAINS*//
	    for(int j=0;j<400;j++){
	    	for(int k=0;k<100;k++){
	    		if(j==monitor_neurons[k]){
	    			if(v_sender[j]>30.0) fprintf(outf2, "1\t");
	    			else fprintf(outf2, "0\t");
	    		}
	    	}
	    }   
	    for(int j=0;j<400;j++){
	    	for(int k=100;k<200;k++){
	    		if(j==monitor_neurons[k]){
	    			if(v_receiver[j]>30.0) fprintf(outf2, "1\t");
	    			else fprintf(outf2, "0\t");
	    		}
	    	}
	    }
	    fprintf(outf2, "\n");
    }  
    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/



    /*******************************************/
    /*Check if a neuron fired (=1) or not (=0)*/
    /**/for(int j=0;j<500;j++){
    /**/  /*Sender*/
    /**/  if(v_sender[j]>30.0) spikes_S[j]=1;
    /**/  else spikes_S[j]=0;
    /**/  /*Receiver*/
    /**/  if(v_receiver[j]>30.0) spikes_R[j]=1;
    /**/  else spikes_R[j]=0;
    /**/}
    /*******************************************/

    /*******************************************/
    /*Create the synapse weight,i.e., conductance X #ofinputs*/
    /*SENDER*/
    /**/for(int j=0;j<500;j++){
    /**/  int nsynE=0;
    /**/  int nsynI=0;
    /**/  int nsynE_ext=0;
    /**/  for(int k=0;k<50;k++){
    /**/    if(neigh_send_internal[j][k]<400) nsynE+=spikes_S[neigh_send_internal[j][k]];
    /**/    else nsynI+=spikes_S[neigh_send_internal[j][k]];
    /**/  }
    /**/  IsynExcSender[j]=gES*nsynE;
    /**/  IsynInhSender[j]=gIS*nsynI;
    /**/  for(int k=0;k<20;k++) nsynE_ext+=spikes_R[neigh_send_external[j][k]];
    /**/  IsynExcSender_ext[j]=gRS*nsynE_ext;
    /**/}
    /**//*RECEIVER*/
    /**/for(int j=0;j<500;j++){
    /**/  int nsynE=0;
    /**/  int nsynI=0;
    /**/  int nsynE_ext=0;
    /**/  for(int k=0;k<50;k++){
    /**/    if(neigh_receiv_internal[j][k]<400) nsynE+=spikes_R[neigh_receiv_internal[j][k]];
    /**/    else nsynI+=spikes_R[neigh_receiv_internal[j][k]];
    /**/  }
    /**/  IsynExcReceiver[j]=gER*nsynE;
    /**/  IsynInhReceiver[j]=gIR*nsynI;
    /**/  for(int k=0;k<20;k++) nsynE_ext+=spikes_S[neigh_receiv_external[j][k]];
    /**/  IsynExcReceiver_ext[j]=gSR*nsynE_ext;
    /**/}
    /*******************************************/

    /*******************************************/
    /*Reset V if the neurons fire - Compute the r fraction for synapses - Update the Neuron Variables*/
    /*SENDER*/
    /**/for(int j=0;j<500;j++){
    /**/  int poissonSpike=0;
    /**/  if(poisson_sender>ran2(&seed)) poissonSpike=1;
    /**/  else poissonSpike=0;
    /**/
    /**/  if(v_sender[j]>30.0){/*if fired reset*/
    /**/    v_sender[j]=c_sender[j];
    /**/    u_sender[j]=u_sender[j]+d_sender[j];
    /**/  }
    /**/
    /**/  rsynE_sender[j]=rsynE_sender[j]+dt*(1.0/5.26)*(-1.0*rsynE_sender[j]+IsynExcSender[j]+IsynExcSender_ext[j]+gEPoisson*poissonSpike);
    /**/  rsynI_sender[j]=rsynI_sender[j]+dt*(1.0/5.6)*(-1.0*rsynI_sender[j]+IsynInhSender[j]);
    /**/          
    /**/  double IDC=0.0;
    /**/  if(j<400) IDC=0.;
    /**/ 
    /**/  v_sender[j]=v_sender[j]+dt*(0.04*v_sender[j]*v_sender[j]+5.0*v_sender[j]+140.-u_sender[j]+rsynE_sender[j]*(0.0-v_sender[j])+rsynI_sender[j]*(-65.0-v_sender[j])+IDC);
    /**/  u_sender[j]=u_sender[j]+dt*(a_sender[j]*(b_sender[j]*v_sender[j]-u_sender[j]));

    if (j == 0){
            Iexc1 = rsynE_sender[j]*(0.0-v_sender[j])+rsynI_sender[j]*(-65.0-v_sender[j])+IDC;
            rexc1 = rsynE_sender[j];
            rinh1 = rsynI_sender[j];
          } else if (j == 450){
            Iinh1 = rsynE_sender[j]*(0.0-v_sender[j])+rsynI_sender[j]*(-65.0-v_sender[j])+IDC;
          }
    /**/}
    /**//*RECEIVER*/
    /**/for(int j=0;j<500;j++){
    /**/  int poissonSpike=0;
    /**/  if(poisson_receiver>ran2(&seed)) poissonSpike=1;
    /**/  else poissonSpike=0;
    /**/
    /**/  if(v_receiver[j]>30.0){/*if fired reset*/
    /**/    v_receiver[j]=c_receiv[j];
    /**/    u_receiver[j]=u_receiver[j]+d_receiv[j];
    /**/  }
    /**/
    /**/  rsynE_receiver[j]=rsynE_receiver[j]+dt*(1.0/5.26)*(-1.0*rsynE_receiver[j]+IsynExcReceiver[j]+IsynExcReceiver_ext[j]+gEPoisson*poissonSpike);
    /**/  rsynI_receiver[j]=rsynI_receiver[j]+dt*(1.0/5.6)*(-1.0*rsynI_receiver[j]+IsynInhReceiver[j]);
    /**/
    /**/  double IDC=0.0;
    /**/  if(j<400) IDC=25.;
    /**/
    /**/  v_receiver[j]=v_receiver[j]+dt*(0.04*v_receiver[j]*v_receiver[j]+5.0*v_receiver[j]+140.-u_receiver[j]+rsynE_receiver[j]*(0.0-v_receiver[j])+rsynI_receiver[j]*(-65.0-v_receiver[j])+IDC);
    /**/  u_receiver[j]=u_receiver[j]+dt*(a_receiv[j]*(b_receiv[j]*v_receiver[j]-u_receiver[j]));

          if (j == 0){
            Iexc2 = rsynE_receiver[j]*(0.0-v_receiver[j])+rsynI_receiver[j]*(-65.0-v_receiver[j])+IDC;
            rexc2 = rsynE_receiver[j];
            rinh2 = rsynI_receiver[j];
          } else if (j == 450){
            Iinh2 = rsynE_receiver[j]*(0.0-v_receiver[j])+rsynI_receiver[j]*(-65.0-v_receiver[j])+IDC;
          }
    /**/}

    
    /*******************************************/
    double tmp[10] = {Iexc2, Iinh2, Iexc1, Iinh1, rexc2, rinh2, rexc1, rinh1, v_receiver[0], v_receiver[400]};
// fprintf(outi, "%f,%f,%f,%f,%f,%f,%f,%f,\n",);
fwrite(tmp, sizeof(double), 10, out_monitor);


  if (t%10000 ==0){
    fprintf(stderr, "%6d done\n", t);
  }  
    
 
    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
    /*CREATE A ~LOCAL FIELD POTENTIAL SIGNAL*/
    /**/double auxmemS=0.0;
    /**/double auxmemR=0.0;
    /**/for(int j=0;j<400;j++){
    /**/  auxmemS+=v_sender[j]/400.0;
    /**/  if(j<400) auxmemR+=v_receiver[j]/400.0;
    /**/}
    /**/memV_Sender[t]=auxmemS;
    /**/memV_Receiver[t]=auxmemR;
    /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

  }/*END LOOPING TIME*/
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/

  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
  /**//*Save data*/
  /**/outf = fopen(datamempotential, "w");
  // /**/fprintf(outf, "#This File Contains the Membrane Potential of each Population, Sender e Receiver(only exc) populations.\n");
  // /**/fprintf(outf, "#Time(ms) <V>Sender <V>Receiver\n");
  /**/for(int t=transient;t<totaltime;t++){
  /**/  if((t)%100==0){
    fprintf(outf, "%.3lf \t %lf \t %lf \n",t*dt,memV_Sender[t],memV_Receiver[t]);
  }

  /**/}
  /**/fclose(outf);
  /**/fclose(outf2);
  fclose(out_monitor);
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/



delete [] memV_Sender;
delete [] memV_Receiver;

return 0;
}






/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
/*RANDOM NUMBER - DO NOT EVEN COME HERE*/
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
  int j;
  long k;
  static long idum2 = 123456789, iy = 0, iv[NTAB];
  float temp;

  if(*idum <= 0)
    {
      if(-(*idum) < 1)   *idum = 1;
      else   *idum = -(*idum);
      idum2 = (*idum);
      for(j = NTAB+7;j >= 0;j--)
	{
	  k = (*idum)/IQ1;
	  *idum = IA1*(*idum-k*IQ1)-k*IR1;
	  if(*idum < 0)   *idum += IM1;
	  if(j < NTAB)   iv[j] = *idum;
	}
      iy = iv[0];
    }

  k = (*idum)/IQ1;
  *idum = IA1*(*idum-k*IQ1)-k*IR1;
  if(*idum < 0)   *idum += IM1;

  k = idum2/IQ2;
  idum2 = IA2*(idum2-k*IQ2)-k*IR2;
  if(idum2 < 0)   idum2 += IM2;

  j = iy/NDIV;
  iy = iv[j]-idum2;
  iv[j] = *idum;
  if(iy < 1)   iy += IMM1;

  if((temp = AM*iy) > RNMX)   return RNMX;
  else   return temp;

}
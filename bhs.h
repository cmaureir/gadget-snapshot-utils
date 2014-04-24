#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include<fstream>
#include<string>
#include <string.h>
#include <sstream>
#include <math.h>

//#include "laod_snapshot.h"
//#include "nr3.h"
//#include "sort.h"

#define pi 3.141592654 
#define GRAVITY 6.672e-8
#define BOLTZMANN 1.3806e-16
#define PROTONMASS 1.6726e-24
#define Msun 1.989e33
#define parsec 3.08567758e18
#define year 31556926 

int allocate_memory(void);


using namespace std; 

int BHS(int i, double Udist, double Uvel,double Umasa, double Mgas_init, int ID_BH1, int ID_BH2)
{	
  
  double  Utime=Udist/Uvel;
  double  G=GRAVITY/pow(Udist,3) * Umasa * pow(Utime,2);
  double  Xh= 0.76;  ///* mass fraction of hydrogen */
  double mu= 4.0/(3*Xh+1) * PROTONMASS; // asumiendo gas neutral
  // por particula  mu= 4.0/(3*Xh+1+4*Xh*P[i].Ne) * PROTONMASS;
  double gamma=5.0/3.0;
    //  double temp_to_u_ratio=mu*((gamma-1)/BOLTZMANN)*PROTOMASS
  double Temp;
  double Uenergy= Umasa * pow(Udist,2) / pow(Utime,2);

  stringstream num;		
  num<<i;
  /*
  string NomPDF="PDF_"+num.str();	
    string NomDENSGAS="DENS_GAS"+num.str();	
    string NomDENSSTAR="DENS_STAR"+num.str();	
  */
        string NomTEMPGAS="RAD_TEMP_TIME_GAS"+num.str();	
	string NomRhoCrit="RhoCrit";	
	/*
  string NomHIST="HIST_RHO_"+num.str();	
  string NomPHASE="RHO_TEMP"+num.str();	

  ofstream PDF_OUT(NomPDF.c_str());
  ofstream DENS_GAS(NomDENSGAS.c_str());
  ofstream DENS_STARS(NomDENSSTAR.c_str());
	*/
  ofstream RAD_TEMP_TIME_GAS(NomTEMPGAS.c_str());
/*      ofstream HIST_OUT(NomHIST.c_str());
      ofstream PHASE_OUT(NomPHASE.c_str());
  */
  ofstream RhoCrit_OUT(NomRhoCrit.c_str());
  FILE *DeltaGas_OUT;
  DeltaGas_OUT=fopen("DeltaGas_Time","a");

  FILE *time_a_OUT;
  time_a_OUT=fopen("DistBHs_Time","a");

  FILE *MBH1_Time;
  MBH1_Time=fopen("MBH1_Time","a");

  FILE *MBH2_Time;
  MBH2_Time=fopen("MBH2_Time","a");


  double Udens=Umasa/(Udist*Udist*Udist);
  double MasaU, DensU, RgasU, RstarU, MstarU;
  //      MasaU=P[i].Mass*Umasa;  //  en unidades de gr
  //     DensU=(P[i].Rho*Udens)/PROTONMASS;  // densidad por unidad de masade proton 1/(cm^3)
  //    RgasU=Rgas[i-1]*Udist; // en unidades de cm
  //GENERANDO CURVA DE DENSIDAD MAXIMA RESUELTA EN SPH //////////////////
  


  double RhoCrit,RhoNumCrit;
  for(double logT=1;logT<=6*log(10);logT=logT+10){
    RhoCrit=exp(logT)*exp(logT)*exp(logT)*(3/(4*pi))*(5*BOLTZMANN/(GRAVITY*PROTONMASS*1.22))*(5*BOLTZMANN/(GRAVITY*PROTONMASS*1.22))*(5*BOLTZMANN/(GRAVITY*PROTONMASS*1.22))*((Ngas/120)*(1/MasaU))*((Ngas/120)*(1/MasaU));
    RhoNumCrit=RhoCrit/PROTONMASS;
    RhoCrit_OUT<<RhoNumCrit<<" "<<exp(logT)<<endl;
    //fprintf(RhoCrit_OUT," %e %d \n",RhoCrit,T);
  }
  
  double Radio=0.0;
  Temp=0.0;
  double Mgas=0.0;
  double MBH=0.0;
  double Pos_BH1[3], Pos_BH2[3];   

  
  // cout<<" Vtot= "<<Vtot<<"   Drho= "<<Drho<<endl;	
  for(int i=0;i<NumPart;i++){ 
    if(P[i].Type==0){
      Temp=(P[i].U*Uenergy/Umasa)*(4.0/(3*Xh+1))*PROTONMASS*((gamma-1)/BOLTZMANN);
      Radio=sqrt(P[i].Pos[0]*P[i].Pos[0]+P[i].Pos[1]*P[i].Pos[1]+P[i].Pos[2]*P[i].Pos[2]);
      RAD_TEMP_TIME_GAS<<Radio<<" "<<Temp<<" "<<header1.time*(Utime/year)<<endl;
      Mgas=Mgas+P[i].Mass;
    }
    
    if(P[i].Type==5 ){  
      MBH=P[i].Mass*Umasa/Msun;
	
      if(Id[i]==ID_BH1){ 
           fprintf(MBH1_Time,"%e %e \n",header1.time*(Utime/year),MBH);
	   Pos_BH1[0]=P[i].Pos[0];
	   Pos_BH1[1]=P[i].Pos[0];
	   Pos_BH1[2]=P[i].Pos[0];

	}
      if(Id[i]==ID_BH2){
	   fprintf(MBH2_Time,"%e %e \n",header1.time*(Utime/year),MBH);
	   Pos_BH2[0]=P[i].Pos[0];
	   Pos_BH2[1]=P[i].Pos[1];
	   Pos_BH2[2]=P[i].Pos[2];

	}
    }
    
  }
  
  
  double r2_BHs=(Udist*Udist/(parsec*parsec))*(Pos_BH1[0]-Pos_BH2[0])*(Pos_BH1[0]-Pos_BH2[0])+(Pos_BH1[1]-Pos_BH2[1])*(Pos_BH1[0]-Pos_BH2[1])+(Pos_BH1[0]-Pos_BH2[2])*(Pos_BH1[0]-Pos_BH2[2]);
  fprintf(time_a_OUT,"%e %e \n",header1.time*(Utime/year),sqrt(r2_BHs)); // distancia en parsecs

  
  double   DeltaGas=Mgas/Mgas_init;
  fprintf(DeltaGas_OUT,"%e %e \n",header1.time*(Utime/year),DeltaGas);
  
  
  RAD_TEMP_TIME_GAS.close();
  fclose(DeltaGas_OUT);
  fclose(time_a_OUT);
  fclose(MBH1_Time);
  fclose(MBH2_Time);
  
  
}


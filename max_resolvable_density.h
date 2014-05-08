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
#define THOMPSON     6.65245e-25
#define C           2.9979e10

int allocate_memory(void);


using namespace std; 

int max_resolvable_density(int i, int j, double Udist, double Uvel,double Umasa)
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

  if(i==j) system("mkdir BHs");	

  stringstream num;		
  num<<i;
        string NomTEMPGAS="RAD_TEMP_TIME_GAS"+num.str();	
	string NomRhoCrit="RhoCrit";	

  ofstream RAD_TEMP_TIME_GAS(NomTEMPGAS.c_str());

  ofstream RhoCrit_OUT(NomRhoCrit.c_str());

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
  double Pos_BH1[3], Pos_BH2[3], Vel_BH1[1], Vel_BH2[3];   

  
  // cout<<" Vtot= "<<Vtot<<"   Drho= "<<Drho<<endl;	
  for(int i=1;i<=NumPart;i++){ 
//    cout<<"NumPat "<<NumPart<<endl; 
    if(P[i].Type==0){
//	cout<<"Type 0 index --> "<<i<<"   NumPart-index --> "<<NumPart-i<<endl; 

     Temp=(P[i].U*Uenergy/Umasa)*(4.0/(3*Xh+1))*PROTONMASS*((gamma-1)/BOLTZMANN);
      Radio=sqrt(P[i].Pos[0]*P[i].Pos[0]+P[i].Pos[1]*P[i].Pos[1]+P[i].Pos[2]*P[i].Pos[2]);
      RAD_TEMP_TIME_GAS<<Radio<<" "<<Temp<<" "<<header1.time*(Utime/year)<<endl;
      Mgas=Mgas+P[i].Mass;
    }	    
  }

  RAD_TEMP_TIME_GAS.close();
}




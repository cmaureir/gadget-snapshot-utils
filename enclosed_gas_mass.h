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

int enclosed_gas_mass(int i, double Utime, double Udist, double Umasa)
{

 //////////////// SORTING PARTICLES BY RADIUS//////////////////////////////////
  //cout<<"Ordenando R"<<endl;
  double Rmax0;      
  NRvector<double> R(NumPart);
  NRvector<double> auxR(NumPart);
  for(long n=1; n<=NumPart; n++){
    R[n-1]=sqrt(P[n].Pos[0]*P[n].Pos[0]+P[n].Pos[1]*P[n].Pos[1]+P[n].Pos[2]*P[n].Pos[2]);
  } 
  Indexx indiceR(R);  
  indiceR.sort(auxR);
  //double Rmax=Rcil[indiceR.indx[NumPart-1]];
  //cout<<"Rmax = "<<Rmax<<endl;
  //if(i==0) Rmax0=Rmax;

 /////////// CALCULING ENCLOSED GAS MASS ///////////////////////////////////////


   FILE *enclosed_gas_mass; // 
   enclosed_gas_mass=fopen("MgasR100pc_MgasR500pc_MgasR1Kpc_Time","a");


  cout<<"Calculing enclosed mass..."<<endl;
 
  double MgasR100pc, MgasR500pc, MgasR1Kpc,Mgas;
  Mgas=MgasR100pc=MgasR500pc=MgasR1Kpc=0.0;
	
  float r_in_parsec;
  for(long i=1;i<=NumPart;i++){  
    if(P[indiceR.indx[i]+1].Type==0){
	Mgas=Mgas+P[indiceR.indx[i]+1].Mass;
       r_in_parsec=R[indiceR.indx[i]]*Udist/parsec;
       if(r_in_parsec>99 && r_in_parsec<101) MgasR100pc=Mgas; 	
       if(r_in_parsec>499 && r_in_parsec<501) MgasR500pc=Mgas; 	
       if(r_in_parsec>999 && r_in_parsec<1001){
	  MgasR1Kpc=Mgas;	
	  fprintf(enclosed_gas_mass, "%e %e %e %e \n",MgasR100pc, MgasR500pc, MgasR1Kpc,header1.time*Utime/year );
          break;
	}  	

    }
}
    


}


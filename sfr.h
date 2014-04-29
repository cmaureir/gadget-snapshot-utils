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


int write_Mstar_dens(int i,double U_Mass, double Utime)
{
  system("mkdir SFR");
  //cout<<"Generando Archivos de Salida...."<<endl;
  stringstream num;
  num<<i;
  stringstream Periodo;
  Periodo<<i;
  
  string NomDens="SFR/Dens_"+Periodo.str();
  ofstream Dens(NomDens.c_str());

  FILE *MstarsT; // 
  MstarsT=fopen("SFR/M_Stars_Time","a");

/////////// DESITY AND TOTAL STAR MASS ////////////////////////////////////////////
  cout<<"Extrayendo densidad..."<<endl;
  double Mstars=0.0;
  float r;
  for(long i=1;i<=NumPart;i++){
    if(P[i].Type==4){
    Mstars=Mstars+P[i].Mass;
    }
    if(P[i].Type==0){
      Dens<<P[i].Rho<<endl;
    }
  }
  double Time_in_years=header1.time*Utime/year;
  double Mstar_in_solar_mass=Mstars*U_Mass/Msun;
  fprintf(MstarsT, "%e %e \n",Time_in_years,Mstar_in_solar_mass);

  Dens.close(); 
  fclose(MstarsT);


}


int SFR(){

   ifstream SFR_out("SFR/M_Stars_Time");
   

   system("wc SFR/M_Stars_Time > SFR/wcSFR");
   system("awk '{print $1}' SFR/wcSFR > SFR/lineas");
   system("rm SFR/wcSFR");
        
    
   ifstream lineas("SFR/lineas", ios::in);
            
   int Nlineas=0;
   lineas >> Nlineas;
                    

   NRvector<double> Mstar(Nlineas);
   NRvector<double> Time(Nlineas);


   int i=0;
   while(i<Nlineas){
     SFR_out >> Time[i] >>Mstar[i]; 
     if(i==1) cout<<Time[i]<<" "<<Mstar[i]<<endl;
     if(i==Nlineas-1) cout<<Time[i]<<" "<<Mstar[i]<<endl;
     i++;
     //cout<<"Reading M_Stars_Time line"<<i<<endl;
   }

   string NomdSfr="SFR/dSfr";
   ofstream dSfr_out(NomdSfr.c_str());

   double dm;
   double dt;
   for(int k=1;k<Nlineas-1;k++){
      cout<<"k "<<k<<endl;
      dm=Mstar[k+1]-Mstar[k-1];
      dt=Time[k+1]-Time[k-1];
      dSfr_out<<Time[k]<<" "<<dm/dt<<endl;
   }

   dSfr_out.close();
   SFR_out.close();
}

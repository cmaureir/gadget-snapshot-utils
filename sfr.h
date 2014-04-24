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


int write_Mstar_dens(int i,double U_Mass)
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

//////////////// ORDENANDO POR R CILINDRICO//////////////////////////////////
//cout<<"Ordenando R"<<endl;
  double Rmax0;
  NRvector<double> Rcil(NumPart);
  NRvector<double> auxRcil(NumPart);

  for(long n=1; n<=NumPart; n++){
    Rcil[n-1]=sqrt(P[n].Pos[0]*P[n].Pos[0]+P[n].Pos[1]*P[n].Pos[1]);
  }

  Indexx indiceR(Rcil);
  indiceR.sort(auxRcil);
  double Rmax=Rcil[indiceR.indx[NumPart-1]];
  //cout<<"Rmax = "<<Rmax<<endl;
  if(i==0) Rmax0=Rmax;

/////////// DENSIDAD Y MASA SOLAR/////////////////////////////////////////////
  cout<<"Extrayendo densidad..."<<endl;
  double Mstars=0.0;
  float r;
  for(long i=1;i<=NumPart;i++){
    if(P[i].Type==4){
    Mstars=Mstars+P[i].Mass*U_Mass;
    }
    if(P[i].Type==0){
      Dens<<P[i].Rho<<endl;
    }
  }

  fprintf(MstarsT, "%f %f \n",header1.time,Mstars );

  Dens.close(); 
  fclose(MstarsT);


}


int SFR(){
   ifstream SFR_out("SFR/M_Stars_Time");
   int Nlineas=0;
   while(!SFR_out.eof()){
    Nlineas++;
    cout<<Nlineas<<endl;
   }
   cout<<Nlineas<<endl;
   NRvector<double> Mstar(Nlineas);
   NRvector<double> Time(Nlineas);


   int i=0;
   while(! SFR_out.eof()){
     SFR_out >> Mstar[i] >>Time[i]; 
     //cout<<Mstar[i]<<" "<<Time[i]<endl;
     i++;
     cout<<"i "<<i<<endl;
   }

   string NomdSfr="SFR/dSfr";
   ofstream dSfr_out(NomdSfr.c_str());

   double dm;
   double dt;
   for(int k=1;k<Nlineas/2.0;k++){
      cout<<"k "<<k<<endl;
      dm=Mstar[k+1]-Mstar[k-1];
      dt=Time[k+1]-Time[k-1];
      dSfr_out<<dm/dt<<" "<<Time[k]<<endl;
   }

   dSfr_out.close();
   SFR_out.close();
}

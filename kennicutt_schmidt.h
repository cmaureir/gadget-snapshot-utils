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

int data_for_keniccutt_schmidt(int i, int DDENS)
{       
  //cout<<"Generando Archivos de Salida...."<<endl;
  stringstream num;
  num<<i;

  
  printf("El tiempo en que se capturo este snapshot es %f \n", Time);
  printf("El numero total de particulas es %i \n", NumPart);
  printf("El numero de particulas de gas es %i \n", header1.npart[0]);
//  printf("La masa de los BHs es %f \n \n", header1.mass[5]);



  //////////////// ORDENANDO PARTICULAS POR R CILINDRICO//////////////////////////////////
  cout<<"Ordenando R"<<endl;
  NRvector<double> Rcil(NumPart);
  NRvector<double> auxRcil(NumPart);

  double Mgas, Mstars;
  Mgas=Mstars=0.0;
    
  cout<<" Ngas = "<<Ngas<<endl;
  for(long n=1; n<=NumPart; n++){
    if(P[n].Type==0) Mgas=Mgas+P[n].Mass;
    if(P[n].Type==4) Mstars=Mstars+P[n].Mass;
    Rcil[n-1]=sqrt(P[n].Pos[0]*P[n].Pos[0]+P[n].Pos[1]*P[n].Pos[1]);
  }
  cout<<"Rcil"<<endl;
 
  Indexx indiceR(Rcil);
  indiceR.sort(auxRcil);
 double Rmax=Rcil[indiceR.indx[NumPart-1]];

  double   Mgas_Rstar,Mstar_Rstar, Mgas_Rgas,Mstar_Rgas;
  ///////////// MASA ENCERRADA ///////////////////////////////////
  //cout<<"Masa Encerrada"<<endl;  

  //  system("rm snaps_cantidades");
     
  FILE *Snaps;
  Snaps=fopen("snaps_cantidades","a");


  double R_Medi, RStar90, R05, R, RGas90,Mgas_R05, Mstar_R05, Mgas_R08, Mstar_R08;
  long double MassGas=0.0;
  long double MassStar=0.0;


  for(long k=0;k<NumPart;k++){
    if(P[indiceR.indx[k]+1].Type==0){
      MassGas=MassGas+P[indiceR.indx[k]+1].Mass;
    }
    if(P[indiceR.indx[k]+1].Type==4){
      MassStar=MassStar+P[indiceR.indx[k]+1].Mass;
    }

    if(MassStar>0.94*Mstars &&MassStar<0.96*Mstars){
      RStar90=sqrt(P[indiceR.indx[k]+1].Pos[0]*P[indiceR.indx[k]+1].Pos[0]+P[indiceR.indx[k]+1].Pos[1]*P[indiceR.indx[k]+1].Pos[1]);
      Mgas_Rstar=MassGas;
      Mstar_Rstar=MassStar;
    }
  if(MassGas>0.94*Mgas &&MassGas < 0.96*Mgas){
      RGas90=sqrt(P[indiceR.indx[k]+1].Pos[0]*P[indiceR.indx[k]+1].Pos[0]+P[indiceR.indx[k]+1].Pos[1]*P[indiceR.indx[k]+1].Pos[1]);
      Mgas_Rgas=MassGas;
      Mstar_Rgas=MassStar;
    }

    R=sqrt(P[indiceR.indx[k]+1].Pos[0]*P[indiceR.indx[k]+1].Pos[0]+P[indiceR.indx[k]+1].Pos[1]*P[indiceR.indx[k]+1].Pos[1]);
    if(R<Rmax*0.52 && R>Rmax*0.48){
      Mgas_R05=MassGas;
      Mstar_R05=MassStar;
    }

    if(R<Rmax*0.82 && R>Rmax*0.78){
      Mgas_R08=MassGas;
      Mstar_R08=MassStar;
    }

  }

  //  double SigmaGas=MassEnc/(2*pi*R90*R90);
  fprintf(Snaps, "%f %f %f %f %f %f %f %f %f %f %f %f\n",Time, RStar90, Mstar_Rstar, Mgas_Rstar, RGas90, Mstar_Rgas, Mgas_Rgas, Rmax, Mgas_R05, Mstar_R05, Mgas_R08, Mstar_R08);

  /* 
  Dopen.close();
  Mass.close();
  Hsml.close();
  Hsmlmin.close();
  */ 
  fclose(Snaps);

}


void MeanFilter(){

  system("wc sfr.txt > wcSFR");
  system("awk '{print $1}' wcSFR > lineas");
  system("rm wcSFR");

  ifstream fileSFR("sfr.txt", ios::in);
  ifstream lineas("lineas", ios::in);
   
  FILE *SFR_MEAN;
  SFR_MEAN=fopen("sfr_mean","a");

  int L;

  lineas >> L;
   
  int Ndat=4;
  NRvector<double> Dat(Ndat+1);
  NRvector<double> auxDat(Ndat+1);
  NRvector<double> SFR(L);
  NRvector<double> SFR2(L);   
  NRvector<double> T(L);
  double a,b,c,sfr,t2;
  
  for(int i=0;i<L;i++){
    fileSFR >> T[i] >> a >> b >> c >> SFR[i];
  }  
    
  for(int i=0;i<L;i++){
    SFR2[i]=SFR[i];
  }
 /////////// PROMEDIO ////////////////////////
  for(int i=(1+(Ndat/2));i<L-(1+(Ndat/2));i++){
    int nn=1+Ndat;
    for(int n=1;n<=Ndat/2;n++){
      if(SFR[i-n]==0) nn--;
      if(SFR[i+n]==0) nn--;
      SFR2[i]=SFR2[i]+SFR[i+n]+SFR[i-n];
    }    
    SFR2[i]=SFR2[i]/(nn);
  }
  /////////////////////////////////////////////

  /////// MEDIA //////////////////////////////
  int Ndat2=4;
  for(int i=0;i<L;i++){
    SFR[i]=SFR2[i];
  }

  for(int i=(1+(Ndat2/2));i<L-(1+(Ndat2/2));i++){
    Dat[0]=SFR[i];
    for(int n=1;n<=Ndat2/2;n++){
      Dat[n+Ndat2/2]=SFR[i+n];
      Dat[n]=SFR[i-n];
    }
    Indexx indiceDat(Dat);
    indiceDat.sort(auxDat);
    SFR2[i]=Dat[indiceDat.indx[1+Ndat2/2]];
  } 
    
  ///////////////////////////////////////////
  for(int i=0;i<L;i++){
    fprintf(SFR_MEAN,"%f %f \n", T[i],SFR2[i]);
  }

  fclose(SFR_MEAN);

}


void Kennicut_Schmidt(double Udist,double Uvel,double Umasa, int Nfin, int Nin, int paso){

   

  system("rm t_dt.txt");
  //  system("rm snaps_cantidades");

  MeanFilter();
  int L=(Nfin-Nin)/paso;
   
  FILE *Kennicutt;
  Kennicutt=fopen("Kennicutt_Schmidt","a");


  FILE *KennU;
  KennU=fopen("Kennicutt_Schmidt_U","a");
  FILE *KennU_Cal;
  KennU_Cal=fopen("Kennicutt_Schmidt_U_Cal","a");

  FILE *Kennicutt_medio;
  Kennicutt_medio=fopen("Kennicutt_Schmidt_medio","a");
  FILE *Kennicutt_medio_U_Cal;
  Kennicutt_medio_U_Cal=fopen("Kennicutt_Schmidt_medio_U_Cal","a");
   
system("sed -n -e '/Step/p' timings.txt > prueba");
  system("awk '{print $4, $6}' prueba > prueba2"); //  T, dT
  system("sed -e '1,2d'  prueba2 > t_dt.txt");
  system("rm prueba prueba2");

  ifstream fileTime("t_dt.txt", ios::in);
  ifstream file("snaps_cantidades", ios::in);
  ifstream fileSFR("sfr_mean", ios::in);


  NRvector<double> SFR(L);


  NRvector<double> T(L);
  NRvector<double> Rstar(L);
  NRvector<double> Mstar_Rstar(L);
  NRvector<double> Mgas_Rstar(L);
  NRvector<double> Rgas(L);
  NRvector<double> Mstar_Rgas(L);
  NRvector<double> Mgas_Rgas(L);
  NRvector<double> Mstar_R05(L);
  NRvector<double> Mgas_R05(L);
  NRvector<double> Mstar_R08(L);
  NRvector<double> Mgas_R08(L);
  NRvector<double> Rmax(L);


  NRvector<double> DT(L);
  double t, t2, dt, sfr;
for(int i=0;i<L;i++){
    file >> T[i]>> Rstar[i] >> Mstar_Rstar[i] >> Mgas_Rstar[i] >> Rgas[i] >> Mstar_Rgas[i] >> Mgas_Rgas[i] >> Rmax[i] >> Mgas_R05[i] >> Mstar_R05[i] >> Mgas_R08[i] >> Mstar_R08[i];
  }

  int k=0;
  cout<<" L = "<<L<<endl;
  do{
    fileTime >> t >> dt;
       cout<<t<<" "<<T[k]<<endl;        
    if(t-dt<=T[k] && t+dt>=T[k]){
      DT[k]=dt;
      cout<<t<<" "<<T[k]<<" "<<dt<<" k="<<k<<endl;
      k++;

    }
  }while(k<L);

  int h=0;
  do{
    fileSFR >> t2 >> sfr ; 
    if(t2-DT[h]<=T[h] && t2+DT[h]>=T[h]){
      SFR[h]=sfr;
      cout<<t2<<" "<<T[h]<<" "<<SFR[h]<<" h="<<h<<endl;
      h++;     
    } 
  }while(h<L);

  cout<<"FIn Comparacion"<<endl;
  double TU, DTU, RstarU, RgasU, Mstar_RstarU, Mstar_RgasU,Mgas_RstarU,Mgas_RgasU,SFRU, Mgas_R05U, Mstar_R05U, Mgas_R08U, Mstar_R08U, RmaxU,D_Mstar_RstarU, D_Mstar_RgasU,D_Mstar_R05U,D_Mstar_R08U;

  NRvector<double> T_medio(L-1);
  NRvector<double> DT_medio(L-1);

  NRvector<double> MstarRstar_medio(L-1);
  NRvector<double> MstarRgas_medio(L-1);
  NRvector<double> MstarR05_medio(L-1);
  NRvector<double> MstarR08_medio(L-1);
    
  NRvector<double> D_MstarRstar_medio(L-1);
  NRvector<double> D_MstarRgas_medio(L-1);
  NRvector<double> D_MstarR05_medio(L-1);
  NRvector<double> D_MstarR08_medio(L-1);

  NRvector<double> MgasRstar_medio(L-1);
  NRvector<double> MgasRgas_medio(L-1);   
  NRvector<double> MgasR05_medio(L-1);
  NRvector<double> MgasR08_medio(L-1);
  NRvector<double> Dens_medio(L-1);   

  NRvector<double> Rstar_medio(L-1);
  NRvector<double> Rgas_medio(L-1);
  NRvector<double> R05_medio(L-1);
  NRvector<double> R08_medio(L-1);
  NRvector<double> Rmax_medio(L-1);

 for(int n=0;n<L-1;n++){
    T_medio[n]=0.5*(T[n+1]+T[n]);
    DT_medio[n]=(T[n+1]-T[n]);
    cout<<n<<endl;
    //    Mstar_medio[n]=

    D_MstarRstar_medio[n]=(Mstar_Rstar[n+1]-Mstar_Rstar[n]);
    D_MstarRgas_medio[n]=(Mstar_Rgas[n+1]-Mstar_Rgas[n]);
    D_MstarR05_medio[n]=(Mstar_R05[n+1]-Mstar_R05[n]);
    D_MstarR08_medio[n]=(Mstar_R08[n+1]-Mstar_R08[n]);

    MstarRstar_medio[n]=0.5*(Mstar_Rstar[n+1]+Mstar_Rstar[n]);
    MstarRgas_medio[n]=0.5*(Mstar_Rgas[n+1]+Mstar_Rgas[n]);
    MstarR05_medio[n]=0.5*(Mstar_R05[n+1]+Mstar_R05[n]);
    MstarR08_medio[n]=0.5*(Mstar_R08[n+1]+Mstar_R08[n]);
    
    MgasRstar_medio[n]=0.5*(Mgas_Rstar[n+1]+Mgas_Rstar[n]);
    MgasRgas_medio[n]=0.5*(Mgas_Rgas[n+1]+Mgas_Rgas[n]);
    MgasR05_medio[n]=0.5*(Mgas_R05[n+1]+Mgas_R05[n]);
    MgasR08_medio[n]=0.5*(Mgas_R08[n+1]+Mgas_R08[n]);

    Rstar_medio[n]=0.5*(Rstar[n+1]+Rstar[n]);
    Rgas_medio[n]=0.5*(Rgas[n+1]+Rgas[n]);
    Rmax_medio[n]=0.5*(Rmax[n+1]+Rmax[n]);

    cout<<n<<"CantidadesMEdias"<<endl;

    //fprintf(Kennicutt_medio, "%e %e %e %e %e %e %e %e %e %e %e %e %e \n",T_medio[n], DT_medio[n],Rstar_medio[n],Rgas_medio[n],Rmax_medio[n],MstarRstar_medio[n],MgasRstar_medio[n],MstarRgas_medio[n],MgasRg
    
    cout<<n<<"Prit1"<<endl; 

    TU=(T_medio[n]*(Udist/Uvel))/31536000; // Yrs
    DTU=(DT_medio[n]*(Udist/Uvel))/31536000;

    
    RstarU=Rstar_medio[n]*Udist/3.08567758e18;           // En unidades de Parsec    
    RgasU=Rgas_medio[n]*Udist/3.08567758e18;           // En unidades de Parsec    

    // SFRU=Mstar_medio[n]*Umasa/1.91891e33;              // Masa de Estrellas formadas en T_medio[n] en unidad de Msolares    
    
  RstarU=Rstar_medio[n]*Udist/3.08567758e18;           // En unidades de Parsec    
    RgasU=Rgas_medio[n]*Udist/3.08567758e18;           // En unidades de Parsec    

    RmaxU=Rmax_medio[n]*Udist/3.08567758e18;           // En unidades de Parsec    
    

    D_Mstar_RstarU=D_MstarRstar_medio[n]*Umasa/1.91891e33; // Masa en Unidad de Msolares
    D_Mstar_RgasU=D_MstarRgas_medio[n]*Umasa/1.91891e33;
    D_Mstar_R05U=D_MstarR05_medio[n]*Umasa/1.91891e33;
    D_Mstar_R08U=D_MstarR08_medio[n]*Umasa/1.91891e33;  

    Mstar_RstarU=MstarRstar_medio[n]*Umasa/1.91891e33; // Masa en Unidad de Msolares
    Mstar_RgasU=MstarRgas_medio[n]*Umasa/1.91891e33;
    Mstar_R05U=MstarR05_medio[n]*Umasa/1.91891e33;
    Mstar_R08U=MstarR08_medio[n]*Umasa/1.91891e33;
    
    Mgas_RstarU=MgasRstar_medio[n]*Umasa/1.91891e33;
    Mgas_RgasU=MgasRgas_medio[n]*Umasa/1.91891e33;
    Mgas_R05U=MgasR05_medio[n]*Umasa/1.91891e33;
    Mgas_R08U=MgasR08_medio[n]*Umasa/1.91891e33;
    
    cout<<n<<"Unidades"<<endl;
    //    double SFRU_DTU_TOT=(D_Mstar/DTU);
    double SFRU_DTU_Rgas=(D_Mstar_RgasU/DTU);
    double SFRU_DTU_Rstar=(D_Mstar_RstarU/DTU);
    double SFRU_DTU_R05=(D_Mstar_R05U/DTU);
    double SFRU_DTU_R08=(D_Mstar_R08U/DTU);

    double SigmaGasRgasB=Mgas_RgasU/(pi*RgasU*RgasU);
    double SigmaGasRstarB=Mgas_RstarU/(pi*RstarU*RstarU);
    double SigmaGasR05B=Mgas_R05U/(pi*(RmaxU*0,5)*(RmaxU*0,5));
    double SigmaGasR08B=Mgas_R08U/(pi*(RmaxU*0,8)*(RmaxU*0,8));  

    double SigmaRgasB=1000000*SFRU_DTU_Rgas/(pi*RgasU*RgasU);
    double SigmaRstarB=1000000*SFRU_DTU_Rstar/(pi*RstarU*RstarU);
    double SigmaR05B=1000000*SFRU_DTU_R05/(pi*(RmaxU*0,5)*(RmaxU*0,5));
    double SigmaR08B=1000000*SFRU_DTU_R08/(pi*(RmaxU*0,8)*(RmaxU*0,8));

    cout<<n<<"Sigmas"<<endl; 
    fprintf(Kennicutt_medio_U_Cal, "%e %e %e %e %e %e %e %e \n",SigmaGasRgasB,SigmaGasRstarB,SigmaGasR05B,SigmaGasR08B,SigmaRgasB,SigmaRstarB,SigmaR05B,SigmaR08B);
    cout<<n<<"fin"<<endl; 
  }




  for(int n=0;n<L;n++){
    cout<<"Ultimo for 1,  n= "<<n<<endl;
    cout<<"T ="<< T[n]<<endl;
    cout<<"DT ="<< DT[n]<<endl;
    cout<<"SFR ="<< SFR[n]<<endl;
    cout<<"Rstar ="<< Rstar[n]<<endl;
    cout<<"Rgas ="<< Rgas[n]<<endl;
    cout<<"Mgas_Rstar ="<< Mgas_Rstar[n]<<endl;
    cout<<"Mstar_Rstar ="<< Mstar_Rstar[n]<<endl;
    cout<<"Mgas_Rgas ="<< Mgas_Rgas[n]<<endl;  
    cout<<"Mstar_Rgas ="<< Mstar_Rgas[n]<<endl;
    
    fprintf(Kennicutt, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f \n",T[n],DT[n],SFR[n],Rstar[n],Mstar_Rstar[n],Mgas_Rstar[n],Rgas[n],Mstar_Rgas[n],Mgas_Rgas[n],Rmax[n],Mgas_R05[n],Mstar_R05[n],Mgas_R08[n],Mstar_R08[n]);
    cout<<"Ultimo for 2,  n= "<<n<<endl;   
    TU=(T[n]*(Udist/Uvel))/31536000; // Yrs
    DTU=(DT[n]*(Udist/Uvel))/31536000;
    SFRU=SFR[n]*Umasa/1.91891e33;    // Masa de Estrellas formadas en T[n] en unidad de Msolares    
    RstarU=Rstar[n]*Udist/3.08567758e18;           // En unidades de Parsec    
    RgasU=Rgas[n]*Udist/3.08567758e18;           // En unidades de Parsec    

    RmaxU=Rmax[n]*Udist/3.08567758e18;           // En unidades de Parsec    
    
    Mstar_RstarU=Mstar_Rstar[n]*Umasa/1.91891e33; // Masa en Unidad de Msolares
    Mstar_RgasU=Mstar_Rgas[n]*Umasa/1.91891e33;
    Mstar_R05U=Mstar_R05[n]*Umasa/1.91891e33;
    Mstar_R08U=Mstar_R08[n]*Umasa/1.91891e33;
    
    Mgas_RstarU=Mgas_Rstar[n]*Umasa/1.91891e33;
    Mgas_RgasU=Mgas_Rgas[n]*Umasa/1.91891e33;
    Mgas_R05U=Mgas_R05[n]*Umasa/1.91891e33;
    Mgas_R08U=Mgas_R08[n]*Umasa/1.91891e33;

    cout<<"Ultimo for 3,  n= "<<n<<endl;

    cout<<"TU ="<< TU<<endl; 
    cout<<"DTU ="<< DTU<<endl; 
    cout<<"SFRU ="<< SFRU<<endl; 
    cout<<"RstarU ="<< RstarU<<endl; 
    cout<<"RgasU ="<< RgasU<<endl; 
    cout<<"Mgas_RstarU ="<< Mgas_RstarU<<endl; 
    cout<<"Mstar_RstarU ="<< Mstar_RstarU<<endl; 
    cout<<"Mgas_RgasU ="<< Mgas_RgasU<<endl;   
    cout<<"Mstar_RgasU ="<< Mstar_RgasU<<endl; 

    double SFRU_DTU=(SFRU/DTU);
    double SigmaRstar=Mgas_RstarU/(pi*RstarU*RstarU);
    double SigmaRgas=Mgas_RgasU/(pi*RgasU*RgasU);

    double SigmaR05=Mgas_R05U/(pi*(RmaxU*0,5)*(RmaxU*0,5));
    double SigmaSFRU_R05=SFRU_DTU/(pi*(RmaxU*0,5)*(RmaxU*0,5));  
      
    double SigmaSFRURstar=SFRU_DTU/(pi*RstarU*RstarU);
    double SigmaSFRURgas=SFRU_DTU/(pi*RgasU*RgasU);
    
    fprintf(KennU, "%e %e %e %e %e %e %e %e %e \n", TU,DTU,SFRU_DTU, RstarU,Mstar_RstarU,Mgas_RstarU, RgasU,Mstar_RgasU,Mgas_RgasU);
    cout<<"Ultimo for 4,  n= "<<n<<endl;
    fprintf(KennU_Cal,"%e %e %e %e %e %e %e %e \n",TU,SFRU_DTU,SigmaRstar,SigmaRgas,SigmaR05,SigmaSFRURstar,SigmaSFRURgas,SigmaSFRU_R05);
    cout<<"Ultimo for,  n= "<<n<<endl; 
  } 

  fclose(Kennicutt);
  fclose(Kennicutt_medio);
  fclose(Kennicutt_medio_U_Cal);
  fclose(KennU);
  fclose(KennU_Cal);
}

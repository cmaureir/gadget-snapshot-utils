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

int BHS(int i, int j, double Udist, double Uvel,double Umasa, double Mgas_init, int ID_BH1, int ID_BH2, int BHS_TYPE)
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
  DeltaGas_OUT=fopen("Time_DeltaGas","a");

  FILE *EnclosedGas_OUT;
  EnclosedGas_OUT=fopen("BHs/Time_GasInsideRBH","a");

  FILE *time_a_OUT;
  time_a_OUT=fopen("BHs/Time_DistBHs","a");

  FILE *MBH1_Time;
  MBH1_Time=fopen("BHs/Time_MBH1","a");

  FILE *MBH2_Time;
  MBH2_Time=fopen("BHs/Time_MBH2","a");


  FILE *VelBH1_Time;
  VelBH1_Time=fopen("BHs/Time_VelBH1","a");

  FILE *VelBH2_Time;
  VelBH2_Time=fopen("BHs/Time_VelBH2","a");

  FILE *PosBH1_Time;
  PosBH1_Time=fopen("BHs/Time_PosBH1","a");

  FILE *PosBH2_Time;
  PosBH2_Time=fopen("BHs/Time_PosBH2","a");


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
    if(P[i].Type==BHS_TYPE ){  
       cout<<"NumPat "<<NumPart<<endl; 

	cout<<"Type BHs index --> "<<i<<"    NumPart-index --> "<<NumPart-i<<endl; 

      MBH=P[i].Mass*Umasa/Msun;
      if(Id[i]==ID_BH1){ 
	cout<<"Firts BH"<<endl;
           fprintf(MBH1_Time,"%e %e \n",header1.time*(Utime/year),MBH);
	cout<<"end of print Mass Firts BH"<<endl;

	   Pos_BH1[0]=P[i].Pos[0];
	   Pos_BH1[1]=P[i].Pos[1];
	   Pos_BH1[2]=P[i].Pos[2];
	cout<<"end of reading positions Firts BH"<<endl;

	   Vel_BH1[0]=P[i].Vel[0];
	   Vel_BH1[1]=P[i].Vel[1];
	   Vel_BH1[2]=P[i].Vel[2];
	cout<<"end of reading velocity Firts BH"<<endl;

	   fprintf(VelBH1_Time,"%e %e %e %e \n",header1.time*(Utime/year),Vel_BH1[0],Vel_BH1[1],Vel_BH1[2]);
	   fprintf(PosBH1_Time,"%e %e %e %e \n",header1.time*(Utime/year),Pos_BH1[0],Pos_BH1[1],Pos_BH1[2]);

	cout<<"end of writing velocity and positions first BH"<<endl;

	}
      if(Id[i]==ID_BH2){
	cout<<"Second BH"<<endl;

	   fprintf(MBH2_Time,"%e %e \n",header1.time*(Utime/year),MBH);
	cout<<"end of print Mass Second BH"<<endl;

	   Pos_BH2[0]=P[i].Pos[0];
	   Pos_BH2[1]=P[i].Pos[1];
	   Pos_BH2[2]=P[i].Pos[2];
	cout<<"end of reading position Second BH"<<endl;


	   Vel_BH2[0]=P[i].Vel[0];
	   Vel_BH2[1]=P[i].Vel[1];
	   Vel_BH2[2]=P[i].Vel[2];

	cout<<"end of reading velocity Second BH"<<endl;

	   fprintf(VelBH2_Time,"%e %e %e %e \n",header1.time*(Utime/year),Vel_BH2[0],Vel_BH2[1],Vel_BH2[2]);
	   fprintf(PosBH2_Time,"%e %e %e %e \n",header1.time*(Utime/year),Pos_BH2[0],Pos_BH2[1],Pos_BH2[2]);

	cout<<"end of writing velocity and positions second BH"<<endl;

	}
    }
    
  }
  
  
  double r2_BHs=(Udist*Udist/(parsec*parsec))*(Pos_BH1[0]-Pos_BH2[0])*(Pos_BH1[0]-Pos_BH2[0])+(Pos_BH1[1]-Pos_BH2[1])*(Pos_BH1[1]-Pos_BH2[1])+(Pos_BH1[2]-Pos_BH2[2])*(Pos_BH1[2]-Pos_BH2[2]);
  fprintf(time_a_OUT,"%e %e \n",header1.time*(Utime/year),sqrt(r2_BHs)); // distancia en parsecs

  
  double   DeltaGas=Mgas/Mgas_init;
  fprintf(DeltaGas_OUT,"%e %e \n",header1.time*(Utime/year),DeltaGas);


///  MASS INSIDE BLACK HOLES ORBIT//////////////////
  double Menclosed=0.0;
	cout<<"Computing Mass Inside Black Holes... "<<endl; 
  double Radio2_parsec=0.0;
  double Radio2=0.0;	
  for(int k=1;k<=NumPart;k++){ 
    if(P[k].Type!=BHS_TYPE){
//      	 cout<<"Particle is not a BH ..."<<endl;		
      Radio2=P[k].Pos[0]*P[k].Pos[0]+P[k].Pos[1]*P[k].Pos[1]+P[k].Pos[2]*P[k].Pos[2];
      Radio2_parsec=(Udist/(parsec))*Radio2;
      if(Radio2_parsec<=r2_BHs){
	 Menclosed=Menclosed+P[k].Mass;	
//	 cout<<"Particle Enclosed by BHs orbit"<<endl;
//	 cout<<"Particle Mass = "<<P[k].Mass<<"   Enclosed Mass ="<<Menclosed<<endl;
	
	} 
     }
   }
  
  fprintf(EnclosedGas_OUT,"%e %e \n",header1.time*(Utime/year),Menclosed);
  

////////////////////////////////////////////////

  RAD_TEMP_TIME_GAS.close();
  fclose(DeltaGas_OUT);
  fclose(time_a_OUT);

  fclose(EnclosedGas_OUT);
  
  fclose(VelBH1_Time);	
  fclose(VelBH2_Time);	

  fclose(MBH1_Time);
  fclose(MBH2_Time);
  
  
}


////// COPUTE THE RATIO BETWEEN THE MASS ACREATION RATE AND THE EDDINGTON MASS ACREATION RATE LIMIT //////////////////////// 

void Medd_frac(){

 ifstream MBH1_MASS("BHs/Time_MBH1");
 ifstream MBH2_MASS("BHs/Time_MBH2");

 system("wc BHs/Time_MBH1 > BHs/wcMBH1");
 system("awk '{print $1}' BHs/wcMBH1 > BHs/lineas");
 system("rm BHs/wcMBH1");
    
 ifstream lines("BHs/lineas", ios::in);
 int Nline=0;
 lines >> Nline;
 system("rm BHs	/lineas");


  double MassBH1[Nline];	
  double MassBH2[Nline];	
  double Time[Nline];

  int i=0;
  while(i<Nline){
     MBH1_MASS >> Time[i] >> MassBH1[i]; 
     MBH2_MASS >> Time[i] >> MassBH2[i]; 
     i++;
   }


 string NomBH1="BHs/Time_MeddFractionBH1";	
 ofstream TIME_MEDD_FRAC_MBH1(NomBH1.c_str());

 string NomBH2="BHs/Time_MeddFractionBH2";	
 ofstream TIME_MEDD_FRAC_MBH2(NomBH2.c_str());

   double dm1, dm2;
   double dt, Mdot1,Mdot2,Medd1,Medd2;
   Mdot1=Mdot2=Medd1=Medd2=dm1=dm2=dt=0.0;

   double Medd_constant=4.3134e-12;

   for(int k=0;k<Nline-1;k++){
      dm1=MassBH1[k+1]-MassBH1[k];
      dm2=MassBH2[k+1]-MassBH2[k];
      dt=Time[k+1]-Time[k];
        
      Medd1=MassBH1[k]*Medd_constant; 
      Medd2=MassBH2[k]*Medd_constant; 

      Mdot1=dm1/dt;        	
      Mdot2=dm2/dt;        	
          
      TIME_MEDD_FRAC_MBH1<<(Time[k+1]+Time[k])/2.0<<" "<<Mdot1/Medd1<<endl;
      TIME_MEDD_FRAC_MBH2<<(Time[k+1]+Time[k])/2.0<<" "<<Mdot2/Medd2<<endl;
   }

   TIME_MEDD_FRAC_MBH1.close();
   TIME_MEDD_FRAC_MBH2.close();
   MBH1_MASS.close();
   MBH2_MASS.close();
   lines.close();


}

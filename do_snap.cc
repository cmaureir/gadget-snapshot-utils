#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include<fstream>
#include<string>
#include <string.h>
#include <sstream>
#include <math.h>

#include "load_snapshot.h"
#include "bhs.h"

#include "nr3.h"
#include "sort.h"

#include "read_pos_vel_u_type.h"
#include "enclosed_gas_mass.h" 
#include "sfr.h"
#include "kennicutt_schmidt.h"

#define pi 3.141592654 
#define GRAVITY 6.672e-8
#define BOLTZMANN 1.3806e-16
#define PROTONMASS 1.6726e-24
#define Msun 1.989e33
#define parsec 3.08567758e18
#define year 31556926 

int allocate_memory(void);


using namespace std; 

int main(int argc, char **argv)
{
  //char *basename;
  char input_fname[200], path[200], *basename;
  basename=argv[1];
  
  int  type, snapshot_number, files;
  sprintf(path, ".");
  //sprintf(basename, "snapshot");
  //  system("mkdir Graficos");
  
  double Udist, Uvel, Umasa;
  //ESTO ES PARA SIMULACION DE GALAXIAS
  
  //Udist=3.78e21;  // en cm
  //Uvel=1.4272e8;  // en cm/s
  //Umasa=1.154e45; // en gramos
  
  int j,DN,N,BHS_TYPE;
  
  
  cout<<"Ingrese el Numero del Pimer Snapshot "<<endl;
  cin >> j;
  cout<<"Ingrese el Numero del Ultimo Snapshot "<<endl;
  cin >> N;
  cout<<"Ingrese cada cuantos snapshots desea analizar"<<endl;
  cin >> DN; 

  cout<<"Ingrese el TYpe de los BHs"<<endl;
  cin >> BHS_TYPE; 
  
  cout<<"Ingrese la Unidad de Distancia del Sistema en cm"<<endl;
  cin>>Udist;
  
  cout<<"Ingrese la Unidad de Velocidad del Sistema en cm/s"<<endl;
  cin>>Uvel;
  
  cout<<"Ingrese la Unidad de Masa del sistema en g"<<endl;
  cin>>Umasa;
  
  double Utime;
  Utime=Udist/Uvel;
	
  double Mgas_init, Mstars_init;
  Mstars_init=Mgas_init=0.0;
  
  int ID_BH[2];
  ID_BH[0]=ID_BH[1]=0;
  int p=0;

  for(int i=j; i<=N ; i=i+DN){
    snapshot_number= i;                    /* number of snapshot */
    files=1;    
    cout<<"Leyendo Snapshot "<<i<<endl;
    sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
    load_snapshot(input_fname, files);
    
    if(i==j){
      
      for(int k=0;k<NumPart;k++){ 
	if(P[k].Type==0){
	  Mgas_init=Mgas_init+P[k].Mass;
	}
	if(P[k].Type==2 || P[k].Type==4){
	  Mstars_init=Mstars_init+P[k].Mass;
          
	}
	if(P[k].Type==BHS_TYPE){
	  ID_BH[p]=Id[k];
	  p++;
	}
	
      }
      
    }
    
    //reordering();  /* call this routine only if your ID's are set properly */


     read_pos_vel_u(i,j);

//    BHS(i,j,Udist,Uvel,Umasa,Mgas_init,ID_BH[0],ID_BH[1],BHS_TYPE); //


//    write_Mstar_dens(i,Umasa,Utime); // for SFR()

//     data_for_keniccutt_schmidt(i,j); // for Kennicut_Schmidt()
//     enclosed_gas_mass(i,Utime, Udist, Umasa); // this routine has a break !! 

  }
//  SFR();
//  Kennicut_Schmidt(Udist,Uvel,Umasa,N,j,DN);
}





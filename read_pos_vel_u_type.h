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

int read_pos_vel_u(int i, int j)
{	
  
  if(i==j) system("mkdir Pos_Vel_Mass_U");	

  stringstream num;		
  num<<i;

  string NAME_Out_Type0="Pos_Vel_Mass_U/Pos_Vel_Mass_U_Type0_"+num.str();
  ofstream Out_Type0(NAME_Out_Type0.c_str());
  
  string NAME_Out_Type1="Pos_Vel_Mass_U/Pos_Vel_Mass_U_Type1_"+num.str();
  ofstream Out_Type1(NAME_Out_Type1.c_str());

  string NAME_Out_Type2="Pos_Vel_Mass_U/Pos_Vel_Mass_U_Type2_"+num.str();	
  ofstream Out_Type2(NAME_Out_Type2.c_str());

  string NAME_Out_Type3="Pos_Vel_Mass_U/Pos_Vel_Mass_U_Type3_"+num.str();	
  ofstream Out_Type3(NAME_Out_Type3.c_str());

  string NAME_Out_Type4="Pos_Vel_Mass_U/Pos_Vel_Mass_U_Type4_"+num.str();	
  ofstream Out_Type4(NAME_Out_Type4.c_str());

  string NAME_Out_Type5="Pos_Vel_Mass_U/Pos_Vel_Mass_U_Type5_"+num.str();	
  ofstream Out_Type5(NAME_Out_Type5.c_str());

  for(int i=1;i<=NumPart;i++){ 
    if(P[i].Type==0){
	Out_Type0<<P[i].Pos[0]<<" "<<P[i].Pos[1]<<" "<<P[i].Pos[2]<<P[i].Vel[0]<<" "<<P[i].Vel[1]<<" "<<P[i].Vel[2]<<" "<<P[i].Mass<<" "<<P[i].U<<endl;
    }
    if(P[i].Type==1){
	Out_Type1<<P[i].Pos[0]<<" "<<P[i].Pos[1]<<" "<<P[i].Pos[2]<<P[i].Vel[0]<<" "<<P[i].Vel[1]<<" "<<P[i].Vel[2]<<" "<<P[i].Mass<<" "<<P[i].U<<endl;
    }
    if(P[i].Type==2){
	Out_Type2<<P[i].Pos[0]<<" "<<P[i].Pos[1]<<" "<<P[i].Pos[2]<<P[i].Vel[0]<<" "<<P[i].Vel[1]<<" "<<P[i].Vel[2]<<" "<<P[i].Mass<<" "<<P[i].U<<endl;
    }
    if(P[i].Type==3){
	Out_Type3<<P[i].Pos[0]<<" "<<P[i].Pos[1]<<" "<<P[i].Pos[2]<<P[i].Vel[0]<<" "<<P[i].Vel[1]<<" "<<P[i].Vel[2]<<" "<<P[i].Mass<<" "<<P[i].U<<endl;
    }
    if(P[i].Type==4){
	Out_Type4<<P[i].Pos[0]<<" "<<P[i].Pos[1]<<" "<<P[i].Pos[2]<<P[i].Vel[0]<<" "<<P[i].Vel[1]<<" "<<P[i].Vel[2]<<" "<<P[i].Mass<<" "<<P[i].U<<endl;
    }
    if(P[i].Type==5){
	Out_Type5<<P[i].Pos[0]<<" "<<P[i].Pos[1]<<" "<<P[i].Pos[2]<<P[i].Vel[0]<<" "<<P[i].Vel[1]<<" "<<P[i].Vel[2]<<" "<<P[i].Mass<<" "<<P[i].U<<endl;
    }

  }

  Out_Type0.close();
  Out_Type1.close();
  Out_Type2.close();
  Out_Type3.close();
  Out_Type4.close();
  Out_Type5.close();

//  system("cat Pos_Vel_Mass_U/Pos_Vel_Mass_U_Type0 Pos_Vel_Mass_U/Pos_Vel_Mass_U_Type1 Pos_Vel_Mass_U/Pos_Vel_Mass_U_Type2 Pos_Vel_Mass_U/Pos_Vel_Mass_U_Type3 Pos_Vel_Mass_U/Pos_Vel_Mass_U_Type4 Pos_Vel_Mass_U/Pos_Vel_Mass_U_Type5 > Pos_Vel_Mass_U/Pos_Vel_Mass_U");

}


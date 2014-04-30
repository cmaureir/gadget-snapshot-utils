#include <iostream> 
#include<fstream>
#include<iomanip>
#include<string>
#include<iomanip>
#include<stdio.h>
#include<stdlib.h>
#include<cstdlib>
#include <string.h>
#include <sstream>
#include<time.h>
#include <math.h>
#include<algorithm>


#include "/escorpio/ldelvall/Tesis/SimGadget/ICs/Clases/nr.h"
#include "/escorpio/ldelvall/Tesis/SimGadget/ICs/Clases/ran0.cpp"
#include "/escorpio/ldelvall/Tesis/SimGadget/ICs/Clases/ran1.cpp"
#include "/escorpio/ldelvall/Tesis/SimGadget/ICs/Clases/ran2.cpp"
#include "/escorpio/ldelvall/Tesis/SimGadget/ICs/Clases/ran3.cpp"



#define pi 3.141592654

using namespace std; 

/////////////////////////////////////////////////////////////////////////////////////////////

//FUNCION GENERADORA DE THETA RANDOM (PARA DISTRIBUCION SOLO DEPENDIENTE DE r)

double theta_ran(int j)
{
  double  b, theta;
  int idum;
  int t0 =int(double(time(NULL))/10000000);
  
  idum=t0*(j+7)+t0;
  idum=NR::ran0(idum)*2000000000;
  idum=-200000000*NR::ran1(idum)+j*200;//
  b = NR::ran2(idum);
  theta=pi*2*b;
  
  return theta;
}

//FUNCION FACTORIAL

long double PITA(double k) 
{ 
  long double f=1;
  double i = k; 
  while ( i!=0) {
    f=f*(k/i+1);
    i--;
  }
  return f; 
} 

//Funcion Pot

double pot(double d, int n){
  return (n<2) ? d: d * pot(d,n-1);
}

///////////////////////////////////////////////////////////////////////
//        Realizacion de distribucion espacial del disco de gas      //
///////////////////////////////////////////////////////////////////////

double z1_ran(int j, double z0)
{
  double  z;
  int idum, t0, t1;
  double f, Fun;
  
  int k=0;
  
  do {
    t0 =int(double(5*j-1+time(NULL))/10000000);
    
    
    idum=t0*(3*j+5*k)+t0;
    idum=-NR::ran0(idum)*2000000000;
    idum=-2000000000*NR::ran1(idum);
    
    z =z0*(2*NR::ran2(idum)-1);
    
    t1 =int(double(3*j-2+time(NULL))/10000000);
    
    idum=t1*(2*j+7*k)+t1;
    idum=-NR::ran0(idum)*2000000000;
    idum=-2000000000*NR::ran1(idum);
    
    f=3*(1/cosh(z0/(z0)))*(1/cosh(z0/(z0)))*NR::ran2(idum);
    
    Fun=(1/cosh(z/(z0)))*(1/cosh(z/(z0)));
    k++;
  }while(f>Fun);
  
  return z;
}

double r_RanVon(int j, double RM, double Rc){
  int CONT=0;
  int idum, i,idumF, idumy;
  double f0, r, fRM, x, y;
  i=0;
  while(CONT==0){ 
    int t0 =int(5000+3*j+i*3+double(time(NULL))/10000000);
    int t1=int(3000+2*j+5*i+double(time(NULL))/10000000);
    int ty=int(4340+3*j+2*i+double(time(NULL))/10000000);
    
    idum = t0*(2*j+j*j+2*i*i+7*i);
    idumy = ty*(3*j+j*j+4*i*i+6*i);
    idumF= t1*(5*j+j*j+3*i*i+4*i); 
    
    idum=-NR::ran0(idum)*2000000000;
    idum=-2000000000*NR::ran1(idum);//
   
    idumy=-NR::ran0(idumy)*2000000000;
    idumy=-2000000000*NR::ran1(idumy);//  

    idumF=-NR::ran0(idumF)*2000000000;
    idumF=-2000000000*NR::ran1(idumF);//  
   
    x=RM*NR::ran2(idum);
    y=sqrt((RM*RM)-x*x)*NR::ran2(idumy);
    
    r=sqrt(x*x+y*y);
    f0=2*NR::ran2(idumF);
    
    if(r<Rc) fRM=1;
    else if(r>=Rc)fRM=Rc/r;
    if(f0<fRM) CONT++;
    else if(f0>=fRM) i++ ;
    //cout<<x<<" "<<y<<endl;
  }
  
  return r;
  
}

double Mplum(double, double, double, double);
double Mini(double, double, double);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
  
  double G = 21.9617; //segun GADGET2 //22.0511; // en (40 Pc)*((5*10^9 Mo)^{-1})*(156.46Km/s)^2
  //-------------------------------------------------------------------------
  // Aca generamos el Disco de Gas	
  //-----------------------------------------------------------------------//

 
  int N;
  
  double z0;//2.5;  // en unidades de 40 [pc]
  
  double Mgas; // Masa del disco de gas en unidades 5*10^9 [Masas Solares]
  double MBHs;//0.5*Mgas; // MBH es la masa que tiene cada aguhero negro, los dos tienen masa 2MBH
  double M_plumm;
  double M_Ini;
  double Vtipe, d,dini, RM, RA, dB, Rc, RP,Renc;
  


  N=1000000; //Antes era 200000

  MBHs=2; //La suma de la masa de los BH's
  RA=2;  // Separacion inicial BH's
  z0=3;  // Ancho del disco

  Mgas=1;
  double Mdentro=Mgas;
  M_Ini=MBHs;
  M_plumm=0.5;//1.0;
  RM=30.0;
  Renc=3.0; // Radio en el cual el disco tiene una masa encerrada M(Renc)=Mgas
  Rc=6.0;
  RP=19.125;
  Vtipe=0;
  d=1; //d=1 potencial activo, d=0 potencial desactivado
  dini=1;// Potencial inicial esfera uniforme
  dB=0;
  
  //Mcore=Mencerrada*Rc*Rc/(Rorb*Rorb) 
  //M=Mcore(1+2*(Rdisk-Rc)/Rdisk)
  

//  cout<<"Ingrese el Numero de particulas de gas"<<endl;
  //cin>>N;
  //cout<<"Ingrese la masa de los BHs"<<endl;
  //cin>>MBHs;
  //cout<<"Ingrese la masa del disco encerrada por la orbita de los BHs (para comparar con disco chico debe ser 0.04)"<<endl;
  //cin>>Mgas;
  //cout<<"Ingresar el radio del disco"<<endl;
  //cin>>RM;
 // cout<<"Ingresar el radio de la orbita de los BHs"<<endl;
  //cin>>RA;
  //cout<<"Ingresar el radio del core del disco"<<endl;
  //cin>>Rc;
  
 // cout<<"Ingresar el ancho del Disco en el Core (Rcore)"<<endl;
  //cin>>z0;
//cout<<"Ingresar 0 para calcular velocidad aproximada (solo masas encerrada) y 1 para velocidad exacta"<<endl;
  //cin>>Vtipe;
 // cout<<"Ingresar 0 para un disco sin potencial y 1 para un disco con potencial"<<endl;
  //cin>>d;
 // cout<<"Ingresar 0 para un disco sin Borde y 1 para un disco con Borde"<<endl;
  //cin>>dB;
  
  double Menc=Mgas;
  cout<<"La masa del core es "<<Mgas*Rc*Rc/(Renc*Renc)<<endl;
  Mgas=Mgas*Rc*Rc/(Renc*Renc);
  double Mc=Mgas;
  double n0=Mc/(pi*(Rc)*(Rc));//Este es la densidad superficial del core
  double cs=pi*z0*G*n0 ; //Factor 4 a ojo PORQUE FUNCIONA??
  double CSS=cs;
cout<<"Para este disco Cs2 = "<<cs<<endl;
  Mgas=Mgas+2*pi*n0*Rc*(RM-Rc);
  cout<<"Mgas = "<<Mgas<<endl;

  //double NRbh=N*(Menc/Mgas); //Numero de particulas dentro del radio de orbita de los BHs
  //  int Nc=NRbh*double(Mc/Menc);
  //cout<<"Nc =  "<<Nc<<endl;
  
///////////NOMBRANDO ARCHIVOS DE SALIDA///////////////////////////////
stringstream cs2;
  cs2<<cs;
  stringstream Zmed;
  Zmed<<z0;
  stringstream pott;
  pott<<d;
  stringstream pottIni;
  pottIni<<dini;
  stringstream ROrb;
  ROrb<<RA;
  stringstream RDISC;
  RDISC<<RM;
  stringstream Dens;
  Dens<<n0;
  stringstream masa;
  masa<<Menc/(MBHs);
  stringstream NumPart;
  NumPart<<N;
  
  string NomPos, NomVel, NomRV;
  if(Vtipe==1){
    NomPos="Pos_Cs"+cs2.str()+"_Z0_"+Zmed.str()+"_MassRatio"+masa.str()+"_Rdisk"+RDISC.str()+"_Rorb"+ROrb.str()+"_PotEsf"+pott.str()+"_PotIni"+pottIni.str()+"_NumPart"+NumPart.str();
    NomVel="Vel_Cs"+cs2.str()+"_Z0_"+Zmed.str()+"_MassRatio"+masa.str()+"_Rdisk"+RDISC.str()+"_Rorb"+ROrb.str()+"_PotEsf"+pott.str()+"_PotIni"+pottIni.str()+"_NumPart"+NumPart.str();
    NomRV="RV_Cs"+cs2.str()+"_Z0_"+Zmed.str()+"_MassRatio"+masa.str()+"_Rdisk"+RDISC.str()+"_Rorb"+ROrb.str()+"_PotEsf"+pott.str()+"_PotIni"+pottIni.str()+"_NumPart"+NumPart.str();
  }
  else if(Vtipe==0){
    NomPos="Pos_Cs"+cs2.str()+"_Z0_"+Zmed.str()+"_MassRatio"+masa.str()+"_Rdisk"+RDISC.str()+"_Rorb"+ROrb.str()+"_PotEsf"+pott.str()+"_PotIni"+pottIni.str()+"_NumPart"+NumPart.str()+"VelAprox";
    NomVel="Vel_Cs"+cs2.str()+"_Z0_"+Zmed.str()+"_MassRatio"+masa.str()+"_Rdisk"+RDISC.str()+"_Rorb"+ROrb.str()+"_PotEsf"+pott.str()+"_PotIni"+pottIni.str()+"_NumPart"+NumPart.str()+"VelAprox";
    NomRV="RV_Cs"+cs2.str()+"_Z0_"+Zmed.str()+"_MassRatio"+masa.str()+"_Rdisk"+RDISC.str()+"_Rorb"+ROrb.str()+"_PotEsf"+pott.str()+"_PotIni"+pottIni.str()+"_NumPart"+NumPart.str()+"VelAprox";       
  }
  ofstream xyz(NomPos.c_str());
  ofstream VxVyVz(NomVel.c_str());
  ofstream RV(NomRV.c_str());

  ////////////////////////////////////////////////////////////////////////
  
  long double x, y, z, r ,Vx, Vy, Vz, theta;
   
  int i=0;
  int k=0;
  
  double dx,dy,dz,DS;
  double D=0.0; 
  long double vc;


  int Nin=0;
  while(i<N){
    
    r=r_RanVon(i+k,RM,Rc);    
    theta=theta_ran(2*i)+2*pi*theta_ran(3*i);
    
    x=r*cos(theta);
    y=r*sin(theta);
    

    if(r<Rc){    
      z=z1_ran(i,z0);
      vc=sqrt(pi*G*n0*r+G*d*Mplum(r,M_plumm,RP,Renc)/r+G*dini*Mini(r,M_Ini,RA)/r);
      if(r<Renc){
	Nin++;
      }
    }
    else if(r>=Rc){
      z=z1_ran(i,z0*(r/Rc));
      vc=sqrt(pi*G*n0*(2*Rc-(Rc*Rc/r))-CSS+G*d*Mplum(r,M_plumm,RP,Renc)/r+G*dini*Mini(r,M_Ini,RA)/r);
    }
    Vx=vc*sin(theta);
    Vy=-vc*cos(theta);
    Vz=0.0;
    
    xyz << x <<" "<< y <<" "<< z <<endl;
    VxVyVz<<Vx<<" "<< Vy <<" "<< Vz <<endl;
    RV<<r<<" "<<vc<<endl;
    //    cout<<r<<" "<<vc<<endl;
    i++;	
    
  }
  cout<<i<<endl;
  cout <<"Num Part Dentro Renc=3 "<<Nin<<endl;
  double mNgas=Mdentro/double(Nin);
  cout <<"Masa por particula Menc/Nin = "<<mNgas<<endl;
  
  ///////////////////////////////////////
  //                BHs               //
  ///////////////////////////////////////
  long double SUM00;
  
  vc=sqrt(pi*G*n0*RA+G*(MBHs/2)/(RA*4)+G*d*Mplum(RA,M_plumm,RP,Renc)/RA);
  
  double xbh1=RA;
  double ybh1=0;
  double zbh1=0;
  
  double Vxbh1=0;
  double Vybh1=-vc;
  double Vzbh1=0;
  
  xyz << xbh1 <<" "<< ybh1 <<" "<< zbh1 <<endl;
  VxVyVz<<Vxbh1<<" "<< Vybh1 <<" "<< Vzbh1 <<endl;
  
  double xbh2=-RA;
  double ybh2=0;
  double zbh2=0;
  
  double Vxbh2=0;
  double Vybh2=vc;
  double Vzbh2=0;
  
  xyz << xbh2 <<" "<< ybh2 <<" "<< zbh2 <<endl;
  VxVyVz<<Vxbh2<<" "<< Vybh2 <<" "<< Vzbh2 <<endl;
  
  
  
  RV.close();
  xyz.close();
  VxVyVz.close();
  
  
}



////////////// Aca definimos la distribucion axisimetrica que representa a los BHs ///////////////////////////

double Mplum(double r, double Mp, double RP, double Renc){
	double M;
	M=Mp*(r*r*r/pot(sqrt(r*r+RP*RP),3))*(pot(sqrt(Renc*Renc+RP*RP),3)/(Renc*Renc*Renc));
	return M;
}

double Mini(double r, double MBH, double RA){
	double M;
	M=MBH*r*r*r/(RA*RA*RA);
	if(r>=RA) M=MBH;
	return M;
}

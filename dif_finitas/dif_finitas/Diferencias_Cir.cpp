#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

int main()
{
	int i,j,k;
	
	double rho1=1.0909;
	double rho2=6.0;
	double rho3=8.2283;
    double rhog=0.0; 
    double z0=1.43;
    
	double dr=0.001*(rho1+rho2+rho3+2*rhog); 
	double dz=0.01*2.0*z0; 
	
	int Vac_points=40;//En cada lado en rho
	
	int Nrho=(int)((rho1+rho2+rho3+2*rhog)/dr+0.5+Vac_points); //+0.5 porque al castear se trunca
	int Nz=(int)(2.0*z0/dz+0.5); 
	
	int j0=(int)(z0/dz+0.5);
	
	//Mem alloc e inicialización del grid y el bool auxiliar - CHECKED
	double **Grid=(double**)malloc(Nrho*sizeof(double*));
	for(i=0;i<Nrho;i++)
	{
		Grid[i]=(double*)malloc(Nz*sizeof(double));
		for(j=0;j<Nz;j++)
			Grid[i][j]=0;
	}
	bool **Electrode=(bool**)malloc(Nrho*sizeof(bool*));
	for(i=0;i<Nrho;i++)
	{
		Electrode[i]=(bool*)malloc(Nz*sizeof(bool));
		for(j=0;j<Nz;j++)
			Electrode[i][j]=false;
	}
	
	double V1=100.;
	double V2=0.;
	double V3=0.;
	
	for(i=0;i<Nrho;i++)
         for(j=0;j<1;j++)//(10*rhog/dr+0.5);j++)
         {
                if(i*dr<=rho1){Electrode[i][j]=true;Grid[i][j]=V1;}
                if(i*dr>=rho1+rhog&&i*dr<=rho2){Electrode[i][j]=true;Grid[i][j]=V2;}
                if(i*dr>=rhog+rho2&&i*dr<=rho3){Electrode[i][j]=true;Grid[i][j]=V3;}
         }
	
	double ar1=(dz*dz)/(2*(dr*dr+dz*dz)), az=ar1*(dr*dr)/(dz*dz);
	double *ar2=(double*)malloc(Nrho*sizeof(double));  //Pointer?
	double a0r=(dz*dz)/(dr*dr+2.*dz*dz), a0z=0.5*(dr*dr)/(dz*dz)*a0r; //Nueva forma de saltar la singularidad (MODIF)
	for(i=1;i<Nrho;i++) ar2[i]=ar1/(2.*i); //MODIF
	
	for(k=0;k<5e4;k++)
	{
		if(k%1000==0)cout<<k<<endl;
		for(i=0;i<Nrho;i++)
			for(j=0;j<Nz;j++)
				if(!Electrode[i][j])
				{
					//Corners
					if(i==0&&j==0)Grid[i][j]=2*a0r*Grid[i+1][j]+a0z*Grid[i][j+1]; //MODIF
					if(i==0&&j==Nz-1)Grid[i][j]=2*a0r*Grid[i+1][j]+a0z*Grid[i][j-1]; //MODIF
					if(i==Nrho-1&&j==0)Grid[i][j]=((ar1-ar2[i])*Grid[i-1][j]+az*Grid[i][j+1])/(1-az);
					if(i==Nrho-1&&j==Nz-1)Grid[i][j]=(ar1-ar2[i])*Grid[i-1][j]+az*Grid[i][j];
					//Sides
					if(i==0&j!=0&&j!=Nz-1)Grid[i][j]=2*a0r*Grid[i+1][j]+a0z*(Grid[i][j+1]+Grid[i][j-1]); //MODIF
					if(i==Nrho-1&j!=0&&j!=Nz-1)Grid[i][j]=(ar1-ar2[i])*Grid[i-1][j]+az*(Grid[i][j+1]+Grid[i][j-1]);
					if(i!=0&&i!=Nrho-1&&j==0)Grid[i][j]=ar1*(Grid[i+1][j]+Grid[i-1][j])+ar2[i]*(Grid[i+1][j]-Grid[i-1][j])+az*Grid[i][j+1];
					if(i!=0&&i!=Nrho-1&&j==Nz-1)Grid[i][j]=ar1*(Grid[i+1][j]+Grid[i-1][j])+ar2[i]*(Grid[i+1][j]-Grid[i-1][j])+az*Grid[i][j-1];
					//Rest
					if(i!=0&&i!=Nrho-1&&j!=0&&j!=Nz-1)Grid[i][j]=ar1*(Grid[i+1][j]+Grid[i-1][j])+ar2[i]*(Grid[i+1][j]-Grid[i-1][j])+az*(Grid[i][j+1]+Grid[i][j-1]);
				}
	}
	
	ofstream out;
	out.open("Test.dat");
	
	for(i=0;i<Nrho;i++)
		for(j=0;j<Nz;j++)
		{
			//if(j==0)out<< i*dr  << "\t" << Grid[i][j] << "\n" ;
			out<< i*dr << "\t" << j*dz << "\t" << Grid[i][j] << "\n" ;
		out<< -i*dr << "\t" << j*dz << "\t" << Grid[i][j] << "\n" ;
		}
	
	out.close();
	 
	//Derivative at 0,0
	
	double E=(-Grid[0][j0+2]+8.0*Grid[0][j0+1]-8*Grid[0][j0-1]+Grid[0][j0-2])/(12.0*dz);
	E*=1000;
	
	double du;
	double At = 187;//40;
	double Z = 75;//20;
	double mass = abs((Z - 1))*9.10938188*1e-31 + Z*1.6726231*1e-27 + (At - Z)*1.67492728*1e-27;
	double T=300;
	double fz=55e3;
	double A = sqrt(2*1.381*1e-23*T/(mass*pow(fz*2*acos(-1.0),2.0)));
	
	double C = 10*1e-12;
	double q = 1.6*1e-19;
	du=A*E*q/(1.0*C*V1);
	At = 40;
	Z = 20;
	mass = abs((Z - 1))*9.10938188*1e-31 + Z*1.6726231*1e-27 + (At - Z)*1.67492728*1e-27;
	
	double b=(q*E*du)/(2.0*mass*2*fz*acos(-1.0)*V1);
	
	cout<<"Du: \t"<<du<<endl;
	cout<<"b: \t"<<b<<endl;
	system("pause");
	
	return 0;
}

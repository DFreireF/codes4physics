#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

int main()
{
	int i,j,k;
	
	double rho0=6.0; //El parámetro que lo define todo
	double dr=0.01*rho0; //Debe ser ~0.01*rho0 o menor
	double dz=0.01*rho0; //Idem
	int Nrho=(int)(1.1*rho0/dr+0.5); //+0.5 porque al castear se trunca
	int Nz=(int)( ( (0.4+10.58/1.0239 )*rho0 ) / dz + 0.5 ); //Idem
	if(Nz%2==0) Nz++; //Para que Nz sea impar y haya siempre un centro
	int j0=Nz/2+1; //Centro en j0*dr
	
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
	
	double Weight=1.25;
	
	/*double Ve=-0.0448963;
	double Vc=-0.0209951;
	double V0=0.0551037;
	double Ve2=-0.0448963;
	double Vc2=-0.0209951;
	double V02=0.0551037;*/
	
	// damos valores a los potenciales en los electrodos. 'c' de corrección
	double Ve1=0.;
	double Ve2=0.;
	double Vc1=20.;
	double Vc2=0.;
	double V0 =0.;
	
	
	//Inicialización
	int imin, imax, jmin, jmax;
	imin=(int)(rho0/dr+0.5);//Radio interno electrodos
	imax=(int)( (1+2*0.0303/1.0239)*rho0 /dr + 0.5 ); //Radio externo electrodos (espesor 2*gap asignado arbitrariamente)
	for(i=imin;i<=imax;i++) //Los <= son a propósito!!
	{
		//E2 (z negativos)
		jmin=(int)( ( 0.2 * rho0 ) / dr + 0.5); //0.2rho
		jmax=(int)( ( ( 0.2+4.327/1.0239 )*rho0 )/dr + 0.5); //0.2rho+ze
		for(j=jmin;j<=jmax;j++)
		{
			Electrode[i][j]=true;
			Grid[i][j]=Ve2;
		}
		
		//C2 (z negativos)
		jmin=(int)( ( ( 0.2+4.3573/1.0239 )*rho0 ) / dr + 0.5 );  //0.2rho+ze+zg
		jmax=(int)( ( ( 0.2+5.1924/1.0239 )*rho0 ) / dr + 0.5 );  //0.2rho+ze+zg+zc
		for(j=jmin;j<=jmax;j++)
		{
			Electrode[i][j]=true;
			Grid[i][j]=Vc2;
		}
		
		//Ring
		jmin=(int)( ( (0.2+5.2227/1.0239 )*rho0 ) / dr + 0.5 ); //0.2rho+ze+zg+zc+zg
		jmax=(int)( ( (0.2+5.3573/1.0239 )*rho0 ) / dr + 0.5 ); //0.2rho+ze+zg+zc+zg+zr
		for(j=jmin;j<=jmax;j++)
		{
			Electrode[i][j]=true;
			Grid[i][j]=V0;
		}
		
		//C1 (z positivos)
		jmin=(int)( ( (0.2+5.3876/1.0239 )*rho0 ) / dr + 0.5 ); //0.2rho+ze+zg+zc+zg+zr+zg
		jmax=(int)( ( (0.2+6.2227/1.0239 )*rho0 ) / dr + 0.5 ); //0.2rho+ze+zg+zc+zg+zr+zg+zc
		for(j=jmin;j<=jmax;j++)
		{
			Electrode[i][j]=true;
			Grid[i][j]=Vc1;
		}
		
		//E1 (z positivos)
		jmin=(int)( ( (0.2+6.2530/1.0239 )*rho0 ) / dr + 0.5 ); //0.2rho+ze+zg+zc+zg+zr+zg+zc+zg
		jmax=(int)( ( (0.2+10.5800/1.0239 )*rho0 ) / dr + 0.5 ); //0.2rho+ze+zg+zc+zg+zr+zg+zc+zg+ze
		for(j=jmin;j<=jmax;j++)
		{
			Electrode[i][j]=true;
			Grid[i][j]=Ve1;
		}
	}
	
	double ar1=(dz*dz)/(2*(dr*dr+dz*dz)), az=ar1*(dr*dr)/(dz*dz);
	double *ar2=(double*)malloc(Nrho*sizeof(double));
	double a0r=(dz*dz)/(dr*dr+2.*dz*dz), a0z=0.5*(dr*dr)/(dz*dz)*a0r; //Nueva forma de saltar la singularidad (MODIF)
	for(i=1;i<Nrho;i++) ar2[i]=ar1/(2.*i); //MODIF
	
	for(k=0;k<1e5;k++)
	{
		if(k%1000==0)cout<<k<<endl; //contador cada 1000 pasos
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
	//escribimos resultados en fichero
	ofstream out;
	out.open("Test.dat");
	
	for(i=0;i<Nrho;i++)
		for(j=0;j<Nz;j++)
		{
			out<< i*dr << "\t" << (j-j0)*dz << "\t" << Grid[i][j] << "\n" ;
			out<< -i*dr << "\t" << (j-j0)*dz << "\t" << Grid[i][j] << "\n" ;
		}
	
	out.close();
	 
	//Derivative at 0,0
	
	/*double E=(-Grid[0][j0+2]+8.0*Grid[0][j0+1]-8*Grid[0][j0-1]+Grid[0][j0-2])/(12.0*dz);
	E*=1000;
	
	double du;
	double At = 187;//40;
	double Z = 75;//20; //renio
	double mass = abs((Z - 1))*9.10938188*1e-31 + Z*1.6726231*1e-27 + (At - Z)*1.67492728*1e-27;
	double T=4.;
	double fz=100*1e3;
	double A = sqrt(2*1.381*1e-23*T/(mass*pow(fz*2*acos(-1.0),2.0)));
	
	double C = 10*1e-12;
	double q = 1.6*1e-19;
	du=A*E*q/(1.0*C*Vc1);
	At = 40;
	Z = 20;
	mass = abs((Z - 1))*9.10938188*1e-31 + Z*1.6726231*1e-27 + (At - Z)*1.67492728*1e-27;
	
	double b=(q*du*E)/(2.0*mass*2*fz*acos(-1.0)*Vc1);
	
	cout<<"Du: \t"<<du<<endl;
	cout<<"b: \t"<<b<<endl;
	system("pause");*/
	
	return 0;
}

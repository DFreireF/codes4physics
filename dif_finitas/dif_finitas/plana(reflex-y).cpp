#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

int main()
{
	int i, j, k, m;
	int Nx, Ny, Nz;
	double d, le, lc, lr, s0, s1;
	double dx=0.01, dy=0.01, dz=0.01;
	double ax=1./(2.*(1+dx*dx/(dy*dy)+dx*dx/(dz*dz))), ay=1./(2.*(1+dy*dy/(dx*dx)+dy*dy/(dz*dz))), az=1./(2.*(1+dz*dz/(dx*dx)+dz*dz/(dy*dy)));
	
	//Mem alloc and init
	double ***Grid=(double***)malloc(Nx*sizeof(double**));
	bool ***Electrode=(bool***)malloc(Nx*sizeof(bool**));
    for(i=0;i<Nx;i++)
    {
        Grid[i]=(double**)malloc(Ny*sizeof(double*));
		Electrode[i]=(bool**)malloc(Ny*sizeof(bool*));
        for(j=0;j<Ny;j++)
        {
            Grid[i][j]=(double*)malloc(Nz*sizeof(double));
			Electrode[i][j]=(bool*)malloc(Nz*sizeof(bool));
            for(k=0;k<Nz;k++)
			{
				Grid[i][j][k]=0.;
				Electrode[i][j][k];
			}
        }

    }
	
	double Ve1 =0.;
	double Vc1 =50.;
	double Vr =0.;
	double Ve2=0.;
	double Vc2=0.;
	double Vg=0.;
	
	//Inicializa electrodos
	int imin, imax, jmax, kmin, kmax;
	jmax=(int)( d/dz + 0.5 ); //0.5 porque al castear se trunca. Electrodos de espesor 2*gap (ojo!: reflexión en y)
	
	
	for(m=0;m<1e5;m++)
	{
		if(m%10000==0)cout<<k<<endl;
		
		for(i=0;i<Nx;i++)
			for(j=0;j<Ny;j++)
				for(k=0;k<Nz;k++)
					if(!Electrode[i][j])
					{
						//Vértices
						if(i==0&&j==0&&k==0) 				Grid[i][j][k]=2.*ax*Grid[i+1][j][k]+2.*ay*Grid[i][j+1][j]+az*Grid[i][j][k+1];	//Reflex x & reflex y & tierra z
						else if(i==0&&j==0&&k==Nz-1)		Grid[i][j][k]=2.*ax*Grid[i+1][j][k]+2.*ay*Grid[i][j+1][j]+az*Grid[i][j][k-1];	//Reflex x & reflex y & tierra z
						else if(i==0&&j==Ny-1&&k==0)		Grid[i][j][k]=2.*ax*Grid[i+1][j][k]+ay*Grid[i][j-1][j]+az*Grid[i][j][k+1];		//Reflex x & tierra y & tierra z
						else if(i==0&&j==Ny-1&&k==Nz-1)		Grid[i][j][k]=2.*ax*Grid[i+1][j][k]+ay*Grid[i][j-1][j]+az*Grid[i][j][k-1];		//Reflex x & tierra y & tierra z
						else if(i==Nx-1&&j==0&&k==0)		Grid[i][j][k]=ax*Grid[i-1][j][k]+2.*ay*Grid[i][j+1][j]+az*Grid[i][j][k+1];		//Tierra x & reflex y & tierra z
						else if(i==Nx-1&&j==0&&k==Nz-1)		Grid[i][j][k]=ax*Grid[i-1][j][k]+2.*ay*Grid[i][j+1][j]+az*Grid[i][j][k-1];		//Tierra x & reflex y & tierra z
						else if(i==Nx-1&&j==Ny-1&&k==0)		Grid[i][j][k]=ax*Grid[i-1][j][k]+ay*Grid[i][j-1][j]+az*Grid[i][j][k+1];			//Tierra x & tierra y & tierra z
						else if(i==Nx-1&&j==Ny-1&&k==Nz-1)	Grid[i][j][k]=ax*Grid[i-1][j][k]+ay*Grid[i][j-1][j]+az*Grid[i][j][k-1];			//Tierra x & tierra y & tierra z
						//Aristas
						else if(i==0&&j==0)			Grid[i][j][k]=2*ax*Grid[i+1][j][k]+2.*ay*Grid[i][j+1][k]+az*(Grid[i][j][k+1]+Grid[i][j][k-1]);	//Reflex x & reflex y
						else if(i==0&&j==Ny-1)		Grid[i][j][k]=2*ax*Grid[i+1][j][k]+ay*Grid[i][j-1][k]+az*(Grid[i][j][k+1]+Grid[i][j][k-1]);		//Reflex x & tierra y
						else if(i==Nx-1&&j==0)		Grid[i][j][k]=ax*Grid[i-1][j][k]+2.*ay*Grid[i][j+1][k]+az*(Grid[i][j][k+1]+Grid[i][j][k-1]);	//Tierra x & reflex y
						else if(i==Nx-1&&j==Ny-1)	Grid[i][j][k]=ax*Grid[i-1][j][k]+ay*Grid[i][j-1][k]+az*(Grid[i][j][k+1]+Grid[i][j][k-1]);		//Reflex x & tierra y
						else if(i==0&&k==0)			Grid[i][j][k]=2*ax*Grid[i+1][j][k]+ay*(Grid[i][j+1][k]+Grid[i][j-1][k])+az*Grid[i][j][k+1];		//Reflex x & tierra z
						else if(i==0&&k==Nz-1)		Grid[i][j][k]=2*ax*Grid[i+1][j][k]+ay*(Grid[i][j+1][k]+Grid[i][j-1][k])+az*Grid[i][j][k-1];		//Reflex x & tierra z
						else if(i==Nx-1&&k==0)		Grid[i][j][k]=ax*Grid[i+1][j][k]+ay*(Grid[i][j+1][k]+Grid[i][j-1][k])+az*Grid[i][j][k+1];		//Tierra x & tierra z
						else if(i==Nx-1&&k==Nz-1)	Grid[i][j][k]=ax*Grid[i+1][j][k]+ay*(Grid[i][j+1][k]+Grid[i][j-1][k])+az*Grid[i][j][k-1];		//Tierra x & tierra z
						else if(j==0&&k==0)			Grid[i][j][k]=ax*(Grid[i+1][j][k]+Grid[i-1][j][k])+2.*ay*Grid[i][j+1][k]+az*Grid[i][j][k+1];	//Reflex y & tierra z
						else if(j==0&&k==Nz-1)		Grid[i][j][k]=ax*(Grid[i+1][j][k]+Grid[i-1][j][k])+2.*ay*Grid[i][j+1][k]+az*Grid[i][j][k-1];	//Reflex y & tierra z
						else if(j==Ny-1&&k==0)		Grid[i][j][k]=ax*(Grid[i+1][j][k]+Grid[i-1][j][k])+ay*Grid[i][j-1][k]+az*Grid[i][j][k+1];		//Tierra y & tierra z
						else if(j==Ny-1&&k==Nz-1)	Grid[i][j][k]=ax*(Grid[i+1][j][k]+Grid[i-1][j][k])+ay*Grid[i][j-1][k]+az*Grid[i][j][k-1];		//Tierra y & tierra z
						//Caras
						else if(i==0)		Grid[i][j][k]=2.*Grid[i+1][j][k]+ay*(Grid[i][j+1][k]+Grid[i][j-1][k])+az*(Grid[i][j][k+1]+Grid[i][j][k-1]);		//Reflex x
						else if(i==Nx-1)	Grid[i][j][k]=Grid[i-1][j][k]+ay*(Grid[i][j+1][k]+Grid[i][j-1][k])+az*(Grid[i][j][k+1]+Grid[i][j][k-1]);		//Tierra x
						else if(j==0)		Grid[i][j][k]=ax*(Grid[i+1][j][k]+Grid[i-1][j][k])+2.*ay*Grid[i][j+1][k]+az*(Grid[i][j][k+1]+Grid[i][j][k-1]);	//Reflex y
						else if(j==Ny-1)	Grid[i][j][k]=ax*(Grid[i+1][j][k]+Grid[i-1][j][k])+ay*Grid[i][j-1][k]+az*(Grid[i][j][k+1]+Grid[i][j][k-1]);		//Tierra y
						else if(k==0)		Grid[i][j][k]=ax*(Grid[i+1][j][k]+Grid[i-1][j][k])+ay*(Grid[i][j+1][k]+Grid[i][j-1][k])+az*Grid[i][j][k+1];		//Tierra z
						else if(k==Nz-1)	Grid[i][j][k]=ax*(Grid[i+1][j][k]+Grid[i-1][j][k])+ay*(Grid[i][j+1][k]+Grid[i][j-1][k])+az*Grid[i][j][k-1];		//Tierra z
						//Resto
						else	Grid[i][j][k]=ax*(Grid[i+1][j][k]+Grid[i-1][j][k])+ay*(Grid[i][j+1][k]+Grid[i][j-1][k])+az*(Grid[i][j][k+1]+Grid[i][j][k-1]);
					}
	}
	ofstream out;out.open("Test.dat");
	
	for(i=0;i<Nrho;i++)
		for(j=0;j<Nz;j++)
			out<<j<<"\t"<<i<<"\t"<<Grid[i][j]<<"\n";
	
	 
	//Derivative at 0,0 [55,275]
	
	double E=(-Grid[55][277]+8.0*Grid[55][276]-8*Grid[55][274]+Grid[55][273])/(12.0*dz);
	E*=1000;
	
	double du;
	double At = 187;//40;
	double Z = 75;//20;
	double mass = abs((Z - 1))*9.10938188*1e-31 + Z*1.6726231*1e-27 + (At - Z)*1.67492728*1e-27;
	double T=300;
	double fz=100*1e3;
	double A = sqrt(2*1.381*1e-23*T/(mass*pow(fz*2*acos(-1.0),2.0)));
	
	double C = 10*1e-12;
	double q = 1.6*1e-19;
	du=A*E*q/(1.0*C*Vc);
	At = 40;
	Z = 20;
	mass = abs((Z - 1))*9.10938188*1e-31 + Z*1.6726231*1e-27 + (At - Z)*1.67492728*1e-27;
	
	double b=(q*du*E)/(2.0*mass*2*fz*acos(-1.0)*Vc);
	
	cout<<"Du: \t"<<du<<endl;
	cout<<"b: \t"<<b<<endl;
	
	system("pause");
	//system("pause");
	return 0;
}

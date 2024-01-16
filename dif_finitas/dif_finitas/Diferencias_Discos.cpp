#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

int main()
{
    int i,j,k;
    
    double rmax=6., zmax=3., dr=0.01, dz=0.01; //mm
    int Nrho=(int)(rmax/dr+1.5), Nz=(int)(zmax/dr+1.5); //Ojo! redondeo (+0.5). El +1 está justificado (hasta z=+1.50 en lugar de +1.49)
    
    //Mem alloc & init
    double **Grid=(double**)malloc(Nrho*sizeof(double*));
    bool **Electrode=(bool**)malloc(Nrho*sizeof(bool*));
    for(i=0;i<Nrho;i++)
    {
        Grid[i]=(double*)malloc(Nz*sizeof(double));
        Electrode[i]=(bool*)malloc(Nz*sizeof(bool));
        for(j=0;j<Nz;j++)
        {
            Grid[i][j]=0;
            Electrode[i][j]=false;
        }    
    }
    
    double res;
    double Weight=1.25;
    
   /* double Ve=12;
    double Vc=9;
    double Vr=6;
    double Vi=3;
    double Vi2=3;
    double Vr2=6;
    double Vc2=9;
    double Ve2=12;*/
    
    double Ve1=0.;
    double Vc1=10.;
    double Vr1=0.;
    double Vi1=0.;
    double Vi2=0.;
    double Vr2=0.;
    double Vc2=0.;
    double Ve2=0.;
    
    double li1=0.20, li2=0.54, lr1=0.55, lr2=0.95, lc1=0.96, lc2=1.49, le1=1.50, le2=3.40;
    double hi=1., hr=1., hc=0.5, he=0.5;
    
    int jmin, jmax, imin, imax;
    //Electrodo i
    jmin=(int)((zmax/2.+hi)/dz+0.5); //Donde empieza el de arriba
    jmax=(int)((zmax/2.-hi)/dz+0.5); //Donde acaba el de abajo
    imin=(int)(li1/dr+0.5);
    imax=(int)(li2/dr+0.5);
    for(i=imin;i<=imax;i++) //Los <= son a propósito
    {
        for(j=0;j<=jmax;j++)
        {
            Grid[i][j]=Vi1;
            Electrode[i][j]=true;
        }
        for(j=jmin;j<Nz;j++)
        {
            Grid[i][j]=Vi2;
            Electrode[i][j]=true;
        }
    }
    //Electrodo r
    jmin=(int)((zmax/2.+hr)/dz+0.5);
    jmax=(int)((zmax/2.-hr)/dz+0.5);
    imin=(int)(lr1/dr+0.5);
    imax=(int)(lr2/dr+0.5);
    for(i=imin;i<=imax;i++)
    {
        for(j=0;j<=jmax;j++)
        {
            Grid[i][j]=Vr1;
            Electrode[i][j]=true;
        }
        for(j=jmin;j<Nz;j++)
        {
            Grid[i][j]=Vr2;
            Electrode[i][j]=true;
        }
    }
    //Electrodo c
    jmin=(int)((zmax/2.+hc)/dz+0.5);
    jmax=(int)((zmax/2.-hc)/dz+0.5);
    imin=(int)(lc1/dr+0.5);
    imax=(int)(lc2/dr+0.5);
    for(i=imin;i<=imax;i++)
    {
        for(j=0;j<=jmax;j++)
        {
            Grid[i][j]=Vc1;
            Electrode[i][j]=true;
        }
        for(j=jmin;j<Nz;j++)
        {
            Grid[i][j]=Vc2;
            Electrode[i][j]=true;
        }
    }
    //Electrodo e
    jmin=(int)((zmax/2.+he)/dz+0.5);
    jmax=(int)((zmax/2.-he)/dz+0.5);
    imin=(int)(le1/dr+0.5);
    imax=(int)(le2/dr+0.5);
    for(i=imin;i<=imax;i++)
    {
        for(j=0;j<=jmax;j++)
        {
            Grid[i][j]=Ve1;
            Electrode[i][j]=true;
        }
        for(j=jmin;j<Nz;j++)
        {
            Grid[i][j]=Ve2;
            Electrode[i][j]=true;
        }
    } 
    
    ofstream out;
    out.open("Test-discos.dat");
    
    for(i=0;i<Nrho;i++)
        for(j=0;j<Nz;j++)
            out << i*dr << "\t" << j*dz << "\t" << Grid[i][j] << "\n";
    for(i=0;i<Nrho;i++)
        for(j=0;j<Nz;j++)
            out << -i*dr << "\t" << j*dz << "\t" << Grid[i][j] << "\n";
    out.close();
    
    double ar1=(dz*dz)/(2*(dr*dr+dz*dz)), az=ar1*(dr*dr)/(dz*dz);
    double *ar2=(double*)malloc(Nrho*sizeof(double));
    double a0r=(dz*dz)/(dr*dr+2.*dz*dz), a0z=0.5*(dr*dr)/(dz*dz)*a0r; //Nueva forma de saltar la singularidad (MODIF)
    for(k=0;k<Nrho;k++) ar2[k]=ar1/(2.*k+0.5);
    //for(k=0;k<Nrho;k++) ar2[k]=ar1/(2.*k+1.);//Wrong
    
    for(k=0;k<1e5;k++)
    {
        if(k%1000==0)cout<<k<<endl;
        for(i=0;i<Nrho;i++)
        for(j=0;j<Nz;j++)
        {
            if(!Electrode[i][j])
            {
               /*
                //Corners
                if(i==0&&j==0)Grid[i][j]=((ar1+ar2[i])*Grid[i+1][j]+az*Grid[i][j+1])/(1-ar1+ar2[i]-az);
                if(i==0&&j==Nz-1)Grid[i][j]=((ar1+ar2[i])*Grid[i+1][j]+az*Grid[i][j-1])/(1-ar1+ar2[i]);
                if(i==Nrho-1&&j==0)Grid[i][j]=((ar1-ar2[i])*Grid[i-1][j]+az*Grid[i][j+1])/(1-az);
                if(i==Nrho-1&&j==Nz-1)Grid[i][j]=(ar1-ar2[i])*Grid[i-1][j]+az*Grid[i][j];
                //Sides
                if(i==0&j!=0&&j!=Nz-1)Grid[i][j]=((ar1+ar2[i])*Grid[i+1][j]+az*(Grid[i][j+1]+Grid[i][j-1]))/(1-ar1+ar2[i]);
                if(i==Nrho-1&j!=0&&j!=Nz-1)Grid[i][j]=(ar1-ar2[i])*Grid[i-1][j]+az*(Grid[i][j+1]+Grid[i][j-1]);
                if(i!=0&&i!=Nrho-1&&j==0)Grid[i][j]=(ar1*(Grid[i+1][j]+Grid[i-1][j])+ar2[i]*(Grid[i+1][j]-Grid[i-1][j])+az*Grid[i][j+1])/(1-az);
                if(i!=0&&i!=Nrho-1&&j==Nz-1)Grid[i][j]=ar1*(Grid[i+1][j]+Grid[i-1][j])+ar2[i]*(Grid[i+1][j]-Grid[i-1][j])+az*Grid[i][j-1];
                //Rest
                if(i!=0&&i!=Nrho-1&&j!=0&&j!=Nz-1)Grid[i][j]=ar1*(Grid[i+1][j]+Grid[i-1][j])+ar2[i]*(Grid[i+1][j]-Grid[i-1][j])+az*(Grid[i][j+1]+Grid[i][j-1]);*/
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
    }
    
    //ofstream out;
    out.open("Test.dat");
    
    for(i=0;i<Nrho;i++)
        for(j=0;j<Nz;j++)
            out<<j<<"\t"<<i<<"\t"<<Grid[i][j]<<"\n";
    for(i=0;i<Nrho;i++)
        for(j=0;j<Nz;j++)
            out<<j<<"\t"<<-i<<"\t"<<Grid[i][j]<<"\n";
    out.close();
    
     
    //Derivative at 0,0 [0,150]
    int j0=(int)(zmax/(2.*dz)+0.5);
	double E=(-Grid[0][j0+2]+8.0*Grid[0][j0+1]-8*Grid[0][j0-1]+Grid[0][j0-2])/(12.0*dz);
    E*=1000;
    
    double du;
    double At = 187;//40;
    double Z = 75;//20;
    double mass = abs((Z - 1))*9.10938188*1e-31 + Z*1.6726231*1e-27 + (At - Z)*1.67492728*1e-27;
    double T=300.;
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
    system("pause");
    
    return 0;
}

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

void WriteRes (double** Grid, double dr, double dz, string name,int Nrho,int Nz);
void PlotRes (string name);

int main()
{
    int i,j,k;
    string name="ringV2";
    double dr=0.01;
    double dz=0.01;
    int Nrho=300;
    int Nz=300;
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
    
    double Ve=10.;
    double Vc=10.;
    double Vr=10.;
    double Vi=10.;//Go big to avoid falldown before reaching the centre
    double Vi2=0.;
    double Vr2=0.;
    double Vc2=0.;
    double Ve2=0.;
    
    /*
    #define li1 1//Distancia electrodo central centro trampa
#define li2 5.8
#define hi  12.4//L en el TFG

#define lr1 6
#define lr2 8
#define hr  9

#define wc 2.0//Anchura electrodo
#define dc 12.5//Posición centro electrodo
#define hc 5

#define we 1.5
#define de 17
#define he 6.5*/
    int L=124;
    
    int li=10;
    int li2=58;
    int hi=150-L;
    
    int lr=60;
    int lr2=80;
    int hr=150-90;
    
    int lc=125-20;
    int lc2=125+20;
    int hc=150-50;
    
    int le=170-15;
    int le2=170+15;
    int he=150-65;
    
    for(i=0;i<Nrho;i++)
    for(j=0;j<Nz;j++)     
    {
         if(i>=li&&i<li2&&j<hi){Grid[i][j]=Vi ;Electrode[i][j]=true;}
         if(i>=lr&&i<lr2&&j<hr){Grid[i][j]=Vr ;Electrode[i][j]=true;}
         if(i>=lc&&i<lc2&&j<hc){Grid[i][j]=Vc ;Electrode[i][j]=true;}
         if(i>=le&&i<le2&&j<he){Grid[i][j]=Ve ;Electrode[i][j]=true;}
         
         if(i>=li&&i<li2&&abs(j-Nz)<hi){Grid[i][j]=Vi2 ;Electrode[i][j]=true;}
         if(i>=lr&&i<lr2&&abs(j-Nz)<hr){Grid[i][j]=Vr2 ;Electrode[i][j]=true;}
         if(i>=lc&&i<lc2&&abs(j-Nz)<hc){Grid[i][j]=Vc2 ;Electrode[i][j]=true;}
         if(i>=le&&i<le2&&abs(j-Nz)<he){Grid[i][j]=Ve2 ;Electrode[i][j]=true;}
         
    }    
    
    ofstream out;
    out.open("Test-ring.dat");
    
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
    	WriteRes(Grid,dr,dz,name,Nrho,Nz);
	PlotRes(name);
    //out;
    out.open("Test.dat");
    
    for(i=0;i<Nrho;i++)
        for(j=0;j<Nz;j++)
            out<<j<<"\t"<<i<<"\t"<<Grid[i][j]<<"\n";
    for(i=0;i<Nrho;i++)
        for(j=0;j<Nz;j++)
            out<<j<<"\t"<<-i<<"\t"<<Grid[i][j]<<"\n";
    out.close();
    
     
    //Derivative at 0,0 [0,150]
    
    double E=(-Grid[0][152]+8.0*Grid[0][151]-8*Grid[0][149]+Grid[0][148])/(12.0*dz);
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
    du=A*E*q/(1.0*C*Ve);
    At = 40;
    Z = 20;
    mass = abs((Z - 1))*9.10938188*1e-31 + Z*1.6726231*1e-27 + (At - Z)*1.67492728*1e-27;
    
    double b=(q*E*du)/(2.0*mass*2*fz*acos(-1.0)*Ve);
    
    cout<<"Du: \t"<<du<<endl;
    cout<<"b: \t"<<b<<endl;
    system("pause");
    
    return 0;
}

void WriteRes (double** Grid, double dr, double dz, string name,int Nrho,int Nz)
{
	ofstream fileres, filez, filer;
	
	fileres.open((name+"-res.dat").c_str());
   	 for(int i=0;i<Nrho;i++)
   	 {
    		for(int j=0;j<Nz;j++)
    		fileres << i*dr << "\t" << j*dz << "\t" << Grid[i][j] << endl;
   		fileres << endl;
   	 }
	fileres.close();
	
	filez.open((name+"-axisz.dat").c_str());
	for(int j=0;j<Nz;j++)
		filez << j*dz << "\t" << Grid[0][j] << endl;
	filez.close();
	
	filer.open((name+"-axisr.dat").c_str());
	for(int i=0;i<Nrho;i++)
		filer << i*dr << "\t" << Grid[i][0] << endl;
	filer.close();
	
	return;
}
void PlotRes (string name)
{
	FILE *gp1;
	
	/*Contour plot: pantalla*/
	gp1=popen(GNUPLOT_PATH,"w");
	fprintf(gp1,"unset key\n");
	fprintf(gp1,"set term wxt\n");
	fprintf(gp1,"set size ratio -1\n");
	fprintf(gp1,"set output\n");
	fprintf(gp1,"set title 'Potencial: resultados'\n");
	fprintf(gp1,"set xlabel 'r'\n");
	fprintf(gp1,"set ylabel 'z'\n");
	fprintf(gp1,"set zlabel 'V'\n");
	fprintf(gp1,"set contour base\n");
	fprintf(gp1,"set pm3d depthorder\n");
	fprintf(gp1,"set view map\n");
	fprintf(gp1,"set cntrparam levels 50\n");
	fprintf(gp1,"splot '");
	fprintf(gp1,(name+"-res.dat").c_str());
	fprintf(gp1,"' w pm3d linewidth 1\n");
	fflush(gp1);
	pclose(gp1);
	/*Contour plot: png*/
	gp1=popen(GNUPLOT_PATH,"w");
	fprintf(gp1,"unset key\n");
	fprintf(gp1,"set term pngcairo size 2200,2000\n");
	fprintf(gp1,"set size ratio -1\n");
	fprintf(gp1,"set output '");
	fprintf(gp1,(name+"-res.png").c_str());
	fprintf(gp1,"'\n");
	fprintf(gp1,"set termoption font ',35'\n");
	fprintf(gp1,"set title 'Potencial: resultados' font ',40'\n");
	fprintf(gp1,"set xlabel 'r' font ',35'\n");
	fprintf(gp1,"set ylabel 'z' font ',35'\n");
	fprintf(gp1,"set zlabel 'V' font ',35'\n");
	fprintf(gp1,"set tics font ',30'\n");
	fprintf(gp1,"set contour base\n");
	fprintf(gp1,"set pm3d depthorder\n");
	fprintf(gp1,"set view map\n");
	fprintf(gp1,"set cntrparam levels 50\n");
	fprintf(gp1,"splot '");
	fprintf(gp1,(name+"-res.dat").c_str());
	fprintf(gp1,"' w pm3d linewidth 1\n");
	fflush(gp1);
	pclose(gp1);
	
	/*Eje z: pantalla*/
	gp1=popen(GNUPLOT_PATH,"w");
	fprintf(gp1,"unset key\n");
	fprintf(gp1,"set term wxt\n");
	fprintf(gp1,"set output\n");
	fprintf(gp1,"set title 'Potencial: eje z'\n");
	fprintf(gp1,"set xlabel 'z'\n");
	fprintf(gp1,"set ylabel 'V'\n");
	fprintf(gp1,"plot '");
	fprintf(gp1,(name+"-axisz.dat").c_str());
	fprintf(gp1,"' w lines\n");
	fflush(gp1);
	pclose(gp1);
	/*Eje z: png*/
	gp1=popen(GNUPLOT_PATH,"w");
	fprintf(gp1,"unset key\n");
	fprintf(gp1,"set term pngcairo size 2200,2000\n");
	fprintf(gp1,"set output '");
	fprintf(gp1,(name+"-axisz.png").c_str());
	fprintf(gp1,"'\n");
	fprintf(gp1,"set termoption font ',35'\n");
	fprintf(gp1,"set title 'Potencial: eje z' font ',40'\n");
	fprintf(gp1,"set xlabel 'z' font ',35'\n");
	fprintf(gp1,"set ylabel 'V' font ',35' offset 1,0\n");
	fprintf(gp1,"set tics font ',30'\n");
	fprintf(gp1,"plot '");
	fprintf(gp1,(name+"-axisz.dat").c_str());
	fprintf(gp1,"' w lines lw 3\n");
	fflush(gp1);
	pclose(gp1);
	
	/*Eje r: pantalla*/
	gp1=popen(GNUPLOT_PATH,"w");
	fprintf(gp1,"unset key\n");
	fprintf(gp1,"set term wxt\n");
	fprintf(gp1,"set output\n");
	fprintf(gp1,"set title 'Potencial: eje r'\n");
	fprintf(gp1,"set xlabel 'r'\n");
	fprintf(gp1,"set ylabel 'V'\n");
	fprintf(gp1,"plot '");
	fprintf(gp1,(name+"-axisr.dat").c_str());
	fprintf(gp1,"' w lines\n");
	fflush(gp1);
	pclose(gp1);
	/*Eje r: png*/
	gp1=popen(GNUPLOT_PATH,"w");
	fprintf(gp1,"unset key\n");
	fprintf(gp1,"set term pngcairo size 2200,2000\n");
	fprintf(gp1,"set output '");
	fprintf(gp1,(name+"-axisr.png").c_str());
	fprintf(gp1,"'\n");
	fprintf(gp1,"set termoption font ',35'\n");
	fprintf(gp1,"set title 'Potencial: eje r' font ',40'\n");
	fprintf(gp1,"set xlabel 'r' font ',35'\n");
	fprintf(gp1,"set ylabel 'V' font ',35' offset 1,0\n");
	fprintf(gp1,"set tics font ',30'\n");
	fprintf(gp1,"plot '");
	fprintf(gp1,(name+"-axisr.dat").c_str());
	fprintf(gp1,"' w lines lw 3\n");
	fflush(gp1);
	pclose(gp1);
	
	return;
}

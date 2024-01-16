#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <fstream>
using namespace std;

#define Nxx 400
#define Nyy 400
#define itmax 20000
#define delta 1e-3
const double pi=acos(-1);

#define GNUPLOT_PATH "/usr/bin/gnuplot -persist"
double Maximo (double d1, double d2);
void WriteInit (double** Grid, double dx, double dy, string name);
void WriteRes (double** Grid, double dx, double dy, string name);
void PlotInit (string name);
void PlotRes (string name);


double g1 (double x, double y)
{
  	double g1;
 	g1=sinh(3*pi*x)*sin(3*pi*y)/10;
  	return g1;
}
double g2 (double x, double y)
{
  	double g2;
 	g2=4*x*y*(x*cos(2*x*y)-y*sin(2*x*y));
  	return g2;
}
double g3 (double x, double y)
{
  	double g3;
 	g3=(x*cos(2*x*y)-y*sin(2*x*y))*exp(x*x-y*y);
  	return g3;
}
double g4 (double x, double y)
{
  	double g4;
 	g4=(x*cos(2*x*y)+y*sin(2*x*y))*exp(-x*x+y*y);
  	return g4;
}

int main(void)
{
	string name="name";
	int i,j,iter;
	int Nx=Nxx;
	int Ny=Nyy;
	double x,y,aux,max,xi;
	double dx=1./(Nx-1); 
	double dy=1./(Ny-1);	
	//Mem alloc e inicialización del grid y el bool auxiliar
	double **Grid=(double**)malloc(Nx*sizeof(double*));
	for(i=0;i<Nx;i++)
	{
		Grid[i]=(double*)malloc(Ny*sizeof(double));
		for(j=0;j<Ny;j++)
			Grid[i][j]=0;
	}
	bool **Electrode=(bool**)malloc(Nx*sizeof(bool*));
	for(i=0;i<Nx;i++)
	{
		Electrode[i]=(bool*)malloc(Ny*sizeof(bool));
		for(j=0;j<Ny;j++)
			Electrode[i][j]=false;
	}

	//Definimos condiciones de contorno************************
	j=0;
	for(i=0;i<Nx;i++)
	{
		x=i*dx;
		y=j*dy;
		Electrode[i][j]=true;
		Grid[i][j]=g2(x,y);
	}
	j=Ny-1;
	for(i=0;i<Nx;i++)
	{
		x=i*dx;
		y=j*dy;
		Electrode[i][j]=true;
		Grid[i][j]=g2(x,y);
	}
	i=0;
	for(j=0;j<Ny;j++) 
	{
		x=i*dx;
		y=j*dy;
		Electrode[i][j]=true;
		Grid[i][j]=g2(x,y);
	}
	i=Nx-1;
	for(j=0;j<Ny;j++)
	{
		x=i*dx;
		y=j*dy;
		Electrode[i][j]=true;
		Grid[i][j]=g2(x,y);
	}
	WriteInit(Grid,dx,dy,name);
	PlotInit(name);
	//Empieza Gauss Seidel*******************************************
	do{ 
		iter++; 
		if(iter%1000==0) cout << iter << " iteraciones calculadas." << endl;
		max=0.;
		for(i=0;i<Nx;i++)
			for(j=0;j<Ny;j++)
				if(!Electrode[i][j])
				{
					aux=Grid[i][j];
					Grid[i][j]=(Grid[i+1][j]+Grid[i-1][j]+Grid[i][j+1]+Grid[i][j-1])/4;
					if(Grid[i][j]!=0) //Si no se hace la comprobacion se puede dividir por cero (error)
					{
						xi=abs((Grid[i][j]-aux)/Grid[i][j]);
						max=Maximo(max,xi);
					}
				}
	}while((max>delta)&&(iter<itmax));
	
	if(max<=delta)
		cout << "Convergencia alcanzada a las " << iter << " iteraciones." << endl << endl;
	else
		cout << "Convergencia no alcanzada en " << iter << " iteraciones." << endl << endl;
		
	WriteRes(Grid,dx,dy,name);
	PlotRes(name);
	//escribimos resultados en fichero
	ofstream out;
	out.open("Test.dat");
	
	for(i=0;i<Nx;i++)
		for(j=0;j<Ny;j++)
			out<< i*dx << "\t" << j*dy << "\t" << Grid[i][j] << "\n" ;
	out.close();
	
	return 0; //Acaba el programa********************************************************************
}



/*Esta funcion escribe en un fichero las condiciones iniciales del potencial en
el formato requerido por gnuplot para hacer un contour plot.*/
void WriteInit (double** Grid,double dx, double dy, string name)
{
	int Nx=Nxx, Ny=Nyy;
	ofstream fileinit;
	
	fileinit.open((name+"-init.dat").c_str());
    	for(int i=0;i<Nx;i++){
    		for(int j=0;j<Ny;j++)
    		{
    		fileinit << i*dx << "\t" << j*dy << "\t" << Grid[i][j] << endl;}
   		fileinit << endl;
    		
    	}
	fileinit.close();
	
	return;
}

/*Esta funcion escribe en un fichero los resultados en el formato requerido por
gnuplot para hacer un contour plot, asi como el potencial a lo largo de los ejes.*/
void WriteRes (double** Grid, double dx, double dy, string name)
{
	int Nx=Nxx, Ny=Nyy;
	ofstream fileres, filez, filer;
	
	fileres.open((name+"-res.dat").c_str());
   	 for(int i=0;i<Nx;i++)
   	 {
    		for(int j=0;j<Ny;j++)
    		fileres << i*dx << "\t" << j*dy << "\t" << Grid[i][j] << endl;
   		fileres << endl;
   	 }
	fileres.close();
	
	filez.open((name+"-axisz.dat").c_str());
	for(int j=0;j<Ny;j++)
		filez << j*dy << "\t" << Grid[0][j] << endl;
	filez.close();
	
	filer.open((name+"-axisr.dat").c_str());
	for(int i=0;i<Nx;i++)
		filer << i*dx << "\t" << Grid[i][0] << endl;
	filer.close();
	
	return;
}

void PlotInit (string name)
{
	FILE *gp1;
	
	/*Contour plot: pantalla*/
	gp1=popen(GNUPLOT_PATH,"w");
	fprintf(gp1,"unset key\n");
	fprintf(gp1,"set term wxt\n");
	fprintf(gp1,"set size square\n");
	fprintf(gp1,"set output\n");
	fprintf(gp1,"set title 'Potencial: condiciones iniciales'\n");
	fprintf(gp1,"set xlabel 'x'\n");
	fprintf(gp1,"set ylabel 'y'\n");
	fprintf(gp1,"set zlabel 'V'\n");
	fprintf(gp1,"set contour base\n");
	fprintf(gp1,"set pm3d depthorder\n");
	fprintf(gp1,"set view map\n");
	fprintf(gp1,"unset contours\n");
	fprintf(gp1,"splot '");
	fprintf(gp1,(name+"-init.dat").c_str());
	fprintf(gp1,"' w pm3d linewidth 1\n");
	fflush(gp1);
	pclose(gp1);
	/*Contour plot: png*/
	gp1=popen(GNUPLOT_PATH,"w");
	fprintf(gp1,"unset key\n");
	fprintf(gp1,"set term pngcairo size 2200,2000\n");
	fprintf(gp1,"set size ratio -1\n");
	fprintf(gp1,"set output '");
	fprintf(gp1,(name+"-init.png").c_str());
	fprintf(gp1,"'\n");
	fprintf(gp1,"set termoption font ',35'\n");
	fprintf(gp1,"set title 'Potencial: condiciones iniciales' font ',40'\n");
	fprintf(gp1,"set xlabel 'x' font ',35'\n");
	fprintf(gp1,"set ylabel 'y' font ',35'\n");
	fprintf(gp1,"set zlabel 'V' font ',35'\n");
	fprintf(gp1,"set tics font ',30'\n");
	fprintf(gp1,"set contour base\n");
	fprintf(gp1,"set pm3d depthorder\n");
	fprintf(gp1,"set view map\n");
	fprintf(gp1,"unset contours\n");
	fprintf(gp1,"splot '");
	fprintf(gp1,(name+"-init.dat").c_str());
	fprintf(gp1,"' w pm3d linewidth 1\n");
	fflush(gp1);
	pclose(gp1);
		
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
	fprintf(gp1,"set xlabel 'x'\n");
	fprintf(gp1,"set ylabel 'y'\n");
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
	fprintf(gp1,"set xlabel 'x' font ',35'\n");
	fprintf(gp1,"set ylabel 'y' font ',35'\n");
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
	
	/*Eje y: pantalla*/
	gp1=popen(GNUPLOT_PATH,"w");
	fprintf(gp1,"unset key\n");
	fprintf(gp1,"set term wxt\n");
	fprintf(gp1,"set output\n");
	fprintf(gp1,"set title 'Potencial: eje y'\n");
	fprintf(gp1,"set xlabel 'y'\n");
	fprintf(gp1,"set ylabel 'V'\n");
	fprintf(gp1,"plot '");
	fprintf(gp1,(name+"-axisz.dat").c_str());
	fprintf(gp1,"' w lines\n");
	fflush(gp1);
	pclose(gp1);
	/*Eje y: png*/
	gp1=popen(GNUPLOT_PATH,"w");
	fprintf(gp1,"unset key\n");
	fprintf(gp1,"set term pngcairo size 2200,2000\n");
	fprintf(gp1,"set output '");
	fprintf(gp1,(name+"-axisz.png").c_str());
	fprintf(gp1,"'\n");
	fprintf(gp1,"set termoption font ',35'\n");
	fprintf(gp1,"set title 'Potencial: eje y' font ',40'\n");
	fprintf(gp1,"set xlabel 'y' font ',35'\n");
	fprintf(gp1,"set ylabel 'V' font ',35' offset 1,0\n");
	fprintf(gp1,"set tics font ',30'\n");
	fprintf(gp1,"plot '");
	fprintf(gp1,(name+"-axisz.dat").c_str());
	fprintf(gp1,"' w lines lw 3\n");
	fflush(gp1);
	pclose(gp1);
	
	/*Eje x: pantalla*/
	gp1=popen(GNUPLOT_PATH,"w");
	fprintf(gp1,"unset key\n");
	fprintf(gp1,"set term wxt\n");
	fprintf(gp1,"set output\n");
	fprintf(gp1,"set title 'Potencial: eje x'\n");
	fprintf(gp1,"set xlabel 'x'\n");
	fprintf(gp1,"set ylabel 'V'\n");
	fprintf(gp1,"plot '");
	fprintf(gp1,(name+"-axisr.dat").c_str());
	fprintf(gp1,"' w lines\n");
	fflush(gp1);
	pclose(gp1);
	/*Eje x: png*/
	gp1=popen(GNUPLOT_PATH,"w");
	fprintf(gp1,"unset key\n");
	fprintf(gp1,"set term pngcairo size 2200,2000\n");
	fprintf(gp1,"set output '");
	fprintf(gp1,(name+"-axisr.png").c_str());
	fprintf(gp1,"'\n");
	fprintf(gp1,"set termoption font ',35'\n");
	fprintf(gp1,"set title 'Potencial: eje x' font ',40'\n");
	fprintf(gp1,"set xlabel 'x' font ',35'\n");
	fprintf(gp1,"set ylabel 'V' font ',35' offset 1,0\n");
	fprintf(gp1,"set tics font ',30'\n");
	fprintf(gp1,"plot '");
	fprintf(gp1,(name+"-axisr.dat").c_str());
	fprintf(gp1,"' w lines lw 3\n");
	fflush(gp1);
	pclose(gp1);
	
	return;
}

double Maximo (double d1, double d2)
{
	if(d1>d2) d2=d1;
	return d2;
}

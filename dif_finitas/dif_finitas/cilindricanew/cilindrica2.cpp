#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

#define GNUPLOT_PATH "/usr/bin/gnuplot -persist"
void WriteInit (double** Grid, double dr, double dz, string name,int Nrho,int Nz);
void WriteRes (double** Grid, double dr, double dz, string name,int Nrho,int Nz);
void PlotInit (string name);
void PlotRes (string name);
//*****************************************************************************
void Fit();
void Parametros (double p[]);
double resta (double Vc1,double Vc2,double beta);
double suma (double Vc1,double Vc2,double beta);
//*****************************************************************************
//void cambio(bool sigue,int aux2,double Vc1,double Vc2,double beta);
double Maximo (double actual, double anterior,double max);

#define itmax 20000 //nº maximo de iteraciones del GAUSS_SEIDEL
#define delta 1e-3 //valor mínimo de convergencia (variación mínima) del GAUSS_SEIDEL

int main()
{
	int i,j,k,aux2;
	double a,b,c,d,deltaa,aux,beta,max=0.,actual=0.,anterior=0.;
	beta=0.2; // El paso en los potenciales
	aux2=-1; // Empezamos disminuyendo Vc
	deltaa=1e-6; // valor máximo de C4=p[3]
	a=b=c=d=0.;
	bool pC4=false;
	string name="cil";
	double p[4]={10.,10.,10.,10.}; //inicializamos a 10, para que sea mayor que deltaa y entre al algoritmo de corrección.
	
	//*********************************************************
	double rho0=6.0; //El parámetro que lo define todo
	double dr=0.01*rho0; //Debe ser ~0.01*rho0 o menor
	double dz=0.01*rho0; //Idem
	int Nrho=(int)(1.1*rho0/dr+0.5); //+0.5 porque al castear se trunca
	int Nz=(int)( ( (0.4+10.58/1.0239 )*rho0 ) / dz + 0.5 ); //Idem
	if(Nz%2==0) Nz++; //Para que Nz sea impar y haya siempre un centro
	int j0=Nz/2+1; //Centro en j0*dr
	
	int imin, imax, jmin, jmax;
	imin=(int)(rho0/dr+0.5);//Radio interno electrodos
	imax=(int)( (1+2*0.0303/1.0239)*rho0 /dr + 0.5 ); //Radio externo electrodos (espesor 2*gap asignado arbitrariamente)
	
	//Mem alloc e inicialización del grid (mallado de los potenciales) y el bool auxiliar - CHECKED
	double **Grid=(double**)malloc(Nrho*sizeof(double*));
	bool **Electrode=(bool**)malloc(Nrho*sizeof(bool*));
	for(i=0;i<Nrho;i++)
	{
		Grid[i]=(double*)malloc(Nz*sizeof(double));
		for(j=0;j<Nz;j++)
			Grid[i][j]=0.;
	}
	for(i=0;i<Nrho;i++)
	{
		Electrode[i]=(bool*)malloc(Nz*sizeof(bool));
		for(j=0;j<Nz;j++)
			Electrode[i][j]=false;
	}
	
	// damos valores a los potenciales en los electrodos. 'c' de corrección, 'e' de endcap y '0' de ring.
	double Ve1=0.;
	double Ve2=0.;
	double Vc1=-1.2;
	double Vc2=-1.2;
	double V0 =10.;
	//definimos las variables que aparecen de la resolución de la ecuación de Laplace en cilíndricas con el método de dif finitas
	double ar1=(dz*dz)/(2*(dr*dr+dz*dz)), az=ar1*(dr*dr)/(dz*dz);
	double *ar2=(double*)malloc(Nrho*sizeof(double));
	double a0r=(dz*dz)/(dr*dr+2.*dz*dz), a0z=0.5*(dr*dr)/(dz*dz)*a0r; //Nueva forma de saltar la singularidad (MODIF)
	//for(i=1;i<Nrho;i++) ar2[i]=ar1/(2.*i); //MODIF ***¿Qué pasa para i=0? 3.51759e+1800 eso pasa. esto se traduce

	for(i=0;i<Nrho;i++) ar2[i]=ar1/(2.*i+0.5); // El fallo anterior no era muy notable, ya que al hacer el gauss seidel tarde 
	//cout<<ar2[0];				// o temprano converge al mismo valor. Tan solo se ven modificados ligeramente los 
						// potenciales alrededor de r=0. Quizás se observan por no dejarlos converger completamente.
						
						
	/****************************************************************************************************************/
	//--------------COMIENZA ALGORITMO DE DEFINICIÓN DE LA GEOMETRÍA DE LA TRAMPA (ASIGNACIÓN DE POTENCIALES Y
	//--------------ELECTRODOS EN LAS POSICIONES DEL MALLADO), ALGORITMO DE GAUSS-SEIDEL Y ALGORITMO DE MODIFICACIÓN
	//--------------DE LOS POTENCIALES DE CORRECCIÓN
	/****************************************************************************************************************/
	
Aqui:	//etiqueta del GOTO utilizado para modificar Vc
	for(i=0;i<Nrho;i++)
	{
		Grid[i]=(double*)malloc(Nz*sizeof(double));
		for(j=0;j<Nz;j++)
			Grid[i][j]=0.;
	}
	
	for(i=0;i<Nrho;i++)
	{
		Electrode[i]=(bool*)malloc(Nz*sizeof(bool));
		for(j=0;j<Nz;j++)
			Electrode[i][j]=false;
	}
	
	Vc1=Vc2; //para ahorrar operaciones en caso de que se quieran modificar los potenciales simétricamente.
	
	//Inicialización
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
		
		//Ring
		jmin=(int)( ( (0.2+5.2227/1.0239 )*rho0 ) / dr + 0.5 ); //0.2rho+ze+zg+zc+zg
		jmax=(int)( ( (0.2+5.3573/1.0239 )*rho0 ) / dr + 0.5 ); //0.2rho+ze+zg+zc+zg+zr
		for(j=jmin;j<=jmax;j++)
		{
			Electrode[i][j]=true;
			Grid[i][j]=V0;
		}
		
		//E1 (z positivos)
		jmin=(int)( ( (0.2+6.2530/1.0239 )*rho0 ) / dr + 0.5 ); //0.2rho+ze+zg+zc+zg+zr+zg+zc+zg
		jmax=(int)( ( (0.2+10.5800/1.0239 )*rho0 ) / dr + 0.5 ); //0.2rho+ze+zg+zc+zg+zr+zg+zc+zg+ze
		for(j=jmin;j<=jmax;j++)
		{
			Electrode[i][j]=true;
			Grid[i][j]=Ve1;
		}
		//C2 (z negativos)
		jmin=(int)( ( ( 0.2+4.3573/1.0239 )*rho0 ) / dr + 0.5 );  //0.2rho+ze+zg
		jmax=(int)( ( ( 0.2+5.1924/1.0239 )*rho0 ) / dr + 0.5 );  //0.2rho+ze+zg+zc
		for(j=jmin;j<=jmax;j++)
		{
			Electrode[i][j]=true;
			Grid[i][j]=Vc2;
		}
		
		//C1 (z positivos)
		jmin=(int)( ( (0.2+5.3876/1.0239 )*rho0 ) / dr + 0.5 ); //0.2rho+ze+zg+zc+zg+zr+zg
		jmax=(int)( ( (0.2+6.2227/1.0239 )*rho0 ) / dr + 0.5 ); //0.2rho+ze+zg+zc+zg+zr+zg+zc
		for(j=jmin;j<=jmax;j++)
		{
			Electrode[i][j]=true;
			Grid[i][j]=Vc1;
		}
	}
	WriteInit(Grid,dr,dz,name,Nrho,Nz);
	PlotInit(name);
	//Termina la inicialización.
	//COMIENZA GAUSS-SEIDEL	
	k=0;
	do{
		k++;
		max=0.;
		if(k%1000==0)cout<<k<<"\t iteraciones"<<endl; //contador cada 1000 pasos
		for(i=0;i<Nrho;i++)
			for(j=0;j<Nz;j++)
				if(!Electrode[i][j])
				{
					//Corners
					if(i==0&&j==0){
						anterior=Grid[i][j];
						Grid[i][j]=2*a0r*Grid[i+1][j]+a0z*Grid[i][j+1];
						actual=Grid[i][j];
						max=Maximo(actual,anterior,max);
						}
					if(i==0&&j==Nz-1){
						anterior=Grid[i][j];
						Grid[i][j]=2*a0r*Grid[i+1][j]+a0z*Grid[i][j-1];
						actual=Grid[i][j];
						max=Maximo(actual,anterior,max);} 
					if(i==Nrho-1&&j==0){
						anterior=Grid[i][j];
						Grid[i][j]=((ar1-ar2[i])*Grid[i-1][j]+az*Grid[i][j+1])/(1-az);
						actual=Grid[i][j];
						max=Maximo(actual,anterior,max);}
					if(i==Nrho-1&&j==Nz-1){
						anterior=Grid[i][j];
						Grid[i][j]=(ar1-ar2[i])*Grid[i-1][j]+az*Grid[i][j];
						actual=Grid[i][j];
						max=Maximo(actual,anterior,max);}
					//Sides
					if(i==0&j!=0&&j!=Nz-1){
						anterior=Grid[i][j];
						Grid[i][j]=2*a0r*Grid[i+1][j]+a0z*(Grid[i][j+1]+Grid[i][j-1]);
						actual=Grid[i][j];
						max=Maximo(actual,anterior,max);} 
					if(i==Nrho-1&j!=0&&j!=Nz-1){
						anterior=Grid[i][j];
						Grid[i][j]=(ar1-ar2[i])*Grid[i-1][j]+az*(Grid[i][j+1]+Grid[i][j-1]);
						actual=Grid[i][j];
						max=Maximo(actual,anterior,max);}
					if(i!=0&&i!=Nrho-1&&j==0){
						anterior=Grid[i][j];
						Grid[i][j]=ar1*(Grid[i+1][j]+Grid[i-1][j])+ar2[i]*(Grid[i+1][j]-Grid[i-1][j])+az*Grid[i][j+1];
						actual=Grid[i][j];
						max=Maximo(actual,anterior,max);}
					if(i!=0&&i!=Nrho-1&&j==Nz-1){
						anterior=Grid[i][j];
						Grid[i][j]=ar1*(Grid[i+1][j]+Grid[i-1][j])+ar2[i]*(Grid[i+1][j]-Grid[i-1][j])+az*Grid[i][j-1];
						actual=Grid[i][j];
						max=Maximo(actual,anterior,max);}
					//Rest
					if(i!=0&&i!=Nrho-1&&j!=0&&j!=Nz-1){
						anterior=Grid[i][j];
						Grid[i][j]=ar1*(Grid[i+1][j]+Grid[i-1][j])+ar2[i]*(Grid[i+1][j]-Grid[i-1][j])+az*(Grid[i][j+1]+Grid[i][j-1]);
						actual=Grid[i][j];
						max=Maximo(actual,anterior,max);}
				}
	}while((max>delta)&&(k<itmax));

	//escribimos resultados en fichero
	WriteRes(Grid,dr,dz,name,Nrho,Nz);
	PlotRes(name);
	ofstream out;
	out.open("Test.dat");
    
   	for(i=0;i<Nrho;i++){
        	for(j=0;j<Nz;j++){
            		out << i*dr << "\t" << (j-j0)*dz << "\t" << Grid[i][j] << "\n";}
            	out<<"\n";}
    	for(i=0;i<Nrho;i++){
        	for(j=0;j<Nz;j++){
            		out << -i*dr << "\t" << (j-j0)*dz << "\t" << Grid[i][j] << "\n";}
            	out<<"\n";}
	out.close();
	
	//*********************************
	
	//sigue=false;
	aux=abs(p[2]);
	Fit();
	Parametros(p);
	if (abs(p[2])<deltaa) pC4=true;
	if (!pC4){
		if(abs(p[2])<aux){
			//sigue=true;
			if(aux2==1){//El Vc1 no lo tocamos para ahorrar cuentas, ya que Vc1=Vc2, a menos que los quieras cambiar independientemente el uno del otro (perderías simetría).
				//Vc1=suma(Vc1,Vc2,beta);
				Vc2=suma(Vc1,Vc2,beta);}
			else if(aux2==-1){
				//Vc1=resta(Vc1,Vc2,beta);
				Vc2=resta(Vc1,Vc2,beta);}
			}
		else{
			if(aux2==1){
				//Vc1=resta(Vc1,Vc2,beta);
				Vc2=resta(Vc1,Vc2,beta);
				aux2=-1;}
			else if(aux2==-1){
				//Vc1=suma(Vc1,Vc2,beta);
				Vc2=suma(Vc1,Vc2,beta);
				aux2=1;}
			}
		goto Aqui;}
		
	delete Grid; //liberamos memoria
	delete Electrode; //liberamos memoria
	
	return 0;
}
/*Esta funcion escribe en un fichero las condiciones iniciales del potencial en
el formato requerido por gnuplot para hacer un contour plot.*/
void WriteInit (double** Grid,double dr, double dz, string name,int Nrho,int Nz)
{
	ofstream fileinit;
	
	fileinit.open((name+"-init.dat").c_str());
    	for(int i=0;i<Nrho;i++){
    		for(int j=0;j<Nz;j++)
    		{
    			fileinit << i*dr << "\t" << j*dz << "\t" << Grid[i][j] << "\n";
    		}
    	fileinit<<"\n";}
	fileinit.close();
	
	return;
}

/*Esta funcion escribe en un fichero los resultados en el formato requerido por
gnuplot para hacer un contour plot, asi como el potencial a lo largo de los ejes.*/
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

void PlotInit (string name)
{
	FILE *gp1;
	
	/*Contour plot: pantalla*/
	gp1=popen(GNUPLOT_PATH,"w");
	fprintf(gp1,"unset key\n");
	fprintf(gp1,"set term wxt\n");
	fprintf(gp1,"set size ratio -1\n");
	fprintf(gp1,"set output\n");
	fprintf(gp1,"set title 'Potencial: condiciones iniciales'\n");
	fprintf(gp1,"set xlabel 'r'\n");
	fprintf(gp1,"set ylabel 'z'\n");
	fprintf(gp1,"set zlabel 'V'\n");
	fprintf(gp1,"set contour base\n");
	fprintf(gp1,"set pm3d depthorder\n");
	fprintf(gp1,"set view map\n");
	fprintf(gp1,"unset contours\n");
	fprintf(gp1,"set cntrparam levels 50\n");
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
	fprintf(gp1,"set xlabel 'r' font ',35'\n");
	fprintf(gp1,"set ylabel 'z' font ',35'\n");
	fprintf(gp1,"set zlabel 'V' font ',35'\n");
	fprintf(gp1,"set tics font ',30'\n");
	fprintf(gp1,"set contour base\n");
	fprintf(gp1,"set pm3d depthorder\n");
	fprintf(gp1,"set view map\n");
	fprintf(gp1,"unset contours\n");
	fprintf(gp1,"set cntrparam levels 50\n");
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
void Fit ()
{
	FILE *gp1;
	gp1=popen(GNUPLOT_PATH,"w");
	fprintf(gp1,"set output\n");
	//y=z ; x=Rho
	fprintf(gp1,"f(x,y)=a+b*(y*y-x*x/2)+c*(y**4-3*x*x*y*y+3*x**4/8)+d*(y**6-15*y**4*x*x/2+45*x**4*y*y/8-5*x**6/16)\n");
	fprintf(gp1,"a=10\n");
	fprintf(gp1,"b=10\n");
	fprintf(gp1,"c=0.00001\n");
	fprintf(gp1,"d=0.00000001\n");
	fprintf(gp1,"fit f(x,y) 'Test.dat' using 1:2:3 via a,b,c,d\n");
	fprintf(gp1,"set print 'paramfit.dat' \n");
	fprintf(gp1,"print a,b,c,d\n");
	//Ya hemos hecho el ajuste, ahora pintamos.
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
	fprintf(gp1,"splot 'Test.dat' w pm3d linewidth 1, f(x,y)\n");
	fflush(gp1);
	pclose(gp1);
	/*Contour plot: png
	gp1=popen(GNUPLOT_PATH,"w");	
	fprintf(gp1,"unset key\n");
	fprintf(gp1,"set term pngcairo size 2200,2000\n");
	fprintf(gp1,"set size ratio -1\n");
	fprintf(gp1,"set output 'fit.png' \n");
	fprintf(gp1,"set title 'Potencial: resultados'\n");
	fprintf(gp1,"set xlabel 'r'\n");
	fprintf(gp1,"set ylabel 'z'\n");
	fprintf(gp1,"set zlabel 'V'\n");
	fprintf(gp1,"set contour base\n");
	fprintf(gp1,"set pm3d depthorder\n");
	fprintf(gp1,"set view map\n");
	fprintf(gp1,"set cntrparam levels 50\n");
	fprintf(gp1,"splot 'Test.dat' w pm3d linewidth 1, f(x,y)\n");	
	fflush(gp1);
	pclose(gp1);*/
	
	return;
}

/*Lee los parametros del ajuste y los almacena en el vector p[].p[2]=C4*/
void Parametros (double p[])
{
   	ifstream file;
	file.open("paramfit.dat");
	file >>p[0]>>p[1]>>p[2]>>p[3];
	file.close();
	
	return;
}
double suma (double Vc1,double Vc2,double beta)
{
	Vc1=Vc1+beta;
	Vc2=Vc2+beta;
	cout<<"Aumento\t"<<Vc2<<"\n";
	return Vc2;
}
double resta (double Vc1,double Vc2,double beta)
{
	Vc1=Vc1-beta;
	Vc2=Vc2-beta;
	cout<<"Disminuyo\t"<<Vc2<<"\n";
	return Vc2;
}
/*void cambio(bool sigue,int aux2,double Vc1,double Vc2,double beta)
{
	double auxxx;
	if(sigue){
		if(aux2==1){//suma
			auxxx=suma(Vc1,Vc2,beta);}
		else if(aux2==-1){
			auxxx=resta(Vc1,Vc2,beta);}
		}
	else{
		if(aux2==1){//suma
			auxxx=resta(Vc1,Vc2,beta);
			aux2=-1;}
		else if(aux2==-1){
			auxxx=suma(Vc1,Vc2,beta);
			aux2=1;}
		}
	return;
}*/
double Maximo (double actual, double anterior,double max)
{
	double epsilon; //error relativo
	if (actual!=0){
		epsilon=(actual-anterior)/actual;
		if(epsilon>max) max=epsilon;}
	
	return max;
}

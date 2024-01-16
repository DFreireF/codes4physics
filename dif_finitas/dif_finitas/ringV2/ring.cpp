#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

#define GNUPLOT_PATH "/usr/bin/gnuplot -persist"
void WriteRes (double** Grid, double dr, double dz, string name,int Nrho,int Nz);
void PlotRes (string name);
void PlotInit (string name);
void Fit();
void Parametros (double p[]);
double resta (double Vc,double Vc2,double beta);
double suma (double Vc,double Vc2,double beta);
double Maximo (double actual, double anterior,double max);
#define itmax 20000 //nº maximo de iteraciones del GAUSS_SEIDEL
#define delta 1e-4 //valor mínimo de convergencia (variación mínima) del GAUSS_SEIDEL


int main()
{
    int i,j,k;
    string name="ringV2";
    double dr=0.01;
    double dz=0.01;
    int Nrho=300;
    int Nz=300;
    double **Grid=(double**)malloc(Nrho*sizeof(double*));
    bool **Electrode=(bool**)malloc(Nrho*sizeof(bool*));
    
	int j0=Nz/2+1; //Centro en j0*dr
	double a,b,c,d,deltaa,aux,beta,max=0.,actual=0.,anterior=0.;
	int aux2=-1; // Empezamos --- Vc (+1 para sumar, -1 para restar)
	beta=4; // El paso en los potenciales
	deltaa=1e-6; // valor máximo de C4=p[3]
	a=b=c=d=0.;
	bool pC4=false;
	double p[4]={10.,10.,10.,10.};
 
    //*********************
    double ar1=(dz*dz)/(2*(dr*dr+dz*dz)), az=ar1*(dr*dr)/(dz*dz);
    double *ar2=(double*)malloc(Nrho*sizeof(double));
    double a0r=(dz*dz)/(dr*dr+2.*dz*dz), a0z=0.5*(dr*dr)/(dz*dz)*a0r; //Nueva forma de saltar la singularidad (MODIF)
    for(k=0;k<Nrho;k++) ar2[k]=ar1/(2.*k+0.5);
    //for(k=0;k<Nrho;k++) ar2[k]=ar1/(2.*k+1.);//Wrong //Yo creo que está igual de bien que con 0.5
    //************************
    
    double Ve=0;
    double Vr=6.5;
    double Vc=5.5;
    double Vi=13;
    double Vi2=13;
    double Vc2=5.5;
    double Vr2=6.5;
    double Ve2=0;
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
    
Aqui:
	Vc=Vc2; //esto lo hago para ahorrarme operaciones en el algoritmo. Modifico los 2 los mismo, no por separado.
	
    for(i=0;i<Nrho;i++)
    {
        Grid[i]=(double*)malloc(Nz*sizeof(double));
        for(j=0;j<Nz;j++)
        Grid[i][j]=0;
    }
    for(i=0;i<Nrho;i++)
    {
        Electrode[i]=(bool*)malloc(Nz*sizeof(bool));
        for(j=0;j<Nz;j++)
        Electrode[i][j]=false;
    }
    for(i=0;i<Nrho;i++){
    	for(j=0;j<Nz;j++){     
        	if(i>=li&&i<li2&&j<hi){Grid[i][j]=Vi; Electrode[i][j]=true;}
         	if(i>=lr&&i<lr2&&j<hr){Grid[i][j]=Vc; Electrode[i][j]=true;}
         	if(i>=lc&&i<lc2&&j<hc){Grid[i][j]=Vr; Electrode[i][j]=true;}
         	if(i>=le&&i<le2&&j<he){Grid[i][j]=Ve; Electrode[i][j]=true;}
         
         	if(i>=li&&i<li2&&abs(j-Nz)<hi){Grid[i][j]=Vi2; Electrode[i][j]=true;}
         	if(i>=lr&&i<lr2&&abs(j-Nz)<hr){Grid[i][j]=Vc2; Electrode[i][j]=true;}
         	if(i>=lc&&i<lc2&&abs(j-Nz)<hc){Grid[i][j]=Vr2; Electrode[i][j]=true;}
         	if(i>=le&&i<le2&&abs(j-Nz)<he){Grid[i][j]=Ve2; Electrode[i][j]=true;}
         }
    }    

    ofstream out;
    out.open("Test-ring.dat");
    
    for(i=0;i<Nrho;i++){
        for(j=0;j<Nz;j++){
            out << i*dr << "\t" << (j-j0)*dz << "\t" << Grid[i][j] << "\n";}
            out<<"\n";}
    for(i=0;i<Nrho;i++){
        for(j=0;j<Nz;j++){
            out << -i*dr << "\t" << (j-j0)*dz << "\t" << Grid[i][j] << "\n";}
            out<<"\n";}
        
    out.close();
    PlotInit(name);
    k=0;
do{
	k++;
	max=0.;    
	if(k%1000==0)cout<<k<<endl;
        for(i=0;i<Nrho;i++){
        	for(j=0;j<Nz;j++){
            	if(!Electrode[i][j]){
                	//Corners
			if(i==0&&j==0){
				anterior=Grid[i][j];
				Grid[i][j]=2*a0r*Grid[i+1][j]+a0z*Grid[i][j+1];
				actual=Grid[i][j];
				max=Maximo(actual,anterior,max);}
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
        	}
        }
    }while((max>delta)&&(k<itmax));
    	//WriteRes(Grid,dr,dz,name,Nrho,Nz);
	//PlotRes(name);
	ofstream out2;
	out2.open("Test.dat");
    
   	for(i=0;i<Nrho;i++){
        	for(j=0;j<Nz;j++){
            		out2 << i*dr << "\t" << (j-j0)*dz << "\t" << Grid[i][j] << "\n";}
            	out2<<"\n";}
    	for(i=0;i<Nrho;i++){
        	for(j=0;j<Nz;j++){
            		out2 << -i*dr << "\t" << (j-j0)*dz << "\t" << Grid[i][j] << "\n";}
            	out2<<"\n";}
	out2.close();
	
    // convergencia voltajes
	aux=abs(p[2]);
	Fit();
	Parametros(p);
	if (abs(p[2])<deltaa) pC4=true;
	if (!pC4){
		if(abs(p[2])<aux){
			if(aux2==1){//suma
				//Vc=suma(Vc,Vc2,beta);
				Vc2=suma(Vc,Vc2,beta);}
			else if(aux2==-1){
				//Vc=resta(Vc,Vc2,beta);
				Vc2=resta(Vc,Vc2,beta);}
			}
		else{
			if(aux2==1){//suma
				//Vc=resta(Vc,Vc2,beta);
				Vc2=resta(Vc,Vc2,beta);
				aux2=-1;}
			else if(aux2==-1){
				//Vc=suma(Vc,Vc2,beta);
				Vc2=suma(Vc,Vc2,beta);
				aux2=1;}
			}
		goto Aqui;}
    
    return 0;
}

void WriteRes (double** Grid, double dr, double dz, string name,int Nrho,int Nz)
{
	ofstream fileres, filez, filer;
	
	fileres.open((name+"-res.dat").c_str());
   	 for(int i=0;i<Nrho;i++)
   	 {
    		for(int j=0;j<Nz;j++){
    			fileres << i*dr << "\t" << j*dz << "\t" << Grid[i][j] << endl;}
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
void PlotInit (string name)
{
	FILE *gp1;
	/*Contour plot: pantalla*/
	gp1=popen(GNUPLOT_PATH,"w");
	fprintf(gp1,"unset key\n");
	fprintf(gp1,"set term wxt\n");
	fprintf(gp1,"set size ratio -1\n");
	//fprintf(gp1,"set output\n");
	fprintf(gp1,"set title 'Potencial: condiciones iniciales'\n");
	fprintf(gp1,"set xlabel 'r'\n");
	fprintf(gp1,"set ylabel 'z'\n");
	fprintf(gp1,"set zlabel 'V'\n");
	fprintf(gp1,"set contour base\n");
	fprintf(gp1,"set pm3d depthorder\n");
	fprintf(gp1,"set view map\n");
	fprintf(gp1,"unset contours\n");
	fprintf(gp1,"set cntrparam levels 50\n");
	fprintf(gp1,"splot 'Test-ring.dat' using 1:2:3 w pm3d linewidth 1\n");
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
	fprintf(gp1,"splot 'Test-ring.dat' w pm3d linewidth 1\n");
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
	fprintf(gp1,"a=0.1\n");
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

/*Lee los parametros del ajuste.*/
void Parametros (double p[])
{
   	ifstream file;
	file.open("paramfit.dat");
	file >>p[0]>>p[1]>>p[2]>>p[3];
	file.close();	
	return;
}
double suma (double Vc,double Vc2,double beta)
{
	Vc=Vc+beta;
	Vc2=Vc2+beta;
	cout<<"Aumento a\t"<<Vc2<<"\n";
	return Vc2;
}
double resta (double Vc,double Vc2,double beta)
{
	Vc=Vc-beta;
	Vc2=Vc2-beta;
	cout<<"Disminuyo a\t"<<Vc2<<"\n";
	return Vc2;
}
double Maximo (double actual, double anterior,double max)
{
	double epsilon; //error relativo
	if (actual!=0){
		epsilon=(actual-anterior)/actual;
		if(epsilon>max) max=epsilon;}
	
	return max;
}

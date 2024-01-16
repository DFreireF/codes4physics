#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#define NELECMAX 5
#define MAXR 10000
#define MAXZ 10000
//#define GNUPLOT_PATH "D:\\gnuplot\\bin\\gnuplot.exe -persist"
#define GNUPLOT_PATH "/usr/bin/gnuplot -persist"
using namespace std;
bool ReadPar (double rint[], double rext[], double z[], double velec[], int& nelec, double& dr, double& dz, 
				double& delta, string& name, double& rmax, double& zmax);
void InitChart (double rint[], double rext[], double z[], double velec[], int nelec, bool electrode[][MAXZ], 
                double V[][MAXZ], double dr, double dz, double rmax, double zmax);
void Solve (bool electrode[][MAXZ], double V[][MAXZ], double a2[], double dr, double dz, double rmax, 
				double zmax, int itmax, double delta);
void WriteInit (double V[][MAXZ], double rmax, double zmax, double dr, double dz, string name);
void WriteRes (double V[][MAXZ], double rmax, double zmax, double dr, double dz, string name);
void PlotInit (string name);
void PlotRes (string name);
double Maximo (double d1, double d2);

double V[MAXR][MAXZ];
bool electrode[MAXR][MAXZ];
double ar2[MAXR];
/*Se declaran como variables globales para evitar desbordar la pila
(solo es critico en el caso de V, que ocupa 8 bits por elemento)*/

int main (void)
{
	double rint[NELECMAX], rext[NELECMAX], z[NELECMAX], velec[NELECMAX];
	bool read;
	int nelec, Nr, Nz, itmax=300000;
	double dr, dz, rmax, zmax, delta;
	string name;
	
	/*rint es el radio interior del electrodo, rext su radio exterior,
	z la distancia desde su borde al plano z=0 y velec su potencial.*/
    
    cout << "Leyendo fichero..." << endl;
	read=ReadPar(rint,rext,z,velec,nelec,dr,dz,delta,name,rmax,zmax);
	
	if(read)
	{
		cout << "Lectura de fichero finalizada. " << endl << endl;
        cout << "Inicializando matriz de potencial..." << endl; 
		InitChart(rint,rext,z,velec,nelec,electrode,V,dr,dz,rmax,zmax);
		cout << "Matriz de potencial inicializada." << endl << endl;
		WriteInit(V,rmax,zmax,dr,dz,name);
		PlotInit(name);
		cout << "Iniciando calculos..." << endl << endl;
		Solve(electrode,V,ar2,dr,dz,rmax,zmax,itmax,delta);
		WriteRes(V,rmax,zmax,dr,dz,name);
		PlotRes(name);
	}
	
	else
        cout << "Error de lectura." << endl;
	
	system("pause");
	
	return 0;
}

/*Lee los parametros desde el fichero por defecto. Devuelve true si se adquieren los
parametros correctamente y false si no.*/
bool ReadPar (double rint[], double rext[], double z[], double velec[], int& nelec, double& dr, double& dz, double& delta, string& name, double& rmax, double& zmax)
{
    ifstream file;
    string aux;
    bool read=false;
    int i=0;
        
	file.open("input.inp");
	
	if(file.is_open())
	{
		read=true;
		file >> name;
		getline(file,aux); //getline se usa para limpiar lineas de comentario de los ficheros
		file >> dr; 
		getline(file,aux);
		file >> dz;
		getline(file,aux);
		file >> rmax; 
		getline(file,aux);
		file >> zmax;
		getline(file,aux);
		file >> delta;
		getline(file,aux);
		file >> nelec;
		getline(file,aux);
		
		
		if(nelec<1||nelec>NELECMAX) read=false;
		
		else
		{
            do
            {
				getline(file,aux);
				file >> rint[i];
				getline(file,aux);			
				file >> rext[i];
				getline(file,aux);
				file >> z[i];
				getline(file,aux);
				file >> velec[i];
				getline(file,aux);
				i++;
            }while(i<nelec&&!file.eof());
		}
		if(file.eof()) read=false;
		file.close();
	}
	
	return read;
}

void InitChart (double rint[], double rext[], double z[], double velec[], int nelec, bool electrode[][MAXZ], 
               double V[][MAXZ], double dr, double dz, double rmax, double zmax)
{
	int imin, imax, jmin, Nr=(int)(rmax/dr-0.5), Nz=(int)(zmax/dz-0.5); //Se suma 0.5 a la expresion exacta (-1) porque (int) trunca, no redondea
    
	//Inicializamos todo a cero...
    for(int i=0;i<=Nr;i++)
		for(int j=0;j<=Nz;j++)
		{
			electrode[i][j]=false;
			V[i][j]=0.;
		}
	
	//... y luego fijamos los electrodos
    for(int k=0;k<nelec;k++)
    {
    	imin=(int)(rint[k]/dr);
    	imax=(int)(rext[k]/dr);
    	jmin=(int)(z[k]/dz);
    	for(int i=imin;i<=imax;i++)
    		for(int j=jmin;j<=Nz;j++)
    		{
    			electrode[i][j]=true;
                V[i][j]=velec[k];
    		}
    }
    
    return;
}

/*La matriz electrode[i][j] toma el valor true si el punto i,j corresponde a un electrodo
(y, por tanto, su potencial no debe ser modificado al calcular, puesto que se trata de una
condicion de contorno). V[i][j] se debe haber inicializado previamente, dando el valor del
potencial del electrodo a aquellos i,j que correspondan a un electrodo, siendo los otros
valores inicializados como corresponda (cero, media, aleatorio, ...).*/
void Solve (bool electrode[][MAXZ], double V[][MAXZ], double a2[], double dr, double dz, double rmax, double zmax, int itmax, double delta)
{
	int iter=0, Nr=(int)(rmax/dr-0.5), Nz=(int)(zmax/dz-0.5);
	double ar1=(dz*dz)/(2*(dr*dr+dz*dz)), az=ar1*(dr*dr)/(dz*dz), aux, max, xi;
	ofstream filetest;
	
	/*Se calculan los ar2[i] para cada i (radio) aparte, para no hacerlo para cada punto en cada iteracion.*/
	for(int i=0;i<=Nr;i++)
		ar2[i]=ar1/(2.*i+1.);
	
	/*Este es el bucle de calculo propiamente dicho*/
	do{ 
		iter++; 
		if(iter%10000==0) cout << iter << " iteraciones calculadas." << endl;
		max=0.;
		
		/*Los calculos se hacen en el orden natural de la matriz (para cada i, j crecientes).*/
		if(!electrode[0][0])
		{
			aux=V[0][0];
			V[0][0]=((ar1+ar2[0])*V[1][0]+az*V[0][1])/(1-ar1+ar2[0]-az);
			if(V[0][0]!=0) //Si no se hace la comprobacion se puede dividir por cero (error)
			{
				xi=abs((V[0][0]-aux)/V[0][0]);
				max=Maximo(max,xi);
			}
		}
		
		for(int j=1;j<Nz;j++)
			if(!electrode[0][j])
			{
				aux=V[0][j];
				V[0][j]=((ar1+ar2[0])*V[1][j]+az*(V[0][j+1]+V[0][j-1]))/(1-ar1+ar2[0]);
				if(V[0][j]!=0)
				{
					xi=abs((V[0][j]-aux)/V[0][j]);
					max=Maximo(max,xi);
				}
			}
		
		if(!electrode[0][Nz])
		{
			aux=V[0][Nz];
			V[0][Nz]=((ar1+ar2[0])*V[1][Nz]+az*V[0][Nz-1])/(1-ar1+ar2[0]);
			if(V[0][Nz]!=0)
			{
				xi=abs((V[0][Nz]-aux)/V[0][Nz]);
				max=Maximo(max,xi);
			}
		}
		
		for(int i=1;i<Nr;i++)
		{
			if(!electrode[i][0])
			{
				aux=V[i][0];
				V[i][0]=(ar1*(V[i+1][0]+V[i-1][0])+ar2[i]*(V[i+1][0]-V[i-1][0])+az*V[i][1])/(1-az);
				if(V[i][0]!=0)
				{
					xi=abs((V[i][0]-aux)/V[i][0]);
					max=Maximo(max,xi);
				}
			}
			
			for(int j=1;j<Nz;j++)
				if(!electrode[i][j])
				{                
					aux=V[i][j];
					V[i][j]=ar1*(V[i+1][j]+V[i-1][j])+ar2[i]*(V[i+1][j]-V[i-1][j])+az*(V[i][j+1]+V[i][j-1]);
					if(V[i][j]!=0)
					{
						xi=abs((V[i][j]-aux)/V[i][j]);
						max=Maximo(max,xi);
					}
				}
				
			if(!electrode[i][Nz])
			{
				aux=V[i][Nz];
				V[i][Nz]=ar1*(V[i+1][Nz]+V[i-1][Nz])+ar2[i]*(V[i+1][Nz]-V[i-1][Nz])+az*V[i][Nz-1];
				if(V[i][Nz]!=0)
				{
					xi=abs((V[i][Nz]-aux)/V[i][Nz]);
					max=Maximo(max,xi);
				}
			}
		}
		
		if(!electrode[Nr][0])
		{
			aux=V[Nr][0];
			V[Nr][0]=((ar1-ar2[Nr])*V[Nr-1][0]+az*V[Nr][1])/(1-az);
			if(V[Nr][0]!=0)
			{
				xi=abs((V[Nr][0]-aux)/V[Nr][0]);
				max=Maximo(max,xi);
			}
		}
		
		for(int j=1;j<Nz;j++)
			if(!electrode[Nr][j])
			{
				aux=V[Nr][j];
				V[Nr][j]=(ar1-ar2[Nr])*V[Nr-1][j]+az*(V[Nr][j+1]+V[Nr][j-1]);
				if(V[Nr][j]!=0)
				{
					xi=abs((V[Nr][j]-aux)/V[Nr][j]);
					max=Maximo(max,xi);
				}
			}
			
		if(!electrode[Nr][Nz])
		{
			aux=V[Nr][Nz];
			V[Nr][Nz]=(ar1-ar2[Nr])*V[Nr-1][Nz]+az*V[Nr][Nz];
			if(V[Nr][Nz]!=0)
			{
				xi=abs((V[Nr][Nz]-aux)/V[Nr][Nz]);
				max=Maximo(max,xi);
			}
		}
		
	}while((max>delta)&&(iter<itmax));
	
	if(max<=delta)
		cout << "Convergencia alcanzada a las " << iter << " iteraciones." << endl << endl;
	else
		cout << "Convergencia no alcanzada en " << iter << " iteraciones." << endl << endl;
    		
	return;
}

/*Esta funcion escribe en un fichero las condiciones iniciales del potencial en
el formato requerido por gnuplot para hacer un contour plot.*/
void WriteInit (double V[][MAXZ], double rmax, double zmax, double dr, double dz, string name)
{
	int Nr=(int)(rmax/dr-0.5), Nz=(int)(zmax/dz-0.5);
	ofstream fileinit;
	
	fileinit.open((name+"-init.dat").c_str());
    for(int i=0;i<=Nr;i++)
    {
    	for(int j=0;j<=Nz;j++)
    		fileinit << (i+0.5)*dr << "\t" << (j+0.5)*dz << "\t" << V[i][j] << endl;
   		fileinit << endl;
    }
	fileinit.close();
	
	return;
}

/*Esta funcion escribe en un fichero los resultados en el formato requerido por
gnuplot para hacer un contour plot, asi como el potencial a lo largo de los ejes.*/
void WriteRes (double V[][MAXZ], double rmax, double zmax, double dr, double dz, string name)
{
	int Nr=(int)(rmax/dr-0.5), Nz=(int)(zmax/dz-0.5);
	ofstream fileres, filez, filer;
	
	fileres.open((name+"-res.dat").c_str());
    for(int i=0;i<=Nr;i++)
    {
    	for(int j=0;j<=Nz;j++)
    		fileres << (i+0.5)*dr << "\t" << (j+0.5)*dz << "\t" << V[i][j] << endl;
   		fileres << endl;
    }
	fileres.close();
	
	filez.open((name+"-axisz.dat").c_str());
	for(int j=0;j<=Nr;j++)
		filez << (j+0.5)*dz << "\t" << V[0][j] << endl;
	filez.close();
	
	filer.open((name+"-axisr.dat").c_str());
	for(int i=0;i<=Nr;i++)
		filer << (i+0.5)*dr << "\t" << V[i][0] << endl;
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

double Maximo (double d1, double d2)
{
	if(d1>d2) d2=d1;
	return d2;
}

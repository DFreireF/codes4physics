# include <iostream>
# include <fstream>
# include <cstdlib>
# include <string>
# include <cmath>
using namespace std ;

const double PI=2.0*atan2(1.0,0.0);

int main () {
	double t,x,dx,dt,tf,courant;
	int Nx,Nt,i,j;
	string numero;
	// Input:
	cout<<" Introduce los parámetros: Nx, Nt, tf:"<< endl;
	cin>>Nx>>Nt>>tf; getline(cin,numero);
	double u[Nx],d2udx2[Nx];
	// Inicialización :
	dx=1.0/(Nx-1);
	dt=tf/(Nt-1);
	courant=dt/(dx*dx);
	cout<<"  dx= "<<dx<<" dt= "<<dt<<"tf= "<<tf<<endl;
	cout<<"  Nx= "<<Nx<<" Nt= "<<Nt<<endl<<" Nº Courant= "<<courant<<endl;
	if (courant>0.5) cout<<" Inestabilidad temporal, disminuir courant\n";
	ofstream fichero ("d.dat");
	fichero.precision(17);
	//Inicialización a t=0
	//u(x,0) = sin(pi*x)
	for (i=0;i<Nx;i++){
		x= i*dx;
		u[i]=sin(PI*x);}
	u[0]=0.0;
	[Nx-1]=0.0;
	for (i=0;i<Nx;i++){
		x=i*dx;
		fichero<<0.0<<" "<<x<<" "<<u[i]<<"\n";}
	fichero<<"\n";
	// Evolución temporal :
	for (j=1;j<Nt;j++){ 
		t=j*dt;
	// Calculammos segunda derivada:
		for (i=1;i<Nx-1;i++) d2udx2[i]=courant*(u[i+1]-2.0*u[i]+u[i-1]);
	// Actualizamos valores :
		for (i=1;i<Nx-1;i++) u[i]+=d2udx2[i];
		for (i=0;i<Nx;i++){
			x=i*dx;
			fichero<<t<<" "<<x<<" "<<u[i]<<"\n";}
		fichero<<"\n";
	}//termina evolución temporal
	fichero.close();
	return 0;
}//termina el programa


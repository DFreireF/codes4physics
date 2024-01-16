# include <iostream>
# include <fstream>
# include <cstdlib>
# include <string>
# include <cmath>
using namespace std ;
const double PI=2.0*atan2(1.0,0.0);

int main () {
	double t,x,dx,dt,tf,courant,prob,r2,x0;
	int Nx,Nt,i,j,nizq,nder;
	string numero;
	// Input:
	cout<<" # Introduce los parámetros: Nx , Nt , tf :"<< endl;
	cin>>Nx>>Nt>>tf; getline(cin,numero);
	double u[Nx],d2udx2[Nx];
	// Inicialización :
	dx=1.0/(Nx-1);
	dt=tf/(Nt-1);
	courant=dt/(dx*dx);
	cout<<" dx= "<<dx<<" dt= "<<dt<<"tf= "<<tf<<endl;
	cout<<" Nx= "<<Nx<<" Nt= "<<Nt<<endl;
	cout<<" Nº Courant= "<<courant<<endl;
	if (courant>0.5) cout<<" Inestabilidad temporal\n";
	ofstream fichero1 ("d2.dat");
	fichero1.precision(17);
	ofstream fichero2 ("d3.dat");
	fichero2.precision(17);
	// Inicialización a t=0
	// u(x,0)=delta(x-x0)=delta(x-(Nx/2-1)*dx)
	for (i=0;i<Nx;i++){
		x=i*dx;
		u[i]=0.0;}
	u[0]= 0.0;
	u[Nx/2-1]=1.0; //aquí colocamos 1 partícula (la única)
	for (i=0;i<Nx;i++){
		x=i*dx;
		fichero1<<0.0<<" "<<x<<" "<<u[i]<<"\n";
		}
	fichero1<< "\n";
	// Evolución temporal:
	for (j=1;j<Nt;j++){ 
		t=j*dt;
		for (i=0;i<Nx;i++){ 
			nder=i+1;
			if (nder>Nx-1) nder=0;
			nizq=i-1;
			if (nizq<0) nizq=Nx-1;
			d2udx2[i]=courant*(u[nizq]-2.0*u[i]+u[nder]);
		}
		prob=0.0; r2=0.0; x0=((Nx/2)-1)*dx;
		for (i=0;i<Nx;i++){
			x=i*dx;
			u[i]+=d2udx2[i];
			prob+=u[i];
			r2+=u[i]*(x-x0)*(x-x0);
			fichero1<<t<<" "<<x<<" "<<u[i]<<"\n";
			}
		fichero1<<"\n";
		fichero2<<"pu"<<t<<" "<<prob<<" "<<r2<<" "<<u[Nx/2-1]<<" "<<u[Nx/4-1]<<" "<<u[0]<<"\n" ;
		}//evolución temporal j
	fichero1.close();
	fichero2.close();
	return 0;
}//programa termina


#include <iostream>
#include <cstdlib>
using namespace std;

double CalcularPi(int n);

int main (void)
{
    int n;
    double PI;
    cout<<"Introduce un valor para calcular el numero PI, a mayor valor, mayor precision: "<<endl<<"------>  ";
    cin>>n;
    PI= CalcularPi(n);
    cout<<endl<<endl<<" PI="<<PI<<endl;
    system("PAUSE");
    return 0;

}

//Función CalcularPi: calcula el número PI a partir de una sucesión de términos
//Parámetros: n
//Devuelve el número PI

double CalcularPi (int n)
{
    int i;
    double suma;
    suma=0.0;
    for (i=0;i<=n;i++)
        {
            if ((i==2) || (i%2==0))
            {
                suma+=(1.0/(2*i+1));
            }
            else   suma-=(1.0/(2*i+1));

        }
    return 4*suma;

}

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

//Funci�n CalcularPi: calcula el n�mero PI a partir de una sucesi�n de t�rminos
//Par�metros: n
//Devuelve el n�mero PI

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

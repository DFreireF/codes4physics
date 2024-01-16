/* Realice un programa que genere 3 numeros enteros aleatorios
y los muestre en pantalla, ordenados de mayor a menor. Use anidamiento
de condicionales, no hay que usar bucles y eso*/

#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main (void)
{
    int n1, n2, n3, mayor, menor, inter;

    srand (time(NULL));

    n1=rand ();
    n2=rand ();
    n3=rand ();
    mayor=0;
    menor=0;
    inter=0;


    if (n1>n2 && n1>n3 && n2>n3)
        {
            mayor=n1;
            inter=n2;
            menor=n3;
        }
    if (n2>n3 && n2>n1 && n1>n3)
        {
            mayor=n2;
            inter=n1;
            menor=n3;
        }
    if (n3>n1 && n3>n2 && n2>n1)
        {
            inter=n2;
            mayor=n3;
            menor=n1;
        }
    if (n1>n2 && n1>n3 && n2<n3)
        {
            mayor=n1;
            inter=n3;
            menor=n2;
        }
    if (n2>n3 && n2>n1 && n1<n3)
        {
            mayor=n2;
            inter=n3;
            menor=n1;
        }
    if (n3>n1 && n3>n2 && n1>n2)
    {
        mayor=n3;
        inter=n1;
        menor=n2;
    }

    cout<<mayor<<">"<<inter<<">"<<menor;

    return 0;

}

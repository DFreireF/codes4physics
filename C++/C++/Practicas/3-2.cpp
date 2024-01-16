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
    mayor=n1;
    menor=n3;


    if (n2>n1 && n2>n3)
    {
        n2=mayor;
    }
    else n3=mayor;

//--------------------------------------------

    if (n1<n2 && n3>n1)
    {
        n1=menor;
    }
    else menor=n2;

//---------------------------------------------

    if (n1>n2 && n1>n3 && n2>n3)
        {
            inter=n2;
    }
    else if (n2>n3 && n2>n1 && n1>n3)
        {
            inter=n1;
        }
    else inter=n3;

    cout<<mayor<<">"<<inter<<">"<<menor;

    return 0;

}

/*Realice un programa que genere 20 números enteros aleatorios comprendidos entre -1000 y 1000 y los muestre en pantalla.
 Además, al final deberá mostrar el valor máximo, el mínimo y el valor medio de los números generados. */

 #include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main (void)
{
    int i,n, mayor, menor, suma, media;
    mayor=-99999;
    menor=999999;
    suma=0;

    srand (time(NULL));
    for (i=0 ; i<=20 ; i++)
    {
       n= rand () % 2001 + -1000;
        if (n>mayor) mayor=n;
        if (n<menor) menor=n;
        suma=suma + n;

    }

    media= suma/20;

    cout<<"Este programa generara 20 numeros aleatorios y mostrara el mayor y menor de \tellos y la media de todos los numeros:\n";
    cout <<endl<< "El numero mayor ingresado es: " <<mayor<<endl;
    cout << "El numero menor ingresado es: " <<menor<<endl;
    cout<<"La media es: "<<media<<endl<<endl;

    system("pause");

    return 0;

}

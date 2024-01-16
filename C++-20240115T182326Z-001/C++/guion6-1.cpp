
#include <iostream>
#include <ctime>

bool es_primo(unsigned int n);
using namespace std;

int main(void)
{
unsigned int n,i;
unsigned long int suma;
clock_t t1, t2;
double numsecs;


//Pedir el número
cout << "Introduzca un numero: ";
cin >> n;

//Realizar la suma acumulada
t1=clock();
suma=0;
for(i=2; i<=n; i++){
if(es_primo(i)) suma=suma+i;}
t2=clock();

//Mostrar el resultado
cout << "Suma de los primos menores o iguales que "; cout << n << ": " << suma << endl;

//Muestro el tiempo invertido en el cálculo
numsecs=(double)(t2-t1)/CLOCKS_PER_SEC;
cout << "Tiempo: " << numsecs << " segundos" << endl;

return 0;
}

//Función que comprueba si un número entero es primo
//Devuelve un booleano: TRUE si es primo, FALSE si no lo es
//Divide el número por todos los enteros menores que él
//Si el resto de alguna de las divisiones es cero
//entonces no es primo

bool es_primo (unsigned int n)
{
unsigned int i;
bool EsPrimo;

EsPrimo=true; for(i=2; i<n; i++)
if(n%i==0) EsPrimo=false;

return EsPrimo;
}

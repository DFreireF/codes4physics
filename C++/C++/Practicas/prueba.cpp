#include <iostream>
#include <fstream>
#include <cmath>
#i nclude <string>
#include <iomanip>

using namespace std;

double CalcularD (int numv[], int nums[], int p)
{
    int i, totalVotos, totalEscanios;
    double suma;
    //Para calcular porcentajes primero calculo totales de votos y de escanios
    totalVotos = totalEscanios = 0;
    for(i=0; i <p; i++)
    {
    totalVotos += numv[i];
    totalEscanios += nums[i];
    }
    //Calcular distorsion, ojo con la división que es real
        suma=0.0;
    for(i=0; i<p; i++)
    suma += fabs(100.0*numv[i]/totalVotos-100.0*nums[i]/totalEscanios);
    return 0.5*suma;
}


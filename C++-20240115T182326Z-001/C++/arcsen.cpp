#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;

double arcsen (float x, int p)
{
    double suma=0.0;

    for(int n=1;n<=p;n++)
    {
        suma+=(((2*n-1)*pow(x, 2*n+1))/((2*n)*(2*n+1)));
    }
    return suma+x;
}
int main(void)
{
    const int pi= 3.141592;
    double x, arcsin =0.0;
    int nveces=0;
    cout<<"Introduzzca un numero para arcsen"<<endl;
    cin>>x;
    cin>>nveces;
    arcsin=arcsen(x,nveces);
    cout<<arcsin<<"           ";
    system("PAUSE");
    return 0;
}

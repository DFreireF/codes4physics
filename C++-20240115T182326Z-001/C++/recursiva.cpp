#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
using namespace std;

int recursiva (int n[])
{
    int l=n.size();
    for (int i=0;i<l;i++)
    {
        n[i]=n[l-i];
    }
    return n;
}
int main(void)
{
    string a, b;
    int l=b.size();
    int nveces=0;
    cout<<"Introduzzca un numero"<<endl;
    cin>>a;
    b= recursiva(a);
    cout<<l<<"  y then "<<b;
    system("PAUSE");
    return 0;
}

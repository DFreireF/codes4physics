#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

int main()
{
    int i,j,k;
    
    double dx=0.01, dy=0.01, dz=0.01;
    double ax=1/(2*(1+dx*dx/(dy*dy)+dx*dx/(dz*dz)));
    double **Grid=(double**)malloc(Nrho*sizeof(double*));
    
    for(i=0;i<Nrho;i++)
    {
        Grid[i]=(double*)malloc(Nz*sizeof(double));
        for(j=0;j<Nz;j++)
        Grid[i][j]=0;
    }
    
    /*
    Reverse ingeniering this thing
    Do[Do[ Loop
    res =  Residue 
    0.25*(MyGrid[[j + 1, i]] + MyGrid[[j - 1, i]] + 
    MyGrid[[j, i + 1]] + MyGrid[[j, i - 1]]) - MyGrid[[j, i]]; 
   MyGrid[[j, i]] = MyGrid[[j, i]] + AA*res; //Original + weight*residue
   error = error + Abs[res], {j, 2, \[Rho]e}], {i, ze + 2, Nz - ze}];*/
    
    double res;
    double Weight=1.25;
    
   /* double Ve=-0.0448963;
    double Vc=-0.0209951;
    double V0=0.0551037;
    double Ve2=-0.0448963;
    double Vc2=-0.0209951;
    double V02=0.0551037;*/
    
    double Ve =1e-22;
    double Vc =50;
    double V0 =1e-22;
    double Ve2=1e-22;
    double Vc2=1e-22;
    double V02=1e-22;
    
    
    for(i=0;i<521;i++)
    {
        j=2;
        if(i<=211)Grid[j][i+15]=Ve;
        if(i>211&&i<=251)Grid[j][i+15]=Vc;
        if(i>251&&i<=259)Grid[j][i+15]=V0;
        if(i==260       )Grid[j][i+15]=V0;
        if(i>260&&i<=268)Grid[j][i+15]=V02;
        if(i>268&&i<=308)Grid[j][i+15]=Vc2;
        if(i>308&&i<=521)Grid[j][i+15]=Ve2;
        j=3;
        if(i<=211)Grid[j][i+15]=Ve;
        if(i>211&&i<=251)Grid[j][i+15]=Vc;
        if(i>251&&i<=259)Grid[j][i+15]=V0;
        if(i==260       )Grid[j][i+15]=V0;
        if(i>260&&i<=268)Grid[j][i+15]=V02;
        if(i>268&&i<=308)Grid[j][i+15]=Vc2;
        if(i>308&&i<=521)Grid[j][i+15]=Ve2;
        j=4;
        if(i<=211)Grid[j][i+15]=Ve;
        if(i>211&&i<=251)Grid[j][i+15]=Vc;
        if(i>251&&i<=259)Grid[j][i+15]=V0;
        if(i==260       )Grid[j][i+15]=V0;
        if(i>260&&i<=268)Grid[j][i+15]=V02;
        if(i>268&&i<=308)Grid[j][i+15]=Vc2;
        if(i>308&&i<=521)Grid[j][i+15]=Ve2;
        
        j=97;
        if(i<=211)Grid[j][i+15]=Ve;
        if(i>211&&i<=251)Grid[j][i+15]=Vc;
        if(i>251&&i<=259)Grid[j][i+15]=V0;
        if(i==260       )Grid[j][i+15]=V0;
        if(i>260&&i<=268)Grid[j][i+15]=V02;
        if(i>268&&i<=308)Grid[j][i+15]=Vc2;
        if(i>308&&i<=521)Grid[j][i+15]=Ve2;
        j=96;
        if(i<=211)Grid[j][i+15]=Ve;
        if(i>211&&i<=251)Grid[j][i+15]=Vc;
        if(i>251&&i<=259)Grid[j][i+15]=V0;
        if(i==260       )Grid[j][i+15]=V0;
        if(i>260&&i<=268)Grid[j][i+15]=V02;
        if(i>268&&i<=308)Grid[j][i+15]=Vc2;
        if(i>308&&i<=521)Grid[j][i+15]=Ve2;
        j=98;
        if(i<=211)Grid[j][i+15]=Ve;
        if(i>211&&i<=251)Grid[j][i+15]=Vc;
        if(i>251&&i<=259)Grid[j][i+15]=V0;
        if(i==260       )Grid[j][i+15]=V0;
        if(i>260&&i<=268)Grid[j][i+15]=V02;
        if(i>268&&i<=308)Grid[j][i+15]=Vc2;
        if(i>308&&i<=521)Grid[j][i+15]=Ve2;
        
        
    }
    for(k=0;k<1e5;k++)
    {
        if(k%10000==0)cout<<k<<endl;
    for(i=0;i<Nrho;i++)
        //that's where the electrodes are
        if(i==2||i==3||i==4||i==98||i==96||i==97)
        {
            for(j=0;j<15;j++)
            {
                if(j==0)res=0.25*(Grid[i+1][j]+Grid[i-1][j]+2.0*Grid[i][j+1])-Grid[i][j];
                else res=0.25*(Grid[i-1][j]+Grid[i+1][j]+Grid[i][j-1]+Grid[i][j+1])-Grid[i][j];
                Grid[i][j]+=res*Weight;
            }
            for(j=536;j<Nz;j++)
            {
                if(j==Nz-1)res=0.25*(Grid[i+1][j]+Grid[i-1][j]+2.0*Grid[i][j-1])-Grid[i][j];
                else res=0.25*(Grid[i-1][j]+Grid[i+1][j]+Grid[i][j-1]+Grid[i][j+1])-Grid[i][j];
                Grid[i][j]+=res*Weight;
            }
        }
        else
        for(j=0;j<Nz;j++)
        {
            //Corners
            if(i==0&&j==0)res=0.25*(2.0*Grid[i+1][j]+2.0*Grid[i][j+1])-Grid[i][j];
            if(i==0&&j==Nz-1)res=0.25*(2.0*Grid[i+1][j]+2.0*Grid[i][j-1])-Grid[i][j];
            if(i==Nrho-1&&j==0)res=0.25*(2.0*Grid[i-1][j]+2.0*Grid[i][j+1])-Grid[i][j];
            if(i==Nrho-1&&j==Nz-1)res=0.25*(2.0*Grid[i-1][j]+2.0*Grid[i][j-1])-Grid[i][j];
            //Sides
            if(i==0&&j!=0&&j!=Nz-1)res=0.25*(2.0*Grid[i+1][j]+Grid[i][j+1]+Grid[i][j-1])-Grid[i][j];
            if(i==Nrho-1&&j!=0&&j!=Nz-1)res=0.25*(2.0*Grid[i-1][j]+Grid[i][j+1]+Grid[i][j-1])-Grid[i][j];
            if(j==0&&i!=0&&i!=Nrho-1)res=0.25*(Grid[i+1][j]+Grid[i-1][j]+Grid[i][j+1])-Grid[i][j];
            if(j==Nz-1&&i!=0&&i!=Nrho-1)res=0.25*(Grid[i+1][j]+Grid[i-1][j]+Grid[i][j-1])-Grid[i][j];
            //Interior
            if(i!=0&&j!=0&&i!=Nrho-1&&j!=Nz-1)res=0.25*(Grid[i-1][j]+Grid[i+1][j]+Grid[i][j-1]+Grid[i][j+1])-Grid[i][j];
            Grid[i][j]+=res*Weight;}
    }
    ofstream out;out.open("Test.dat");
    
    for(i=0;i<Nrho;i++)
        for(j=0;j<Nz;j++)
            out<<j<<"\t"<<i<<"\t"<<Grid[i][j]<<"\n";
    
     
    //Derivative at 0,0 [55,275]
    
    double E=(-Grid[55][277]+8.0*Grid[55][276]-8*Grid[55][274]+Grid[55][273])/(12.0*dz);
    E*=1000;
    
    double du;
    double At = 187;//40;
    double Z = 75;//20;
    double mass = abs((Z - 1))*9.10938188*1e-31 + Z*1.6726231*1e-27 + (At - Z)*1.67492728*1e-27;
    double T=300;
    double fz=100*1e3;
    double A = sqrt(2*1.381*1e-23*T/(mass*pow(fz*2*acos(-1.0),2.0)));
    
    double C = 10*1e-12;
    double q = 1.6*1e-19;
    du=A*E*q/(1.0*C*Vc);
    At = 40;
    Z = 20;
    mass = abs((Z - 1))*9.10938188*1e-31 + Z*1.6726231*1e-27 + (At - Z)*1.67492728*1e-27;
    
    double b=(q*du*E)/(2.0*mass*2*fz*acos(-1.0)*Vc);
    
    cout<<"Du: \t"<<du<<endl;
    cout<<"b: \t"<<b<<endl;
    
    system("pause");
    //system("pause");
    return 0;
}

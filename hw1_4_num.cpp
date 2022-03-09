/* 
    第四题Numerov代码(C++)
*/
#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
using namespace std;

void output(double **array,int m,int n,string filepath)
{   
    ofstream outFile;
    outFile.open(filepath);
    for (int i=0;i<m;i++)
    {
        for (int j=0;j<n;j++)
        {
            outFile << ((double *)array+n*i)[j] << ", ";
        }
        outFile << "\n";
    }
}

double integral(int Nx, double *PHI, double dx)
{
     double s=0;
    for (int i=1;i<=Nx-1;i++)
    {
        s=s+PHI[i]*PHI[i];
    }
    return s*dx;
}

int numerov(int Nx, double *X, double *PHI, double *V, double E, double dx, double eps)
{   
     double f[Nx+1];
    for (int i=0;i<=Nx;i++)
    {
        f[i]=1+(2*(E - V[i]))/12*dx*dx;
    }
    PHI[0]=0;
    PHI[1]=0.01;
    for (int i=1;i<Nx;i++)
    {
        PHI[i+1]= ((12 - 10*f[i])*PHI[i] - f[i-1]*PHI[i-1])/f[i+1];
    }
    double c=1/pow(integral(Nx,PHI,dx),0.5);
    for (int i=0;i<=Nx;i++)
    {
        PHI[i]=PHI[i]*c;
    }

    if (PHI[Nx]>eps)
    {
        return 1;
    }
    else if (PHI[Nx]<-eps)
    {
        return -1;
    }
    else
    {
        return 0;
    }
}

double bi_solve(int Nx, double *X, double *PHI, double *V, double E_low, double E_high, double dx, double eps,int low,int high)
{   
    double E_mid=(E_low+E_high)/2;
    int mid=numerov(Nx,X,PHI,V,E_mid,dx,eps);
    if (mid==0)
    {
        return E_mid;
    }
    else if (mid==high)
    {
        return bi_solve(Nx,X,PHI,V,E_low,E_mid,dx,eps,low,mid);
    }
    else
    {
        return bi_solve(Nx,X,PHI,V,E_mid,E_high,dx,eps,mid,high);
    }
}

int main()
{   
    //parameter
    int Nx=2e2;
     double E;
     double max_step=0.1;
     double dx=10.0/Nx;
     double eps=1e-6;
    //init array
     double X[Nx+1];
     double PHI[Nx+1];
     double V[Nx+1];
    for (int i=0;i<=Nx;i++)
    {
        X[i]=i*dx-5;
        V[i]=0.5*(((X[i]*X[i]-1)*(X[i]*X[i]-1)/4)-X[i]*X[i]);//调整V0
    }

    
    //loop for E
    int task=1;
    double data_E[task];
    double data_PHI[task][Nx+1];
    E=-5;
    int count=0;
    int last;
    int now;
    while (last=numerov(Nx,X,PHI,V,E,dx,eps)==0)//make sure the first E is not a solution
    {
        E=E-max_step;
    }
    while (count <task)
    {   
        E=E+max_step;
        int result=numerov(Nx,X,PHI,V,E,dx,eps);
        if (last*result<0)// come to bisection method
        {
            data_E[count]=bi_solve(Nx,X,PHI,V,E-max_step,E,dx,eps,last,-last);
            for(int i=0;i<=Nx;i++)
            {
                data_PHI[count][i]=PHI[i];
            }
            E=data_E[count];
            count++;
        }
        else if (result==0)
        {
            data_E[count]=E;
            for(int i=0;i<=Nx;i++)
            {
                data_PHI[count][i]=PHI[i];
            }
            count++;
        }
        last=result;
        
    }
    for (int i=0;i<task;i++)
    {
        cout << setprecision(9) <<data_E[i] << endl;
    }
    output((double **)data_PHI,task,Nx+1,"./hw1_4_phi.csv");
    system("pause");
}


#include <iostream>
#include <math.h>
#include <windows.h>

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
    //std::cout << c <<std::endl;
    for (int i=0;i<=Nx;i++)
    {
        PHI[i]=PHI[i]*c;
    }
    std::cout << PHI[Nx] <<std::endl;
    //Sleep(100);

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
    std::cout << "b" <<std::endl;
    std::cout <<"E"<<E_low<<" "<<E_high<<std::endl;
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
    int Nx=1e3;
     double E;
     double max_step=0.1;
     double dx=1.0/Nx;
     double eps=1e-3;
    //init array
     double X[Nx+1];
     double PHI[Nx+1];
     double V[Nx+1];
    for (int i=0;i++;i<=Nx)
    {
        X[i]=i*dx;
        V[i]=0;
    }

    
    //loop for E
    int task=10;
    double data_E[task];
    double data_PHI[task][Nx+1];
    E=-1;
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
    std::cout << data_E[i] << std::endl;
    }
system("pause");
}

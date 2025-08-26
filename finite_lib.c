#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "finite_lib.h"

double Sat_Quan(double x, double Max_V, double Min_V, int m, int k)
{
    double val;
    int p,q;
    
    val = x;
    val = val > Max_V ? Max_V : val;
    val = val < Min_V ? Min_V : val;

    q = Quantize(val, m, k);
    Tru2intS(q,m,&p);
    return (double) p / pow(2, k);
}

int Quantize(double x, int m, int k)
{
    int t, t1;
    t=1<<k;
    t=(int)(x*t+(x>0?0.5 : -0.5));
    t1=(1<<m)-1;
    t1=t1&t;
    return t1;
}

void Tru2intS(int x, int m, int *dx)
{
    int t1=x;
    int t2=-1;
    t1=(x>>(m-1));
    if (t1==0)(*dx)=x;
    else{
        t1=(1<<m)-1;
        t2=t2^t1;
        *dx=t2|x;
    }
}
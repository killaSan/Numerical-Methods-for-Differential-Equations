#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
#define PI 3.14159265

double* progonka(double* a, double* b, double* c, double* f, int n)
{
    double *alpha = new double[n];
    double *beta = new double[n];
    double *sol = new double[n];
    n--; // since we start from x0 (not x1)
    alpha[0] =  b[0]/ c[0];
    beta[0]  =  f[0] / c[0];
    for (int i = 1; i < n; i++)
    {
        alpha[i] = b[i]/(c[i] - a[i]*alpha[i-1]);
        beta[i] = (-a[i]*beta[i-1] + f[i])/(c[i] - a[i]*alpha[i-1]);
    }
    beta[n] = (-a[n]*beta[n-1] + f[n])/(c[n] - a[n]*alpha[n-1]);
    sol[n] = beta[n];

    for (int i = n - 1 ; i >= 0; i--)
    {
        sol[i] = beta[i] - alpha[i]* sol[i+1];
    }

    delete []alpha;
    delete []beta;
    return sol;
}

/*
BegI and EndI are the  beginning and ending of the interval.
interal[N] is the division of the interval [1,-1] on N parts with distance h.
a[N], b[N], c[N] and f[N] are the coeficients with which the the Tridiagonal matrix algorithm is going to
be called.
*/

double* yh(int begI, int endI, int N, double h)
{
    double *interval = new double[N];
    double *sol = new double[N];
    interval[0] = begI;
    interval[N-1] = endI;

    for (int i = 1 ;i<N-1; i++)
    {
        interval[i] = interval[i-1]+h;
    }
    double *f = new double[N];
    for (int i = 0;i<N-1;i++)
    {
       f[i] = (1+interval[i])*sin(interval[i]) + (1-interval[i])*cos(interval[i]);
    }
    f[0] += 2/h; // because of the approximation
    f[N-1] = 1; // because of the approximation

    double *a = new double[N];
    a[0] = a[N-1] = 0;
    for (int i = 1;i < N - 1;i++)
    {
        a[i] = 1/(h*h);
    }

    double *b = new double[N];
    b[N-1] = 0;
    for(int i = 0;i < N-1;i++)
    {
        b[i] = 1/(h*h);
    }

    double *c = new double[N];
    for (int i =0;i < N - 1;i++)
    {
        c[i] = - (2/(h*h)) - interval[i];
    }
    c[N-1] = -1;
    sol = progonka(a,b,c,f,N);
    delete []f;
    delete []interval;
    delete []a;
    delete []b;
    delete []c;
    return sol;
}

int main()
{
    ofstream file;
    file.open("solution.txt");
    double begI = 0;
    double endI = 2*PI;
    int N;
    cout << "N= ";
    cin >> N;
    double h1 = (endI-begI)/N;
    double h2 = (endI-begI)/(2*N);
    double h4 = (endI-begI)/(4*N);
    double *f = new double[N];
    double *s = new double[2*N];
    double *t = new double[4*N];
    f = yh(begI,endI, N, h1);
    s = yh(begI,endI, 2*N, h2);
    t = yh(begI,endI, 4*N, h4);
    double *news = new double[N];
    double *newt = new double[N];
    int j = 0;
    for (int i = 0;i <N;i ++)
    {
        news[i] = s[j];
        j+=2;
    }
    news[N-1] = s[2*N-1];
    j = 0;

    for (int i = 0; i <N;i ++)
    {
        newt[i] = t[j];
        j+=4;
    }
    newt[N-1] = t[4*N-1];
        for (int i =0;i<N;i++)
    {
        file << f[i] << '\t';
    }
    file << endl;
        for (int i =0;i<N;i++)
    {
        file << news[i] << '\t';
    }
    file << endl;
        for (int i =0;i<N;i++)
    {
        file << newt[i] << '\t';
    }
    double *alfa = new double[N];
    file << endl;
    for (int i =0;i< N;i++)
    {
        alfa[i] =  log(abs((f[i] - news[i])/(news[i]-newt[i])))/log(2);
    }
    for (int i = 0;i<N;i++)
    {
        cout << alfa[i] << '\t';
        file << alfa[i] << endl;
    }
    file.close();
    delete []f;
    delete []s;
    delete []t;
    delete []news;
    delete []newt;
    delete []alfa;
    return 0;
}

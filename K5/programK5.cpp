#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;
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

void printMatrix(double **m , int rows, int cols)
{
    for (int i = 0;i < rows; i ++)
    {
        for (int j = 0;j < cols; j++)
        {
            cout << m[i][j] << '\t';
        }
        cout << endl;
    }
}

void deleteMatrix(double **m, int rows)
{
    for (int i = 0;i < rows; i++)
    {
        delete []m[i];
    }
    delete []m;
}

double *getRow(double **m, int rows, int cols, int n)
{
    double *res = new double[cols];
    for (int i=0;i<cols; i++)
    {
        res[i] = m[n][i];
    }
    return res;
}
double **explicitS(double **m, int rows,int cols, double tao, double h,double *x, double *t)// явна схема
{
    for (int i = 0 ; i < rows ; i ++) // u(x,0) = 0;
    {
        m[i][0] = 0;
    }
    for(int i = 0; i < cols; i++) // u(0,t) = 0;
    {
        m[0][i] = 0;
    }
    double *f = new double[cols];
    for (int i = 0;i < cols;i++)
    {
        f[i] = x[i];
    }

    for (int j = 1; j < rows; j++)
    {
        for (int i = 1; i < cols - 1; i++)
        {
            double a = tao/(h*h);
            double b = (h*h - (2*tao))/(h*h);
            m[j][i] = a*(m[j-1][i+1] + m[j-1][i-1]) + b*m[j-1][i] + tao* f[i];
        }
        double a = (4*h*tao)/(2*h*tao + 2*tao - h*h);
        double b = (h*h*tao)/(2*h*tao + 2*tao - h*h);
        double c = (h*h)/(2*h*tao + 2*tao - h*h);
        double d = (2*tao)/(2*h*tao + 2*tao - h*h);
        m[j][cols - 1] = t[j]*a - f[cols-1]*b - m[j-1][cols - 1]*c + d*m[j][cols - 2];
    }
    printMatrix(m,rows,cols);
    return m;
}

double **implicit(double **m, int rows,int cols, double tao, double h,double *x, double *t) // неявна схема
{
    double *f = new double[cols];
    for (int i = 0;i < cols;i++)
    {
        f[i] = x[i];
    }
    for (int i = 1; i < rows ; i++)
    {
        double *prevLevel = new double[cols];
        prevLevel = getRow(m,rows,cols,i-1); // getting the y-s from j-1-th level
        double *A = new double[cols];
        A[0] = 0;
        A[cols - 1] = 1/h;
        for (int k = 1; k<cols - 1; k++)
        {
            A[k] = 1/(h*h);
        }
        double *B = new double[cols];
        B[0] = 0;
        B[cols-1] = 0;
        for (int k = 1;k<cols-1; k++)
        {
            B[k] = 1/(h*h);
        }
        double *C = new double[cols];
        C[0]= -1;
        C[cols-1] = -(1 + (1/h) - (h/(2*tao)));
        for (int k = 1; k < cols-1;k++)
        {
            C[k] = - ((1/tao) + (2/(h*h)));
        }
        double *F = new double[cols];
        F[0] = 0;
        F[cols - 1] = - 2*t[i] + ((h*f[cols - 1])/2) + ((h*prevLevel[cols-1])/(2*tao)) ;
        for (int k = 1; k< cols-1;k++)
        {
            F[k] = - f[k] - (prevLevel[k]/tao);
        }
        double *res = new double[cols];
        res = progonka(A,B,C,F,cols);
        for (int w = 0; w < cols;w++)
        {
            m[i][w] = res[w];
        }
        delete []prevLevel;
        delete []A;
        delete []B;
        delete []C;
        delete []F;
        delete []res;
    }
    return m;
}


int main()
{
int n, m;
    cout << "n: ";
    cin >> n; // cols
    cout << "m: ";
    cin >> m; // rows
    double beginX = 0;
    double endX = 1;
    double beginT = 0;
    double endT = 2;
    double h = (endX - beginX)/n;
    double tao = (endT - beginT)/m;
    cout << h << endl << tao << endl;
    double *x = new double[n+1];
    x[0] = beginX;
    for (int i = 1;i < n + 1;i++)
    {
        x[i] = h + x[i-1];
    }
    double * t = new double[m+1];
    t[0]= beginT;
    for (int i = 1; i < m+1;i++)
    {
        t[i] = tao + t[i-1];
    }
    int minXT = min(n+1,m+1);
    double **exact = new double* [m+1];
    for (int i = 0; i < m+1;i++)
    {
        exact[i] = new double[n + 1];
    }
    for (int i = 0 ;i < m+1;i++)
    {
        for (int j = 0;j < n+1 ;j++)
        {
            exact[i][j] = t[i]*x[j];
        }
    }
    cout << endl;
    double ** matrix = new double* [m+1];
    for (int i = 0;i < m+1; i++) // n+1
    {
        matrix[i] = new double[n+1]; // m+1
    }
    matrix = implicit(matrix, m+1,n+1,tao, h, x,t);
    cout << endl;
    for (int i = 0;i<m+1;i++)
    {
        for (int j = 0;j <n+1; j++)
        {
          //  cout << fixed;
            cout << abs(exact[i][j] - matrix[i][j]) << '\t';
        }
        cout << endl;
    }

    deleteMatrix(matrix,m+1);
    deleteMatrix(exact, m+1);
    delete []x;
    delete []t;
}

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;

const double eps = pow(10, -8);
double **initialMat(double **m, int n)
{
    double x1 = 0;
    double x2 = 0;
    for (int i = 0 ;i < n;i++)
    {
        for (int j = 0; j < n;j++)
        {
            m[i][j] = 0;
        }
    }
    return m;
}

double **createMatrix(int n)
{
    double **res = new double*[n];
    for (int i = 0; i< n;i++)
    {
        res[i] = new double[n];
    }
    return res;
}
void printMatrix(double **m, int rows, int cols)
{
    for (int i = 0; i < rows; i ++)
    {
        for (int j = 0; j < cols; j++)
        {
            cout << m[i][j] << '\t';
        }
        cout << endl;
    }
}

void deleteMatrix(double **m, int rows)
{
    for (int i = 0; i < rows; i++)
    {
        delete []m[i];
    }
    delete []m;
}

double **yakobi(double **m, int n, double h)
{
    double x1 = 0;
    double x2 = 0;
    double ** f = createMatrix(n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            f[i][j] = (1 - x1)*(1 - x1) * (1 - x2)*(1 - x2) + x1 * x2;
          // f[i][j] = x1 + cos(3*x2);
        //  f[i][j] = x1*x1 + 2*x2;
            x1 +=h;
        }
        x2 +=h;
        x1 = 0;
    }
    for (int i = 1; i < n - 1; i++)
    {
        for (int j = 1; j < n-1; j++)
        {
            bool k = 1;
            while (k)
            {
                double old = m[i][j];
                m[i][j] = (m[i-1][j] + m[i+1][j] + m[i][j-1] + m[i][j+1] + h*h*f[i][j])/4;
                if (((abs (m[i][j] - old))/(abs(m[i][j])) - eps) < 0.00001)
                {
                    k = 0;
                }
            }
        }
    }
    deleteMatrix(f,n);
    return m;
}


int main()
{
    ofstream file ;
    file.open( "res.txt");
    if(!file)
    {
        return -1;
    }
    int n;
    cin >> n;
    double *part = new double[n+1];
    part[0] = 0;
    double h = 1/double(n);
    for (int i = 1;i < n+1;i++)
    {
        part[i] = part[i-1] + h;
    }
    double **matrix = createMatrix(n+1);
    matrix = initialMat(matrix,n+1);
    double **matrix2 = createMatrix(2*n + 1);
    matrix2 = initialMat(matrix2,2*n+1);
    double **matrix4 = createMatrix(4*n + 1);
    matrix4 = initialMat(matrix4,4*n+1);
    double **alpha = createMatrix(n+1);
    matrix = yakobi(matrix,n+1,h);
    matrix2 = yakobi(matrix2, 2*n+1, h/2);
    matrix4 = yakobi(matrix4, 4*n+1, h/4);
    for (int i = 0 ; i<n+1; i++)
    {
        file << "Level: " << i <<" = " << part[i] << endl;
        for (int j = 0; j<n+1; j++)
        {
            alpha[i][j] = log(abs((matrix[i][j] - matrix2[2*i][2*j])/(matrix2[2*i][2*j] - matrix4[4*i][4*j]) ))/ log(2);
            file << alpha[i][j] << '\t';
        }
        file << endl;
    }

    deleteMatrix(matrix,n+1);
    deleteMatrix(matrix2,2*n+1);
    deleteMatrix(matrix4,4*n+1);
    deleteMatrix(alpha,n+1);
    delete []part;
    file.close();
    return 0;
}

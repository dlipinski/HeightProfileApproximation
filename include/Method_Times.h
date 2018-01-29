#ifndef METHOD_TIMES_H
#include <iostream>
#include <ctime>
#include <stdio.h>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <time.h>
#include "MyMatrix.h"
using namespace std;
class Method_Times
{
    public:
        Method_Times();

        //-------------------------------------------------------------------partial Gauss
        double Gauss_partial_time() {
            int n = get_rows();

            clock_t begin = clock();
            for (int i = 0; i < n; i++) {
                // znajd wiersz z maksymalnym elementem
                double maxEl = abs(matrix[i][i]);
                int maxRow = i;
                for (int k = i+1; k < n; k++) {
                    if (abs(matrix[k][i]) > maxEl) {
                        maxEl = abs(matrix[k][i]);
                        maxRow = k;
                    }
                }
                // zamieñ maksymalny wiersz z obecnym
                for (int k = i; k < n+1; k++) {
                    double pom = matrix[maxRow][k];
                    matrix[maxRow][k] = matrix[i][k];
                    matrix[i][k] = pom;
                }
                // wyprowad zera przed obecnym wierszem
                for (int k = i+1; k < n; k++) {
                    double c = -matrix[k][i] / matrix[i][i];
                    for (int j = i; j < n+1; j++) {
                        if (i == j) {
                            if(matrix[k][j] > 0)
                                matrix[k][j] = 0;
                        } else {
                            matrix[k][j] += c * matrix[i][j];
                        }
                    }
                }
            }
            // rozwi¹¿ Ax = B za pomoc¹ powsta³ej macierzy trójk¹tnej
            vector<double> x(n);
            for (int i=n-1; i>=0; i--) {
                x[i] = matrix[i][n] / matrix[i][i];
                for (int k=i-1;k>=0; k--) {
                    matrix[k][n] -= matrix[k][i] * x[i];
                }
            }
            clock_t end = clock();
            return ());
        }
//-------------------------------------------------------------------partial Gauss optimized
        vector<double> Gauss_partial_optimized() {
            int n = get_rows();
            for (int i = 0; i < n; i++) {
                // znajd wiersz z maksymalnym elementem
                double maxEl = abs(matrix[i][i]);
                int maxRow = i;
                for (int k = i+1; k < n; k++) {
                    if (abs(matrix[k][i]) > maxEl) {
                        maxEl = abs(matrix[k][i]);
                        maxRow = k;
                    }
                }
                // zamieñ maksymalny wiersz z obecnym
                for (int k = i; k < n+1; k++) {
                    double pom = matrix[maxRow][k];
                    matrix[maxRow][k] = matrix[i][k];
                    matrix[i][k] = pom;
                }
                // wyprowad zera przed obecnym wierszem
                for (int k = i+1; k < n; k++) {
                    double c = -matrix[k][i] / matrix[i][i];
                    for (int j = i; j < n+1; j++) {
                        if (i == j) {
                            if(matrix[k][j] > 0  && matrix[k][j] !=0)
                                matrix[k][j] = 0;
                        } else {
                            matrix[k][j] += c * matrix[i][j];
                        }
                    }
                }
            }
            // rozwi¹¿ Ax = B za pomoc¹ powsta³ej macierzy trójk¹tnej
            vector<double> x(n);
            for (int i=n-1; i>=0; i--) {
                x[i] = matrix[i][n] / matrix[i][i];
                for (int k=i-1;k>=0; k--) {
                    matrix[k][n] -= matrix[k][i] * x[i];
                }
            }
            return x;
        }

//-----------------------------------------------------------Jacobi
        vector<double> Jacobi(int iter,double initial){
            //Create,resize and fill matrix A
            vector < vector < double > > A;
            A.resize(width-1);
            for(int i=0;i<width-1;i++){
                A[i].resize(width-1,0);
                for(int j=0;j<width-1;j++)
                    A[i][j]=matrix[i][j];
            }
            //Create,resize and fill vector B
            vector <double > B;
            B.resize(width-1,0);
            for(int i=0;i<width-1;i++)
                B[i]=matrix[i][width-1];
            //Create matrix X (answer)
            vector <double > X;       // <--- x(k+1)
            X.resize(width-1,initial);
            //Create matrix X (previous answer)
            vector <double > _X;      // <--- x(k)
            _X.resize(width-1,0);
            double sum=0;

            //Start counting
            for(int it=0;it<iter;it++){
                for(int i=0;i<width-1;i++){
                    sum=0;
                    for(int j=0;j<width-1;j++){
                        if(j!=i)
                            sum+=A[i][j] * _X[j]; //sum(AijXj)
                    }
                    sum=-sum;// -sum(AijXj)
                    sum+=B[i]; // -sum(AijXj)+Bi
                    sum/=A[i][i]; // ( -sum(AijXj)+Bi ) / Aii
                    X[i]=sum; // X =
                }
                _X=X;
            }
            return X;
        }
//-------------------------------------------------------------------Gauss Seidel
        vector<double> Siedel(int iter,double initial){
            //Create,resize and fill matrix A
            vector < vector < double > > A;
            A.resize(width-1);
            for(int i=0;i<width-1;i++){
                A[i].resize(width-1,0);
                for(int j=0;j<width-1;j++)
                    A[i][j]=matrix[i][j];
            }
            //Create,resize and fill vector B
            vector <double > B;
            B.resize(width-1,0);
            for(int i=0;i<width-1;i++)
                B[i]=matrix[i][width-1];
            //Create matrix X (answer)
            vector <double > X;            // <--- x(k+1)
            X.resize(width-1,initial);
            //Create matrix _X (previous answer)
            vector <double > _X;           // <--- x(k)
            _X.resize(width-1,initial);
            double sum1=0,sum2=0,sum=0;

            //Start counting
            for(int it=0;it<iter;it++){
                for(int i=0;i<width-1;i++){
                    sum=sum1=sum2=0;
                    for(int j=0;j<i-1;j++){
                            sum1+=A[i][j] * X[j];
                    }
                    sum1=-sum1;
                    for(int j=i+1;j<width-1;j++){
                            sum2+=A[i][j] * _X[j];
                    }
                    sum2=-sum2;
                    sum=sum1+sum2+B[i];
                    sum/=A[i][i];
                    X[i]=sum;
                }
                _X=X;
            }
            return X;
        }


};

#endif // METHOD_TIMES_H

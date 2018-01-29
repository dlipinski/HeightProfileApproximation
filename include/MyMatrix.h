#ifndef MYMATRIX_H
#define MYMATRIX_H
#include <iostream>
#include <ctime>
#include <stdio.h>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <fstream>
#include "Eigen/Dense"

/*
    Author: Dawid Lipiński
*/

using namespace std;


class MyMatrix {
 private:
  vector<vector <double> > matrix;
  int height;
  int width;

 double abs(double x) {
        return (x >= 0) ? x : -x;
    }

 public:

    MyMatrix(int _rows, int _cols) {
        matrix.resize(_rows);
        double  _initial=0;
        for (int i=0; i<matrix.size(); i++) {
        matrix[i].resize(_cols, _initial);
        }
        height = _rows;
        width = _cols;
    }

    MyMatrix& operator=(const vector< vector <double> > tab) {
        for(int i=0; i< height;i++)
            for(int j=0; j<width;j++)
                matrix[i][j]=tab[i][j];
        return *this;
    }

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//                                   METODY OBLICZENIOWE
//////////////////////////////////////////////////////////////////////////////////////////////////


//-------------------------------------------------------------------partial Gauss
        double Gauss_partial_time() {
            int n = width-1;
            clock_t begin = clock();
            for (int i = 0; i < n; i++) {
                double maxEl = abs(matrix[i][i]);
                int maxRow = i;
                for (int k = i+1; k < n; k++) {
                    if (abs(matrix[k][i]) > maxEl) {
                        maxEl = abs(matrix[k][i]);
                        maxRow = k;
                    }
                }
                for (int k = i; k < n+1; k++) {
                    double pom = matrix[maxRow][k];
                    matrix[maxRow][k] = matrix[i][k];
                    matrix[i][k] = pom;
                }
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
            vector<double> x(n);
            for (int i=n-1; i>=0; i--) {
                x[i] = matrix[i][n] / matrix[i][i];
                for (int k=i-1;k>=0; k--) {
                    matrix[k][n] -= matrix[k][i] * x[i];
                }
            }
            clock_t end = clock();
            return ( double(end - begin) / CLOCKS_PER_SEC);
        }
//-------------------------------------------------------------------partial Gauss optimized
        double Gauss_partial_optimized_time() {
            int n = width-1;
            clock_t begin = clock();
            for (int i = 0; i < n; i++) {
                double maxEl = abs(matrix[i][i]);
                int maxRow = i;
                for (int k = i+1; k < n; k++) {
                    if (abs(matrix[k][i]) > maxEl) {
                        maxEl = abs(matrix[k][i]);
                        maxRow = k;
                    }
                }
                for (int k = i; k < n+1; k++) {
                    double pom = matrix[maxRow][k];
                    matrix[maxRow][k] = matrix[i][k];
                    matrix[i][k] = pom;
                }
                for (int k = i+1; k < n; k++) {
                    if(matrix[k][i] !=0){  //OPTYMALIZACJA
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
                if((i%100)==0) cout << i << endl;
            }
            vector<double> x(n);
            for (int i=n-1; i>=0; i--) {
                x[i] = matrix[i][n] / matrix[i][i];
                for (int k=i-1;k>=0; k--) {
                    matrix[k][n] -= matrix[k][i] * x[i];
                }
            }
            clock_t end = clock();
            return (double(end - begin) / CLOCKS_PER_SEC);
        }

//-------------------------------------------------------------------Gauss Seidel
        double Siedel_time(){
            int cols = width-1;
            int rows = cols;
            vector < vector < double > > A;
            A.resize(rows);
            for(int i=0;i<cols;i++){
                A[i].resize(cols,0);
                for(int j=0;j<cols;j++)
                    A[i][j]=matrix[i][j];
            }
            vector <double > B;
            B.resize(cols,0);
            for(int i=0;i<cols;i++)
                B[i]=matrix[i][cols];
            vector <long double > X;            // <--- x(k+1)
            X.resize(cols,0);
            vector <long double > _X;           // <--- x(k)
            _X.resize(cols,0);
            long double sum1=0,sum2=0;
            double eps = 0.0000000001;
            int condition_counter = 0;
            bool condition_fulfilled=false;
            int check =0;
            clock_t begin = clock();
            while(true){
                check ++;
                _X=X;
                for(int i=0;i<cols;i++){
                    sum1=0;sum2=0;
                    for(int j=0;j<i-1;j++)
                            sum1+=A[i][j] * X[j];
                    for(int j=i;j<cols;j++)
                            sum2+=A[i][j] * _X[j];
                    X[i]=(-sum1-sum2+B[i])/A[i][i];
                    if( abs(X[i]-_X[i])<=eps ) condition_counter++;
                }
                if(check < 1020 && check > 1000) cout<< X[1] <<"  "<<_X[1] <<endl;
                if(condition_counter >= cols)
                    break;
                else
                    condition_counter = 0;
            }
            clock_t end = clock();
            return (double(end - begin) / CLOCKS_PER_SEC);
        }
//-------------------------------------------------------------------Gauss Seidel Eigen
        double Siedel_Eigen_time(){
            int cols = width-1;
            int rows = cols;

            Eigen::SparseMatrix<double> A(rows,rows);
            for(int i=0;i<rows;i++)
                for(int j=0;j<rows;j++)
                    if(matrix[i][j]!=0)
                        A.insert(i,j)=matrix[i][j];


            Eigen::VectorXd A_Diag(rows);
            A_Diag = A.diagonal();

            Eigen::VectorXd B(rows);
            for(int i=0;i<rows;i++)
                if(matrix[i][rows]!=0)
                    B[i]=matrix[i][rows];


            Eigen::VectorXd X(rows);//<-- k+1
            Eigen::VectorXd _X(rows);//<-- k


            long double sum1=0,sum2=0;
            double eps = 0.0000000001;

            clock_t begin = clock();

            while(true){

                for(int i=0;i<cols;i++){
                    sum1=0;sum2=0;
                    for (Eigen::SparseMatrix<double>::InnerIterator it(A,i); it; ++it){
                        if(it.index()<i)
                            sum1+=it.value() * X[it.index()];
                        if(it.index()>i)
                            sum2+=it.value() * _X[it.index()];
                    }
                    X[i]=(-sum1-sum2+B[i])/A_Diag[i];
                }
                if((X-_X).cwiseAbs().maxCoeff()<eps)
                    break;
                _X=X;

            }

            clock_t end = clock();
            return (double(end - begin) / CLOCKS_PER_SEC);
        }

//-------------------------------------------------------------------Eigen LU
    double Eigen_time(){
        Eigen::SparseMatrix<double> A(width-1,width-1);
        for(int i=0;i<width-1;i++){
            for(int j=0;j<width-1;j++)
                A.insert(i,j) = matrix[i][j];
        }

        Eigen::SparseLU <Eigen::SparseMatrix<double> > solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        Eigen::SparseVector<double> B(width-1);
        for(int i=0;i<width-1;i++)
            B.insert(i) = matrix[i][width-1];

        Eigen::SparseVector<double> ans(width-1);
        clock_t begin = clock();
        ans = solver.solve(B);
        clock_t end = clock();

        return (double(end - begin) / CLOCKS_PER_SEC);
    }

};


vector<double> Eigen_Solve(vector<vector <double> > _A, vector<double> _B){
        Eigen::MatrixXd A(_A.size(),_A.size());
        for(int i=0;i<_A.size();i++)
            for(int j=0;j<_B.size();j++)
                A(i,j)=_A[i][j];
        Eigen::VectorXd B(_B.size());
        for(int i=0;i<_B.size();i++)
            B(i)=_B[i];
        Eigen::VectorXd X=A.partialPivLu().solve(B);
        vector <double> ans(_B.size());
        for(int i=0;i<_B.size();i++) ans[i]=X(i);
        return ans;
    }

#endif // MYMATRIX_H

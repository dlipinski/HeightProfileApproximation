#ifndef APPROXIMATION_H
#define APPROXIMATION_H
#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include "MyMatrix.h"


using namespace std;
double power(double a, int b);

class Approximation
{
    public:
        vector <double> X;
        vector <double> Y;
        int size=0;
        vector <double> polynominal;
        int degree=0;

        Approximation();

        Approximation(vector <double> _X, vector <double> _Y, int _degree){
            X=_X;
            Y=_Y;
            size = _X.size();
            degree=_degree+1;
            polynominal.resize(degree,0);
        };

        void create_polynominal(){
            int s_size = (degree * 2) -1;

            vector <double> s(s_size,0);
            for(int i=0;i<s_size;i++){
                double sum=0;
                for(int j=0; j<size;j++)
                    sum+=power(X[j],i);
                s[i]=sum;
            }

            vector <double> t(degree,0);
            for(int i=0;i<degree;i++){
                double sum=0;
                for(int j=0; j<size; j++)
                    sum+= (Y[j] * power(X[j],i));
                t[i]=sum;
            }

            vector<vector<double> > A; A.resize(degree); for(int i=0; i< degree;i++) A[i].resize(degree,0);
            for(int i=0; i<degree; i++)
                for(int j=0;j<degree;j++)
                    A[i][j]=s[j+i];
            polynominal = Eigen_Solve(A,t);
        }

        void print_polynominal(){
            cout << "y(x) = ";
            for(int i=0;i<degree-1;i++)
                cout <<"("<< polynominal[i] << ")x^"<<i<<" +";
            cout <<"("<< polynominal[degree-1] <<")x^"<<degree-1<< endl;
        }

        string to_string(string arg){
            ostringstream oss;
            oss << arg << ": ";
            oss << "{ y(x) = ";
            for(int i=degree-1;i>0;i--)
                oss <<"("<< polynominal[i] << ") x ^{"<<i<<"} +";
            oss <<"("<< polynominal[0] <<")";
            oss << " }"<< endl;
            string var = oss.str();
            return var;
        }

        double solve(double arg){
            double sum=polynominal[0];
            for(int i=1;i<degree;i++)
                sum+= polynominal[i]*power(arg,i);
            return sum;
        }

};

double power(double a, int b){
    if(b==0)
        return 1;
    double ans=a;
    for(int i=0;i<b-1;i++)
        ans*=a;
    return ans;
}
#endif // APPROXIMATION_H

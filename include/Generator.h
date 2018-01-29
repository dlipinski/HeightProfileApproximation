#ifndef GENERATOR_H
#define GENERATOR_H
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <cstdlib>
#define def_limit 100000
#define prec 5
#include <limits>
#include <unistd.h>
#include <vector>
#include "State.h"
/*
    Author: Dawid Lipiński
*/
int indexOf(vector<int> p, int n){
    for(int i = 0; i < p.size(); i++){
        if(p[i]==n) return i;
    }
    return -1;
}

void print(vector< vector<double> > A, vector<double> B) {
    int n = A.size();
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            (A[i][j] == 0 ) ? printf("%-2d\t", A[i][j]) : printf("%-2.2f\t", A[i][j]);
		}
    }
}

void fillBoard(vector<int> &p, int n){
    int number = -n;
    for(int i = 0; i < 2*n+1; i++){
        p[i] = number++;
    }
}

int mod(int a, int b) { return (a % b + b) % b; }

vector < vector <double > > generate_matrix(int board_size){

    vector < vector <double> > ans;

    ans.resize(board_size);
    for(int i=0;i<board_size;i++)
        ans[i].resize(board_size+1,0);

    for(int i=0;i<board_size;i++)
        for(int j=0;j<board_size;j++)
            if(i==j)
                ans[i][j]=1;

    for (int i=0;i<board_size;i+=(board_size/3)){
        cout << i << endl;
        ans[i][board_size]=0.5;
    }

    return ans;
}




#endif // GENERATOR_H

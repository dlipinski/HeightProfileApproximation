#include <iostream>
#include <ctime>
#include "MyMatrix.h"
#include "Generator.h"
#include "Approximation.h"
#include "MyFunctions.h"
#include "Data.h"
#include <map>
/*
    Author: Dawid Lipiński
    Created: 11.01.2018
    Modified: 14.01.2018
    Description: Program wykonuje pomiary czasów wykonywania się obliczeń macierzy metodami:
                    -Gaussa
                    -Gaussa zoptymalizowanego
                    -Siedla
                    -SpraseMatrix
                 Nastepnie za pomocą aproksymacji średniokwadratowej dyskretnej oraz danych zebranych wcześniej, tworzy wielomiany aproksymacyjne.
                 Umożliwiają one obliczenie czasu wykonywania działania daną metodą na podstawie danego rozmiaru macierzy.
*/


typedef std::numeric_limits< double > dbl;

using namespace std;


int main()
{
    Data data;
    string types[3] = {"flat","one_in_middle","two_on_edges"};
    int rarities[5]  = {1,2,10,20,50};
    double indexes[5]  = {0,0.25,0.5,0.75,1};
    //WPYW LICZBY PUNKTÓW WĘZŁOWYCH NA WYNIKI
    for(int type=0;type<3;type++){
        stringstream ss;
        ss <<"errors/"<< types[type] << ".csv";
        string errors_file_name = ss.str();
        stringstream sss;
        sss <<"functions/"<< types[type] << ".csv";
        cout << sss.str()<<endl;
        string functions_file_name = sss.str();
        ofstream errors;
        ofstream functions;
        errors.open (errors_file_name);
        functions.open (functions_file_name);
        errors << "rarity;index;error"<<endl;
        functions <<"real;1;2;10;20;50"<<endl;

        for(int rarity=0;rarity<5;rarity++){
            for(int index=0;index<5;index++){
                double first = data.get_first_raw(types[type],indexes[index]);
                double real_second = data.get_real_second(types[type],first);
                double calc_second = data.lagrange_part(rarities[rarity],types[type],first);
                double error = abs(real_second - calc_second);
                errors << rarities[rarity] <<";"<<indexes[index]<<";"<<error<<endl;
            }
        }
         for(double a=0;a<1;a+=0.01){
                double first = data.get_first_raw(types[type],a);
                double real_second = data.get_real_second(types[type],first);
                functions<<real_second;
                for(int rarity=0;rarity<5;rarity++){
                    double calc_second = data.lagrange_part(rarities[rarity],types[type],first);
                    functions <<";"<<calc_second;
                }
                functions<<endl;
            }

        errors.close();
        functions.close();
    }

    return 0;

}




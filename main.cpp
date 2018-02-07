#include <iostream>
#include <ctime>
#include <map>
#include "Data.h"

/*
    Author: Dawid Lipiński
    Created: 01.02.2018
    Modified: 07.02.2018
    Description:
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
        ss <<"errors/lagrange_"<< types[type] << ".csv";
        string l_errors_file_name = ss.str();

        ss.str(string());
        ss <<"functions/lagrange_"<< types[type] << ".csv";
        string l_functions_file_name = ss.str();

        ss.str(string());
        ss <<"errors/splines_"<< types[type] << ".csv";
        string s_errors_file_name = ss.str();

        ss.str(string());
        ss <<"functions/splines_"<< types[type] << ".csv";
        string s_functions_file_name = ss.str();

        ofstream l_errors;
        ofstream l_functions;
        ofstream s_errors;
        ofstream s_functions;

        l_errors.open (l_errors_file_name);
        l_functions.open (l_functions_file_name);
        l_errors << "rarity;index;error"<<endl;
        l_functions <<"real;1;2;10;20;50"<<endl;

        s_errors.open (s_errors_file_name);
        s_functions.open (s_functions_file_name);
        s_errors << "rarity;index;error"<<endl;
        s_functions <<"real;1;2;10;20;50"<<endl;

        for(int rarity=0;rarity<5;rarity++){
            for(int index=0;index<5;index++){
                double first = data.get_first_raw(types[type],indexes[index]);
                double real_second = data.get_real_second(types[type],first);

                double l_calc_second = data.lagrange_part(rarities[rarity],types[type],first);
                double s_calc_second = data.splines_part(rarities[rarity],types[type],first);

                l_errors << rarities[rarity] <<";"<<indexes[index]<<";"<<abs(real_second - l_calc_second)<<endl;
                s_errors << rarities[rarity] <<";"<<indexes[index]<<";"<<abs(real_second - s_calc_second)<<endl;
            }
        }
         for(double a=0;a<1;a+=0.01){
                double first = data.get_first_raw(types[type],a);
                double real_second = data.get_real_second(types[type],first);

                l_functions<<real_second;
                s_functions<<real_second;

                for(int rarity=0;rarity<5;rarity++){

                    l_functions <<";"<<data.lagrange_part(rarities[rarity],types[type],first);
                    s_functions <<";"<<data.splines_part(rarities[rarity],types[type],first);
                }
                l_functions<<endl;
                s_functions<<endl;
            }

        l_errors.close();
        l_functions.close();
        s_errors.close();
        s_functions.close();

    }

    return 0;

}




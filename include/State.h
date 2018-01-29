#ifndef STATE_H
#define STATE_H
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <set>
#include <map>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
using namespace std;

class State
{
    public:
        int player, x, y, wsp, g1, g2,s1,s2, win_loose=-1;
        double wage;
        map<int, bool> mushrooms;

        State();

        State(int p, int x_, int y_, int wsp_, int g_1, int g_2, map<int,bool> gameField,int n,double _wage){
            player = p;
            x = x_;
            y = y_;
            wsp = wsp_;
            s1=-n;
            s2=n;
            mushrooms = gameField;
            g1 = g_1;
            g2 = g_2;
            wage=_wage;
        }

        void print(){
            string result;

            result = "("+string_me(wage)+")" + "P" + to_string(player) + "(" + to_string(x) + ", " + to_string(y) + ")";

        }
        string string_me(double a){
            ostringstream strs;
            strs << a;
            string str = strs.str();


            return str;
        }
        string int_to_string(int a){
            ostringstream strs;
            strs << a;
            string str = strs.str();


            return str;
        }
        string to_stringg(int _zero_one){
            string result;
            if(_zero_one==2)
                result = "("+string_me(wage)+")" + "P" + to_string(player) + "(" + to_string(x) + ", " + to_string(y) + ")" ;
            else
                result = "("+string_me(wage)+")" + int_to_string(_zero_one);
            return result;
        }

        string to_stringasdg(){
            string result;
            result = "("+string_me(wage)+")" + "P" + to_string(player) + "(" + to_string(x) + ", " + to_string(y) + ")" + "[" + to_string(g1) + "," + to_string(g2) + "]" ;


            return result;
        }
        string toString(){
            string result;

            result = to_string(player) + to_string(x) + to_string(y);
            for(int i = s1; i<= s2; i++){
                if(mushrooms[i])
                    result += to_string(i);
            }
            return result;
        }

        int status(){
            if(x == 0){
                if(g1 > g2){
                    return 1;
                } else if(g1 < g2){
                    return 0;
                } else {
                    return 1;
                }
            }
            else if(y == 0){
                if(g1 > g2){
                    return 1;
                } else if(g1 < g2){
                    return 0;
                } else {
                    return 0;
                }
            }
            else{
                return -1;
            }
        }



};


#endif // STATE_H

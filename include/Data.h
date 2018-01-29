#ifndef DATA_H
#define DATA_H
#include <iostream>
#include <map>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
using namespace std;

vector<string> split(const string &s, char delim);

class Data
{
    public:
        map<double,double> flat;
        map<double,double> one_in_middle;
        map<double,double> two_on_edges;

        //Kontruktor wczytuje dane z plików
        Data(){
            string line;
            ifstream file;
            file.open("data/flat.csv");
            for( string line; getline( file, line ); )
            {
               vector<string> datas = split(line,',');
               flat[stod(datas[0])] = stod(datas[1] );
            }
            file.close();

            file.open("data/one_in_middle.csv");
            for( string line; getline( file, line ); )
            {
               vector<string> datas = split(line,',');
               one_in_middle[stod(datas[0])] = stod(datas[1] );
            }
            file.close();

            file.open("data/two_on_edges.csv");
            for( string line; getline( file, line ); )
            {
               vector<string> datas = split(line,',');
               two_on_edges[stod(datas[0])] = stod(datas[1] );
            }
            file.close();
        }

        //Zwraca przerzedzoną mapę, biorąc co number element głównej
        map<double,double> submap(string data_type, int number){
            map<double,double> answer;
                int i=0;
                if(data_type =="flat"){
                    for (map<double, double>::iterator it = flat.begin(); it != flat.end(); it++ )
                        {
                            if(i++ % number == 0)
                                answer[it->first]=it->second;
                        }
                }
                if(data_type =="one_in_middle"){
                    for (map<double, double>::iterator it = one_in_middle.begin(); it != one_in_middle.end(); it++ )
                        {
                            if(i++ % number == 0)
                                answer[it->first]=it->second;
                        }
                }
                if(data_type =="two_on_edges"){
                    for (map<double, double>::iterator it = two_on_edges.begin(); it != two_on_edges.end(); it++ )
                        {
                            if(i++ % number == 0)
                                answer[it->first]=it->second;
                        }
                }
            return answer;
        }

        //Dzieli stringa jak split
        vector<string> split(const string &s, char delim) {
            stringstream ss(s);
            string item;
            vector<string> tokens;
            while (getline(ss, item, delim)) {
                tokens.push_back(item);
            }
            return tokens;
        }

        //Zwraca odległość znajdująca się w pewnej części map
        double get_first_raw(string data_type,double a){
            int b = flat.size();
            int i=0;
            int c = a * b;
            if(data_type=="flat")
                for (map<double, double>::iterator it = flat.begin(); it != flat.end(); it++ )
                        {
                            if(i++  > c)
                                return it->first;
                        }
            if(data_type=="one_in_middle")
                for (map<double, double>::iterator it = one_in_middle.begin(); it != one_in_middle.end(); it++ )
                        {
                            if(i++  >c)
                                return it->first;
                        }
            if(data_type=="two_on_edges")
                for (map<double, double>::iterator it = two_on_edges.begin(); it != two_on_edges.end(); it++ )
                        {
                            if(i++  > c)
                                return it->first;
                        }
            return 0;
        }

        //Zwraca prawdziwą wysokośc dla danej odległości
        double get_real_second(string data_type,double first){
                if(data_type=="flat"){
                    return flat[first];
                }
                if(data_type=="one_in_middle"){
                    return one_in_middle[first];
                }
                if(data_type=="two_on_edges"){
                    return two_on_edges[first];
                }
                return 0;
        }

        //Zwraca obliczoną wysokośc lagrangiem na danej odległości
        double lagrange_part(int part,string data_type, double first){
                double result = 0; // Initialize result
                map<double,double> to_calc;
                int m = 0;
                to_calc = submap(data_type,part);
                for (map<double, double>::iterator i = to_calc.begin(); i != to_calc.end(); i++ )
                    {
                        double term = i->second;
                        for (map<double, double>::iterator j = to_calc.begin(); j != to_calc.end(); j++ )
                        {
                            if (j!=i)
                                term = term*(first - j->first)/double(i->first - j->first);
                        }
                        result += term;
                    }
            return result;
        }
};

#endif // DATA_H

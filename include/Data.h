#ifndef DATA_H
#define DATA_H
#include <iostream>
#include <map>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "Eigen/Dense"

using namespace std;

vector<string> split(const string &s, char delim);

vector<double> Eigen_Solve_me(vector<vector <double> > _A, vector<double> _B){
        Eigen::MatrixXd A(_A.size(),_A.size());
        for(unsigned int i=0;i<_A.size();i++)
            for(unsigned int j=0;j<_B.size();j++)
                A(i,j)=_A[i][j];
        Eigen::VectorXd B(_B.size());
        for(unsigned int i=0;i<_B.size();i++)
            B(i)=_B[i];
        Eigen::VectorXd X=A.partialPivLu().solve(B);
        vector <double> ans(_B.size());
        for(unsigned int i=0;i<_B.size();i++) ans[i]=X(i);
        return ans;
    }
double power(double a, int b){
    if(b==0)
        return 1;
    double ans=a;
    for(int i=0;i<b-1;i++)
        ans*=a;
    return ans;
}

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
        double lagrange_part(int part,string data_type, double arg){
                double result = 0; // Initialize result
                map<double,double> to_calc;
                to_calc = submap(data_type,part);
                for (map<double, double>::iterator i = to_calc.begin(); i != to_calc.end(); i++ )
                    {
                        double term = i->second;
                        for (map<double, double>::iterator j = to_calc.begin(); j != to_calc.end(); j++ )
                        {
                            if (j!=i)
                                term *=(arg - j->first)/double(i->first - j->first);
                        }
                        result += term;
                    }
            return result;
        }

        double splines_part(int part, string data_type, double arg){
                double result = 0; // Initialize result
                map<double,double> to_calc;
                to_calc = submap(data_type,part);
                int n = to_calc.size();
                int it =0;
                for (map<double, double>::iterator i = to_calc.begin(); next(i,1) != to_calc.end(); i++ )
                {
                    if (i->first == arg)
                        return i->second;
                    else
                        if ((i->first < arg) && (next(i,1)->first > arg)){

                        double argi = i->first;
                        double argi_n = next(i,1)->first;
                        double vali = i->second;
                        double vali_n = next(i,1)->second;

                        double argg = argi_n - argi;
                        double Vi = 0;
                        double Vi_n = 0;

                        vector < vector <double> > V;
                        vector <double> TEMP;

                        vector <double> BS;
                        BS.resize(n,0);
                        V.resize(n);
                        for(unsigned int a=0;a<V.size();a++)
                            V[a].resize(n,0);
                        TEMP.resize(n,0);


                        V[0][0]=2;
                        V[0][1]=0;
                        V[0][n]=0;

                        int temp=0;
                        for (map<double, double>::iterator j = next(to_calc.begin(),1); next(j,1) != to_calc.end(); j++ ){

                            double _arggj = j->first - prev(j)->first;
                            double _arggj_n = next(j,1)->first - j->first;

                            double _valj = j->second;
                            double _valj_n = next(j,1)->second;
                            double _valj_p = prev(j) ->second;

                            V[temp][temp-1] = _arggj/ (_arggj + _arggj_n);
                            V[temp][temp]=2.0;
                            V[temp][temp+1] = _arggj_n / (_arggj + _arggj_n);
                            BS[temp] = 6/(_arggj + _arggj_n) * (((_valj_n-_valj)/(_arggj_n))-((_valj-_valj_p)/(_arggj)));
                            temp++;

                        }

                        V[n-1][n-2]=0;
                        V[n-1][n-1]=0;
                        V[n-1][n]=0;

                        TEMP = Eigen_Solve_me(V,BS);

                        Vi = TEMP[it];
                        Vi_n = TEMP[it+1];

                        result = Vi * (  (pow(argi_n-arg,3) / (6*argg))  )
                                    + Vi_n * (  pow(arg-argi,3) / (6*argg)  )
                                    + (((vali_n-vali)/argg) - ((argg*(Vi_n-Vi)/6)) * (arg - argi))
                                    + vali - Vi * ((arg*arg)/6);
                        return result;
                    }
                    it++;

                }
                return 0;
        }
};



#endif // DATA_H

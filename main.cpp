//#pragma GCC optimize(3,"inline","Ofast")
//version 0.7.2
#include <iostream>
#include "computational.h"
#include<chrono>
#include<string>
using namespace std;
using namespace malg;
using namespace chrono;
int main()
{
    Matrix<double> A({{0, 2}, {2, 1}});
    function<double(double,double)> g = [](double x, double y) { return y; };
    function<Matrix<double>(double,Matrix<double>)> f= [&](double x, Matrix<double> y) 
    { return A*y+Matrix<double>({{3*x},{0}},2,1); };
    double x0=0;
    Matrix<double> y0({{1},{1}},2,1);
    double x =0.1;
    cout<<DSolve(f, {x0,y0},x,"modified Euler")<<endl;
    cout << DSolve(g, {0.0, 1.0}, 1.0) << endl;
    system("pause");
    return 0;
}
                  
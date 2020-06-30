//#pragma GCC optimize(3,"inline","Ofast")
//version 0.7.2
#include <iostream>
#include "computational.h"
#include<chrono>
#include<ctime>
#include<string>
using namespace std;
using namespace malg;
using namespace chrono;
int main()
{
    /*Matrix<double> A({{0, 2}, {2, 1}});
    function<double(double,double)> g = [](double x, double y) { return y; };
    function<Matrix<double>(double,Matrix<double>)> f= [&](double x, Matrix<double> y) 
    { return A*y+Matrix<double>({{3*x},{0}},2,1); };
    double x0=0;
    Matrix<double> y0({{1},{1}},2,1);
    double x =0.1;
    cout<<DSolve(f, {x0,y0},x,"modified Euler")<<endl;
    cout << DSolve(g, {0.0, 1.0}, 1.0) << endl;*/

    /*auto f = [](Matrix<double> x) { return x(0,0) * x(0,0) + exp(x(1,0)); };
    cout << D<double>(f, 0)(Matrix<double>({{5}, {1}},2,1))<<endl;
    cout << D<double>(f, 1)(Matrix<double>({{10}, {1}},2,1))<<endl;*/

    auto f1=[](Matrix<double> x) { return x(0) * x(0) -x(1)-1; };
    auto f2=[](Matrix<double> x) { return x(0) * x(0) -4*x(0)+x(1)*x(1)-x(1)+3.25; };
    Matrix<double> x0({{0}, {0}}, 2, 1);
    clock_t st, et;
    st = clock();
    for (int i = 0; i < 1000;++i)
        FindRoot<double>({f1, f2}, x0, "Newton");
    et = clock();
    cout << "Newton " << et - st << "ms" << endl;
    st = clock();
    for (int i = 0; i < 1000;++i)
        FindRoot<double>({f1, f2}, x0, "Broyden");
    et = clock();
    cout << "Broyden " << et - st << "ms" << endl;
    system("pause");
    return 0;
}
                            
//#pragma GCC optimize(3,"inline","Ofast")
//version 0.7.5
#include <iostream>
//#include "algebra.h"
#include "computational.h"
#include<ctime>
#include <cmath>
using namespace std;
using namespace malg;
int main()
{
    //auto f = [](double t, double y) { return 1-t+4*y; };
    //cout << DSolve<double, double>(f, {0, 1}, 0.5,10,"backward Euler") << endl;
    Matrix<double> A({{1,0},{0,1e-10}},"");
    cout << Cond(A, 2) << endl;
    system("pause");
    return 0;
}
                                    
                    
                               
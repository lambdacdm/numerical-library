//#pragma GCC optimize(3,"inline","Ofast")
//version 0.7.4
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
    Matrix<double> A({{12, -51, 4}, {6, 167, -68}, {-4, 24, -41}});
    cout << A << endl;
    auto Slist = SchurDecomposition(A);
    cout << Slist[0] << endl;
    cout << Slist[1] << endl;
    cout << Slist[0]*Slist[1]*Slist[2] << endl;
    clock_t st, et;
    st = clock();
    cout << EigenValueMatrix(A)<< endl;
    et = clock();
    cout << et - st << "ms" << endl;
    system("pause");
    return 0;
}
                                    
                    
                         
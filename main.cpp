//#pragma GCC optimize(3,"inline","Ofast")
//version 0.7.3
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
    //Matrix<double> A({{2.3, 4.5}, {6.7, -1.2}});
    //auto A = 2.0*Eye<double>(10);
    cout << A << endl;
    clock_t st, et;
    st = clock();
    cout << EigenValueMatrix(A) << endl;
    et = clock();
    cout << et - st << "ms" << endl;
    system("pause");
    return 0;
}
                                    
                    
                
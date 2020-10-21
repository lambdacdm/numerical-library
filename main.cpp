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
    Matrix<double> A({{-1, -1}, {2, 2}}, "");
    cout << EigenValueMatrix(A) << endl;
    system("pause");
    return 0;
}
                                    
                    
                               
//#pragma GCC optimize(3,"inline","Ofast")
//version 0.7.5
#include <iostream>
#include "computational.h"
#include <cmath>
using namespace std;
using namespace malg;
int main()
{
    Plot<double>([](double x) { return sin(x); },{0,2*Pi<double>});
    //绘制sinx在x属于[0,2pi]的图像
    system("pause");
    return 0;
}
                                    
                    
                                  
//#pragma GCC optimize(3,"inline","Ofast")
//version 0.7.1
#include <iostream>
#include "computational.h"
#include<chrono>
#include<string>
using namespace std;
using namespace malg;
using namespace chrono;
int main()
{
    Polynomial<int> f({0,1,2,3});
    Polynomial<int> g({5,3,7,9});
    cout << f * g << endl;
    system("pause");
    return 0;
}
                  
//#pragma GCC optimize(3,"inline","Ofast")
//version 0.7.1
#include <iostream>
#include "computational.h"
#include<chrono>
#include<string>
using namespace std;
using namespace malg;
using namespace chrono;
void test(int n)
{
    auto A=Generate<double>(n,n);
    auto B = Generate<double>(n,n);
    clock_t st,et;
    st = clock();
    A *B;
    et = clock();
    cout << "单线程 " << et - st << "ms" << endl;
    st = clock();
    ThreadMatProd(A, B);
    et = clock();
    cout << "多线程 " << et - st << "ms" << endl;
}
void exam(int n)
{
    srand(time(0));
    for(int i=1;i<=n;++i)
    {
        auto A=Generate<double>(i,i);
        auto B=Generate<double>(i,i);
        auto C_1 = A * B;
        auto C_2=ThreadMatProd(A,B);
        if(C_1==C_2)
            cout << "通过" << endl;
        else
        {
            cout << "不通过" << endl;
            cout << A << endl;
            cout << B << endl;
            cout << C_1 << endl;
            cout << C_2 << endl;
            break;
        }
    }
}
int main()
{
    srand(time(0));
    int n = 800;
    auto a=Generate<double>(n,n);
    auto e = Eye<double>(n);
    auto st = system_clock::now();
    LinearSolve(a, e, "LUP");
    auto et = system_clock::now();
    auto dur = duration<double>(et - st);
    cout << dur.count()<< endl;
    system("pause");
    return 0;
}
                
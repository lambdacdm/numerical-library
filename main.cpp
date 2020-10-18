//#pragma GCC optimize(3,"inline","Ofast")
//version 0.7.4
#include <iostream>
#include "computational.h"
#include<ctime>
using namespace std;
using namespace malg;
int main()
{
    int n = 4096;
    Matrix<double> A = Generate<double>(n, n);
    cout << "order=" << n << endl;
    int start, end;

    start = clock();
    auto LU = LUDecomposition(A);
    end = clock();
    auto L = LU[0];
    auto U = LU[1];
    cout << "Ordinary LU time=" << double(end - start) / 1000 << "s" << endl;
    cout << "Ordinary LU error=" << Norm(A-L*U, 1)<<endl;

    start = clock();
    LU = LUDivideConquer(A);
    end = clock();
    L = LU[0];
    U = LU[1];
    cout << "Divide Conquer LU time=" << double(end - start) / 1000 << "s" << endl;
    cout << "Divide Conquer LU error=" << Norm(A-L*U, 1)<<endl;

    start = clock();
    LU = LUPDecomposition(A);
    end = clock();
    L = LU[0];
    U = LU[1];
    auto P = LU[2];
    cout << "LUP time=" << double(end - start) / 1000 << "s" << endl;
    cout << "LUP error=" << Norm(A-Transpose(P)*L*U, 1)<<endl;

    system("pause");
    return 0;
}
                                    
                    
                               
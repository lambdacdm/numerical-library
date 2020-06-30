#ifndef COMPUTATIONAL_H
#define COMPUTATIONAL_H
#include<iostream>
#include<vector>
#include<cmath>
#include<iomanip>
#include <typeinfo>
#include<functional>
#include<cstdlib>
#include<complex>
#include<algorithm>
#include<deque>
#include<thread>
#include<execution>
#include<utility>
using std::cerr;
using std::complex;
using std::cout;
using std::deque;
using std::get;
using std::ios;
using std::ostream;
using std::pair;
using std::ref;
using std::string;
using std::swap;
using std::thread;
using std::vector;

namespace malg{

//-----常量部分-----
template<class DC> const DC Pi = 3.1415926535897932;
template <class DC>
const std::function<DC(DC)> Identity =[](DC x) { return x; };
template<class DC>
const std::function<DC(DC, DC)> Multiply = [](DC x, DC y) { return x * y; };

//-----声明部分-----

//杂例I
template <class DM>
double Norm(const vector<DM> &, double);
template<class DJ> DJ Abs(DJ);
template <class DJ> vector<DJ> Range(DJ, DJ,DJ);
template <class DJ> vector<DJ> Range(DJ, DJ);
template <class DJ> vector<DJ> Rand(int);
template <class DJ> DJ Limit(std::function<DJ(DJ)>, DJ);
template <class DJ> DJ N(DJ,int);
inline int Sig(){return rand()%2?1:-1;}//random sign(+/-)

//大整数类
class BigInt
{
    public:
        BigInt cutzero();
        BigInt();
        BigInt(string);
        BigInt(int);
        BigInt(bool);

        friend BigInt operator-(const BigInt &);
        friend bool operator==(const BigInt&, const BigInt&);
        friend bool operator!=(const BigInt &, const BigInt &);
        friend bool operator<(const BigInt&, const BigInt&);
        friend bool operator>(const BigInt&, const BigInt&);
        friend bool operator<=(const BigInt&, const BigInt&);
        friend bool operator>=(const BigInt&, const BigInt&);
        friend ostream &operator << (ostream &, const BigInt&);
        friend BigInt operator>>(const BigInt&, int);
        friend BigInt operator+(const BigInt &,const BigInt&);
        friend BigInt operator-(const BigInt &,const BigInt&);
        friend BigInt operator*(const BigInt &,const BigInt&);
        friend BigInt operator/(const BigInt&, const BigInt&);
        friend BigInt operator^(const BigInt&, int);
        BigInt &operator+=(const BigInt &);
        BigInt &operator-=(const BigInt &);

        explicit operator int();
        void Reassign(string );
        void _digitChange(string);
        void _signChange(bool);
        friend string Get_digit(const BigInt&);
        friend string GetString(const BigInt &);
        friend bool Sign(const BigInt&);
        friend unsigned Length(const BigInt &);

    private:
        string _digit;
        bool _sign;
};
vector<long long> CompressBit(const string&);
string BitToString(const vector<long long> &);

class HighPrecision: public BigInt
{
    public:
        HighPrecision CutTail();
        HighPrecision();
        HighPrecision(string);
        HighPrecision(const char[]);
        explicit HighPrecision(const BigInt&);
        HighPrecision(double);
        HighPrecision(string, int);
        HighPrecision(double, int);

        friend ostream &operator<<(ostream &, const HighPrecision &);
        friend HighPrecision operator<<(const HighPrecision &,int);
        friend HighPrecision operator+(const HighPrecision &, const HighPrecision &);
        friend HighPrecision operator- (const HighPrecision &);
        friend HighPrecision operator-(const HighPrecision &, const HighPrecision &);
        friend HighPrecision operator*(const HighPrecision &, const HighPrecision &);
        friend HighPrecision operator/(const HighPrecision &, const HighPrecision &);

        friend unsigned DecimalLength(const HighPrecision &);
        friend BigInt IntegerPart(const HighPrecision &);
        friend HighPrecision DecimalPart(const HighPrecision &);
        friend int Order(const HighPrecision &);
        friend string TopKDigit(const HighPrecision &,unsigned);
        friend unsigned TotalLength(const HighPrecision &);
        friend unsigned SignificantLength(const HighPrecision &);
        friend HighPrecision SignificantFigure(const HighPrecision &, unsigned);

    private:
        string _decimal;
};

//多项式类
template<class DC>
class Polynomial{
public:
    Polynomial();
    Polynomial(const vector<DC>&);
    Polynomial(const vector<DC>&, int);
    Polynomial(DC);

    Polynomial<DC> &operator=(const Polynomial<DC>&);
    DC operator()(DC);
    template <class DD> friend Polynomial<DD> operator+(const Polynomial<DD> &,const Polynomial<DD> &);
    template <class DD> friend Polynomial<DD> operator-(const Polynomial<DD> &,const Polynomial<DD> &);
    template <class DD> friend Polynomial<DD> operator+(DD, const Polynomial<DD>&);
    template <class DD> friend Polynomial<DD> operator*(DD, const Polynomial<DD>&);
    template <class DD> friend Polynomial<DD> operator*(const Polynomial<DD>&,const Polynomial<DD>&);
    template <class DD> friend Polynomial<DD> operator/(const Polynomial<DD>&,const Polynomial<DD>&);
    template <class DD> friend Polynomial<DD> operator^(const Polynomial<DD>&,int);
    template<class DD> friend ostream &operator << (ostream &,const Polynomial<DD> &);
    template <class DD> friend Polynomial<DD> operator>>(const Polynomial<DD>&, int);
    Polynomial<DC> &operator+=(const Polynomial<DC>&);
    Polynomial<DC> &operator-=(const Polynomial<DC>&);

    template <class DD> friend DD Get(const Polynomial<DD> &, DD);
    template <class DD> friend vector<DD> GetCoef(const Polynomial<DD> &);
    template <class DD> friend void Disp(const Polynomial<DD>&);
    template <class DD> friend int Deg(const Polynomial<DD>&);
    template <class DD> friend Polynomial<DD> SubPoly(const Polynomial<DD>&, int, int);
    template <class DD> friend Polynomial<DD> D(const Polynomial<DD>&);
    template <class DD> friend Polynomial<DD> Integrate(const Polynomial<DD> &);
    template <class DD> friend DD Integrate(const Polynomial<DD> &,DD,DD);
    template <class DD> friend Polynomial<DD> Schmidt(int,DD,DD);
    template <class DD> friend Polynomial<DD> Legendre(int n);

private:
    int poly_deg;
    vector<DC> poly_coeff;
};
Polynomial<double> operator*(const Polynomial<double>&,const Polynomial<double>&);
Polynomial<float> operator*(const Polynomial<float>&,const Polynomial<float>&);
Polynomial<complex<double>> operator*(const Polynomial<complex<double>> &, const Polynomial<complex<double>> &);
Polynomial<complex<float>> operator*(const Polynomial<complex<float>> &, const Polynomial<complex<float>> &);
Polynomial<unsigned> operator*(const Polynomial<unsigned> &, const Polynomial<unsigned> &);
Polynomial<unsigned long long> operator*(const Polynomial<unsigned long long> &, const Polynomial<unsigned long long> &);


//矩阵类
template<class DM>
class Matrix
{
public:
    Matrix();
    Matrix(DM);
    Matrix(int,int);
    Matrix(const vector<vector<DM>>&, int, int);
    Matrix(const vector<vector<DM>>&);
    Matrix(const vector<DM> &);

    Matrix<DM> &operator=(const Matrix<DM> &);
    DM &operator()(int,int);
    const vector<DM> &operator[](int);
    template <class DB> friend bool operator==(const Matrix<DB> &, const Matrix<DB> &);
    template<class DB> friend ostream &operator<< (ostream &,const Matrix<DB> &);
    template<class DB> friend Matrix<DB> operator+(const Matrix<DB> &,const Matrix<DB>&);
    template<class DB> friend Matrix<DB> operator-(const Matrix<DB> &,const Matrix<DB>&);
    template <class DB>friend Matrix<DB> operator-(const Matrix<DB> &);
    template<class DB> friend Matrix<DB> operator*(const Matrix<DB> &,const Matrix<DB>&);
    template<class DB> friend Matrix<DB> operator/(const Matrix<DB> &,const Matrix<DB>&);
    template<class DB> friend Matrix<DB> operator^(const Matrix<DB> &,int);
    template <class DB> friend Matrix<DB> operator+(DB, const Matrix<DB>&);
    template <class DB> friend Matrix<DB> operator+(const Matrix<DB>&,DB);
    template <class DB> friend Matrix<DB> operator-(DB, const Matrix<DB>&);
    template <class DB> friend Matrix<DB> operator-(const Matrix<DB>&,DB);
    template <class DB> friend Matrix<DB> operator*(DB, const Matrix<DB>&);
    template <class DB> friend Matrix<DB> operator*(const Matrix<DB>&,DB);
    template<class DB> friend Matrix<DB> operator/(const Matrix<DB> &,DB);
    Matrix<DM> &operator+=(const Matrix<DM> &);
    Matrix<DM> &operator-=(const Matrix<DM> &);

    template <class DB> friend DB Get(const Matrix<DB> &, int, int);
    template<class DB> friend void Clear(Matrix<DB>&);
    template<class DB> friend void Disp(const Matrix<DB>&);
    template <class DB> friend Matrix<DB> Generate(int, int);
    template <class DB> friend Matrix<DB> Act(std::function<DB(DB)>, const Matrix<DB> &);
    template <class DB> friend Matrix<DB> Act(std::function<DB(DB, DB)>, const Matrix<DB> &, const Matrix<DB> &);
    template <class DB> friend Matrix<DB> Table(std::function<DB(DB)>, const vector<DB>&);
    template<class DB> friend Matrix<DB> Transpose(const Matrix<DB>&); 
    template<class DB> friend Matrix<DB> Diagonal(const Matrix<DB>&);
    template <class DB> friend Matrix<DB> DiagonalMatrix(const vector<DB> &);
    template<class DB> friend int RowSize(const Matrix<DB>&);
    template<class DB> friend int ColumnSize(const Matrix<DB>&);
    template <class DB> friend void Resize(Matrix<DB> &,int,int);
    template <class DB> friend Matrix<DB> SubRow(const Matrix<DB> &, int, int);
    template <class DB> friend Matrix<DB> SubMatrix(const Matrix<DB> &, int, int, int, int);
    template <class DB> friend void DeleteRow(Matrix<DB> &,int n);
    template <class DB> friend void ReplaceRow(Matrix<DB> &, int n,const Matrix<DB>&);
    //解线性方程组
    template<class DB> friend Matrix<DB> RowCat(const Matrix<DB>&, const Matrix<DB>&);
    template <class DB> friend Matrix<DB> RowCat(const vector<Matrix<DB>>&);
    template <class DB> friend Matrix<DB> ColumnCat(const Matrix<DB> &, const Matrix<DB> &);
    template<class DB> friend void SwapRow(Matrix<DB>&,int,int);
    template <class DB> friend Matrix<DB> USplit(const Matrix<DB>&);
    template <class DB> friend vector<Matrix<DB>> GaussianElimination(const Matrix<DB>&,const Matrix<DB>&);
    template <class DB> friend Matrix<DB> LUCompactDecomposition(const Matrix<DB>&);
    template <class DB> friend Matrix<DB> LUPCompactDecomposition(const Matrix<DB>&, vector<int> &, int&, DB&);
    template <class DB> friend Matrix<DB> GaussianElimination(const Matrix<DB>&);
    template <class DB> friend int MatrixRank(const Matrix<DB>&);
    template <class DB> friend DB Det(const Matrix<DB>&);
    template <class DB> friend vector<Matrix<DB>> LUDecomposition(const Matrix<DB>&);
    template <class DB> friend vector<Matrix<DB>> LUPDecomposition(const Matrix<DB> &);
    template <class DB> friend vector<vector<DB>> Tridiagonal(const Matrix<DB>&);
    template <class DB> friend vector<vector<DB>> TridiagonalSeparation(const vector<vector<DB>>&);
    template <class DB> friend Matrix<DB> CholeskyCompactDecomposition(const Matrix<DB>&);
    template <class DB> friend vector<Matrix<DB>> CholeskyDecomposition(const Matrix<DB> &);
    template <class DB> friend vector<Matrix<DB>> JacobiIterationMatrix(const Matrix<DB>&, const Matrix<DB>&);
    template <class DB> friend vector<Matrix<DB>> GSIterationMatrix(const Matrix<DB>&, const Matrix<DB>&);
    template <class DB> friend Matrix<DB> LUSolve(const Matrix<DB>&,const Matrix<DB>&);
    template <class DB> friend Matrix<DB> LUPSolve(const Matrix<DB>&,const vector<int>&,const Matrix<DB>&);
    template <class DB> friend Matrix<DB> LSolve(const Matrix<DB>&,const Matrix<DB>&);
    template <class DB> friend Matrix<DB> USolve(const Matrix<DB>&,const Matrix<DB>&);
    template <class DB> friend Matrix<DB> TridiagonalSolve(const vector<vector<DB>>&, const Matrix<DB>&);
    template <class DB> friend Matrix<DB> Iteration(const Matrix<DB>&, const Matrix<DB>&);
    template<class DB> friend Matrix<DB> LinearSolve(const Matrix<DB>&,const Matrix<DB>&,string,double);
    template<class DB> friend Matrix<DB> LinearSolve(const Matrix<DB>&,const Matrix<DB>&,string);
    template<class DB> friend Matrix<DB> LinearSolve(const Matrix<DB>&,const Matrix<DB>&);
    //求逆
    template<class DB> friend Matrix<DB> Eye(int);
    template<class DB> friend Matrix<DB> Inverse(const Matrix<DB>&);
    //范数、条件数
    template <class DB> friend DB Trace(const Matrix<DB>&);
    template <class DB> friend double Norm(const Matrix<DB>&, double);
    template <class DB> friend double Norm(const Matrix<DB>&, string);
    template <class DB> friend double Cond(const Matrix<DB>&, double);
    template <class DB> friend double Cond(const Matrix<DB>&, string);
    //插值
    template <class DB> friend vector<std::function<DB(DB)>> LagrangeBasis(const Matrix<DB> &);
    template <class DB> friend DB LBasisToInterp(const vector<std::function<DB(DB)>>&,const Matrix<DB> &,DB);
    template <class DB> friend Matrix<DB> DifferenceQuotientTable(const Matrix<DB> &);
    template <class DB> friend Matrix<DB> DividedDifferenceTable(const Matrix<DB> &);
    template <class DB> friend DB DifferenceQuotient(const Matrix<DB> &);
    template <class DB> friend DB DifferenceQuotient(std::function<DB(DB)>, const Matrix<DB> &);
    template <class DB> friend DB DQTableToInterp(const Matrix<DB> &, const Matrix<DB> &,DB);
    template <class DB> friend DB DDTableToInterp(const Matrix<DB> &, const Matrix<DB> &,DB);
    template <class DB> friend DB _2P3DHermiteInterp(const Matrix<DB> &,int i,int j,DB);
    template <class DB> friend DB _3P3DHermiteInterp(const Matrix<DB> &,DB);
    template <class DE> friend Matrix<DE> SplineSlope(const Matrix<DE> &);
    template <class DE> friend DE Interpolation(const Matrix<DE> &,DE,string);
    template <class DE> friend DE Interpolation(const Matrix<DE> &,DE);
    template <class DE> friend Matrix<DE> Interpolation(const Matrix<DE> &A, const Matrix<DE> &x, string);
    template <class DE> friend std::function<DE(DE)> Interpolation(const Matrix<DE> &,string);
    //拟合与逼近
    template <class DB> friend Matrix<DB> PseudoInverse(const Matrix<DB>&);
    template <class DB> friend Matrix<DB> LeastSquares(const Matrix<DB>&, const Matrix<DB>&);
    template <class DB> friend Polynomial<DB> PolyFit(const Matrix<DB> &, const vector<Polynomial<DB>> &);
    template <class DB> friend Polynomial<DB> PolyFit(const Matrix<DB>&,int);
    template <class DB> friend Polynomial<DB> PolyFit(std::function<DB(DB)>,const vector<Polynomial<DB>>&,DB,DB);
private:
    vector<vector<DM>> value;
    int row_num;
    int column_num;
};
template <class DB> const Matrix<DB> EmptyMatrix = Matrix<DB>(0, 0);

//积分
template <class DI> DI Integrate(std::function<DI(DI)>,DI,DI,string);
template <class DI> DI Integrate(std::function<DI(DI)>,DI,DI);

//微分
template <class DD> DD D(std::function<DD(DD)>, DD);

//一元函数求根
template <class DD> DD Iteration(std::function<DD(DD)>, DD);
template<class DD> DD FindRoot(std::function<DD(DD)>, DD,DD);
template <class DD> DD FindRoot(std::function<DD(DD)>, DD,string);
template <class DD> DD FindRoot(std::function<DD(DD)>, DD);
template <class DD> DD FindRoot(std::function<DD(DD)>, DD,string,int);
template <class DD> DD Sqrt(DD);
template <> double Sqrt(double);

//杂例II
BigInt Factorial(int);
BigInt Fibonacci_r(int);
BigInt Fibonacci_p(int);
BigInt Fibonacci_m(int);
BigInt Fibonacci(int);


//-----定义部分-----


//快速傅里叶变换
inline int Rev(int k,int m)
{
    int s = 0;
    while(m>0)
    {
        --m;
        s += ((k%2) << m);
        k = (k >> 1);
    }
    return s;
}
template<class DC> vector<DC> bit_reverse_copy(const vector<DC> &a)
{
    int n = a.size();
    int m = log2(n);
    vector<DC> A(n);
    for(int k=0;k<n;++k)
        A[Rev(k,m)] = a[k];
    return A;
}
template<class DC> vector<complex<DC>> iterative_FFT(const vector<complex<DC>>&a,bool _jud)
{
    vector<complex<DC>> A = bit_reverse_copy(a);
    int n = a.size();
    int m= 1;
    complex<DC> omega_m;
    complex<DC> omega;
    for (int s = 1; s <= log2(n);++s)
    {
        m= m<< 1;
        if(_jud)
            omega_m=complex<DC>(cos(2*Pi<DC>/m),sin(2*Pi<DC>/m));
        else
            omega_m=complex<DC>(cos(2*Pi<DC>/m),0-sin(2*Pi<DC>/m));
        for (int k = 0; k < n;k=k+m)
        {
            omega =complex<DC>(1,0);
            int halfm = m >> 1;
            for (int j = 0; j < halfm;++j)
            {
                complex<DC> t=omega*A[k+j+halfm];
                complex<DC> u = A[k + j];
                A[k+j]=u+t;
                A[k + j + halfm] = u - t;
                omega = omega * omega_m;
            }
        }
    }
    return A;
}
template<class DC> vector<complex<DC>> recursive_FFT(const vector<complex<DC>>& a,bool _jud)
{
    int n = a.size();
    if(n==1)
       return a;
    int halfn = (n >> 1);
    complex<DC> omega_n;
    if(_jud)
        omega_n=complex<DC>(cos(2 * Pi<DC>/n),sin(2*Pi<DC>/n));
    else
        omega_n=complex<DC>(cos(2 * Pi<DC>/n),sin(0-2*Pi<DC>/n));
    complex<DC> omega{1, 0};
    vector<complex<DC>> a_even(halfn);
    vector<complex<DC>> a_odd(halfn);
    for (int i = 0;i<halfn;++i)
    {
        int twicei = i << 1;
        a_even[i] = a[twicei];
        a_odd[i] = a[twicei+1];
    }
    vector<complex<DC>> y_even= recursive_FFT(a_even,_jud);
    vector<complex<DC>> y_odd = recursive_FFT(a_odd,_jud);
    vector<complex<DC>> y(n);
    for (int k = 0; k < halfn;++k)
    {
        complex<DC> temp = omega * y_odd[k];
        y[k] = y_even[k] + temp;
        y[k + halfn] = y_even[k] - temp;
        omega = omega * omega_n;
    }
    return y;
}
template<class DC> vector<complex<DC>> DFT(const vector<complex<DC>> &a)
{
    return iterative_FFT(a, true);
}
template<class DC> vector<DC> Product(const vector<DC>& a,const vector<DC> &b)
{
    int n = a.size()<b.size()?a.size():b.size();
    vector<DC> c(n);
    for (int i = 0; i < n;++i)
        c[i] = a[i] * b[i];
    return c;
}
template<class DC> vector<complex<DC>> IDFT(const vector<complex<DC>>& a)
{
    unsigned n = a.size();
    vector<complex<DC>> r = iterative_FFT(a, false);
    for (unsigned i=0;i<n;++i)
        r[i] = r[i] /DC(n);
    return r;
}
template<class DC> vector<complex<DC>> FFT(const vector<complex<DC>> &x,const vector<complex<DC>> &y)
{
    int n = x.size();
    vector<complex<DC>> a=x;
    vector<complex<DC>> b=y;
    a.resize(n << 1);
    b.resize(n << 1);
    auto r=IDFT(Product(DFT(a), DFT(b)));
    r.pop_back();
    return r;
}
template<class DC> vector<DC> RealFFT(const vector<DC> &x,const vector<DC> &y)
{
    int n = x.size();
    vector<complex<DC>> a(n);
    vector<complex<DC>> b(n);
    for (int i = 0; i < n;++i)
    {
        a[i] = complex<DC>(x[i], 0);
        b[i] = complex<DC>(y[i], 0);
    }
    auto c = FFT(a, b);
    int m = c.size();
    vector<DC> r(m);
    for (int i = 0;i<m;++i)
        r[i]=c[i].real();
    return r;
}
template<class DC> vector<vector<DC>> _AddZero(const vector<DC> &x,const vector<DC> &y)
{
    unsigned n = x.size() > y.size() ? x.size() : y.size();
    unsigned nearest = 1 << unsigned(log2(n));
    if(n!=nearest)
        nearest = nearest << 1;
    vector<DC> a=x;
    vector<DC> b=y;
    a.resize(nearest);
    b.resize(nearest);
    return {a, b};
}
template<class DC> vector<complex<DC>> Convolution(const vector<complex<DC>>&x,const vector<complex<DC>>& y)
{
    unsigned n = x.size() > y.size() ? x.size() : y.size();
    unsigned nearest = 1 << unsigned(log2(n));
    if(n!=nearest)
        nearest = nearest << 1;
    vector<complex<DC>> a=x;
    vector<complex<DC>> b=y;
    a.resize(nearest);
    b.resize(nearest);
    auto r=FFT(a, b);
    r.resize(x.size() + y.size() - 1);
    return r;
}
template<class DC> vector<DC> RealConvolution(const vector<DC>&x,const vector<DC>& y)
{
    vector<vector<DC>> ab = _AddZero(x, y);
    auto r=RealFFT(ab[0], ab[1]);
    r.resize(x.size() + y.size() - 1);
    return r;
}
template<class DC> DC ModPower(DC a,DC n,DC mod)
{
    DC m=a;
    DC b=1;
    while(n>=1)
    {
        if(n&1)
        {
            b =1ll*m*b%mod;
        }
        n=n>>1;
        m=1ll*m*m%mod;
    }
    return b;
}
template<class DC> vector<DC> ModProduct(const vector<DC>& a,const vector<DC> &b,DC mod)
{
    int n = a.size()<b.size()?a.size():b.size();
    vector<DC> c(n);
    for (int i = 0; i < n;++i)
        c[i] = 1ll*a[i]*b[i]%mod;
    return c;
}
template<class DC> vector<DC> iterative_NTT(const vector<DC>&a,bool _jud)
{
    const DC g = 3;
    const DC gi = 332748118;
    const DC mod = 998244353;
    vector<DC> A= bit_reverse_copy(a);
    unsigned n = a.size();
    DC omega, omega_m;
    unsigned m, halfm;
    for (unsigned s = 1; s <= unsigned(log2(n));++s)
    {
        m= 1<<s;
        if(_jud)
            omega_m=ModPower(g,(mod-1)/DC(m),mod);
        else
            omega_m=ModPower(gi,(mod-1)/DC(m),mod);
        for (unsigned k = 0; k < n;k=k+m)
        {
            omega =1;
            halfm = m >> 1;
            for (unsigned j = 0; j < halfm;++j)
            {
                DC t=1ll*omega*A[k+j+halfm]%mod;
                DC u = A[k + j];
                A[k+j]=(u+t)%mod;
                A[k + j + halfm] =(mod+u - t)%mod;
                omega =1ll*omega * omega_m%mod;
            }
        }
    }
    if(!_jud)
    {
        DC n_inverse=ModPower(DC(n),mod-2,mod);
        for (unsigned i = 0; i < n;++i)
            A[i]=1ll*A[i]*n_inverse%mod;//这里要乘长度的逆元
    }
    return A;
}
template<class DC> vector<DC> NTT(const vector<DC> &x, const vector<DC> &y)
{
    unsigned n=x.size();
    vector<DC> a = x;
    vector<DC> b = y;
    a.resize(n << 1);
    b.resize(n << 1);
    const DC mod=998244353;
    auto r=iterative_NTT(ModProduct(iterative_NTT(a, true), iterative_NTT(b, true),mod), false);
    r.pop_back();
    return r;
}
template<class DC> vector<DC> IntConvolution(const vector<DC> &x,const vector<DC> &y)
{
    vector<vector<DC>> ab = _AddZero(x, y);
    auto r=NTT(ab[0], ab[1]);
    r.resize(x.size()+y.size()-1);
    return r;
}
template<class DC> vector<DC> CarryBit(const vector<DC>& a,DC bit)
{
    vector<DC> r=a;
    unsigned n = r.size();
    for (unsigned i = 0; i <n-1;++i)
    {
        if(r[i]>=bit)
        {
            r[i + 1] += r[i] / bit;
            r[i] = r[i] % bit;
        }
    }
    while(r[r.size()-1]>=bit)
    {
        r.emplace_back(r[r.size()-1]/bit);
        r[r.size() - 2] = r[r.size()-2] % bit;
    }
    return r;
}

//杂例I
template<class DM> double Norm(const vector<DM> &a, double p)
{
    if(p==INFINITY)
    {
        double max = 0;
        for (long long unsigned int i = 0; i < a.size();++i)
        {
            if(Abs(a[i])>max)
                max = Abs(a[i]);
        }
        return max;
    }
    double s = 0;
    for (long long unsigned int i = 0; i < a.size();++i)
    {
        s += pow(Abs(a[i]), p);
    }
    s = pow(s, 1 / p);
    return s;
}
template<class DJ> DJ Abs(DJ a)
{
    if(a>0)
        return a;
    return -a;
}
template <class DJ> vector<DJ> Range(DJ a, DJ b,DJ c)
{
    vector<DJ> r;
    DJ s = a;
    while(s<=b)
    {
        r.emplace_back(s);
        s = s + DJ(c);
    }
    return r;
}
template <class DJ> vector<DJ> Range(DJ a, DJ b)
{
    return Range(a, b, DJ(1));
}
template <class DJ> vector<DJ> Rand(int a)
{
    vector<DJ> r(a);
    for (int i = 0;i<a;++i)
        r[i]=std::rand()*Sig();
    return r;
}
template <class DJ> DJ Limit(std::function<DJ(DJ)> f, DJ x0)
{
    const DJ h = 1e-7;
    return 0.5 * (f(x0+h)+f(x0-h));
}
template <class DJ> DJ N(DJ a,int n)
{
    cout << std::setprecision(n) << a << '\n';
    return a;
}

//多项式类
const Polynomial<double> X({0, 1}, 1);
template<class DC> Polynomial<DC>::Polynomial()
{
    poly_deg = 0;
    poly_coeff.push_back(0);
}
template<class DC> Polynomial<DC>::Polynomial(const vector<DC> &coeff)
{
    poly_deg = coeff.size() - 1;
    poly_coeff = coeff;
}
template<class DC> Polynomial<DC>::Polynomial(const vector<DC> &s,int deg)
{
    poly_deg = deg;
    poly_coeff = s;
}
template<class DC> Polynomial<DC>::Polynomial(DC r)
{
    poly_deg = 0;
    poly_coeff.push_back(r);
}
template<class DC> Polynomial<DC> &Polynomial<DC>::operator=(const Polynomial<DC> &f)
{
    poly_deg = f.poly_deg;
    poly_coeff=f.poly_coeff;
    return *this;
}
template<class DC> DC Polynomial<DC>::operator()(DC t)
{
    DC s=poly_coeff[poly_deg];
    for (int i = poly_deg-1; i >= 0;--i)
        s = s * t + poly_coeff[i];
    return s;
}
template<class DC> Polynomial<DC> operator+(const Polynomial<DC> &a,const Polynomial<DC> &b)
{
    bool here = true;
    int degmax = a.poly_deg;
    if(b.poly_deg>a.poly_deg)
    {
        degmax=b.poly_deg;
        here = false;
    }
    vector<DC> r(degmax + 1);
    Polynomial<DC> c;
    c.poly_deg = degmax;
    c.poly_coeff = r;
    if(here)
    {
        for (int i = 0; i <=b.poly_deg;++i)
        {
            c.poly_coeff[i] = a.poly_coeff[i] + b.poly_coeff[i];
        }
        for (int i = b.poly_deg + 1;i<=c.poly_deg;++i)
        {
            c.poly_coeff[i] = a.poly_coeff[i];
        }
        return c;
    }
    for (int i = 0; i <=a.poly_deg;++i)
    {
        c.poly_coeff[i] = a.poly_coeff[i] + b.poly_coeff[i];
    }
    for (int i = a.poly_deg + 1;i<=c.poly_deg;++i)
    {
        c.poly_coeff[i] = b.poly_coeff[i];
    }
    return c;
}
template<class DC> Polynomial<DC> operator-(const Polynomial<DC> &a,const Polynomial<DC> &b)
{
    bool here = true;
    int degmax = a.poly_deg;
    if(b.poly_deg>a.poly_deg)
    {
        degmax=b.poly_deg;
        here = false;
    }
    vector<DC> r(degmax + 1);
    Polynomial<DC> c;
    c.poly_deg = degmax;
    c.poly_coeff = r;
    if(here)
    {
        for (int i = 0; i <=b.poly_deg;++i)
        {
            c.poly_coeff[i] = a.poly_coeff[i] - b.poly_coeff[i];
        }
        for (int i = b.poly_deg + 1;i<=c.poly_deg;++i)
        {
            c.poly_coeff[i] = a.poly_coeff[i];
        }
        return c;
    }
    for (int i = 0; i <=a.poly_deg;++i)
    {
        c.poly_coeff[i] = a.poly_coeff[i] - b.poly_coeff[i];
    }
    for (int i = a.poly_deg + 1;i<=c.poly_deg;++i)
    {
        c.poly_coeff[i] = DC(0)-b.poly_coeff[i];
    }
    return c;
}
template <class DD> Polynomial<DD> operator+(DD a,const Polynomial<DD> &f)
{
    Polynomial<DD> p = f;
    p.poly_coeff[0] += a;
    return p;
}
template <class DD> Polynomial<DD> operator*(DD a,const Polynomial<DD> &f)
{
    vector<DD> r;
    for (int i = 0; i <= f.poly_deg;++i)
    {
        r.push_back(a * f.poly_coeff[i]);
    }
    Polynomial<DD> p;
    p.poly_deg = f.poly_deg;
    p.poly_coeff = r;
    return p;
}
template <class DD> Polynomial<DD> operator*(const Polynomial<DD> &a, const Polynomial<DD> &b)
{
    int n = a.poly_deg > b.poly_deg ? a.poly_deg : b.poly_deg;
    if (a.poly_deg == 0)
        return a.poly_coeff[0] * b;
    if(b.poly_deg==0)
        return b.poly_coeff[0] * a;
    if(n==1)
    {
        vector<DD> r(3);
        r[0]=a.poly_coeff[0]*b.poly_coeff[0];
        r[2]=a.poly_coeff[1]*b.poly_coeff[1];
        r[1] = (a.poly_coeff[0] + a.poly_coeff[1]) * (b.poly_coeff[0] + b.poly_coeff[1]) - r[0] - r[2];
        Polynomial<DD> f;
        f.poly_deg=2;
        f.poly_coeff = r;
        return f;
    }
    Polynomial<DD> p=a;
    Polynomial<DD> q=b;
    if(a.poly_deg>b.poly_deg)
    {
        q.poly_coeff.insert(q.poly_coeff.end(), a.poly_deg - b.poly_deg, 0);
        q.poly_deg = n;
    }
    if(b.poly_deg>a.poly_deg)
    {
        p.poly_coeff.insert(p.poly_coeff.end(), b.poly_deg - a.poly_deg, 0);
        p.poly_deg = n;
    }
    Polynomial<DD> p0 = SubPoly(p, 0, n / 2);
    Polynomial<DD> p1 =SubPoly(p, n/2+1, n);
    Polynomial<DD> q0=SubPoly(q,0,n/2);
    Polynomial<DD> q1 = SubPoly(q, n/2+1, n);
    Polynomial<DD> r0 = p0 * q0;
    Polynomial<DD> r1=p1*q1;
    Polynomial<DD> r2=(p0+p1)*(q0+q1);
    return r0 + ((r2 - r0 - r1) >> (n/2+1)) + (r1 >> ((n/2+1)<<1));
}
template <class DD> Polynomial<DD> operator/(const Polynomial<DD>& f,const Polynomial<DD>& g)
{
    cerr<<"尚待开发"<<'\n';
    return f;
}
template<class DC> Polynomial<DC> operator^(const Polynomial<DC> &a,int n)
{
    if(n<0)
    {
        cerr << "错误：多项式不能有负的幂次。" << '\n';
        return a;
    }
    Polynomial<DC> m=a;
    Polynomial<DC> b({1},0);
    while(n>=1)
    {
        if(n&1)
        {
            b = m * b;
        }
        n=n>>1;
        m = m * m;
    }
    return b;
}
template<class DC> ostream & operator<<(ostream & os, const Polynomial<DC> &f)
{
    os << f.poly_coeff[0];
    for (int i =1;i<=f.poly_deg;++i)
    {
        if(f.poly_coeff[i]==0)
            continue;
        if(f.poly_coeff[i]>0)
        {
            os << '+';
        }
        os << f.poly_coeff[i] << "x" ;
        if(i!=1)
        {
            os << '^'<<i;
        }
    }
    return os;
}
template<class DC> Polynomial<DC> &Polynomial<DC>::operator+=(const Polynomial<DC> &b)
{
    *this = (*this) + b;
    return *this;
}
template<class DC> Polynomial<DC> &Polynomial<DC>::operator-=(const Polynomial<DC> &b)
{
    *this = (*this)- b;
    return *this;
}
template <class DD> Polynomial<DD> operator>>(const Polynomial<DD> &f, int a)
{
    vector<DD> r(f.poly_deg+1+a);
    for (int i = 0; i <= f.poly_deg;++i)
    {
        r[i + a] = f.poly_coeff[i];
    }
    Polynomial<DD> p;
    p.poly_deg = f.poly_deg + a;
    p.poly_coeff = r;
    return p;
}
template <class DD> DD Get(const Polynomial<DD> &f, DD t)
{
    DD s=f.poly_coeff[f.poly_deg];
    for (int i = f.poly_deg-1; i >= 0;--i)
        s = s * t + f.poly_coeff[i];
    return s;
}
template <class DD> vector<DD> GetCoef(const Polynomial<DD> &f)
{
    return f.poly_coeff;
}
template<class DC> void Disp(const Polynomial<DC> &f)
{
    cout << f<<'\n'<<'\n';
}
template<class DC> int Deg(const Polynomial<DC> &f)
{
    return f.poly_deg;
}
template <class DC> Polynomial<DC> SubPoly(const Polynomial<DC> &f, int a, int b)
{
    if(a<0)
        a = 0;
    if(b>f.poly_deg)
        b = f.poly_deg;
    vector<DC> r;
    for (int i = a; i <= b;++i)
    {
        r.push_back(f.poly_coeff[i]);
    }
    Polynomial<DC> p;
    p.poly_deg=b-a;
    p.poly_coeff = r;
    return p;
}
template <class DC> Polynomial<DC> D(const Polynomial<DC> &f)
{
    if(f.poly_deg==0)
    {
        Polynomial<DC> d({0},0);
        return d;
    }
    int n=f.poly_deg;
    vector<DC> r;
    for (int i = 1; i <= n;++i)
    {
        r.push_back(i * f.poly_coeff[i]);
    }
    Polynomial<DC> d;
    d.poly_deg = n - 1;
    d.poly_coeff = r;
    return d;
}
template <class DD> Polynomial<DD> Integrate(const Polynomial<DD> &f)
{
    if(f.poly_deg==0 && f.poly_coeff[0]==0)
    {
        return Polynomial<DD>({0}, 1);
    }
    int deg = f.poly_deg;
    vector<DD> r(deg+2);
    for (int i = 1; i <= deg + 1;++i)
    {
        r[i] = f.poly_coeff[i - 1] / i;
    }
    return Polynomial<DD>(r, deg + 1);
}
template <class DD> DD Integrate(const Polynomial<DD> &f,DD a,DD b)
{
    Polynomial<DD> g=Integrate(f);
    return g(b) - g(a);
}
template <class DD> Polynomial<DD> Schmidt(int n,DD a,DD b)
{
    if(n<0)
    {
        return Polynomial<DD>({0}, 0);
    }
    vector<Polynomial<DD>> phi(n+1);
    vector<DD> innerproduct(n);
    phi[0] = Polynomial<DD>({1}, 0);
    for (int k = 1; k <= n;++k)
    {
        vector<DD> r(k+1);
        r[k] = DD(1);
        Polynomial<DD> f(r, k);
        phi[k] = f;
        innerproduct[k - 1] = Integrate(phi[k - 1] * phi[k - 1], a, b);
        for (int j = 0; j < k;++j)
        {
            phi[k] = phi[k] - (Integrate(phi[j] * f, a, b) / innerproduct[j])*phi[j];
        }
    }
    return phi[n];
}
template<class DD> Polynomial<DD> Legendre(int n)
{
    if(n<0)
        return Polynomial<DD>({0}, 0);
    vector<Polynomial<DD>> p({
        Polynomial<DD>({1}, 0),
        Polynomial<DD>({0, 1}, 1),
        Polynomial<DD>({-0.5, 0, 1.5}, 2),
        Polynomial<DD>({0, -1.5, 0, 2.5}, 3),
        Polynomial<DD>({0.375, 0, -3.75, 0, 4.375}, 4)
    });
    if(n<=4)
        return p[n];
    Polynomial<DD> pre1 = p[4];
    Polynomial<DD> pre2 = p[3];
    Polynomial<DD> now;
    for (int k = 5; k <= n;++k)
    {
        now = (DD(2 *k- 1) /k) *(pre1>>1) - (DD(k - 1) /k) * pre2;
        pre2 = pre1;
        pre1 = now;
    }
    return now;
}

//矩阵类
template<class DB> Matrix<DB>::Matrix()
{
    row_num=1;
    column_num = 1;
    value = {{0}};
}
template<class DB> Matrix<DB>::Matrix(DB a)
{
    row_num = 1;
    column_num = 1;
    value ={{a}};
}
template<class DB> Matrix<DB>::Matrix(int r,int c)
{
    row_num=r;
    column_num=c;
    vector<DB> rr(c);
    for(int i=0;i<r;++i)
    {
        value.push_back(rr);
    }
}
template<class DB> Matrix<DB>::Matrix (const vector<vector<DB>> &a,int r,int c)
{
    row_num=r;
    column_num = c;
    value = a;
}
template<class DB> Matrix<DB>::Matrix (const vector<vector<DB>> &a)
{
    row_num = a.size();
    column_num = a[0].size();
    value = a;
}
template<class DB> Matrix<DB>::Matrix (const vector<DB> &a)
{
    row_num = 1;
    column_num=a.size();
    value = {a};
}
template<class DB> Matrix<DB>& Matrix<DB>::operator=(const Matrix<DB> &b)
{
    row_num = b.row_num;
    column_num = b.column_num;
    value = b.value;
    return *this;
}
template<class DB> DB& Matrix<DB>::operator()(int r,int c)
{
    if(r>=row_num || c>=column_num)
    {
        cerr << "错误：矩阵下标越界" << '\n';
        return value[0][0];
    }
    return value[r][c];
}
template<class DB> const vector<DB>& Matrix<DB>::operator[](int r)
{
    if(r>=row_num)
    {
        cerr << "错误：矩阵下标越界" << '\n';
        return value[0];
    }
    return value[r];
}
template <class DB> bool operator==(const Matrix<DB> &A, const Matrix<DB> &B)
{
    int r=A.row_num;
    int c=A.column_num;
    if(r!=B.row_num ||c!=B.column_num)
        return false;
    for (int i = 0; i < r;++i)
        if(A.value[i]!=B.value[i])
            return false;
    return true;
}
template<class DB> ostream &operator<< (ostream &os,const Matrix<DB> &a)
{
    for (int i = 0;i<a.row_num;++i)
    {
        os << '[' << a.value[i][0];
        for (int j = 1;j<a.column_num;++j)
        {
            os <<", "<<a.value[i][j];
        }
        os << "]" << '\n';
    }
    return os;
}
template<class DB> Matrix<DB> operator+(const Matrix<DB> &a,const Matrix<DB>&b)
{
    if(a.row_num!=b.row_num ||a.column_num!=b.column_num)
    {
        cerr << "错误：规格不同的矩阵不能相加。" << '\n';
        return a;
    }
    int r = a.row_num;
    int c = a.column_num;
    Matrix<DB> d(r,c);
    for (int i = 0; i < r;++i)
    {
        for (int j = 0; j < c;++j)
        {
            d.value[i][j] = a.value[i][j] + b.value[i][j];
        }
    }
    return d;
}
template<class DB> Matrix<DB> operator-(const Matrix<DB> &a,const Matrix<DB>&b)
{
    if(a.row_num!=b.row_num ||a.column_num!=b.column_num)
    {
        cerr << "错误：规格不同的矩阵不能相减。" << '\n';
        return a;
    }
    int r = a.row_num;
    int c = a.column_num;
    Matrix<DB> d(r,c);
    for (int i = 0; i < r;++i)
    {
        for (int j = 0; j < c;++j)
        {
            d.value[i][j] = a.value[i][j] - b.value[i][j];
        }
    }
    return d;
}
template <class DB> Matrix<DB> operator-(const Matrix<DB> &b)
{
    int r = b.row_num;
    int c = b.column_num;
    Matrix<DB> d(r,c);
    for (int i = 0; i < r;++i)
    {
        for (int j = 0; j < c;++j)
        {
            d.value[i][j] = DB(0)-b.value[i][j];
        }
    }
    return d;
}
template<class DB> Matrix<DB> operator*(const Matrix<DB> &a,const Matrix<DB>&b)
{
    if(a.column_num!=b.row_num)
    {
        cerr << "错误：矩阵相乘时，前一个矩阵的列数必须等于后一个矩阵的行数。" << '\n';
        return a;
    }
    int p=a.row_num;
    int q=a.column_num;
    int r = b.column_num;
    Matrix<DB> d(p, r);
    int kernel =4;
    int pend = p / kernel * kernel;
    int rend = r / kernel * kernel;
    for (int i = 0; i < pend; i+=kernel) 
        for (int k = 0; k < q; ++k) 
        {
            for (int j= 0; j< rend; j+=kernel) 
                for (int iker = 0; iker < kernel;++iker)
                    for (int jker = 0; jker < kernel; ++jker)
                        d.value[i + iker][j + jker] += a.value[i + iker][k] * b.value[k][j + jker]; 
            for (int iker = 0; iker < kernel;++iker)
                for (int jker= rend; jker <r; ++jker)
                    d.value[i + iker][jker] += a.value[i + iker][k] * b.value[k][jker]; 
        }                                  
    for (int k = 0; k < q; ++k) 
    {
        for (int j= 0; j< rend; j+=kernel) 
            for (int iker =pend; iker < p;++iker)
                for (int jker = 0; jker < kernel; ++jker)
                    d.value[iker][j + jker] += a.value[iker][k] * b.value[k][j + jker]; 
        for (int iker = pend; iker < p;++iker)
            for (int jker= rend; jker <r; ++jker)
                d.value[iker][jker] += a.value[iker][k] * b.value[k][jker]; 
    }  
    return d;
}
template <class DC> Matrix<DC> operator+(DC r, const Matrix<DC> &a)
{
    Matrix<DC> d(a.row_num,a.column_num);
    for (int i = 0; i < a.row_num;++i)
    {
        for (int j = 0; j < a.column_num;++j)
        {
            d.value[i][j] = r+a.value[i][j];
        }
    }
    return d;
}
template <class DC> Matrix<DC> operator+(const Matrix<DC> &a,DC r)
{
    return r + a;
}
template <class DC> Matrix<DC> operator-(DC r,const Matrix<DC> &a)
{
    Matrix<DC> d(a.row_num,a.column_num);
    for (int i = 0; i < a.row_num;++i)
    {
        for (int j = 0; j < a.column_num;++j)
        {
            d.value[i][j] = r-a.value[i][j];
        }
    }
    return d;
}
template <class DC> Matrix<DC> operator-(const Matrix<DC> &a,DC r)
{
    Matrix<DC> d(a.row_num,a.column_num);
    for (int i = 0; i < a.row_num;++i)
    {
        for (int j = 0; j < a.column_num;++j)
        {
            d.value[i][j] = a.value[i][j]-r;
        }
    }
    return d;
}
template <class DC> Matrix<DC> operator*(DC r, const Matrix<DC> &a)
{
    Matrix<DC> d(a.row_num,a.column_num);
    for (int i = 0; i < a.row_num;++i)
    {
        for (int j = 0; j < a.column_num;++j)
        {
            d.value[i][j] = r * a.value[i][j];
        }
    }
    return d;
}
template <class DC> Matrix<DC> operator*(const Matrix<DC> &a,DC r)
{
    return r * a;
}
template<class DB> Matrix<DB> operator/(const Matrix<DB> &a,DB r)
{
    Matrix<DB> d(a.row_num,a.column_num);
    for (int i = 0; i < a.row_num;++i)
    {
        for (int j = 0; j < a.column_num;++j)
        {
            d.value[i][j] = a.value[i][j]/r;
        }
    }
    return d;
}
template<class DB> Matrix<DB> &Matrix<DB>::operator+=(const Matrix<DB> &b)
{
    *this = *this + b;
    return *this;
}
template<class DB> Matrix<DB> &Matrix<DB>::operator-=(const Matrix<DB> &b)
{
    *this = *this - b;
    return *this;
}
template <class DB> DB Get(const Matrix<DB> &A, int r, int c)
{
    return A.value[r][c];
}
template<class DB> void Clear(Matrix<DB> &a)
{
    a.value.clear();
}
template<class DB> void Disp(const Matrix<DB> &a)
{
    cout << a << '\n';
}
template <class DB> Matrix<DB> Generate(int r, int c)
{
    Matrix<DB> G(r,c);
    for(int i=0;i<r;++i)
        G.value[i] = Rand<DB>(c);
    return G;
}
template <class DB> Matrix<DB> Act(std::function<DB(DB)> f, const Matrix<DB> &A)
{
    Matrix<DB> r(A.row_num,A.column_num);
    for (int i = 0; i < A.row_num;++i)
        for (int j = 0; j < A.column_num;++j)
            r.value[i][j] = f(A.value[i][j]);
    return r;
}
template <class DB> 
Matrix<DB> Act(std::function<DB(DB, DB)> f, const Matrix<DB> &A, const Matrix<DB> &B)
{
    int r=A.row_num;
    int c=A.column_num;
    if(r!=B.row_num || c!=B.column_num)
    {
        cerr<<"错误：两个矩阵规格不同，无法作用"<<'\n';
        return A;
    }
    Matrix<DB> R(r,c);
    for (int i = 0; i < r;++i)
        for (int j = 0; j <c;++j)
            R.value[i][j] = f(A.value[i][j],B.value[i][j]);
    return R;
}
template <class DB> Matrix<DB> Table(std::function<DB(DB)> f,const vector<DB>& r)
{
    vector<DB> vt;
    for(auto i:r)
        vt.push_back(f(i));
    Matrix<DB> M({vt});
    return M;
}
template<class DB> Matrix<DB> Transpose(const Matrix<DB> &a)
{
    Matrix<DB> T(a.column_num, a.row_num);
    for (int i = 0; i < a.row_num;++i)
    {
        for (int j =0; j<a.column_num;++j)
            T.value[j][i] = a.value[i][j];
    }
    return T;
}
template<class DB> Matrix<DB> Diagonal(const Matrix<DB> &a)
{
    if(a.row_num!=a.column_num)
    {
        cerr << "错误：非方阵不能取其对角元素" << '\n';
        return a;
    }
    int n = a.row_num;
    Matrix<DB> d(n,1);
    for (int i = 0; i < n;++i)
        d.value[i][0] = a.value[i][i];
    return d;
}
template<class DB> Matrix<DB> DiagonalMatrix(const vector<DB> &d)
{
    int n = d.size();
    Matrix<DB> A(n,n);
    for(int i=0;i<n;++i)
        A.value[i][i]=d[i];
    return A;
}
template<class DB> int RowSize(const Matrix<DB> &a)
{
    return a.row_num;
}
template<class DB> int ColumnSize(const Matrix<DB> &a)
{
    return a.column_num;
}
template <class DB> void Resize(Matrix<DB> &A,int r,int c)
{
    if(r<A.row_num)
        A.value.resize(r);
    if(r>A.row_num)
        for (int i = A.row_num; i < r;++i)
            A.value.emplace_back(vector<DB>(c));
    if(c==A.column_num)
        return;
    for (auto item : A.value)
        item.resize(c);
    A.row_num=r;
    A.column_num=c;
}
template <class DB> Matrix<DB> SubRow(const Matrix<DB> &A, int a, int b)
{
    Matrix<DB> r(b-a+1,A.column_num);
    for (int i = a; i <= b;++i)
        r.value[i-a]=A.value[i];
    return r;
}
template <class DB> Matrix<DB> SubMatrix(const Matrix<DB> &A, int lur, int luc, int rdr, int rdc)
{
    //左上点的坐标：(lur,luc)，右下点的坐标：(rdr,rdc)
    Matrix<DB> r(rdr - lur + 1, rdc - luc + 1);
    for (int i = lur; i <= rdr;++i)
    {
        for (int j = luc; j <= rdc;++j)
        {
            r.value[i - lur][j - luc] = A.value[i][j];
        }
    }
    return r;
}
template <class DB> void DeleteRow(Matrix<DB> &A,int n)
{
    if(n<0||n>=A.row_num)
    {
        cerr<<"错误：矩阵行数越界"<<'\n';
        return ;
    }
    --A.row_num;
    A.value.erase(A.value.begin() + n);
}
template <class DB> void ReplaceRow(Matrix<DB> &A,int n,const Matrix<DB> &B)
{
    if(n<0||n>=A.row_num)
    {
        cerr<<"错误：矩阵行数越界"<<'\n';
        return ;
    }
    if(B.row_num==0)
    {
        DeleteRow(A, n);
        return;
    }
    A.value[n] = B.value[0];
}
template<class DB> void _MatProd(const Matrix<DB> &A,const Matrix<DB> &B,Matrix<DB> &C)
{
    C = A * B;
}
template<class DB> Matrix<DB> ThreadMatProd(const Matrix<DB> &A,const Matrix<DB> &B)
{
    int p=RowSize(A);
    int q=RowSize(B);
    int r=ColumnSize(B);
    int halfp=p/2;
    int halfq=q/2;
    int halfr=r/2;
    Matrix<DB> C(p, r);
    const Matrix<DB> &ALU=SubMatrix(A,0,0,halfp-1,halfq-1);
    const Matrix<DB> &ARU = SubMatrix(A, 0,halfq, halfp- 1, q - 1);
    const Matrix<DB> &ALD=SubMatrix(A,halfp,0,p-1,halfq-1);
    const Matrix<DB> &ARD = SubMatrix(A, halfp, halfq, p - 1, q- 1);
    const Matrix<DB> &BLU=SubMatrix(B,0,0,halfq-1,halfr-1);
    const Matrix<DB> &BRU = SubMatrix(B, 0,halfr, halfq - 1, r- 1);
    const Matrix<DB> &BLD=SubMatrix(B,halfq,0,q-1,halfr-1);
    const Matrix<DB> &BRD = SubMatrix(B,halfq, halfr, q- 1,r- 1);
    Matrix<DB> LULU, RULD, LURU, RURD, LDLU, RDLD, LDRU, RDRD;
    thread t1(_MatProd<DB>,ref(ALU),ref(BLU),ref(LULU));
    thread t2(_MatProd<DB>,ref(ARU),ref(BLD),ref(RULD));
    thread t3(_MatProd<DB>,ref(ALU),ref(BRU),ref(LURU));
    thread t4(_MatProd<DB>,ref(ARU),ref(BRD),ref(RURD));
    thread t5(_MatProd<DB>,ref(ALD),ref(BLU),ref(LDLU));
    thread t6(_MatProd<DB>,ref(ARD),ref(BLD),ref(RDLD));
    thread t7(_MatProd<DB>,ref(ALD),ref(BRU),ref(LDRU));
    thread t8(_MatProd<DB>,ref(ARD),ref(BRD),ref(RDRD));
    t1.join();t2.join();t3.join();t4.join();
    t5.join();t6.join();t7.join();t8.join();
    return ColumnCat(RowCat(LULU+RULD,LURU+RURD),RowCat(LDLU+RDLD,LDRU+RDRD));
}

//解线性方程组
template<class DB> Matrix<DB> RowCat(const Matrix<DB> &a,const Matrix<DB> &b)
{
    if(a.row_num!=b.row_num)
    {
        cerr << "错误：行数不相等的两个矩阵无法行连接。" << '\n';
        return a;
    }
    Matrix<DB> c(a.row_num, a.column_num+b.column_num);
    for (int i = 0; i < a.row_num;++i)
    {
        for (int j = 0;j<a.column_num;++j)
        {
            c.value[i][j] = a.value[i][j];
        }
        for (int j = a.column_num; j < a.column_num + b.column_num;++j)
        {
            c.value[i][j] = b.value[i][j - a.column_num];
        }
    }
    return c;
}
template <class DB> Matrix<DB> RowCat(const vector<Matrix<DB>> &a)
{
    Matrix<DB> r=a[0];
    for (long long unsigned int i = 1;i<a.size();++i)
    {
        r=RowCat(r, a[i]);
    }
    return r;
}
template <class DB> Matrix<DB> ColumnCat(const Matrix<DB> &a, const Matrix<DB> &b)
{
    if(a.column_num!=b.column_num)
    {
        cerr << "错误：列数不相等的两个矩阵无法列连接。" << '\n';
        return a;
    }
    vector<vector<DB>> r = a.value;
    for (int i = 0; i < b.row_num;++i)
        r.emplace_back(b.value[i]);
    return Matrix<DB>(r);
}
template<class DB> void SwapRow(Matrix<DB> &m,int a,int b)
{
    if(a==b)
    {
        return;
    }
    vector<DB> temp;
    temp=m.value[a];
    m.value[a]=m.value[b];
    m.value[b] = temp;
}
template <class DB> Matrix<DB> USplit(const Matrix<DB> &A)
{
    if(A.row_num!=A.column_num)
    {
        cerr << "错误：只有方阵才能进行上三角部分分割。" << '\n';
        return A;
    }
    int n=A.row_num;
    Matrix<DB> U(n,n);
    for (int i = 0; i < n;++i)
    {
        for (int j = i + 1; j < n;++j)
        {
            U.value[i][j] = -A.value[i][j];
        }
    }
    return U;
}
template <class DB> vector<Matrix<DB>> GaussianElimination(const Matrix<DB> &A, const Matrix<DB> &b)
{
    Matrix<DB> GA = A;
    Matrix<DB> Gb = b;
    DB epsilon = 1e-14;
    int n=GA.row_num;
    int m = Gb.column_num;
    DB max = 0;
    int maxposition = 0;
    DB pivot=1;//主元
    DB multiplicator=1;//乘数
    for(int i=0;i<n;++i)
    {
        max = Abs(GA.value[i][i]);
        maxposition = i;
        for (int j = i+1; j < n;++j)
        {
            if(Abs(GA.value[j][i])>max)
            {
                max = Abs(GA.value[j][i]);
                maxposition = j;
            }
        }
        if(max<=epsilon)
        {
            return {GA,Gb};
        }
        SwapRow(GA, i, maxposition);
        SwapRow(Gb, i, maxposition);
        pivot = GA.value[i][i];
        for(int k=i+1;k<n;++k)
        {
            multiplicator = GA.value[k][i]/pivot;
            for (int j = i; j < n;++j)
            {
                GA.value[k][j] = GA.value[k][j] - multiplicator * GA.value[i][j];
            }
            for (int j = 0; j < m;++j)
            {
                Gb.value[k][j] = Gb.value[k][j] - multiplicator * Gb.value[i][j];
            }
        }
    }
    return {GA,Gb};
}
template <class DB> Matrix<DB> LUCompactDecomposition(const Matrix<DB> &A)
{
    if(A.row_num!=A.column_num)
    {
        cerr << "错误：只有方阵才能进行LU分解。" << '\n';
        return A;
    }
    int n = A.row_num;
    Matrix<DB> LU = A;
    DB s = 0;
    for (int i = 0; i < n;++i)
    {
        if(LU.value[i][i]==0)
        {
            cerr << "错误：有一个顺序主子式等于0，不能LU分解。" << '\n';
            return A;
        }
        for (int j = i; j < n;++j)
        {
            s = 0;
            for (int k = 0; k < i;++k)
            {
                s += LU.value[i][k] * LU.value[k][j];
            }
            LU.value[i][j] = LU.value[i][j] - s;
        }
        for (int t = i + 1; t < n;++t)
        {
            s = 0;
            for (int k = 0; k < i;++k)
            {
                s += LU.value[t][k] * LU.value[k][i];
            }
            LU.value[t][i] = (LU.value[t][i] - s) / LU.value[i][i];
        }
    }
    return LU;
}
template <class DB> Matrix<DB> LUPCompactDecomposition(const Matrix<DB> &A, vector<int> &pi, int &rank, DB &det)
{
    Matrix<DB> LU;
    if(A.row_num>A.column_num)
        LU = Transpose(A);
    else
        LU = A;
    DB epsilon = 1e-14;
    int n=LU.row_num;
    int m = LU.column_num;
    DB max = 0;
    int maxposition = 0;
    DB pivot=1;//主元
    DB multiplicator=1;//乘数
    det = 1;
    for (int i = 0; i < n;++i)
        pi.push_back(i);
    for(int i=0;i<n;++i)
    {
        max = Abs(LU.value[i][i]);
        maxposition = i;
        for (int j = i+1; j < n;++j)
        {
            if(Abs(LU.value[j][i])>max)
            {
                max = Abs(LU.value[j][i]);
                maxposition = j;
            }
        }
        if(max<=epsilon)
        {
            rank = i;
            det = 0;
            return LU;
        }
        SwapRow(LU, i, maxposition);
        swap(pi[i], pi[maxposition]);
        pivot = LU.value[i][i];
        det =det*pivot;
        for(int k=i+1;k<n;++k)
        {
            multiplicator = LU.value[k][i]/pivot;
            for (int j = i; j < m;++j)
            {
                LU.value[k][j] = LU.value[k][j] - multiplicator * LU.value[i][j];
            }
            LU.value[k][i] = multiplicator;
        }
    }
    rank = n;
    return LU;
}
template <class DB> Matrix<DB> GaussianElimination(const Matrix<DB> &A)
{
    vector<int> pi;
    int rank = 0;
    DB det = 1;
    Matrix<DB> LU=LUPCompactDecomposition(A, pi, rank, det);
    int n=LU.row_num;
    for (int i = 0;i<n;++i)
    {
        for (int j = 0; j < i;++j)
        {
            LU.value[i][j] = 0;
        }
    }
    return LU;
}
template <class DB> int MatrixRank(const Matrix<DB> &A)
{
    vector<int> pi;
    int rank = 0;
    DB det = 1;
    LUPCompactDecomposition(A, pi, rank, det);
    return rank;
}
template <class DB> DB Det(const Matrix<DB> &A)
{
    vector<int> pi;
    int rank = 0;
    DB det = 1;
    LUPCompactDecomposition(A, pi, rank, det);
    return det;
}
template <class DB> vector<Matrix<DB>> LUDecomposition(const Matrix<DB> &A)
{
    Matrix<DB> LU=LUCompactDecomposition(A);
    int n = LU.row_num;
    Matrix<DB> L = LU;
    Matrix<DB> U = LU;
    for (int i = 0; i < n;++i)
    {
        L.value[i][i] = 1;
        for (int j = 0; j < i;++j)
        {
            U.value[i][j] = 0;
        }
        for (int j = i + 1; j < n; ++j)
        {
            L.value[i][j] = 0;
        }
    }
    return {L, U};
}
template <class DB> vector<Matrix<DB>> LUPDecomposition(const Matrix<DB> &A)
{
    int _rank;
    DB det;
    vector<int> pi;
    const Matrix<DB> &LUP=LUPCompactDecomposition(A, pi, _rank, det);
    int n=A.column_num;
    Matrix<DB> L(n,n);
    Matrix<DB> U(n,n);
    Matrix<DB> P(n,n);
    for (int i = 0; i < n;++i)
    {
        P.value[i][pi[i]] = 1;
        L.value[i][i]=1;
        for(int j=0;j<i;++j)
            L.value[i][j] = LUP.value[i][j];
        for (int j = i;j<n;++j)
            U.value[i][j]=LUP.value[i][j];
    }
    return {L, U, P};
}
template <class DB> vector<vector<DB>> Tridiagonal(const Matrix<DB> &A)
{
    int n = A.row_num;
    vector<DB> a, b, c;
    a.push_back(0);
    for (int i = 0; i < n-1;++i)
    {
        b.push_back(A.value[i][i]);
        c.push_back(A.value[i][i + 1]);
        a.push_back(A.value[i+1][i]);
    }
    b.push_back(A.value[n - 1][n - 1]);
    c.push_back(0);
    return {a, b, c};
}
template <class DB> vector<vector<DB>> TridiagonalSeparation(const vector<vector<DB>> &T)
{
    int n = T[0].size();
    vector<DB> l(n);
    vector<DB> u(n);
    l[0] = 0;
    u[0] = T[1][0];
    for (int i = 1; i < n;++i)
    {
        l[i] = T[0][i] / u[i - 1];
        u[i] = T[1][i] - l[i] * T[2][i - 1];
    }
    return {l, u,T[2]};
}
template <class DB> Matrix<DB> CholeskyCompactDecomposition(const Matrix<DB> &A)
{
    if(A.row_num!=A.column_num)
    {
        cerr << "错误：只有方阵才能进行Cholesky分解。" << '\n';
        return A;
    }
    int n = A.row_num;
    Matrix<DB> L(n,n);
    DB s = 0;
    for (int j = 0; j < n;++j)
    {
        s = 0;
        for (int k = 0; k < j;++k)
        {
            s += L.value[j][k]*L.value[j][k];
        }
        if (A.value[j][j] - s <= 0)
        {
            cerr << "错误：这不是一个正定矩阵，不能Cholesky分解" << '\n';
            return A;
        }
        L.value[j][j] = Sqrt(A.value[j][j] - s);
        for (int i = j+1; i < n;++i)
        {
            s=0;
            for (int k = 0; k <j ;++k)
            {
                s += L.value[i][k] * L.value[j][k];
            }
            L.value[i][j] = (A.value[i][j] -s) / L.value[j][j];
        }    
    }
    return L;
}
template <class DB> vector<Matrix<DB>> CholeskyDecomposition(const Matrix<DB> &A)
{
    const Matrix<DB> &L=CholeskyCompactDecomposition(A);
    return {L,Transpose(L)};
}
template <class DB> vector<Matrix<DB>> JacobiIterationMatrix(const Matrix<DB> &A,const Matrix<DB>&b)
{
    int n = A.row_num;
    int m = b.column_num;
    Matrix<DB> B(n,n);
    Matrix<DB> f(n, m);
    Matrix<DB> x(n, m);
    for (int i = 0; i < n;++i)
    {
        if(A.value[i][i]==0)
        {
            cerr << "错误：Jacobi迭代法需要系数矩阵的主对角元不为0" << '\n';
            return {A,b};
        }
        for (int j = 0;j<i;++j)
        {
            B.value[i][j] = -A.value[i][j] / A.value[i][i];
        }
        B.value[i][i]=0;
        for (int j = i + 1; j < n;++j)
        {
            B.value[i][j] = -A.value[i][j] / A.value[i][i];
        }
    }
    for (int j = 0; j < m;++j)
    {
        for (int i = 0; i < n;++i)
        {
            f.value[i][j] = b.value[i][j] / A.value[i][i];
        }
    }
    return {B, f};
}
template <class DB> vector<Matrix<DB>> GSIterationMatrix(const Matrix<DB> &A,const Matrix<DB> &b)
{
    return {LSolve(A, USplit(A)), LSolve(A, b)};
}
template <class DB> Matrix<DB> LUSolve(const Matrix<DB> &LU,const Matrix<DB> &b)
{
    int n = LU.row_num;
    int m = b.column_num;
    Matrix<DB> y(n, m);
    Matrix<DB> x(n,m);
    DB s = 0;
    for (int j = 0; j < m;++j)
    {
        for (int i = 0; i < n;++i)
        {
            s = 0;
            for (int k = 0; k < i;++k)
            {
                s += LU.value[i][k] *y.value[k][j];
            }
            y.value[i][j] = (b.value[i][j] - s);
        }
        for (int i = n-1; i>=0;--i)
        {
            s = 0;
            for (int k =i+1 ; k <n;++k)
            {
                s += LU.value[i][k] * x.value[k][j];
            }
            x.value[i][j] = (y.value[i][j] - s) / LU.value[i][i];
        }
    }
    return x;
}
template <class DB> Matrix<DB> LUPSolve(const Matrix<DB> &LU, const vector<int> &pi,const Matrix<DB> &b)
{
    int n = LU.row_num;
    int m = b.column_num;
    Matrix<DB> y(n, m);
    Matrix<DB> x(n,m);
    DB s = 0;
    for (int j = 0; j < m;++j)
    {
        for (int i = 0; i < n;++i)
        {
            s = 0;
            for (int k = 0; k < i;++k)
            {
                s += LU.value[i][k] *y.value[k][j];
            }
            y.value[i][j] = (b.value[pi[i]][j] - s);
        }
        for (int i = n-1; i>=0;--i)
        {
            s = 0;
            for (int k =i+1 ; k <n;++k)
            {
                s += LU.value[i][k] * x.value[k][j];
            }
            x.value[i][j] = (y.value[i][j] - s) / LU.value[i][i];
        }
    }
    return x;
}
template<class DB> Matrix<DB> LSolve(const Matrix<DB> &L,const Matrix<DB> &b)
{
    int n = L.row_num;
    int m = b.column_num;
    Matrix<DB> y(n, m);
    DB s = 0;
    for (int j = 0; j < m;++j)
    {
        for (int i = 0; i < n;++i)
        {
            s = 0;
            for (int k = 0; k < i;++k)
            {
                s += L.value[i][k] *y.value[k][j];
            }
            y.value[i][j] = (b.value[i][j] - s)/L.value[i][i];
        }
        
    }
    return y;
}
template<class DB> Matrix<DB> USolve(const Matrix<DB> &U,const Matrix<DB> &y)
{
    int n = U.row_num;
    int m = y.column_num;
    Matrix<DB> x(n,m);
    DB s = 0;
    for (int j = 0; j < m;++j)
    {
        for (int i = n-1; i>=0;--i)
        {
            s = 0;
            for (int k =i+1 ; k <n;++k)
            {
                s += U.value[i][k] * x.value[k][j];
            }
            x.value[i][j] = (y.value[i][j] - s) / U.value[i][i];
        }
    }
    return x;
}
template <class DB> Matrix<DB> TridiagonalSolve(const vector<vector<DB>> &T, const Matrix<DB> &f)
{
    const vector<vector<DB>> &LUC = TridiagonalSeparation(T);
    int n = f.row_num;
    int m = f.column_num;
    Matrix<DB> y(n,m);
    Matrix<DB> x(n, m);
    for (int j = 0; j < m;++j)
    {
        y.value[0][j] = f.value[0][j];
        for (int i = 1; i < n;++i)
        {
            y.value[i][j]=f.value[i][j] - LUC[0][i] * y.value[i - 1][j];
        }
        x.value[n - 1][j] = y.value[n - 1][j]/ LUC[1][n - 1];
        for (int i = n - 2; i >= 0;--i)
        {
            x.value[i][j] = (y.value[i][j]-LUC[2][i]*x.value[i+1][j]) / LUC[1][i];
        }
    }
    return x;
}
template <class DB> Matrix<DB> Iteration(const Matrix<DB> &B,const Matrix<DB> &f)
{
    Matrix<DB> x = f;
    int times = 10;//迭代次数
    for (int i = 1; i < times;++i)
    {
        x=B*x+f;
    }
    return x;
}
template<class DB> Matrix<DB> LinearSolve(const Matrix<DB> &A ,const Matrix<DB> &b,string str,double w)
{
    if(A.row_num!=b.row_num)
    {
        cerr << "错误：系数矩阵与列向量的行数必须相等。" << '\n';
        return b;
    }
    if(A.row_num>A.column_num)
    {
        cerr<<"错误：系数矩阵的行数大于列数，可能是超定的方程组。"<<'\n';
        return b;
    }
    if(A.row_num<A.column_num)
    {
        cerr<<"错误：系数矩阵的行数小于列数，可能是欠定的方程组。"<<'\n';
        return b;
    }
    if(str=="SOR")
    {
        if(w<=0 || w>=2)
        {
            cerr<<"错误：参数w不在(0,2)范围内，这会导致SOR方法不收敛。" <<'\n';
            return b;
        }
        int n = A.row_num;
        int m = b.column_num;
        int times=10;//迭代次数
        int k = 0;
        Matrix<DB> x(n,m);
        DB part = 0;
        DB s = 0;
        vector<DB> rightcoeff;
        for (int i = 0; i < n;++i)
        {
            if(A.value[i][i]==0)
            {
                cerr << "错误：SOR迭代法要求主对角元素全不为0" << '\n';
                return b;
            }
            rightcoeff.push_back(w / A.value[i][i]);
        }
        while(k<times)
        {
            for (int j = 0; j < m;++j)
            {
                for (int i = 0; i < n;++i)
                {
                    part = (1 - w) * x.value[i][j];
                    s = b.value[i][j];
                    for (int r = 0; r < i;++r)
                    {
                        s -= A.value[i][r] * x.value[r][j];
                    }
                    for (int r = i + 1; r < n;++r)
                    {
                        s -= A.value[i][r] * x.value[r][j];
                    }
                    s *=rightcoeff[i]; 
                    x.value[i][j] = part + s;
                }
            }
            ++k;
        }
        return x;
    }
    cerr<<"错误：未定义该方法。"<<'\n';
    return b;
}
template<class DB> Matrix<DB> LinearSolve(const Matrix<DB> &A,const Matrix<DB> &b,string str)
{
    if(A.row_num!=b.row_num)
    {
        cerr << "错误：系数矩阵与列向量的行数必须相等。" << '\n';
        return b;
    }
    if(A.row_num>A.column_num)
    {
        cerr<<"错误：系数矩阵的行数大于列数，可能是超定的方程组。"<<'\n';
        return b;
    }
    if(A.row_num<A.column_num)
    {
        cerr<<"错误：系数矩阵的行数小于列数，可能是欠定的方程组。"<<'\n';
        return b;
    }
    if(str=="Gauss")
    {
        vector<Matrix<DB>> Ab = GaussianElimination(A, b);
        return USolve(Ab[0], Ab[1]);
    } 
    if(str=="LU")
    {
        return LUSolve(LUCompactDecomposition(A), b);
    }
    if(str=="LUP")
    {
        vector<int> pi;
        int rank;
        DB det;
        Matrix<DB> LU=LUPCompactDecomposition(A, pi, rank, det);
        if(rank<A.row_num)
        {
            cerr << "错误：遇到奇异矩阵，无法求解。" << '\n';
            return b;
        }
        return LUPSolve(LU, pi, b);
    }
    if(str=="Thomas"|| str=="chase" || str=="TDMA")
    {
        return TridiagonalSolve(Tridiagonal(A), b);
    }
    if(str=="Cholesky"|| str=="squareroot")
    {
        Matrix<DB> L = CholeskyCompactDecomposition(A);
        return USolve(Transpose(L), LSolve(L, b));
    }
    if(str=="Jacobi")
    {
        /*vector<Matrix<DB>> Bf = JacobiIterationMatrix(A, b);
        return Iteration(Bf[0],Bf[1]);*/
        int n = A.row_num;
        int m = b.column_num;
        int times = 10;//迭代次数
        int k = 0;
        Matrix<DB> xpre(n, m);
        Matrix<DB> x(n, m);
        DB s = 0;
        for (int i = 0; i < n;++i)
        {
            if(A.value[i][i]==0)
            {
                cerr << "错误：Jacobi迭代法要求主对角元素全不为0" << '\n';
                return b;
            }
        }
        while(k<times)
        {
            for (int j = 0; j < m;++j)
            {
                for (int i = 0; i < n;++i)
                {
                    s= b.value[i][j];
                    for (int r = 0; r < i;++r)
                    {
                        s-= A.value[i][r] * xpre.value[r][j];
                    }
                    for (int r = i + 1; r < n;++r)
                    {
                        s-= A.value[i][r] * xpre.value[r][j];
                    }
                    x.value[i][j] = s/ A.value[i][i];
                }
            }
            xpre = x;
            ++k;
        }
        return x;
    }
    if(str=="Gauss-Seidel"|| str=="G-S" || str=="GS")
    {
        int n = A.row_num;
        int m = b.column_num;
        int times = 10;//迭代次数
        int k = 0;
        Matrix<DB> x(n, m);
        DB s = 0;
        for (int i = 0; i < n;++i)
        {
            if(A.value[i][i]==0)
            {
                cerr << "错误：Gauss-Seidel迭代法要求主对角元素全不为0" << '\n';
                return b;
            }
        }
        while(k<times)
        {
            for (int j = 0; j < m;++j)
            {
                for (int i = 0; i < n;++i)
                {
                    s= b.value[i][j];
                    for (int r = 0; r < i;++r)
                    {
                        s-= A.value[i][r] * x.value[r][j];
                    }
                    for (int r = i + 1; r < n;++r)
                    {
                        s-= A.value[i][r] * x.value[r][j];
                    }
                    x.value[i][j] = s/ A.value[i][i];
                }
            }
            ++k;
        }
        return x;
    }
    cerr<<"错误：未定义该方法。"<<'\n';
    return b;
}
template<class DB> Matrix<DB> LinearSolve(const Matrix<DB> &a,const Matrix<DB> &b)
{
    return LinearSolve(a, b, "LUP");
}

//求逆
template<class DB> Matrix<DB> Eye(int n)
{
    if(n<=0)
    {
        cerr << "错误：矩阵的规模必须是正数。" << '\n';
        Matrix<DB> e(1, 1);
        return e;
    }
    Matrix<DB> e(n, n);
    for (int i = 0; i < n;++i)
    {
        e.value[i][i] = 1;
    }
    return e;
}
template<class DB> Matrix<DB> Inverse(const Matrix<DB> &a)
{
    if(a.row_num!=a.column_num)
    {
        cerr << "错误：只有方阵才具有逆矩阵" << '\n';
        return a;
    }
    return LinearSolve(a, Eye<DB>(a.row_num));
}
template<class DB> Matrix<DB> operator/(const Matrix<DB> &a,const Matrix<DB>&b)
{
    if(b.column_num==1 && b.row_num==1)
        return a / Get(b, 0, 0);
    return a * Inverse(b);
}
template<class DB> Matrix<DB> operator^(const Matrix<DB> &a,int n)
{
    if(a.column_num!=a.row_num)
    {
        cerr << "错误：只有方阵才能进行幂运算。" << '\n';
        return a;
    }
    Matrix<DB> m = a;
    Matrix<DB> b = Eye<DB>(a.row_num);
    if(n<0)
    {
        //m = Inverse(m);
        //n = -n;
        cerr << "施工：不打算实现负数次方幂了。每定义一个类都要重载一大堆运算符，累死" << '\n';
    }
    while(n>=1)
    {
        if(n&1)
        {
            b = m * b;
        }
        n=n>>1;
        m = m * m;
    }
    return b;
}
//范数、条件数
template <class DB> DB Trace(const Matrix<DB> &a)
{
    if(a.row_num!=a.column_num)
    {
        cerr << "只有方阵才有迹！" << '\n';
        return 0;
    }
    DB s = 0;
    for (int i = 0; i < a.row_num;++i)
    {
        s += a.value[i][i];
    }
    return s;
}
template <class DB> double Norm(const Matrix<DB> &a, double p)
{
    if(p==INFINITY)
    {
        double max=0;
        double s= 0;
        for (int i = 0; i < a.row_num;++i)
        {
            s = 0;
            for (int j = 0; j < a.column_num;++j)
            {
                s += Abs(a.value[i][j]);
            }
            if(s>max)
                max = s;
        }
        return max;
    }
    if(p==1)
    {
        double max=0;
        double s= 0;
        for (int j = 0; j < a.column_num;++j)
        {
            s = 0;
            for (int i = 0; i < a.row_num;++i)
            {
                s += Abs(a.value[i][j]);
            }
            if(s>max)
                max = s;
        }
        return max;
    }
    cerr << "错误：待开发。" << '\n';
    return 0;
}
template <class DB> double Norm(const Matrix<DB> &a, string str)
{
    if(str=="Frobenius")
    {
        return Sqrt(Trace(a * Transpose(a)));
    }
    cerr << "错误：矩阵范数没有该方法。" << '\n';
    return 0;
}
template <class DB> double Cond(const Matrix<DB> &a, double p)
{
    return Norm(a, p) * Norm(Inverse(a), p);
}
template <class DB> double Cond(const Matrix<DB> &a, string str)
{
    return Norm(a, str) * Norm(Inverse(a), str);
}

//插值
template <class DB> vector<std::function<DB(DB)>> LagrangeBasis(const Matrix<DB> &A)
{
    int n = A.row_num;
    vector<std::function<DB(DB)>> r;
    std::function<DB(DB)> f;
    for (int i = 0; i <n;++i)
    {
        f = [=](DB x) {
            DB s1=DB(1);
            DB s2=DB(1);
            for (int j = 0; j < i;++j)
            {
                s1=s1*(x-A.value[j][0]);
                s2=s2*(A.value[i][0] - A.value[j][0]);
            }
            for (int j = i + 1;j<n;++j)
            {
                s1=s1*(x-A.value[j][0]);
                s2=s2*(A.value[i][0] - A.value[j][0]);
            }
            return s1 / s2;
        };
        r.push_back(f);
    }
    return r;
}
template <class DB> DB LBasisToInterp(const vector<std::function<DB(DB)>> &lb,const Matrix<DB> &A,DB x)
{
    int n = A.row_num;
    DB s = lb[0](x) * A.value[0][1];
    for (int i = 1; i < n;++i)
    {
        s += lb[i](x) * A.value[i][1];
    }
    return s;
}
template <class DB> Matrix<DB> DifferenceQuotientTable(const Matrix<DB> & A)
{
    int n = A.row_num;
    Matrix<DB> table(n,n);
    for (int i = 0; i < n;++i)
        table.value[i][0] = A.value[i][1];
    for (int j = 1; j <n;++j)
    {
        for (int i = j; i < n;++i)
        {
            table.value[i][j] =(table.value[i][j-1]-table.value[i-1][j-1])/(A.value[i][0]-A.value[i-j][0]);
        }
    }
    return table;
}
template <class DB> Matrix<DB> DividedDifferenceTable(const Matrix<DB> &A)
{
    int data_num= A.row_num;
    int der_deg = A.column_num-1;
    int n = data_num * der_deg;
    int factorial = 1;
    vector<DB> z(n);
    for (int i = 0; i <n;++i)
        z[i] = A.value[i / der_deg][0];
    Matrix<DB> table(n, n);
    for (int j = 0; j <der_deg;++j)
    {
        for (int i = j; i < n;++i)
        {
            if((i%der_deg)<j)
                table.value[i][j]=(table.value[i][j-1]-table.value[i-1][j-1])/(z[i]-z[i-j]);
            else
                table.value[i][j] = A.value[i / der_deg][j + 1]/factorial;
        }
        factorial = factorial * (j + 1);
    }
    for(int j=der_deg;j<n;++j)
        for (int i = j; i < n;++i)
            table.value[i][j]=(table.value[i][j-1]-table.value[i-1][j-1])/(z[i]-z[i-j]);
    return table;
}
template <class DB> DB DifferenceQuotient(const Matrix<DB> &A)
{
    int n = A.row_num;
    return DifferenceQuotientTable(A).value[n-1][n-1];
}
template <class DB> DB DifferenceQuotient(std::function<DB(DB)> f, const Matrix<DB> &x)
{
    int n = x.row_num;
    Matrix<DB> A(n,2);
    for (int i = 0; i < n;++i)
    {
        A.value[i][0] = x.value[i][0];
        A.value[i][1] = f(x.value[i][0]);
    }
    return DifferenceQuotient(A);
}
template <class DB> DB DQTableToInterp(const Matrix<DB> &A, const Matrix<DB> &table,DB x)
{
    int n = A.row_num;
    DB s = table.value[n - 1][n - 1];
    for (int i = n - 2; i >= 0;--i)
    {
        s = table.value[i][i] + (x -A.value[i][0]) * s;
    }
    return s;
}
template <class DB> DB DDTableToInterp(const Matrix<DB> &A, const Matrix<DB> &table,DB x)
{
    int data_num = A.row_num;
    int der_deg=A.column_num-1;
    int n = data_num * der_deg;
    vector<DB> nest(n);
    nest[0] = DB(1);
    for(int i=0;i<n-1;++i)
        nest[i+1] = nest[i] * (x - A.value[i/der_deg][0]);
    DB s=DB(0);
    for (int i = 0;i<n;++i)
        s += nest[i] * table.value[i][i];
    return s;
}
template <class DE> DE _2P3DHermiteInterp(const Matrix<DE> &M,int i,int j,DE x)
{
    DE l0x = (x-M.value[j][0]) / (M.value[i][0] - M.value[j][0]);
    DE l1x = (x - M.value[i][0]) / (M.value[j][0] - M.value[i][0]);
    DE l0x2 = l0x * l0x;
    DE l1x2 = l1x * l1x;
    return M.value[i][1] * (1 + 2 * l1x) * l0x2 + M.value[j][1] * (1 + 2 * l0x) * l1x2 +
           M.value[i][2] * (x - M.value[i][0]) * l0x2 + M.value[j][2] * (x - M.value[j][0]) *l1x2;
}
template <class DE> DE _3P3DHermiteInterp(const Matrix<DE> &M,DE x)
{
    std::function<DE(DE)> L2 = [=](DE t) {
        return DQTableToInterp(M, DifferenceQuotientTable(M),t);
    };
    DE k = (M.value[1][2]-D(L2,M.value[1][0]))/
           ((M.value[1][0]-M.value[0][0])*(M.value[1][0]-M.value[2][0]));
    return L2(x) + k * (x - M.value[0][0]) * (x - M.value[1][0]) * (x - M.value[2][0]);
}
template <class DE> Matrix<DE> SplineSlope(const Matrix<DE> &A)
{
    int n = A.row_num;
    vector<DE> h(n);
    vector<DE> lambda(n);
    vector<DE> mu(n);
    Matrix<DE> g(n, 1);
    h[0] = A.value[1][0] - A.value[0][0];
    lambda[0]=DE(0);
    mu[0] = DE(1);
    g.value[0][0] = 3 * (A.value[1][1]-A.value[0][1]) / h[0];
    for (int k = 1; k<= n - 2;++k)
    {
        h[k]=A.value[k+1][0]-A.value[k][0];
        lambda[k] = h[k] / (h[k] + h[k - 1]);
        mu[k]=1-lambda[k];
        g.value[k][0]=3*(lambda[k]*(A.value[k][1]-A.value[k-1][1])/h[k-1]
        +mu[k]*(A.value[k+1][1]-A.value[k][1])/h[k]);
    }
    lambda[n-1]=DE(1);
    mu[n-1]=DE(0);
    g.value[n - 1][0] = 3 * (A.value[n - 1][1] - A.value[n - 2][1]) / h[n - 2];
    Matrix<DE> m=TridiagonalSolve({lambda, vector<DE>(n, DE(2)), mu}, g);
    return m;
}
template <class DE> DE Interpolation(const Matrix<DE> &A,DE x,string str)
{
    if(A.row_num<2)
    {
        cerr << "错误：至少需要两个点才能插值" << '\n';
        return x;
    }
    if(A.column_num<2)
    {
        cerr<<"错误：传入的数据矩阵至少两列"<<'\n';
        return x;
    }
    if(str=="Lagrange")
    {
        return LBasisToInterp(LagrangeBasis(A), A, x);
    }
    if(str=="Newton")
    {
        return DQTableToInterp(A, DifferenceQuotientTable(A),x);
    }
    if(str=="Hermite")
    {
        return DDTableToInterp(A, DividedDifferenceTable(A), x);
    }
    if(str=="non-standard Hermite")
    {
        return _3P3DHermiteInterp(A, x);
    }
    if(str=="linear" || str=="piecewise Lagrange")
    {
        int n = A.row_num;
        int k= 1;
        while(k<n)
        {
            if(x<A.value[k][0])
            {
                return A.value[k - 1][1] * (x - A.value[k][0]) / (A.value[k - 1][0] - A.value[k][0]) + 
                       A.value[k][1] * (x - A.value[k-1][0]) / (A.value[k][0] - A.value[k - 1][0]);
            }
            ++k;
        }
        return A.value[n-2][1]*(x-A.value[n-2][0])/(A.value[n-2][0]-A.value[n-1][0])+
               A.value[n-1][1]*(x-A.value[n-2][0])/(A.value[n-1][0]-A.value[n-2][0]);
    }
    if(str=="pchip" || str=="piecewise Hermite")
    {
        int n = A.row_num;
        int i = 1;
        while(i<n)
        {
            if(x<A.value[i][0])
                return _2P3DHermiteInterp(A,i-1,i, x);
            ++i;
        }
        return _2P3DHermiteInterp(A,n-2,n-1, x);
    }
    if(str=="spline")
    {
        int n = A.row_num;
        int k= 1;
        const Matrix<DE> &m = SplineSlope(A);
        Matrix<DE> B(2, 3);
        while(k<n)
        {
            if(x<A.value[k][0])
            {
                for (int i = 0; i < 2;++i)
                {
                    B.value[i][2] = m.value[k - 1 + i][0];
                    for (int j = 0; j < 2;++j)
                        B.value[i][j] = A.value[k - 1 + i][j];
                }
                return _2P3DHermiteInterp(B, 0, 1, x);
            }               
            ++k;
        }
        for (int i = 0; i < 2;++i)
        {
            B.value[i][2] = m.value[n-2+i][0];
            for (int j = 0; j < 2;++j)
                B.value[i][j] = A.value[n-2+ i][j];
        }
        return _2P3DHermiteInterp(B, 0, 1, x);
    }
    cerr<<"错误：没有定义该方法"<<'\n';
    return x;
}
template <class DE> DE Interpolation(const Matrix<DE> &A,DE x)
{
    if(A.column_num<=2)
        return Interpolation(A,x,"linear");
    return Interpolation(A, x, "piecewise Hermite");
}
template <class DE> Matrix<DE> Interpolation(const Matrix<DE> &A,const Matrix<DE> &x,string str)
{
    if(A.row_num<2)
    {
        cerr << "错误：至少需要两个点才能插值" << '\n';
        return x;
    }
    if(A.column_num<2)
    {
        cerr<<"错误：传入的数据矩阵至少两列"<<'\n';
        return x;
    }
    if(str=="Lagrange")
    {
        vector<std::function<DE(DE)>> L = LagrangeBasis(A);
        Matrix<DE> s(x.row_num, 1);
        for (int i = 0; i < x.row_num;++i)
            s.value[i][0] = LBasisToInterp(L, A, x.value[i][0]);
        return s;
    }
    if(str=="Newton")
    {
        Matrix<DE> D = DifferenceQuotientTable(A);
        Matrix<DE> s(x.row_num, 1);
        for (int i = 0; i < x.row_num;++i)
            s.value[i][0] = DQTableToInterp(A, D, x.value[i][0]);
        return s;
    }
    if(str=="Hermite")
    {
        Matrix<DE> D = DividedDifferenceTable(A);
        Matrix<DE> s(x.row_num, 1);
        for (int i = 0; i < x.row_num;++i)
            s.value[i][0] = DDTableToInterp(A, D, x.value[i][0]);
        return s;
    }
    if(str=="non-standard Hermite")
    {
        Matrix<DE> s(x.row_num,1);
        for (int i = 0; i < x.row_num;++i)
            s.value[i][0] = _3P3DHermiteInterp(A, x.value[i][0]);
        return s;
    }
    if(str=="linear" || str=="piecewise Lagrange")
    {
        int n = A.row_num;
        Matrix<DE> s(x.row_num, 1);
        int k= 1;
        int j = 0;
        while(j<x.row_num && k<n)
        {
            if(x.value[j][0]<A.value[k][0])
            {
                s.value[j][0]=A.value[k - 1][1] * (x.value[j][0]- A.value[k][0]) / (A.value[k - 1][0] - A.value[k][0]) + 
                       A.value[k][1] * (x.value[j][0] - A.value[k-1][0]) / (A.value[k][0] - A.value[k - 1][0]);
                ++j;
                continue;
            }
            ++k;
        }
        for (int t = j; t < x.row_num;++t)
        {
            s.value[t][0] =A.value[n-2][1]*(x.value[t][0]-A.value[n-2][0])/(A.value[n-2][0]-A.value[n-1][0])+
               A.value[n-1][1]*(x.value[t][0]-A.value[n-2][0])/(A.value[n-1][0]-A.value[n-2][0]);
        }
        return s;
    }
    if(str=="pchip" || str=="piecewise Hermite")
    {
        int n = A.row_num;
        Matrix<DE> s(x.row_num, 1);
        int i=1;
        int j = 0;
        while(j<x.row_num && i<n)
        {
            if(x.value[j][0]<A.value[i][0])
            {
                s.value[j][0] = _2P3DHermiteInterp(A,i-1,i, x.value[j][0]);
                ++j;
                continue;
            }
            ++i;
        }
        for (int t = j; t < x.row_num;++t)
            s.value[t][0] = _2P3DHermiteInterp(A,n-2,n-1, x.value[t][0]);
        return s;
    }
    if(str=="spline")
    {
        int n = A.row_num;
        Matrix<DE> s(x.row_num, 1);
        int k= 1;
        int j = 0;
        const Matrix<DE> &m = SplineSlope(A);
        Matrix<DE> B(2, 3);
        while(j<x.row_num && k<n)
        {
            if(x.value[j][0]<A.value[k][0])
            {
                for (int i = 0; i < 2;++i)
                {
                    B.value[i][2] = m.value[k - 1 + i][0];
                    for (int j = 0; j < 2;++j)
                        B.value[i][j] = A.value[k - 1 + i][j];
                }
                s.value[j][0] = _2P3DHermiteInterp(B, 0, 1, x.value[j][0]);
                ++j;
                continue;
            }
            ++k;
        }
        for (int t = j; t < x.row_num;++t)
        {
            for (int i = 0; i < 2;++i)
            {
                B.value[i][2] = m.value[n-2+i][0];
                for (int j = 0; j < 2;++j)
                    B.value[i][j] = A.value[n-2+ i][j];
            }
            s.value[t][0]=_2P3DHermiteInterp(B, 0, 1, x.value[t][0]);
        }
        return s;
    }
    cerr<<"错误：没有定义该方法"<<'\n';
    return x;
}
template <class DE> Matrix<DE> Interpolation(const Matrix<DE> &A, const Matrix<DE> &x)
{
    if(ColumnSize(A)<=2)
        return Interpolation(A,x,"linear");
    return Interpolation(A, x, "piecewise Hermite");
}
template <class DE> std::function<DE(DE)> Interpolation(const Matrix<DE> &A,string str)
{
    return [=](DE x) { return Interpolation(A, x, str); };
}

//积分
template<class DI> DI Integrate(std::function<DI(DI)> f,DI a,DI b,string str,int n)
{
    if(str=="Newton-Cotes"|| str=="N-C")
    {
        if(n<=0)
        {
            cerr<<"错误：Newton-Cotes法要求最后一个参数大于0"<<'\n';
            return (b - a) * f(a);
        }
        if(n>6)
        {
            cerr << "警告：插入过多的点可能会造成Newton-Cotes法数值不稳定。为了确保稳定性，已减少插入的点。" << '\n';
            return Integrate(f, a, b, "Newton-Cotes", 6);
        }
        const vector<vector<DI>> CotesTable({
            {0.5,0.5},
            {1.0/6,4.0/6,1.0/6},
            {0.125,0.375,0.375,0.125},
            {7.0/90,32.0/90,12.0/90,32.0/90,7.0/90},
            {19.0/288,75.0/288,50.0/288,50.0/288,75.0/288,19.0/288},
            {41.0/840,216.0/840,27.0/840,272.0/840,27.0/840,216.0/840,41.0/840}
        }); 
        DI s=0;
        DI divpoint = a;
        DI interval = b - a;
        DI h = interval / n;
        for (int k = 0; k <= n;++k)
        {
            s=s+CotesTable[n-1][k]*f(divpoint);
            divpoint = divpoint + h;
        }
        return interval * s;
    }
    if(str=="Gauss-Legendre" || str=="G-L" || str=="Gauss")
    {
        if(n<0)
        {
            cerr << "错误：Gauss-Legendre法要求至少插入一个点" << '\n';
            return (b - a) * f((a + b) / 2);
        }
        if(n>4)
        {
            cerr<<"警告：不支持插入过多点的Gauss-Legendre法。已自动减少插入的点。"<<'\n';
            return Integrate(f, a, b, "Gauss", 4);
        }
        const vector<vector<DI>> GaussPoint({
            {0},
            {0.5773502692,-0.5773502692},
            {0.7745966692,-0.7745966692,0},
            {0.8611363116,-0.8611363116,0.3399810436,-0.3399810436},
            {0.9061798459,-0.9061798459,0.5384693101,-0.5384693101,0}
        });
        const vector<vector<DI>> GaussCoef({
            {2},
            {1,1},
            {5.0/9,5.0/9,8.0/9},
            {0.3478548451,0.3478548451,0.6521451549,0.6521451549},
            {0.2369268851,0.2369268851,0.4786286705,0.4786286705,0.5688888889}
        });
        DI interval=b-a;
        DI halfinterval=interval/2;
        DI midpoint = (a + b) / 2;
        DI s = 0;
        for (int k = 0; k <=n;++k)
        {
            s = s + GaussCoef[n][k] * f(halfinterval*GaussPoint[n][k]+midpoint);
        }
        return halfinterval * s;
    }
    cerr << "错误：没有定义这种方法" << '\n';
    return f(a) * (b - a);
}
template<class DI> DI CompoundIntegrate(std::function<DI(DI)> f,DI a,DI b,string str,int n,int m)
{
    if(m<0)
    {
        cerr<<"错误：需要划分至少一个区间"<<'\n';
        return f(a) * (b - a);
    }
    DI s = 0;
    DI interval = b - a;
    DI h = interval / m;
    DI leftdiv = a;
    DI rightdiv = a + h;
    for (int i = 0; i < m;++i)
    {
        s = s + Integrate(f, leftdiv, rightdiv, str, n);
        leftdiv = rightdiv;
        rightdiv = rightdiv + h;
    }
    return s;
}
template <class DI> DI Integrate(std::function<DI(DI)> f,DI a,DI b,string str)
{
    if(str=="Romberg")
    {
        const DI epsilon = 1e-14;
        vector<DI> pre(1);
        vector<DI> now(1, (b - a) / 2 * (f(a) + f(b)));
        int k= 0;
        int times_max=20;
        DI accu = 0;
        int powof2 = 1;
        int powof4 = 1;
        do{
            ++k;
            pre = now;
            accu = 0;
            powof2 =(powof2<<1);
            for (int j = 0; j <(powof2>>1);++j)
            {
                accu = accu + f(a + (2 * j + 1) * (b - a) / powof2);
            }
            now[0] = 0.5 * pre[0] + (b - a) / powof2 * accu;
            powof4 = 4;
            for (int m= 1; m< k;++m)
            {
                now[m] = 1 / DI(powof4 - 1) * (powof4 * now[m - 1] - pre[m - 1]);
                powof4 = powof4 * 4;
            }
            now.emplace_back(1 / DI(powof4 - 1)*(powof4 * now[k-1] - pre[k-1]));
        }while(k<times_max && Abs(now[k]-pre[k-1])>=epsilon);
        if(k>=times_max)
            cerr << "警告：未在指定次数内达到预期精度"<<'\n';
        return now[k];
    }
    if(str=="trapz")
    {
        return Integrate(f, a, b, "Newton-Cotes", 1);
    }
    if(str=="compound trapz")
    {
        return CompoundIntegrate(f, a, b, "Newton-Cotes", 1, 100);
    }
    if(str=="Simpson")
    {
        return Integrate(f, a, b, "Newton-Cotes", 2);
    }
    if(str=="compound Simpson")
    {
        return CompoundIntegrate(f, a, b, "Newton-Cotes", 2, 100);
    }
    if(str=="Simpson3/8")
    {
        return Integrate(f, a, b, "Newton-Cotes", 3);
    }
    if(str=="compound Simpson3/8")
    {
        return CompoundIntegrate(f, a, b, "Newton-Cotes", 3, 100);
    }
    if(str=="Cotes"||str=="Milne")
    {
        return Integrate(f, a, b, "Newton-Cotes", 4);
    }
    if(str=="compound Cotes" || str=="compound Milne")
    {
        return CompoundIntegrate(f, a, b, "Newton-Cotes", 4, 100);
    }
    cerr << "错误：没有定义这种方法" << '\n';
    return f(a) * (b - a);
}
template <class DI> DI Integrate(std::function<DI(DI)> f,DI a,DI b)
{
    return Integrate(f, a,b,"Romberg");
}
template<class DI> DI CompoundIntegrate(std::function<DI(DI)>f ,DI a,DI b,string str,int m)
{
    if(m<0)
    {
        cerr<<"错误：需要划分至少一个区间"<<'\n';
        return f(a) * (b - a);
    }
    DI s = 0;
    DI interval = b - a;
    DI h = interval / m;
    DI leftdiv = a;
    DI rightdiv = a + h;
    for (int i = 0; i < m;++i)
    {
        s = s + Integrate(f, leftdiv, rightdiv, str);
        leftdiv = rightdiv;
        rightdiv = rightdiv + h;
    }
    return s;
}

//微分
template <class DD> DD D(std::function<DD(DD)> f,DD x,string str)
{
    if(str=="center")
    {
        DD h = 1e-7;
        return (f(x + h) - f(x - h)) / (2 * h);
    }
    if(str=="forward")
    {
        DD h = 1e-7;
        return (f(x + h) - f(x))/h;
    }
    if(str=="backward")
    {
        DD h = 1e-7;
        return (f(x) - f(x-h))/h;
    }
    cerr << "错误：没有定义该方法。" << '\n';
    return x;
}
template <class DD> DD D(std::function<DD(DD)> f, DD x)
{
    return D(f, x, "center");
}
template<class DD> std::function<DD(DD)> D(std::function<DD(DD)> f)
{
    return [=](DD x) { return D(f, x); };
}
template<class DD> std::function<DD(DD)> D(int n,std::function<DD(DD)> f)
{
    if(n==0)
        return f;
    return D(n - 1, D(f));
}
template<class DD> DD D(int n,std::function<DD(DD)>f, DD x)
{
    return D(n, f)(x);
}

//一元函数求根
template <class DD> DD Iteration(std::function<DD(DD)> phi, DD a)
{
    int times = 100;//迭代次数
    int k = 0;
    DD pre=a;
    DD x = phi(a);
    while (Abs(x-pre)>1e-14 && k < times)
    {
        pre = x;
        x = phi(x);
        ++k;
    } 
    if(k>=times)
    {
        cerr<<"警告：未在指定次数内达到预期精度。"<<'\n';
        cerr << "可能是初值的选择使迭代法不收敛，或者迭代函数的选择使收敛速度很慢。" << '\n';
    }
    return x;
}
template<class DD> DD FindRoot(std::function<DD(DD)> f, DD a,DD b,string str)
{
    if(str=="bisection")
    {
        if(a>=b)
        {
            cerr << "错误：区间左端值需小于右端。" << '\n';
            return a;
        }
        if(f(a)*f(b)>0)
        {
            cerr << "错误：区间两端函数值同号，无法使用二分法。" << '\n';
            return a;
        }
        int times = 30;//迭代次数
        DD _left = a;
        DD _right=b;
        int n = 1;
        DD x = (_left + _right) / 2;
        while(n<=times)
        {
            if(f(x)==0) return x;
            if(f(_left)*f(x)<0)
            {
                _right=x;
            }
            else
            {
                _left = x;
            }
            x = (_left + _right) / 2;
            ++n;
        }
        return x;
    }
    if(str=="secant")
    {
        DD epsilon = 1e-14;
        if(Abs(a-b)<=epsilon)
        {
            cerr << "错误：弦截法的前两个初值选取的过于接近。" << '\n';
            return a;
        }
        int times = 100;//迭代次数
        int k = 0;
        DD pre = a;
        DD now = b;
        DD x =0;
        while (Abs(now-pre)>epsilon && k < times)
        {
            x = now - f(now) / (f(now) - f(pre)) * (now - pre);
            pre = now;
            now = x;
            ++k;
        } 
        if(k>=times)
        {
            cerr<<"警告：未在指定次数内达到预期精度。"<<'\n';
            cerr << "可能是初值的选择使迭代法不收敛，或者迭代函数的选择使收敛速度很慢。" << '\n';
        }
        return x;
    }
    cerr << "错误：没有定义该方法。" << '\n';
    return a;
}
template <class DD> DD FindRoot(std::function<DD(DD)> f, DD a,string str)
{
    if(str=="Newton")
    {
        std::function<DD(DD)> phi = [&](DD x) { return x-f(x) / D(f, x); };
        return Iteration(phi, a);
    }
    if(str=="multiplicity")
    {
        std::function<DD(DD)> F = [&](DD x) { return f(x) / D(f, x); };
        std::function<DD(DD)> phi = [&](DD x) { return x - F(x) / D(F, x); };
        return Iteration(phi, a);
    }
    if(str=="simplified Newton")
    {
        DD g = D(f, a);
        std::function<DD(DD)> phi = [&](DD x) { return x - f(x) / g; };
        return Iteration(phi, a);
    }
    if(str=="downhill")
    {
        DD lambda = 1;//下山因子
        DD b = 0;
        int k = 0;
        int times = 30;//下山次数
        for (k = 0; k <times;++k)
        {
            b = a - lambda * f(a) / D(f, a);
            if(Abs(f(b))<Abs(f(a)))
                break;
            lambda = lambda / 2;
        }
        if(k==times)
        {
            cerr << "警告：超出下山次数，未找到合适的下山因子。" << '\n';
        } 
        std::function<DD(DD)> phi=[&](DD x) { return x - lambda * f(x) / D(f, x); };
        return Iteration(phi, a);
    }
    cerr << "错误：没有定义该方法" << '\n';
    return a;
}
template <class DD> DD FindRoot(std::function<DD(DD)> f, DD a)
{
    return FindRoot(f, a, "Newton");
}
template <class DD> DD FindRoot(std::function<DD(DD)> f, DD a,string str,int m)
{
    if(str=="multiplicity")
    {
        std::function<DD(DD)> phi=[&](DD x) { return x - m* f(x) / D(f, x); };
        return Iteration(phi, a);
    }
    cerr << "错误：没有定义该方法" << '\n';
    return a;
}
template<class DD> DD Sqrt(DD c)
{
    std::function<DD(DD)> phi = [&](DD x) { return (x + c / x) / 2; };
    return Iteration(phi,c);
}

//拟合与逼近
template<class DF> DF InnerProductC(const Matrix<DF> &A,const Matrix<DF> &B)
{
    return Get(Transpose(A)*B,0,0);
}
template<class DF> DF InnerProductR(const Matrix<DF> &A,const Matrix<DF> &B)
{
    return Get(A*Transpose(B),0,0);
}
template<class DI,class DF,class DG> DI InnerProduct(DF f,DG g,DI a,DI b)
{
    std::function<DI(DI)> h = [&](DI x) { return f(x) * g(x); };
    return Integrate(h, a, b);
}
template<class DF> Matrix<DF> FindFit(const Matrix<DF> &data,const vector<std::function<DF(DF)>> &phi)
{
    int m=RowSize(data);
    int n = phi.size();
    Matrix<DF> f(1,m);
    vector<Matrix<DF>> Phi(n,Matrix<DF>(1,m));
    Matrix<DF> Gram(n,n);
    Matrix<DF> B(n, 1);
    for(int j=0;j<m;++j)
    {
        for (int i = 0; i < n;++i)
            Phi[i](0,j)= phi[i](Get(data,j,0));   
        f(0, j) = Get(data, j, 1);  
    }    
    for (int i = 0; i < n;++i)
    {
        for (int j =i; j < n;++j)
            Gram(i, j) = InnerProductR(Phi[i], Phi[j]);
        B(i, 0) = InnerProductR(Phi[i], f);
    }
    for (int i = 0; i < n;++i)
        for (int j = 0; j < i;++j)
            Gram(i, j) = Gram(j, i);
    return Transpose(LinearSolve(Gram,B));
}
template<class DF> Matrix<DF> FindFit(std::function<DF(DF)> f,const vector<std::function<DF(DF)>> &phi,DF a,DF b)
{
    int n = phi.size();
    Matrix<DF> Gram(n,n);
    Matrix<DF> B(n, 1);
    for(int i=0;i<n;++i)
    {
        for (int j = i;j<n;++j)
            Gram(i, j) = InnerProduct(phi[i],phi[j],a,b);
        B(i, 0) = InnerProduct(phi[i], f, a, b);
    }
    for (int i = 0; i < n;++i)
        for (int j = 0; j < i;++j)
            Gram(i, j) = Gram(j, i);
    return Transpose(LinearSolve(Gram,B));
}
template<class DF> DF Fit(const Matrix<DF> &data,const vector<std::function<DF(DF)>> &phi,DF x)
{
    int n = phi.size();
    const Matrix<DF> &Coef=FindFit(data,phi);
    DF s = 0;
    for (int i = 0;i<n;++i)
        s = s + Get(Coef,0,i)* phi[i](x);
    return s;
}
template<class DF> DF Fit(std::function<DF(DF)> f,const vector<std::function<DF(DF)>> &phi,DF a,DF b,DF x)
{
    int n = phi.size();
    const Matrix<DF> &Coef=FindFit(f,phi,a,b);
    DF s = 0;
    for(int i=0;i<n;++i)
        s=s+Get(Coef,0,i)*phi[i](x);
    return s;
}
template<class DF> std::function<DF(DF)> Fit(const Matrix<DF> &data,const vector<std::function<DF(DF)>> &phi)
{
    return [=](double x) { return Fit(data, phi, x); };
} 
template<class DF> std::function<DF(DF)> Fit(std::function<DF(DF)> f,const vector<std::function<DF(DF)>> &phi,DF a,DF b)
{
    return [=](double x) { return Fit(f, phi, a, b, x); };
} 
template <class DB> Matrix<DB> PseudoInverse(const Matrix<DB> &A)
{
    Matrix<DB> T = Transpose(A);
    return Inverse(T * A) * T;
}
template <class DB> Matrix<DB> LeastSquares(const Matrix<DB> &A, const Matrix<DB> &y)
{
    Matrix<DB> T = Transpose(A);
    return LinearSolve(T * A, T * y);
}
template <class DB> Polynomial<DB> PolyFit(const Matrix<DB> &data, const vector<Polynomial<DB>> &phi)
{
    int m = data.row_num;
    int n = phi.size();
    Matrix<DB> A(m,n);
    Matrix<DB> y(m, 1);
    for (int i = 0;i<m;++i)
    {
        for (int j = 0;j<n;++j)
        {
            A.value[i][j] =Get(phi[j],data.value[i][0]);
        }
        y.value[i][0] = data.value[i][1];
    }
    const Matrix<DB> &c = LeastSquares(A, y);
    vector<DB> r(n);
    for (int i = 0; i <n;++i)
        r[i]=c.value[i][0];
    return Polynomial<DB>(r);
}
template <class DB> Polynomial<DB> PolyFit(const Matrix<DB> &xy,int deg)
{
    if(xy.column_num!=2)
    {
        cerr << "错误：传入数据矩阵的列数必须为2" << '\n';
        return X;
    }
    int n = xy.row_num;
    Matrix<DB> A(n, deg + 1);
    Matrix<DB> b(n, 1);
    for (int i = 0; i < n;++i)
    {
        b.value[i][0] = xy.value[i][1];
        A.value[i][0] = 1;
        for(int j=1;j<=deg;++j)
        {
            A.value[i][j] = A.value[i][j - 1] * xy.value[i][0];
        }
    }
    const Matrix<DB> &c=LeastSquares(A, b);
    vector<DB> r;
    for (int i = 0; i <= deg;++i)
        r.push_back(c.value[i][0]);
    Polynomial<DB> f(r);
    return f;
}
template <class DB> Polynomial<DB> PolyFit(std::function<DB(DB)> f,const vector<Polynomial<DB>>&phi,DB a,DB b)
{
    int n = phi.size();
    Matrix<DB> Gram(n, n);
    Matrix<DB> B(n, 1);
    for(int i=0;i<n;++i)
    {
        for (int j = i;j<n;++j)
            Gram(i, j) = InnerProduct(phi[i],phi[j],a,b);
        B(i, 0) = InnerProduct(phi[i], f, a, b);
    }
    for (int i = 0; i < n;++i)
        for (int j = 0; j < i;++j)
            Gram(i, j) = Gram(j, i);
    const Matrix<DB> &c = LinearSolve(Gram, B);
    vector<DB> r(n);
    for (int i = 0; i <n;++i)
        r[i]=c.value[i][0];
    return Polynomial<DB>(r);

}

//常微分方程
template<class DS,class DO,class DP> DO Iteration(const DP &phi,const pair<DS,DO> &initial,DS x,int n)
{
    DS h= (x - get<0>(initial))/n;//步长
    DS xk = get<0>(initial);
    DO yk= get<1>(initial);
    for (int i = 0; i < n;++i)
    {
        phi(xk, yk, h);
        xk = xk + h;
    }
    return yk;
}
template<class DS,class DO> DO RKSolve(std::function<DO(DS,DO)> f,const pair<DS,DO>& initial,DS x,int n)
{
    auto phi = [=](DS &xk, DO &yk, DS &h) {
            DS halfh = h / 2.0;
            DS xk_plus_halfh = xk + halfh;
            DO k1=f(xk,yk);
            DO k2 =f(xk_plus_halfh,yk+halfh*k1);
            DO k3 =f(xk_plus_halfh,yk+halfh*k2);
            DO k4=f(xk+h,yk+h*k3);
            yk = yk + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        };
    return Iteration(phi, initial, x, n);
}
template<class DS,class DO,class DP,class DF> 
DO LinearMultiStep(const DP &phi,const DF &f,const pair<DS,DO> &initial,DS x,int n,int startstep)
{
    DS h = (x - get<0>(initial))/n;
    std::deque<DS> xk={get<0>(initial)};
    std::deque<DO> yk={get<1>(initial)};
    std::deque<DO> fk={f(xk[0],yk[0])};
    for (int i =1; i<startstep;++i)
    {
        yk.push_back(RKSolve(f,{xk.back(),yk.back()},xk.back()+h,1));
        xk.push_back(xk.back()+h);
        fk.push_back(f(xk.back(), yk.back()));
    }
    DO y0adpt = yk.back();
    for (int i = startstep; i <=n;++i)
    {
        phi(xk, yk, fk,h,y0adpt);
        xk.push_back(xk.back() + h);
        fk.push_back(f(xk.back(), yk.back()));
        xk.pop_front();
        yk.pop_front();
        fk.pop_front();
    }
    return yk.back();
}
template<class DS,class DO> DO DSolve(std::function<DO(DS,DO)> f,const pair<DS,DO>& initial,DS x,string str)
{
    if(str=="Euler"|| str=="forward Euler" || str=="explicit Euler")
    {
        auto phi = [=](DS &xk, DO &yk, DS &h) 
        {   yk = yk + h * f(xk, yk);
        };
        return Iteration(phi, initial, x, 100);
    }
    if(str=="backward Euler" || str=="implicit Euler")
    {
        auto phi = [=](DS &xk, DO &yk, DS &h) {
            DS x_plus_h = x + h;
            DO y0=yk+h*f(xk,yk);
            DO y1=yk+h*f(x_plus_h,y0);
            yk = yk + h * f(x_plus_h, y1);
        };
        return Iteration(phi, initial, x, 100);
    }
    if(str=="trapz" || str=="trapezoidal")
    {
        auto phi = [=](DS &xk, DO &yk, DS &h) {
            DO fk=f(xk,yk);
            DS halfh = h / 2.0;
            DS x_plus_h = xk + h;
            DO y0 = yk + h * fk;
            DO y1=yk+halfh*(fk+f(x_plus_h,y0));
            yk = yk + halfh * (fk + f(x_plus_h, y1));
        };
        return Iteration(phi, initial, x, 100);
    }
    if(str=="modified Euler")
    {
        auto phi = [=](DS &xk, DO &yk, DS &h) {
            DO k1=f(xk,yk);
            DO k2=f(xk+h,yk+h*k1);
            yk=yk+h/2.0*(k1+k2);
        };
        return Iteration(phi, initial, x, 100);
    }
    if(str=="midpoint")
    {
        auto phi = [=](DS &xk, DO &yk, DS &h) {
            DS halfh = h / 2.0;
            DO k1=f(xk,yk);
            DO k2=f(xk+halfh,yk+halfh*k1);
            yk=yk+h*k2;
        };
        return Iteration(phi, initial, x, 100);
    }
    if(str=="RK4" || str=="Runge-Kutta")
    {
        return RKSolve(f, initial, x, 100);
    }
    if(str=="Adams")
    {
        auto phi = [=](deque<DS> &xk, deque<DO> &yk, deque<DO> &fk, DS &h,DO &y0adpt) {
            yk.push_back(yk[3]+h/24.0*(55.0*fk[3]-59.0*fk[2]+37.0*fk[1]-9.0*fk[0]));
        };
        return LinearMultiStep(phi, f, initial, x, 100, 4);
    }
    if(str=="Milne")
    {
        auto phi = [=](deque<DS> &xk, deque<DO> &yk, deque<DO> &fk, DS &h,DO &y0adpt) {
            yk.push_back(yk[0]+4.0*h/3.0*(2.0*fk[3]-fk[2]+2.0*fk[1]));
        };
        return LinearMultiStep(phi, f, initial, x, 100, 4);
    }
    if(str=="adaptive Adams")
    {
        auto phi = [=](deque<DS> &xk, deque<DO> &yk, deque<DO> &fk, DS &h,DO &y0adpt) {
            DO y0=yk[3]+h/24.0*(55.0*fk[3]-59.0*fk[2]+37.0*fk[1]-9.0*fk[0]);
            DO y1 = y0 + 251.0 * (yk[3] - y0adpt)/270.0;
            DO y2 = yk[3] + h / 24.0 * (9.0 * f(xk[3] + h, y1) + 19.0 * fk[3] - 5.0* fk[2] + fk[1]);
            yk.push_back(y2-19.0*(y2-y0)/270.0);
            y0adpt = y0;
        };
        return LinearMultiStep(phi, f, initial, x, 100, 4);
    }
    if(str=="adaptive Milne-Hamming")
    {
        auto phi = [=](deque<DS> &xk, deque<DO> &yk, deque<DO> &fk, DS &h,DO &y0adpt) {
            DO y0=yk[0]+4.0*h/3.0*(2.0*fk[3]-fk[2]+2.0*fk[1]);
            DO y1=y0 + 112.0*(yk[3]-y0adpt)/121.0;
            DO y2 = (9.0 * yk[3] - yk[1]) / 8.0 + 3.0 * h / 8.0 * (f(xk[3] + h, y1)+ 2.0 * fk[3] - fk[2]);
            yk.push_back(y2 - 9.0 * (y2 - y0) / 121.0);
            y0adpt = y0;
        };
        return LinearMultiStep(phi, f, initial, x, 100, 4);
    }
    cerr<<"错误：没有定义该方法"<<'\n';
    return x;
}
template<class DS,class DO> DO DSolve(std::function<DO(DS,DO)> f,const pair<DS,DO>& initial,DS x)
{
    return DSolve(f, initial, x, "adaptive Milne-Hamming");
}
template<class DS,class DO> std::function<DO(DS)> DSolve(std::function<DO(DS,DO)> f,const pair<DS,DO> &initial,string str)
{
    return [=](DS x) { return DSolve(f, initial, x, str); };
}
template<class DS,class DO> std::function<DO(DS)> DSolve(std::function<DO(DS,DO)> f,const pair<DS,DO> &initial)
{
    return [=](DS x) { return DSolve(f, initial, x); };
}

}
#endif
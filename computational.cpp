#include "computational.h"
#include<algorithm>
#include<string>
#include<cstdint>
#include<chrono>
using std::to_string;

namespace malg{

//杂例I


//大整数类
BigInt BigInt::cutzero()
{
    int i = _digit.size()-1;
    while(i>=0 && _digit[i]=='0')
    {
        _digit.pop_back();
        --i;
    }
    if(i<0)
        _digit = "0";
    return *this;
}
BigInt::BigInt()
{
    _digit= "0";
    _sign = true;
}
BigInt::BigInt(const string &a)
{
    _digit = a;
    reverse(_digit.begin(), _digit.end());
    if(_digit.back()=='-')
    {
        _sign = false;
        _digit.pop_back();
    }
    else
        _sign = true;
    cutzero();
}
BigInt::BigInt(const char * a)
{
    _digit = a;
    reverse(_digit.begin(), _digit.end());
    if(_digit.back()=='-')
    {
        _sign = false;
        _digit.pop_back();
    }
    else
        _sign = true;
    cutzero();
}
BigInt::BigInt(int a)
{
    _digit=to_string(a);
    reverse(_digit.begin(), _digit.end());
    if(a<0)
    {
        _sign = false;
        _digit.pop_back();
    }
    else
        _sign = true;
}
BigInt::BigInt(bool flag)
{
    _digit = "1";
    if(flag)   
        _sign = true;
    else
        _sign = false;
}
BigInt::BigInt(const string &a, bool flag)
{
    _digit = a;
    _sign = flag;
    cutzero();
}
BigInt operator-(const BigInt &a)
{
    BigInt b = a;
    b._sign = !a._sign;
    return b;
}
bool operator==(const BigInt &a, const BigInt &b)
{
    return (a._digit == b._digit) && (a._sign==b._sign);
}
bool operator!=(const BigInt &a, const BigInt &b)
{
    return !(a == b);
}
bool operator<(const BigInt &a,const BigInt &b)
{
    if(a._sign!=b._sign)
        return b._sign;
    bool sign = a._sign;
    if(a._digit.size()<b._digit.size())
        return sign;
    if(a._digit.size()>b._digit.size())
        return !sign;
    int n = a._digit.size();
    for (int i = n - 1; i >= 0;--i)
    {
        if(a._digit[i]-'0'<b._digit[i ]-'0')
            return sign;
        if(a._digit[i]-'0'>b._digit[i]-'0')
            return !sign;
    }
    return !sign;
}
bool operator>(const BigInt &a,const BigInt &b)
{
    return b < a;
}
bool operator<=(const BigInt &a, const BigInt &b)
{
    return (a == b) || (a < b);
}
bool operator>=(const BigInt &a, const BigInt &b)
{
    return (a == b) || (a > b);
}
ostream & operator<<(ostream & os,const BigInt &a)
{
    os << GetString(a);
    return os;
}
BigInt operator>>(const BigInt &a, int b)
{
    BigInt c = a;
    c._digit.insert(0,b, '0');
    return c;
}
BigInt operator+(const BigInt &a,const BigInt &b)
{
    if(a._sign==false && b._sign==true)
        return b - (-a);
    if(a._sign==true && b._sign==false)
        return a - (-b);
    int m = a._digit.size();
    int n = b._digit.size();
    if(n>m)
        return b + a;
    int flag = 0;
    int s = 0;
    int i = 0;
    BigInt r;
    r._digit = "";
    for (; i < n;++i)
    {
        s = a._digit[i]-'0' + b._digit[i]-'0' + flag;
        r._digit.push_back((s % 10)+'0');
        if(s>=10)
            flag = 1;
        else
            flag = 0;
    }
    for (; i < m;++i)
    {
        s =a._digit[i]-'0' + flag;
        r._digit.push_back((s % 10)+'0');
        if(s>=10)
            flag = 1;
        else
            flag = 0;
    }
    if (flag)
        r._digit.push_back(flag + '0');
    r.cutzero();
    r._sign = a._sign;
    return r;
}
BigInt operator-(const BigInt &a,const BigInt &b)
{
    if(a._sign!=b._sign)
        return a+(-b);
    if(a._sign==false)
        return -((-a)-(-b));
    if(a<b)
        return -(b - a);
    int m = a._digit.size();
    int n = b._digit.size();
    int flag = 0;
    int s = 0;
    int i = 0;
    BigInt r;
    r._digit = "";
    for (; i < n;++i)
    {
        s = a._digit[i]-'0' - (b._digit[i]-'0')+10-flag;
        r._digit.push_back((s % 10)+'0');
        if(s<10)
            flag = 1;
        else
            flag = 0;
    }
    for (; i < m;++i)
    {
        s =a._digit[i]-'0'+10- flag;
        r._digit.push_back((s % 10)+'0');
        if(s<10)
            flag = 1;
        else
            flag = 0;
    }
    if (flag)
        r._digit.push_back(flag + '0');
    r.cutzero();
    return r;
}
BigInt Product_DivideConquer(const BigInt &x,const BigInt &y)
{
    BigInt p = x;
    BigInt q = y;
    if(p._digit.size()<q._digit.size())
        p._digit.append(q._digit.size()-p._digit.size(),'0');
    else if(p._digit.size()>q._digit.size())
        q._digit.append(p._digit.size() - q._digit.size(), '0');
    int m=p._digit.size();
    int n =q._digit.size();
    int k = (m >> 1);
    int u = (n >> 1);
    if(m<=8 && n<=8)
    {
        BigInt r;
        string a = p._digit;
        string b = q._digit;
        reverse(a.begin(),a.end());
        reverse(b.begin(), b.end());
        b=to_string(stoull(a) * stoull(b));
        reverse(b.begin(),b.end());
        r._digit = b;
        return r;
    }
    BigInt p0, p1, q0, q1;
    p0._digit=p._digit.substr(0,k);
    p1._digit=p._digit.substr(k, m - k);
    q0._digit=q._digit.substr(0,u);
    q1._digit=q._digit.substr(u,n-u);
    BigInt r0= Product_DivideConquer(p0,q0);
    BigInt r1 = Product_DivideConquer(p1,q1);
    BigInt r2 = Product_DivideConquer((p0 + p1), (q0 + q1));
    BigInt r=r0 + ((r2 - r0 - r1) >> k) + (r1 >> (k << 1));
    r.cutzero();
    return r;
}
vector<unsigned long long> CompressBit(const string &a_str)
{
    /*uint32_t n = a_str.size();
    uint32_t aq = n/d;
    uint8_t ar = n%d;
    vector<unsigned long long> a(aq);
    string temp;
    for (uint32_t i = 0; i < a.size();++i)
    {
        temp=a_str.substr(i*d, d);
        reverse(temp.begin(), temp.end());
        a[i] =stoull(temp);
    }  
    if(ar)
    {
        temp=a_str.substr(n-ar, ar);
        reverse(temp.begin(),temp.end());
        a.push_back(stoull(temp));
    }*/
    vector<unsigned long long> a(a_str.size());
    std::transform(a_str.begin(), a_str.end(), a.begin(), [](char c)
                   { return c - '0'; });
    return a;
}
BigInt Product_NTT(const BigInt &x,const BigInt &y)
{
    /*if(d>4) d = 4; 
    const vector<uint32_t> mapping({1, 10, 100, 1000, 10000});*/
    auto rr =IntConvolution(CompressBit(x._digit), CompressBit(y._digit));
    auto r = CarryBit<unsigned long long>(rr, 10);
    return BigInt(BitToString(r), (x._sign == y._sign));
}
BigInt operator*(const BigInt &x,const BigInt &y)
{
    uint32_t n=std::max((x._digit).size(),(y._digit).size());
    /*if (n<=32)
        return Product_NTT(x, y, uint8_t(4));
        else if (n<=4096)
            return Product_NTT(x, y, uint8_t(3));     
            else if (n<=262144)
                return Product_NTT(x, y, uint8_t(2));*/
    if (n<=15)
        return Product_DivideConquer(x, y);
    return Product_NTT(x, y);
}
BigInt operator/(const BigInt &p,const BigInt &q)
{
    if(q==BigInt(0))
    {
        cerr << "错误：不能除以0." << '\n';
        return p;
    }
    //未开发
    return p;
}
BigInt operator^(const BigInt &a, int n)
{
    if(n<0)
    {
        cerr << "错误：尚不支持大整数的负数次幂。" << '\n';
        return a;
    }
    BigInt m=a;
    BigInt b("1");
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
BigInt operator^(const BigInt &a,unsigned long long n)
{
    if(n<0)
    {
        cerr << "错误：尚不支持大整数的负数次幂。" << '\n';
        return a;
    }
    BigInt m=a;
    BigInt b("1");
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
BigInt &BigInt::operator+=(const BigInt &b)
{
    *this = *this + b;
    return *this;
}
BigInt &BigInt::operator-=(const BigInt &b)
{
    *this = *this -b;
    return *this;
}
BigInt::operator int()
{
    string _str = _digit;
    reverse(_str.begin(),_str.end());
    if(_sign)
        return stoi(_str);
    return -stoi(_str);
}
void BigInt::Reassign(string str)
{
    *this = BigInt(str);
}
void BigInt::_digitChange(string str)
{
    _digit = str;
}
void BigInt::_signChange(bool sign)
{
    _sign = sign;
}
string Get_digit(const BigInt&a)
{
    return a._digit;
}
string GetString(const BigInt&a)
{
    string _string=a._digit;
    if(a._sign==false)
        _string.push_back('-');
    reverse(_string.begin(), _string.end());
    return _string;
}
bool PositivityTest(const BigInt &a)
{
    return a._sign;
}
unsigned Length(const BigInt &a)
{
    if(a==BigInt(0))
        return 0;
    return a._digit.size();
}

HighPrecision HighPrecision::CutTail()
{
    const string &str = _decimal;
    int d =-1;
    int n = str.size();
    for (int i = 0; i < n;++i)
    {
        if(str[i]!='0')
        {
            d=i;
            break;
        }
    }
    if(d>0)
        _decimal.erase(0,d);
    else if(d<0)
        _decimal = "";
    return *this;
}
HighPrecision::HighPrecision()
{
    _decimal = "";
}
HighPrecision::HighPrecision(string str)
{
    unsigned d = str.size();
    for (unsigned i = 0; i <str.size();++i)
        if(str[i]=='.')
        {
            d = i;
            break;
        }
    Reassign(str.substr(0,d));
    if(d>=str.size()-1)   
        _decimal = "";
    else
    {
        _decimal = str.substr(d + 1);
        reverse(_decimal.begin(), _decimal.end());
    }
    CutTail();
}
HighPrecision::HighPrecision(const char a[])
{
    new (this) HighPrecision(string(a));
}
HighPrecision::HighPrecision(const BigInt &a)
{
    _digitChange(Get_digit(a));
    _signChange(PositivityTest(a));
    _decimal = "";
}
HighPrecision::HighPrecision(double a)
{
    new (this) HighPrecision(to_string(a));
}
HighPrecision::HighPrecision(string a,int n)
{
   *this=(HighPrecision(a))<<n;
}
HighPrecision::HighPrecision(double a,int n)
{
    new (this) HighPrecision(to_string(a), n);
}
ostream &operator<<(ostream &os, const HighPrecision &a)
{ 
    string result = GetString(a);
    if(a._decimal.size()>0)
    {
        string _string = a._decimal;
        reverse(_string.begin(), _string.end());
        result.append(".");
        result.append(_string);
    }
    os << result;
    return os;
}
HighPrecision operator<<(const HighPrecision &a,int n)
{
    if(n==0)
        return a;
    HighPrecision t = a;
    if(n>0)
    {
        const string &cat = t._decimal + Get_digit(t);
        int len = t._decimal.size();
        if(n>len)
        {
            string zeros;
            zeros.insert(zeros.end(), n - len, '0');
            t._decimal = "";
            t._digitChange(zeros + cat);
        }
        else if(n==len)
        {
            t._decimal ="";
            t._digitChange(cat);
        }
        else
        {
            t._decimal = cat.substr(0, len - n);
            t._digitChange(cat.substr(len - n));
        }
        t.cutzero();
        return t;
    }
    int len= Get_digit(t).size();
    string cat=t._decimal+Get_digit(t);
    cat.insert(cat.end(),-n,'0');
    t._digitChange(cat.substr(cat.size() - len));
    t._decimal=cat.substr(0, cat.size() - len);
    t.cutzero();
    t.CutTail();
    return t;
}
HighPrecision operator+(const HighPrecision &a, const HighPrecision &b)
{
    string str1=a._decimal;
    string str2=b._decimal;
    if(str1.size()<str2.size())
        str1.insert(str1.begin(),str2.size()-str1.size(),'0');
    else if(str2.size()<str1.size())
        str2.insert(str2.begin(),str1.size()-str2.size(),'0');
    unsigned dsize= str1.size();
    str1.append(Get_digit(a));
    str2.append(Get_digit(b));
    BigInt c, d;
    c._digitChange(str1);
    d._digitChange(str2);
    c._signChange(PositivityTest(a));
    d._signChange(PositivityTest(b));
    HighPrecision r(c+d);
    const string &str = Get_digit(r);
    if(dsize<str.size())
    {
        r._decimal = str.substr(0,dsize);
        r._digitChange(str.substr(dsize));
    }
    else
    {
        r._decimal = str;
        r._decimal.insert(r._decimal.end(), dsize - str.size(), '0');    
        r._digitChange("0");
    }
    r.CutTail();
    return r;
}
HighPrecision operator-(const HighPrecision &a)
{
    HighPrecision b = a;
    b._signChange(!PositivityTest(a));
    return b;
}
HighPrecision operator-(const HighPrecision &a,const HighPrecision &b)
{
    return a +(-b);
}
HighPrecision operator*(const HighPrecision &a, const HighPrecision &b)
{
    const string &str1 = a._decimal + Get_digit(a);
    const string &str2 = b._decimal + Get_digit(b);
    BigInt c, d;
    c._digitChange(str1);
    d._digitChange(str2);
    c._signChange(PositivityTest(a));
    d._signChange(PositivityTest(b));
    HighPrecision r(c*d);
    const string &str=Get_digit(r);
    unsigned dsize=a._decimal.size() + b._decimal.size();
    if(dsize<str.size())
    {
        r._decimal = str.substr(0,dsize);
        r._digitChange(str.substr(dsize));
    }
    else
    {
        r._decimal = str;
        r._decimal.insert(r._decimal.end(), dsize - str.size(), '0');    
        r._digitChange("0");
    }
    r.CutTail();
    return r;
}
HighPrecision operator/(const HighPrecision &a, const HighPrecision &b)
{
    if(BigInt(b)==0 && b._decimal.size()==0)
    {
        cerr << "错误：除以0" << '\n';
        return HighPrecision("0");
    }  
    const HighPrecision two("2");
    const unsigned times = std::max(double(4),std::log2(SignificantLength(b)));
    int order_b = Order(b);
    const HighPrecision B = b <<(-order_b);
    cout << "B=" << B << std::endl;
    HighPrecision x(100.0/stod(TopKDigit(b,3)));
    unsigned digitcontrol = 2;
    if(x._decimal.size()>digitcontrol)
        x._decimal.erase(0,x._decimal.size()-digitcontrol);
    for (unsigned i = 0; i < times;++i)
    {
        x=x*(two-x*B);
    }
    x=(a*x)<<(-order_b);
    unsigned sflx = SignificantLength(x);
    if(sflx<=1)
        return x;
    x = SignificantFigure(x, sflx >> 1);
    x.CutTail();
    return x;
}
unsigned DecimalLength(const HighPrecision &a)
{
    return a._decimal.size();
}
BigInt IntegerPart(const HighPrecision &a)
{
    BigInt r;
    r._digitChange(Get_digit(a));
    r._signChange(PositivityTest(a));
    return r;
}
HighPrecision DecimalPart(const HighPrecision &a)
{
    HighPrecision r=a;
    r._digitChange("0");
    return r;
}
int Order(const HighPrecision &a)
{
    const BigInt &intpart=BigInt(a);
    if(intpart!=BigInt("0"))
        return Length(a)-1;
    return -(a._decimal.size());
}
string TopKDigit(const HighPrecision &a,unsigned k)
{
    string result;
    if(BigInt(a)==BigInt(0) && a._decimal.size()==0)
    {
        result.insert(result.end(), k, '0');
        return result;
    }
    if(!PositivityTest(a))
        result.push_back('-');
    const string &cat = a._decimal + Get_digit(a);
    int n = cat.size();
    unsigned cnt = 0;
    for (int i = n - 1; i >= 0;--i)
    {
        if(cat[i]!='0' && cnt<k)
        {
            result.push_back(cat[i]);
            ++cnt;
        }
    }
    if(cnt<k)
    {
        result.insert(result.end(), k - cnt, '0');
    }
    return result;
}
unsigned TotalLength(const HighPrecision &a)
{
    if(BigInt(a)==BigInt(0) && a._decimal.size()==0)
        return 0;
    return Length(a) + a._decimal.size();
}
unsigned SignificantLength(const HighPrecision & a)
{
    const string &str = a._decimal + Get_digit(a);
    int n = str.size();
    int lft=-1, rght=-1;
    for (int i = 0; i < n;++i)
        if(str[i]!='0')
        {
            lft = i;
            break;
        }
    if(lft<0)
        return 0;
    for (int i = n-1; i>=0;--i)
        if(str[i]!='0')
        {
            rght = i;
            break;
        }
    return rght - lft + 1;
}
HighPrecision SignificantFigure(const HighPrecision &a, unsigned n)
{
    string str=a._decimal+Get_digit(a);
    int len= str.size();
    unsigned d=len;
    for (int i = len - 1; i >= 0;--i)
        if(str[i]!='0')
        {
            d=i;
            break;
        }
    if(d<0)
        return HighPrecision("0");
    if(d<n)
        return a;
    for (int i = d-n; i >= 0;--i)
        str[i] = '0';
    int decilen = a._decimal.size();
    HighPrecision r;
    r._decimal=str.substr(0,decilen);
    r._digitChange(str.substr(decilen));
    r._signChange(PositivityTest(a));
    r.CutTail();
    return r;
}

//杂例II
BigInt Factorial(int n)
{
    if(n<0)
    {
        cerr<<"错误：负数没有阶乘。"<<'\n';
        return 0;
    }
    BigInt s("1");
    for (int i = 1; i <=n;++i)
    {
        s =s* BigInt(i);
    }
    return s;
}

BigInt Fibonacci_r(int n)
{
    if(n==0)
        return 0;
    if(n==1)
        return 1;
    BigInt pre2(0);
    BigInt pre1(1);
    BigInt s(1);
    for (int i = 2; i <=n;++i)
    {
        s = pre2 + pre1;
        pre2 = pre1;
        pre1 = s;
    }
    return s;
}
BigInt Fibonacci_p(int n)
{
    Matrix<BigInt> A({{1, 1}, {1, 0}}, 2, 2);
    Matrix<BigInt> b({{1},{0}},2, 1);
    return Get((A ^ (n - 1)) ,0, 0);
}
BigInt Fibonacci_m(int n,vector<BigInt> &g)
{
    if(g[n]>0)
        return g[n];
    if(n%2)
        return g[n]=Fibonacci_m(n/2,g)*Fibonacci_m(n/2,g)+Fibonacci_m(n/2+1,g)*Fibonacci_m(n/2+1,g);
    return g[n] = (Fibonacci_m(n / 2 - 1, g) + Fibonacci_m(n / 2 + 1, g)) * Fibonacci_m(n / 2, g);
}
BigInt Fibonacci(int n)
{
    if(n<255)
        return Fibonacci_r(n);
    vector<BigInt> g({0, 1, 1,2,3});
    g.resize(n+1);
    return Fibonacci_m(n, g);
}


//多项式
Polynomial<double> operator*(const Polynomial<double>& f,const Polynomial<double>& g)
{
    return Polynomial<double>(RealConvolution(GetCoef(f), GetCoef(g)));
}
Polynomial<float> operator*(const Polynomial<float>& f,const Polynomial<float>& g)
{
    return Polynomial<float> (RealConvolution(GetCoef(f), GetCoef(g)));
}
Polynomial<complex<double>> operator*(const Polynomial<complex<double>>& f,const Polynomial<complex<double>>& g)
{
    return Polynomial<complex<double>>(Convolution(GetCoef(f), GetCoef(g)));
}
Polynomial<complex<float>> operator*(const Polynomial<complex<float>> &f, const Polynomial<complex<float>> &g)
{
    return Polynomial<complex<float>> (Convolution(GetCoef(f), GetCoef(g)));
}
/*Polynomial<unsigned> operator*(const Polynomial<unsigned> &f, const Polynomial<unsigned> &g)
{
    return Polynomial<unsigned>(IntConvolution(GetCoef(f), GetCoef(g)));
}
Polynomial<unsigned long long> operator*(const Polynomial<unsigned long long> &f, const Polynomial<unsigned long long> &g)
{
    return Polynomial<unsigned long long>(IntConvolution(GetCoef(f), GetCoef(g)));
}*/ /*test*/


//一元函数求根
template<>
double Sqrt<double>(double c)
{
    if(c<0)
    {
        cerr << "错误：根号下不能是负数。" << '\n';
        return NAN;
    }
    double pre=0;
    double x = 0.5*c;
    while (fabs(pre-x)>1e-14)
    {
        pre = x;
        x = 0.5 * (x + c / x);
    } 
    return x;
}

}
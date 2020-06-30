# numerical-library
## Overview 概况
这是一个C++实现的简单的数值库，集合了基础的矩阵运算、方程求解、 插值拟合、数值微积分、微分方程求解等操作。
## Algorithms Used 所用算法
对于下面列出的各种功能，均使用多种算法予以实现。下面列出的算法只是相应功能调用时的默认算法。
* 大整数乘法 —— NTT
* 高精度除法 —— 牛顿迭代法
* 多项式乘法 —— FFT
* 多项式求值 —— 秦九韶算法
* 线性方程组求解 —— LUP分解
* 代数方程(组)求解 —— 牛顿迭代法
* 插值 —— 分段线性插值
* 带导数插值 —— 分段两点三次埃尔米特插值
* 拟合 —— 最小二乘法
* 函数逼近 —— 最小二乘法
* 定积分 —— Romberg算法
* 求导数 —— 中心差商法
* 常微分方程(组)求解 —— 预测-校正的Milne-Hamming公式
## Development Environment 开发环境
C++14 gcc10.1.1
## Change Log 变更日志
目前已发布的最新版本是v0.7.2，详见release
* v0.7.1 首个发布的版本
* v0.7.2 加入常微分方程组求解与代数方程组求解

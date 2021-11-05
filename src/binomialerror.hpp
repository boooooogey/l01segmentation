#pragma once
#include "function.hpp"
#include <limits>

class BinomialError : Function
{
    private:
        //a log(x) + b log(1 - x) + c
        double a;
        double b;
        double c;
    public:
        static double rangeinf;
        static double rangeninf;
        static double domaininf;
        static double domainninf;
        BinomialError();
        BinomialError(const BinomialError& other);
        void set(const double* y, const double* w, const int& i);
        void set(const double& t);
        BinomialError& operator=(const BinomialError& other);
        BinomialError& operator+=(const BinomialError& rhs);
        double operator()(const double& x);
        bool max(double& xprime, double& yprime);
        void solve(const double& t, double& left, double& right, bool& leftexists, bool& rightexists);
        static void converter(double* x, const int& n);
};

BinomialError operator+(BinomialError lhs, const BinomialError& rhs);

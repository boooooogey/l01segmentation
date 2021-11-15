#pragma once
#include "function.hpp"
#include <limits>

class ExponentialError : Function
{
    private:
        //a x + b exp(x) + c
        double a;
        double b;
        double c;
    public:
        static double rangeinf;
        static double rangeninf;
        static double domaininf;
        static double domainninf;
        ExponentialError();
        ExponentialError(const ExponentialError& other);
        ExponentialError(const double* y, const double* w, const int& i);
        ExponentialError(const double& t);
        void set(const double* y, const double* w, const int& i);
        void set(const double& t);
        ExponentialError& operator=(const ExponentialError& other);
        ExponentialError& operator+=(const ExponentialError& rhs);
        double operator()(const double& x);
        bool max(double& xprime, double& yprime);
        void solve(const double& t, double& left, double& right, bool& leftexists, bool& rightexists);
        static void converter(double* x, const int& n);
};

ExponentialError operator+(ExponentialError lhs, const ExponentialError& rhs);

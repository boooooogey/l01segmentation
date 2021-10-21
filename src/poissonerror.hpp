#pragma once
#include "function.hpp"
#include <limits>

class PoissonError : Function
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
        PoissonError();
        PoissonError(const PoissonError& other);
        void set(const double& y, const double& w);
        void set(const double& t);
        PoissonError& operator=(const PoissonError& other);
        PoissonError& operator+=(const PoissonError& rhs);
        double operator()(const double& x);
        bool max(double& xprime, double& yprime);
        void solve(const double& t, double& left, double& right, bool& leftexists, bool& rightexists);
        static void converter(double* x, const int& n);
};

PoissonError operator+(PoissonError lhs, const PoissonError& rhs);

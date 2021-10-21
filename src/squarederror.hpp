#pragma once
#include "function.hpp"
#include <limits>

class SquaredError : Function
{
    private:
        //a x^2 + b x + c
        double a;
        double b;
        double c;
    public:
        static double rangeinf;
        static double rangeninf;
        static double domaininf;
        static double domainninf;
        SquaredError();
        SquaredError(const SquaredError& other);
        void set(const double& y, const double& w);
        void set(const double& t);
        SquaredError& operator=(const SquaredError& other);
        SquaredError& operator+=(const SquaredError& rhs);
        double operator()(const double& x);
        bool max(double& xprime, double& yprime);
        void solve(const double& t, double& left, double& right, bool& leftexists, bool& rightexists);
        static void converter(double* x, const int& n);
};

SquaredError operator+(SquaredError lhs, const SquaredError& rhs);

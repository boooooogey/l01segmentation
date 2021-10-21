#pragma once

class Function{
    public:
        virtual void set(const double& y, const double& w) = 0;
        virtual void set(const double& t) = 0;
        virtual double operator()(const double& x) = 0;
        virtual bool max(double& xprime, double& yprime) = 0;
        virtual void solve(const double& t, double& left, double& right, bool& leftexits, bool& rightexists) = 0;
};

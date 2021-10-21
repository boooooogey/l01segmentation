#include "poissonerror.hpp"
#include <math.h>
#include <iterator>
#include <limits>
#include <boost/math/special_functions/lambert_w.hpp>

double PoissonError::rangeinf = std::numeric_limits<double>::infinity();
double PoissonError::rangeninf = -std::numeric_limits<double>::infinity();
double PoissonError::domaininf = std::numeric_limits<double>::infinity();
double PoissonError::domainninf = -std::numeric_limits<double>::infinity();

PoissonError::PoissonError(){
    a = 0;
    b = 0;
    c = 0;
}

PoissonError::PoissonError(const PoissonError& other){
    a = other.a;
    b = other.b;
    c = other.c;
}

void PoissonError::set(const double& y, const double& w){
    a = y * w;
    b = -w;
    c = 0;
}

void PoissonError::set(const double& t){
    a = 0;
    b = 0;
    c = t;
}

PoissonError& PoissonError::operator=(const PoissonError& rhs){
    this->a = rhs.a;
    this->b = rhs.b;
    this->c = rhs.c;
    return *this;
}

PoissonError& PoissonError::operator+=(const PoissonError& rhs){
    a += rhs.a;
    b += rhs.b;
    c += rhs.c;
    return *this;
}

PoissonError operator+(PoissonError lhs, const PoissonError& rhs){
    lhs += rhs;
    return lhs;
}

double PoissonError::operator()(const double& x){
    return a * x + b * exp(x) + c;
}

bool PoissonError::max(double& xprime, double& yprime){
    if (b == 0){
        if (a > 0){
            xprime = domaininf;
        }
        else if (a < 0){
            xprime = domainninf;
        }
        else{
            return false;
        }
    }
    else if (a == 0){
        if (b > 0){
            xprime = domaininf;
        }
        else if (b < 0){
            xprime = domainninf;
        }
        else{
            return false;
        }
    }
    else if (-a/b > 0){
        xprime = log(-a/b);
    }
    else{
        if (b < 0  && a < 0){
            xprime = domainninf;
        }
        else{
            xprime = domaininf;
        }
    }
    yprime = a * xprime + b * exp(xprime) + c;
    return true;
}

void PoissonError::solve(const double& t, double& left, double& right, bool& leftexists, bool& rightexists){

    double infmin = std::numeric_limits<double>::min();
    using boost::math::lambert_w0;
    using boost::math::lambert_wm1;
    double a = this->a, b=this->b, c = t - this->c; 
    if( a != 0 && b != 0){
        double exin = b * exp(c/a) / a;
        if( -exp(-1) <= exin && exin <= 0 ){
            if (abs(exin) < infmin) exin = 0;
            left = ( c - a * lambert_wm1(exin))/a;
            right = ( c - a * lambert_w0(exin))/a;
            if(left > right){
                double tmp = left;
                left = right;
                right = tmp;
            }
            leftexists = rightexists = true;
        }
        else{
            left = right = 0;
            leftexists = rightexists = false;
        }
    }
    else if( a != 0 && b == 0){
        left = right = c/a;
        //leftexists = rightexists = false;
        leftexists = false;
        rightexists = true;
    }
    else if(a == 0 && b != 0 && c != 0){
        left = right = log(c/b);
        //leftexists = rightexists = true;
        leftexists = true;
        rightexists = false;
    }
    else{
        left = right = 0;
        leftexists = rightexists = false;
    }
}

void PoissonError::converter(double* x, const int& n){
    for(int i = 0; i < n; i++){
        x[i] = exp(x[i]);
    }
}

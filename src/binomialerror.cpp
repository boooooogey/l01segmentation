#include "binomialerror.hpp"
#include <math.h>
#include <iterator>
#include <limits>
#include <boost/math/special_functions/lambert_w.hpp>

void objective(const double& x, const double& a, const double& b, const double& c, double& obj){
    obj = a*x - (a + b) * log(exp(x) + 1) + c; //-b * x - (a + b) * log(1 + exp(-x)) + c;
}

void derivative(const double& x, const double& a, const double& b, double& deriv){
    deriv = (a + b) / (exp(x) + 1) - b;//(-b + a * exp(-x)) / (1 + exp(-x));
}

void newton(double& x, const double& a, const double& b, const double& c){
    double fx, pfx;
    double x0 = x;
    double tol = 1e-6;
    int maxiter = 10;
    while( true ){
        objective(x, a, b, c, fx);
        derivative(x, a, b, pfx);
        x = x - fx/pfx;
        if ( abs(x - x0) <= tol || maxiter == 0 ) break;
        else x0 = x;
        maxiter--;
    }
}

void transformx(double& x){
    x = -log(1 - x) + log(x);
}

void transformy(double& y){
    y = 1/(1+exp(-y));
}

double BinomialError::rangeinf = std::numeric_limits<double>::infinity();
double BinomialError::rangeninf = -std::numeric_limits<double>::infinity();
double BinomialError::domaininf = std::numeric_limits<double>::infinity();
double BinomialError::domainninf = -std::numeric_limits<double>::infinity();

BinomialError::BinomialError(){
    a = 0;
    b = 0;
    c = 0;
}

BinomialError::BinomialError(const BinomialError& other){
    a = other.a;
    b = other.b;
    c = other.c;
}

void BinomialError::set(const double* y, const double* w, const int& i){
    a = y[2*i] * w[i];
    b = (y[2*i+1] - y[2*i]) * w[i];
    c = 0;
}

void BinomialError::set(const double& t){
    a = 0;
    b = 0;
    c = t;
}

BinomialError& BinomialError::operator=(const BinomialError& rhs){
    this->a = rhs.a;
    this->b = rhs.b;
    this->c = rhs.c;
    return *this;
}

BinomialError& BinomialError::operator+=(const BinomialError& rhs){
    a += rhs.a;
    b += rhs.b;
    c += rhs.c;
    return *this;
}

BinomialError operator+(BinomialError lhs, const BinomialError& rhs){
    lhs += rhs;
    return lhs;
}

double BinomialError::operator()(const double& x){
    double logterm;
    if( x > 15 ) logterm = x;
    else logterm = log(exp(x) + 1);
    return a * x - (a + b) * logterm + c;//a * log(x) + b * log(1-x) + c;
}

bool BinomialError::max(double& xprime, double& yprime){
    if(a == 0 && b == 0){
        yprime = c;
        return false;
    }
    else if(a == 0){ 
        xprime = domainninf;
        yprime = c;
    }
    else if(b == 0){
        xprime = domaininf;
        yprime = c;
    }
    else{
        xprime = log(a/b);
        yprime = this->operator()(xprime);
    }
    return true;
}

void BinomialError::solve(const double& t, double& left, double& right, bool& leftexists, bool& rightexists){
    double constant = c - t;
    using boost::math::lambert_w0;
    double tol = 1e-16;
    if( this->operator()(log(a/b))-t <= 0 ){
        leftexists = false;
        rightexists = false;
    }
    else if( a == 0 ){
        leftexists = false;
        rightexists = true;
        if(constant / b > 15) right = constant/b;
        else right = log(exp(constant/b) - 1);
    }
    else if( b == 0 ){
        leftexists = true;
        rightexists = false;
        if(constant / a > 15) left = constant/a;
        else left = log(exp(constant/a) - 1);
    }
    else if( a == b){
        leftexists = true;
        rightexists = true;
        if(constant / a > 15){
            left = -constant;
            right = constant;
        }
        else{
            double delta = sqrt(1 - 4 * exp(-constant/a));
            left = (1 - delta) / 2;
            transformx(left);
            right = (1 + delta) / 2;
            transformx(right);
        }
    }
    else{
        leftexists = true;
        rightexists = true;
        double left_delta = - b/a * exp(-constant/a); if (abs(left_delta) < tol) left_delta = 0;
        double right_delta = - a/b * exp(-constant/b); if (abs(right_delta) < tol) right_delta = 0;
        if(left_delta != 0){
            left = - a/b * lambert_w0(left_delta);
            transformx(left);
            newton(left, a, b, constant);
        }
        else{
            double k = a * log(a) + b * log(b) - (a + b) * log(a + b);
            double yk = log(a) - log(b);
            double delta = sqrt( 2 * (a + b) * (k + constant) / a / b);
            left = yk-delta; 
            //left = tol;
            newton(left, a, b, constant);
        }
        if(right_delta != 0){
            right = 1 + b / a * lambert_w0(right_delta);
            transformx(right);
            newton(right, a, b, constant);
        }
        else{
            double k = a * log(a) + b * log(b) - (a + b) * log(a + b);
            double yk = log(a) - log(b);
            double delta = sqrt( 2 * (a + b) * (k + constant) / a / b);
            right = yk+delta;
            //right = 1 - tol;
            newton(right, a, b, constant);
        }
        //double k = a * log(a) + b * log(b) - (a + b) * log(a + b);
        //double yk = log(a) - log(b);
        //double delta = sqrt( 2 * (a + b) * (k + constant) / a / b);
        //left = yk - delta;
        //newton(left, a, b, constant);
        //left = 1/(1 + exp(-left));
        //right = yk + delta;
        //newton(right, a, b, constant);
        //right = 1/(1 + exp(-right));
    }
}

void BinomialError::converter(double* x, const int& n){
    for(int i = 0; i < n; i++){
        transformy(x[i]);
    }
}

#include "squarederror.hpp"
#include <math.h>
#include <iterator>
#include <limits>

double SquaredError::rangeinf = std::numeric_limits<double>::infinity();
double SquaredError::rangeninf = -std::numeric_limits<double>::infinity();
double SquaredError::domaininf = std::numeric_limits<double>::infinity();
double SquaredError::domainninf = -std::numeric_limits<double>::infinity();

SquaredError::SquaredError(){
    a = 0;
    b = 0;
    c = 0;
}

SquaredError::SquaredError(const SquaredError& other){
    a = other.a;
    b = other.b;
    c = other.c;
}

void SquaredError::set(const double& y, const double& w){
    a = -w;
    b = 2 * y * w;
    c = 0;
}

void SquaredError::set(const double& t){
    a = 0;
    b = 0;
    c = t;
}

SquaredError& SquaredError::operator=(const SquaredError& rhs){
    this->a = rhs.a;
    this->b = rhs.b;
    this->c = rhs.c;
    return *this;
}

SquaredError& SquaredError::operator+=(const SquaredError& rhs){
    a += rhs.a;
    b += rhs.b;
    c += rhs.c;
    return *this;
}

SquaredError operator+(SquaredError lhs, const SquaredError& rhs){
    lhs += rhs;
    return lhs;
}

double SquaredError::operator()(const double& x){
    return (a * x + b) * x + c;
}

bool SquaredError::max(double& xprime, double& yprime){
    xprime = -b/(2 * a);
    yprime = (a * xprime + b) * xprime + c;
    return true;
}

void SquaredError::solve(const double& t, double& left, double& right, bool& leftexists, bool& rightexists){
    double b, c;
    if (this->a == 0){
        left = right = 0;
        leftexists = rightexists = false;
    }
    else{
        b = this->b / this->a;
        c = (this->c - t) / this->a;
        double delta = b * b - 4 * c;
        if (delta == 0){
            left = right = -b / 2;
            leftexists = rightexists = true;
        }
        else if (delta < 0){
            left = right = 0;
            leftexists = rightexists = false;
        }
        else{
            delta = sqrt(delta);
            left = (-b + delta) / 2;
            right = (-b - delta) / 2;
            leftexists = rightexists = true;
            if (left > right){
                double tmp = left;
                left = right;
                right = tmp;
            }
        }
    }
}

void SquaredError::converter(double* x, const int& n){
}

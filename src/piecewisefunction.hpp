#pragma once
#include "range.hpp"
#include <iterator>
#include <limits>
#include <vector>

template <class T>
class PiecewiseFunction
{
    private:
        std::vector<double> knots;
        std::vector<T> pieces;
        int length;
        int capacity;
    public:
        PiecewiseFunction();
        PiecewiseFunction(const T* pieces, const double* knots, const int& length);
        PiecewiseFunction(const PiecewiseFunction& other);
        void resize();
        void reset();
        void append(const T& piece, const double& knot);
        void append(const double* y, const double* w, const double& knot, const int& i);
        void append(const double& t, const double& knot);
        void index(const int& i, T& piece, double& knot);
        int len();
        PiecewiseFunction& operator=(const PiecewiseFunction<T>& other);
        PiecewiseFunction& operator+=(const T& rhs);
        double operator()(const double& x);
        void max(double& xprime, double& yprime);
        void flood(const double& threshold, PiecewiseFunction<T>& out, Range& ranges);
        ~PiecewiseFunction();
};

template <class T>
PiecewiseFunction<T>::PiecewiseFunction(){
    this->pieces.reserve(10);

    this->knots.reserve(10);
}

template <class T>
PiecewiseFunction<T>::PiecewiseFunction(const T* pieces, const double* knots, const int& length){
    this->pieces.resize(length);
    this->knots.resize(length);

    std::copy(pieces, pieces+length, this->pieces.begin());

    std::copy(knots, knots+length, this->knots.begin());

}

template <class T>
PiecewiseFunction<T>::PiecewiseFunction(const PiecewiseFunction& other){
    if(this != &other){
        length = other.length;
        capacity = other.capacity;
        pieces = other.pieces;
        knots = other.knots;
    }
}

template <class T>
void PiecewiseFunction<T>::reset(){
    pieces.resize(0);
    knots.resize(0);
}

template <class T>
void PiecewiseFunction<T>::append(const T& piece, const double& knot){
    pieces.push_back(piece);
    knots.push_back(knot);
}

template <class T>
void PiecewiseFunction<T>::append(const double* y, const double* w, const double& knot, const int& i){
    pieces.push_back(T(y, w, i));
    knots.push_back(knot);
}

template <class T>
void PiecewiseFunction<T>::append(const double& t, const double& knot){
    pieces.push_back(T(t));
    knots.push_back(knot);
}

template <class T>
void PiecewiseFunction<T>::index(const int& i, T& piece, double& knot){
    if(length-1 < i)
        return;
    else{
        piece = pieces[i];
        knot = knots[i];
    }
}

template <class T>
int PiecewiseFunction<T>::len(){
    return pieces.size();
}

template <class T>
PiecewiseFunction<T>& PiecewiseFunction<T>::operator=(const PiecewiseFunction<T>& other){
    pieces.resize(other.pieces.size());
    knots.resize(other.knots.size());
    std::copy(other.pieces.begin(), other.pieces.end(), pieces.begin());
    std::copy(other.knots.begin(), other.knots.end(), knots.begin());
    return *this;
}

template <class T>
PiecewiseFunction<T>& PiecewiseFunction<T>::operator+=(const T& rhs){
    for(int i = 0; i < this->pieces.size(); i++){
        this->pieces[i] += rhs;
    }
    return *this;
}

template <class T>
PiecewiseFunction<T> operator+(PiecewiseFunction<T> lhs, const T& rhs){
    lhs += rhs;
    return lhs;
}

template <class T>
PiecewiseFunction<T> operator+(const T& lhs, PiecewiseFunction<T> rhs){
    rhs += lhs;
    return rhs;
}

template <class T>
double PiecewiseFunction<T>::operator()(const double& x){
    if(pieces.size() == 1){
        return pieces[0](x);
    } 
    int k = 0;
    for(int i = 1; i < pieces.size(); i++){
        if (x > knots[i]){
            k++;
        }
        else{
            break;
        }
    }
    return pieces[k](x);
}

template <class T>
void PiecewiseFunction<T>::max(double& xprime, double& yprime){
    yprime = -std::numeric_limits<double>::infinity();
    xprime = 0;
    double ycurr = 0;
    double xcurr = 0;
    for(int i = 0; i < pieces.size(); i++){
        pieces[i].max(xcurr, ycurr);
        if(i != length-1 && xcurr > knots[i+1]){
            xcurr = knots[i+1];
            ycurr = pieces[i](xcurr);
        }
        else if(xcurr < knots[i]){
            xcurr = knots[i];
            ycurr = pieces[i](xcurr);
        }
        if(yprime < ycurr){
            yprime = ycurr;
            xprime = xcurr;
        }
    }
}
        
template <class T>
void PiecewiseFunction<T>::flood(const double& threshold, PiecewiseFunction<T>& out, Range& ranges){
    bool underwater = true;
    bool leftexists, rightexists;
    double left, right;
    double first, last;

    out.reset();

    if(pieces[0](knots[0]) < threshold){
        first = T::domainninf;
        out.append(threshold, knots[0]);
    }
    else{
        out.append(pieces[0],knots[0]);
        underwater = false;
    }

    for(int i = 0; i < pieces.size()-1; i++){
        pieces[i].solve(threshold, left, right, leftexists, rightexists);
        if(underwater && leftexists && left < knots[i+1] && left > knots[i]){
            out.append(pieces[i], left);
            last = left;
            ranges.add(first, last);
            underwater = false;
        }
        else if ( i != 0 && !underwater){
            out.append(pieces[i], knots[i]);
        }
        if(!underwater && rightexists && right < knots[i+1] && right > knots[i]){
            out.append(threshold, right);
            first = right;
            underwater = true;
        }
    }
    if(underwater){
        last = T::domaininf;
        ranges.add(first,last);
    }
    out.append(T::rangeninf, T::domaininf);
}

template <class T>
PiecewiseFunction<T>::~PiecewiseFunction(){
}


#pragma once
#include "range.hpp"
#include <iterator>
#include <limits>

template <class T>
class PiecewiseFunction
{
    private:
        double* knots;
        T* pieces;
        int length;
        int capacity;
    public:
        PiecewiseFunction();
        PiecewiseFunction(const T* pieces, const double* knots, const int& length);
        PiecewiseFunction(const PiecewiseFunction& other);
        void resize();
        void reset();
        void append(const T& piece, const double& knot);
        void append(const double& y, const double& w, const double& knot);
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
    capacity = 10;

    this->pieces = new T[capacity];

    this->knots = new double[capacity];

    this->length = 0;
}

template <class T>
PiecewiseFunction<T>::PiecewiseFunction(const T* pieces, const double* knots, const int& length){
    capacity = length * 2;

    this->pieces = new T[capacity];
    std::copy(pieces, pieces+length, this->pieces);

    this->knots = new double[capacity];
    std::copy(knots, knots+length, this->knots);

    this->length = length;
}

template <class T>
PiecewiseFunction<T>::PiecewiseFunction(const PiecewiseFunction& other){
    if(this != &other){
        length = other.length;
        capacity = other.capacity;

        pieces = new T[capacity];
        std::copy(other.pieces, other.pieces + length, pieces);

        knots = new double[capacity];
        std::copy(other.knots, other.knots + length, knots);
    }
}

template <class T>
void PiecewiseFunction<T>::resize(){
    T* tmppieces = pieces;
    pieces = new T[capacity];
    std::copy(tmppieces, tmppieces + length, pieces);
    delete[] tmppieces;

    double* tmpknots = knots;
    knots = new double[capacity];
    std::copy(tmpknots, tmpknots + length, knots);
    delete[] tmpknots;
}

template <class T>
void PiecewiseFunction<T>::reset(){
    length = 0;
}

template <class T>
void PiecewiseFunction<T>::append(const T& piece, const double& knot){
    if (length == capacity){
        capacity = capacity * 2;
        resize();
    }
    pieces[length] = piece;
    knots[length]  = knot;
    length++;
}

template <class T>
void PiecewiseFunction<T>::append(const double& y, const double& w, const double& knot){
    if (length == capacity){
        capacity = capacity * 2;
        resize();
    }
    pieces[length].set(y, w);
    knots[length]  = knot;
    length++;
}

template <class T>
void PiecewiseFunction<T>::append(const double& t, const double& knot){
    if (length == capacity){
        capacity = capacity * 2;
        resize();
    }
    pieces[length].set(t);
    knots[length]  = knot;
    length++;
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
    return length;
}

template <class T>
PiecewiseFunction<T>& PiecewiseFunction<T>::operator=(const PiecewiseFunction<T>& other){
    if(this == &other)
        return *this;
    if(capacity > other.length)
        length = other.length;
    else{
        delete[] pieces;
        delete[] knots;
        length = other.length;
        capacity = other.capacity;
        pieces = new T[capacity];
        knots = new double[capacity];
    }
    std::copy(other.pieces, other.pieces + length, pieces);
    std::copy(other.knots, other.knots + length, knots);
    return *this;
}

template <class T>
PiecewiseFunction<T>& PiecewiseFunction<T>::operator+=(const T& rhs){
    for(int i = 0; i < length; i++){
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
    if(length == 1){
        return pieces[0](x);
    } 
    int k = 0;
    for(int i = 1; i < length; i++){
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
    for(int i = 0; i < length; i++){
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

    for(int i = 0; i < length-1; i++){
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
    delete[] pieces;
    delete[] knots;
}


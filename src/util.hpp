#pragma once
#include "piecewisefunction.hpp"
#include "range.hpp"

void backtrace(const RangeList& ranges, const double* xprimes, double* solution);
void backtrace(const RangeList& ranges, const double* xprimes, int& k, int* start, int* end, double* value);
void backtrace(const RangeList* ranges, const int& N, const double* xprimes, double* solution);
void backtrace(const RangeList* ranges, const int& N, const double* xprimes, int* start, int* end, double* value);
int backtrace(const RangeList& ranges, const double& xprimes);

template <typename T>
void approximate(const int& n, const double* y, const double* lambda, const double* w, double* solution){
    PiecewiseFunction<T> f;
    PiecewiseFunction<T> g;

    f.append(y, w, T::domainninf, 0);
    f.append(T::rangeninf, T::domaininf);

    RangeList ranges(n-1);
    double* xprimes = new double[n];

    double yprime;

    T tmp;
    
    for(int i = 1; i < n; i++){
        f.max(xprimes[i-1], yprime);
        f.flood(yprime - lambda[i-1], g, ranges[i-1]);
        tmp.set(y, w, i);
        f = g + tmp;
    }
    f.max(xprimes[n-1], yprime);
    backtrace(ranges, xprimes, solution); 
    T::converter(solution, n);

    delete[] xprimes;
}

template <typename T>
void approximate(const int& n, const double* y, const double* lambda, const double* w, int& k, int* start, int* end, double* value){
    PiecewiseFunction<T> f;
    PiecewiseFunction<T> g;

    f.append(y, w, T::domainninf, 0);
    f.append(T::rangeninf, T::domaininf);

    RangeList ranges(n-1);
    double* xprimes = new double[n];

    double yprime;

    T tmp;
    
    for(int i = 1; i < n; i++){
        f.max(xprimes[i-1], yprime);
        f.flood(yprime - lambda[i-1], g, ranges[i-1]);
        tmp.set(y, w, i);
        f = g + tmp;
    }
    f.max(xprimes[n-1], yprime);
    backtrace(ranges, xprimes, k, start, end, value); 
    T::converter(value, k);

    delete[] xprimes;
}

template <typename T>
void approximate(const int& n, const double* y, const int& N, const double* w, double* solution){

    PiecewiseFunction<T>* fs = new PiecewiseFunction<T>[N];
    PiecewiseFunction<T>* gs = new PiecewiseFunction<T>[N];
    for(int i = 0; i < N; i++){
        fs[i].append(y, w, T::domainninf, 0);
        fs[i].append(T::rangeninf, T::domaininf);
    }

    double* xprimes = new double[n*N];

    RangeList* ranges = new RangeList[N];
    for(int i = 0; i < N; i++){
        ranges[i].resize(n-1);
    }

    double yprime;
    T tmp;

    for(int i = 1; i < n; i++){
        for(int j = N-1; j >= 0; j--){
            yprime = T::rangeninf;

            if( j <= i && j != 0){
                fs[j-1].max(xprimes[(i-1)*N+j],yprime);
            }

            if(j == i){
                gs[j].reset();
                gs[j].append(yprime, T::domainninf);
                gs[j].append(T::rangeninf, T::domaininf);
            }
            else{
                fs[j].flood(yprime, gs[j], ranges[j][i-1]);
            }
            tmp.set(y, w, i);
            fs[j] = gs[j] + tmp;
        }
    }
    fs[N-1].max(xprimes[n*N-1], yprime);
    backtrace(ranges, N, xprimes, solution);
    T::converter(solution, n);
    delete[] fs;
    delete[] gs;
    delete[] xprimes;
    delete[] ranges;
}

template <typename T>
void approximate(const int& n, const double* y, const int& N, const double* w, int* start, int* end, double* value){

    PiecewiseFunction<T>* fs = new PiecewiseFunction<T>[N];
    PiecewiseFunction<T>* gs = new PiecewiseFunction<T>[N];
    for(int i = 0; i < N; i++){
        fs[i].append(y, w, T::domainninf, 0);
        fs[i].append(T::rangeninf, T::domaininf);
    }

    double* xprimes = new double[n*N];

    RangeList* ranges = new RangeList[N];
    for(int i = 0; i < N; i++){
        ranges[i].resize(n-1);
    }

    double yprime;
    T tmp;

    for(int i = 1; i < n; i++){
        for(int j = N-1; j >= 0; j--){
            yprime = T::rangeninf;

            if( j <= i && j != 0){
                fs[j-1].max(xprimes[(i-1)*N+j],yprime);
            }

            if(j == i){
                gs[j].reset();
                gs[j].append(yprime, T::domainninf);
                gs[j].append(T::rangeninf, T::domaininf);
            }
            else{
                fs[j].flood(yprime, gs[j], ranges[j][i-1]);
            }
            tmp.set(y, w, i);
            fs[j] = gs[j] + tmp;
        }
    }
    fs[N-1].max(xprimes[n*N-1], yprime);
    backtrace(ranges, N, xprimes, start, end, value);
    T::converter(value, N);
    delete[] fs;
    delete[] gs;
    delete[] xprimes;
    delete[] ranges;
}

template <typename T>
int findbreakpoint(const int& n, const double* y, const double* w){

    PiecewiseFunction<T>* fs = new PiecewiseFunction<T>[2];
    PiecewiseFunction<T>* gs = new PiecewiseFunction<T>[2];
    for(int i = 0; i < 2; i++){
        fs[i].append(y, w, T::domainninf, 0);
        fs[i].append(T::rangeninf, T::domaininf);
    }

    RangeList ranges(n-1);

    double xprime, yprime;
    T tmp;

    for(int i = 1; i < n; i++){
        for(int j = 1; j >= 0; j--){
            yprime = T::rangeninf;

            if( j <= i && j != 0){
                fs[j-1].max(xprime,yprime);
            }

            if(j == i){
                gs[j].reset();
                gs[j].append(yprime, T::domainninf);
                gs[j].append(T::rangeninf, T::domaininf);
            }
            else if(yprime != T::rangeninf){
                fs[j].flood(yprime, gs[j], ranges[i-1]);
            }
            else{
                gs[j] = fs[j];
            }
            tmp.set(y, w, i);
            fs[j] = gs[j] + tmp;
        }
    }
    fs[1].max(xprime, yprime);
    int bp = backtrace(ranges, xprime);
    delete[] fs;
    delete[] gs;
    return bp;
}

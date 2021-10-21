#include "util.hpp"
#include "range.hpp"
#include <iterator>

void backtrace(const RangeList& ranges, const double* xprimes, double* solution){
    int n = ranges.len() + 1;
    double first, last;
    solution[n-1] = xprimes[n-1];
    for(int i = n-2; i >= 0; i--){
        solution[i] = solution[i+1];
        for(int j = 0; j < ranges.len(i); j++){
            ranges.index(i, j, first, last);
            if(first <= solution[i+1] && solution[i+1] <= last){
                solution[i] = xprimes[i];
                break;
            }
        }
    }
}

void backtrace(const RangeList& ranges, const double* xprimes, int& k, int* start, int* end, double* value){
    int n = ranges.len() + 1;
    double first, last;
    k = 0;
    value[k] = xprimes[n-1];
    end[k] = n;
    for(int i = n-2; i >= 0; i--){
        for(int j = 0; j < ranges.len(i); j++){
            ranges.index(i, j, first, last);
            if(first <= value[k] && value[k] <= last){
                start[k] = i+1;
                k++;
                value[k] = xprimes[i];
                end[k] = i+1;
                break;
            }
        }
    }
    start[k] = 0;
    k++;
}

void backtrace(const RangeList* ranges, const int& N, const double* xprimes, double* solution){
    int n = ranges[0].len() + 1;
    int T = N - 1;
    double first, last;
    solution[n-1] = xprimes[n*N-1];
    for(int i = n-2; i >= 0; i--){
        solution[i] = solution[i+1];
        for(int j = 0; j < ranges[T].len(i); j++){
            ranges[T].index(i, j, first, last);
            if(first <= solution[i+1] && solution[i+1] <= last){
                solution[i] = xprimes[i*N+T];
                T--;
                break;
            }
        }
    }
}

void backtrace(const RangeList* ranges, const int& N, const double* xprimes, int* start, int* end, double* value){
    int n = ranges[0].len() + 1;
    int T = N - 1;
    double first, last;
    int k = 0;
    value[k] = xprimes[n*N-1];
    end[k] = n;
    for(int i = n-2; i >= 0; i--){
        for(int j = 0; j < ranges[T].len(i); j++){
            ranges[T].index(i, j, first, last);
            if(first <= value[k] && value[k] <= last){
                start[k] = i+1;
                k++;
                value[k] = xprimes[i*N+T];
                end[k] = i+1;
                T--;
                break;
            }
        }
    }
    start[k] = 0;
    k++;
}

int backtrace(const RangeList& ranges, const double& xprime){
    int n = ranges.len() + 1;
    double first, last;
    for(int i = n-2; i >= 0; i--){
        for(int j = 0; j < ranges.len(i); j++){
            ranges.index(i, j, first, last);
            if(first <= xprime && xprime <= last){
                return i;
            }
        }
    }
    return -1;
}


#include "util.h"
#include <limits>

double inf = std::numeric_limits<double>::max();
double neginf = -inf;

bool addRange(double * ranges, const double & val, int & ii, const int & max){
    if (ii+1 > max){
        return false;
    }
    ranges[ii] = val;
    ii += 1;
    return true;
}

void backtrace(const double * bstar, const double * ranges, const int * range_inds, const int & range_inds_len, double * sol){
    sol[range_inds_len] = bstar[range_inds_len];
    int s, e;
    for (int i = range_inds_len-1; i > -1; i--){
        s = range_inds[i];
        e = range_inds[i+1];
        sol[i] = sol[i+1];
        for (int j = s; j < e; j++){
            if (ranges[2*j] <= sol[i+1] && ranges[2*j+1] >= sol[i+1]){
                sol[i] = bstar[i];
                break;
            }
        }
    }
}

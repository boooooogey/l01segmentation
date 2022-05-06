#include <Rcpp.h>
#include <math.h>
#include <limits>

using namespace Rcpp;

//[[Rcpp::export]]
List binVector(NumericVector y, int bin_size){
    int n = y.size();
    int binned_size = (int) ceil(n / (double) bin_size);
    NumericVector binned = NumericVector(binned_size);
    IntegerVector start = IntegerVector(binned_size);
    IntegerVector end = IntegerVector(binned_size);
    int counter = bin_size;
    double aggre = 0;
    int binned_index = 0;
    start[0] = 0;
    end[binned_size-1] = n;
    for(int i = 0; i < n; i++){
        aggre += y[i];
        counter--;
        if(counter == 0){
            binned[binned_index] = aggre / bin_size;
            end[binned_index] = i+1;
            binned_index++;
            start[binned_index] = i+1;
            counter = bin_size;
            aggre = 0.0;
        }
    }
    if(counter != bin_size){
        binned[binned_index] = aggre / (bin_size - counter);
    }
    return List::create(Named("values") = binned, Named("start") = start, Named("end") = end);
} 

#include <RcppArmadillo.h>
#include <algorithm>
#include <chrono>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::mat Xhat(const arma::vec & w){
    int n = w.size() + 1;
    arma::mat x(n-1,n);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n-1; j++){
            if (i <= j)
                x(j,i) = ((j+1)/(double)n-1)*w(j);
            else
                x(j,i) = (j+1)/(double)n*w(j);
        }
    }
    return x;
}

arma::mat dotX(const arma::mat & R, const arma::vec & w){
    int p = R.n_rows, n = R.n_cols+1;
    arma::mat out(p,n);
    out.col(0).fill(0);
    for(int i = 1; i < n; i ++) out.col(i) = out.col(i-1) + w(i-1) * R.col(i-1);
    return out;
}

void dotXT(const arma::mat & R, const arma::vec & w, arma::mat & out){
    int p = R.n_rows;
    int n = R.n_cols;
    if (w.size() != n-1)
        stop("Number of weights should be equal to n-1 (n is the number of columns of matrix R).");
    if (out.n_rows != p || out.n_cols != n-1)
        stop("One or both output dimensions are wrong.");

    arma::mat r(p,n);
    r.col(0) = R.col(0);
    for(int i=1; i < n; i++){
        r.col(i) = r.col(i-1) + R.col(i);
    }

    for(int i=0; i < n-1; i++){
        out.col(i) = w(i) * ((i+1) * r.col(n-1) / n - r.col(i));
    }
}

arma::mat dotXT(const arma::mat & R, const arma::vec & w){
    int p = R.n_rows;
    int n = R.n_cols;
    if (w.size() != n-1)
        stop("Number of weights should be equal to n-1 (n is the number of columns of matrix R).");

    arma::mat out(p,n-1);

    arma::mat r(p,n);
    r.col(0) = R.col(0);
    for(int i=1; i < n; i++){
        r.col(i) = r.col(i-1) + R.col(i);
    }

    for(int i=0; i < n-1; i++){
        out.col(i) = w(i) * ((i+1) * r.col(n-1) / n - r.col(i));
    }

    return out;
}

void XdotXT(const arma::uvec & ii1, const arma::uvec & ii2, const arma::vec & w, arma::mat & out){
    int r = ii1.size();
    int c = ii2.size();
    int m = std::max(r,c);
    if (w.size() < m)
        stop("Number of weights should not be less than the lenght of the index vectors.");
    
    if (out.n_rows != r || out.n_cols != c)
        stop("One or both output dimensions are wrong.");

    int n = w.size() + 1;
    for(int i = 0; i < c; i++){
        for(int j = 0; j < r; j++){
            out(j,i) = w(ii1(j)) * w(ii2(i)) * std::min(ii2(i)+1,ii1(j)+1) * (n - std::max(ii2(i)+1,ii1(j)+1)) / n;
        }
    }
}

void XdotXiT(const int & i, const arma::vec & w, arma::mat & out){
    
    if (out.n_rows != w.size() || out.n_cols != 1)
        stop("One or both output dimensions are wrong.");

    int n = w.size() + 1;
    for(int j = 0; j < w.size(); j++){
        out(j,0) = w(j) * w(i) * std::min(i+1,j+1) * (n - std::max(i+1,j+1)) / n;
    }
}

double rowXhatnorm(const double & i, const double & w, const double & n){
    return i * (n - i) * w * w / n;
}

void dotXXT(const arma::mat & R, const arma::vec & w, arma::mat & out){
    int p = R.n_rows;
    int n = R.n_cols;

    if (w.size() != n)
        stop("Number of weights should be equal to n (n is the number of columns of matrix R).");

    if (out.n_rows != p || out.n_cols != n)
        stop("One or both output dimensions are wrong.");

    for(int i = 0; i < n; i++){
        out.col(i) = w(i) * R.col(i);
    }
    arma::vec S = out*arma::linspace<arma::vec>(1, n, n)/(n+1);
    arma::mat T(p,n);
    T.col(n-1) = out.col(n-1);
    for(int i = n-2; i >= 0; i--){
        T.col(i) = out.col(i) + T.col(i+1);
    }
    out.col(0) = S - T.col(0);
    for(int i = 1; i < n; i++){
        out.col(i) = out.col(i-1) + S - T.col(i);
    }
    out.each_row() %= -w.t();
}

void update(const arma::mat & S, const double & lambda, const double & gammai, arma::mat & out, const int & i){
    out.col(i) = std::max(1-lambda/norm(S.col(i)),0.0) * S.col(i) / gammai;
}

void removeb0(const arma::mat & B, arma::uvec & A, int & Alength){
    arma::uvec ii = find(sum(B % B,0) != 0);
    Alength = ii.size();
    A.head_rows(Alength) = A.rows(ii);
}

bool addbn0(arma::uvec & A, int & Alength, const arma::mat & S, const double & lambda){
    arma::rowvec norms = sum(S % S,0);
    norms.elem(A.head_rows(Alength)).fill(0);
    arma::uword i = norms.index_max();
    if (norms(i) < (lambda*lambda))
        return true;
    A(Alength) = i;
    Alength = Alength+1;
    return false;
}

//[[Rcpp::export]]
List blockcoordinatedescent(const arma::mat & Yhat, const double & lambda, const arma::vec & w, const double mintimer = 5, const double tol = 1e-6){
    int p = Yhat.n_rows; int n = Yhat.n_cols;

    arma::mat B(p,n-1,arma::fill::zeros);
    arma::vec maskb(p,arma::fill::zeros);
    arma::vec xxit(n-1,arma::fill::zeros);
    arma::mat C(p,n-1,arma::fill::zeros); 
    arma::mat S(p,n-1,arma::fill::zeros);
    arma::mat BXXT(p,n-1,arma::fill::zeros);
    arma::uvec A(n-1,arma::fill::zeros);
    int Ait = 0;
    int Alength = 0;
    int i = 0;
    bool converged;

    double exctime = 60*mintimer;
    bool timesup = false;
    auto start = std::chrono::system_clock::now();
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed;

    dotXT(Yhat,w,C);
    while(true){
        Ait = 0;
        converged = true;
        // Run on the active set until convergence
        while(Alength != 0){
            i = A(Ait);
            maskb = B.col(i);
            B.col(i).fill(0);
            XdotXiT(i, w, xxit);
            S.col(i) = C.col(i) - B.cols(A.head_rows(Alength)) * xxit.rows(A.head_rows(Alength)); //just select non 0 from B
            update(S,lambda,rowXhatnorm(i+1,w(i),n),B,i);
            //Rcout << "i = " << i << " diff = " << norm(maskb - B.col(i)) << std::endl;
            if(norm(maskb - B.col(i)) > tol) converged = false;
            Ait++;
            if(Ait == Alength){
                if(converged){
                    break;
                }
                else {
                    converged = true;
                    Ait = 0;
                }
            }
        }
        // Converged
        // Remove the groups that are inactive
        if (Alength != 0){
            removeb0(B.cols(A.head_rows(Alength)), A, Alength);
        }
        //Rcout << "(removed) A =" << std::endl;
        //Rcout <<  A.head_rows(Alength) << std::endl;
        // Check to see if KKT satisfied
        dotXXT(B, w, BXXT); //also B has zero here
        S = C - BXXT;
        converged = addbn0(A,Alength,S,lambda);
        //Rcout << "(added) A =" << std::endl;
        //Rcout <<  A.head_rows(Alength) << std::endl;
        if (converged) break;
        end = std::chrono::system_clock::now();
        elapsed = end-start;
        if (exctime <= elapsed.count()){
            warning("Execution time limit has been reached. Stopping execution before convergence.");
            break;
        }
    }
    arma::uvec Asorted = arma::sort(A.head_rows(Alength));
    arma::mat Bsub = B.cols(Asorted);
    return List::create(Named("B") = Bsub, Named("ii") = IntegerVector(Asorted.begin(),Asorted.begin()+Alength)+1,Named("converged") = converged);
}


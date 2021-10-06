#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double g(double theta, int x1, int x2){
    double res = 10*theta*theta - (30+x2)*theta+(x1+x2) ;
    return res ;
}

// [[Rcpp::export]]
double gprime(double theta, int x1, int x2){
    double res = 20*theta-(30+x2) ;
    return res ;
}


// [[Rcpp::export]]
double mix_mle(int x1, int x2, double tol=1e-6){
    double theta_current=(x1+2*x2)/40 ;
    double theta_new = 0 ;
    while(true){
        theta_new = theta_current  - g(theta_current, x1,x2) / gprime(theta_current, x1,x2) ;
        if(abs(g(theta_new, x1, x2)) < tol){
            break ;
        }
        if(theta_new > 1 || theta_new < 0 ){
            Rcpp::Rcout << "Theta is out of parameter space" << std::endl ;
            break ;
        }
        theta_current = theta_new ;
    }
    return theta_new ;
}

// [[Rcpp::export]]
double cauchy_lprime(double theta, arma::vec x){
    int n = x.n_elem ;
    double sum = 0 ;
    for(int i = 0; i<n ; ++i){
        sum = sum + (2*(x(i)-theta)) / (1+(x(i)-theta)*(x(i)-theta) ) ;
    }
    return sum ;
}

// [[Rcpp::export]]
double cauchy_ldprime(double theta, arma::vec x){
    int n = x.n_elem ;
    double sum = 0 ;
    for(int i = 0; i<n ; ++i){
        sum = sum +  (2*( (x(i)-theta)*(x(i)-theta)-1) ) / ((1+(x(i)-theta)*(x(i)-theta))*(1+(x(i)-theta)*(x(i)-theta))) ;
    }
    return sum ;
}

// [[Rcpp::export]]
double cauchy_mle(arma::vec x, std::string method="Newton", double tol=1e-6){
    if(method=="Newton"){
        double theta_current = arma::median(x) ;
        double theta_new = 0 ;
        int iter=0 ;
        while(true){
            theta_new = theta_current  - cauchy_lprime(theta_current, x) / cauchy_ldprime(theta_current, x) ;
            if(abs(cauchy_lprime(theta_new, x)) < tol ){
                break ;
            }
            theta_current = theta_new ;
            iter +=1 ;
            if( iter > 1e+3){
                Rcpp::Rcout << "Too slow for iteration to convege" << std::endl ;
                break ;
            }
        }
        return theta_new ;
    }
    else if(method=="Dekker"){
        double a = - (abs(arma::median(x))+3) ;
        double theta_past = a ;
        double theta_current = abs(arma::median(x))+3 ;
        double theta_new =0 ;
        double proposal ;
        double middle ;
        double ratio ;
        int iter1=0 ; 
        while(cauchy_lprime(theta_past, x)*cauchy_lprime(theta_current, x) > 0 ){
            Rcpp::Rcout << "The initial values do not satisfy starting condition of bisection method. Shifting it a little..." << std::endl ;
            theta_past = theta_past + arma::randn() ;
            theta_current = theta_current + arma::randn() ;
            a = theta_past ;
            iter1 +=1 ;
            if( iter1 > 1e+3){
                Rcpp::Rcout << "Initial value shifting is too complicated" << std::endl ;
                break ;
            }
        }
        int iter2 =0 ;
        while(true){
            ratio = (cauchy_lprime(theta_current, x)-cauchy_lprime(theta_past, x)) / (theta_current-theta_past) ;
            middle = (a + theta_current) / 2 ;
            proposal = cauchy_lprime(theta_current, x)-cauchy_lprime(theta_past, x) !=0 ? theta_current-cauchy_lprime(theta_current, x) /ratio : middle ;
            if(middle < theta_current){
                theta_new = proposal<= theta_current && proposal >=middle ? proposal : middle ;
            }
            else{
                theta_new = proposal<= middle && proposal >= theta_current ? proposal : middle ;
            }
            if(abs(cauchy_lprime(theta_new,x))<tol ){
                break ;
            }

            if(cauchy_lprime(theta_new, x)*cauchy_lprime(a, x) >0 ){
                a = theta_current ;
            }
            if(abs(cauchy_lprime(theta_new,x)) > abs(cauchy_lprime(a, x)) ) {
                std::swap(a, theta_new) ;
            }

            theta_past = theta_current ;
            theta_current = theta_new ;

            iter2 += 1 ;
            if( iter2 > 1e+3){
                Rcpp::Rcout << "Too slow for iteration to converge" << std::endl ;
                break ;
            }
        }
        return theta_new ;
    }
    else{
        Rcpp::Rcout << "Default method is Newton. The other choice is Dekker" << std::endl ;
        return 0 ;
    }
}

// [[Rcpp::export]]
arma::mat MLEwithCRLB(int N=10000, int n=100, std::string method="Newton", double tol=1e-6){
    arma::mat X(N, 2) ;
    for(int i=0 ; i<N ; ++i){
        arma::vec x = Rcpp::rcauchy( n, 0.0, 1.0) ;
        double mle = cauchy_mle(x, method, tol) ;
        if( abs(mle) > 1e+03){
            Rcpp::Rcout << "Abnormal value of MLE is yielded" << std::endl ;
            while(abs(mle) > 1e+03){
                x = Rcpp::rcauchy( n, 0.0, 1.0) ;
                mle = cauchy_mle(x, method, tol) ;
                if( abs(mle) > 1e+03){
                    Rcpp::Rcout << "Abnormal value of MLE is yielded" << std::endl ;
                }
            }
        }
        double obs_info_bound = - 1 / cauchy_ldprime(mle, x) ;
        X.row(i)(0) = mle ;
        X.row(i)(1) = obs_info_bound ;
    }
    return X ;
}



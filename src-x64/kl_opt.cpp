#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(Rcpp)]]

using namespace arma;


// [[Rcpp::export]]
Rcpp::NumericVector calc_vec(Rcpp::NumericVector v1, Rcpp::NumericVector v2){
  Rcpp::NumericVector out;

  //out = sum(v1 * ((v1) - (v2)));
  out = (v1 - v2);
  return(out);
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix KL(const Rcpp::NumericMatrix & x, const Rcpp::NumericMatrix & logx) {
  unsigned int outcols = x.ncol(), i = 0, j = 0;
  double d1, d2;
  Rcpp::NumericMatrix out(outcols, outcols);

  for (j = 0; j < outcols - 1; j++) {
    Rcpp::NumericVector v1 = x.column(j);
    Rcpp::NumericVector lv1 = logx.column(j);
    for (i = 1; i < outcols; i++) {
      Rcpp::NumericVector v2 = x.column(i);
      Rcpp::NumericVector lv2 = logx.column(i);

      //Rcpp::Rcout << "The length is " << v2 << std::endl;
      d1 = sum(v1 * (lv1 - lv2));//calc_vec(v1, v2);//sum(v1 * (log(v1)-log(v2)));
      d2 = sum(v2 * (lv2 - lv1));//calc_vec(v2, v1);//sum(v2 * (log(v2)-log(v1)));
      out(j, i) = d1;
      out(i, j) = d2;
    }
  }

  return out;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix ED2(const Rcpp::NumericMatrix & x) {
  unsigned int outcols = x.ncol(), i = 0, j = 0;
  double d;
  Rcpp::NumericMatrix out(outcols, outcols);

  for (j = 0; j < outcols - 1; j++) {
    Rcpp::NumericVector v1 = x.column(j);
    for (i = j + 1; i < outcols; i++) {
      d = sqrt(sum(pow(v1 - x.column(i), 2.0)));
      out(i, j) = d;
      out(j, i) = d;
    }
  }

  return out;
}

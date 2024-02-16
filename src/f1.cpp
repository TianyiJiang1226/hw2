#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double firstl(NumericVector x, double theta) {
  int n = x.size();
  double s = 0;
  for (int i = 0; i < n; i++) {
    s += (x[i] - theta)/(pow((x[i]-theta),2)+1);
  }
  s = s*2;
  return s;
}

// [[Rcpp::export]]
double secondl(NumericVector x, double theta) {
  int n = x.size();
  double s = 0;
  for (int i = 0; i < n; i++) {
    s += (pow((x[i] - theta),2)-1)/(pow((pow((x[i]-theta),2)+1),2));
  }
  s = s*2;
  return s;
}

// [[Rcpp::export]]
NumericVector bisection(NumericVector x, double bound){
  double theta0 = median(x);
  double theta1 = min(x);
  double theta2 = max(x);
  double f0 = firstl(x,theta0);
  double f1 = firstl(x,theta1);
  int n = 1;
  while (abs(f1) > bound){
    n += 1;
    if(f0 * f1 < 0){
      theta2 = theta1;
      theta1 = (theta0+theta2)/2;
    }
    else{
      theta0 = theta1;
      theta1 = (theta0+theta2)/2;
    }
    f0 = firstl(x,theta0);
    f1 = firstl(x,theta1);
  }
  NumericVector res = NumericVector::create(_["theta1"] = theta1, _["n"] = n);
  return res;
}

// [[Rcpp::export]]
NumericVector newton(NumericVector x, double bound){
  double theta = median(x);
  double f1 = firstl(x,theta);
  double f2 = secondl(x,theta);
  int n = 1;
  while (abs(f1) > bound){
    n += 1;
    theta = theta - f1/f2;
    f1 = firstl(x,theta);
    f2 = secondl(x,theta);
  }
  NumericVector res = NumericVector::create(_["theta"] = theta, _["n"] = n);
  return res;
}

// [[Rcpp::export]]
double fish(NumericVector x, double theta) {
  int n = x.size();
  double s = 0;
  for (int i = 0; i < n; i++) {
    s += (x[i] - theta)/(pow((x[i]-theta),2)+1)*(x[i] - theta)/(pow((x[i]-theta),2)+1);
  }
  s = s*4;
  return s;
}

// [[Rcpp::export]]
NumericVector fisher(NumericVector x, double bound){
  double theta = median(x);
  double f1 = firstl(x,theta);
  int n = 1;
  while (abs(f1) > bound){
    n += 1;
    double fishn = fish(x,theta);
    theta = theta + f1/fishn;
    f1 = firstl(x,theta);
  }
  NumericVector res = NumericVector::create(_["theta"] = theta, _["n"] = n);
  return res;
}

// [[Rcpp::export]]
NumericVector secant(NumericVector x, double bound){
  double theta0 = median(x);
  double theta1 = theta0 + 1;
  double middletheta = theta1;
  double f0 = firstl(x,theta0);
  double f1 = firstl(x,theta1);
  int n = 1;
  while (abs(f1) > bound){
    n += 1;
    middletheta = theta1;
    theta1 = theta1 - f1*(theta1 - theta0)/(f1 - f0);
    theta0 = middletheta;
    f0 = firstl(x,theta0);
    f1 = firstl(x,theta1);
  }
  NumericVector res = NumericVector::create(_["theta1"] = theta1, _["n"] = n);
  return res;
}

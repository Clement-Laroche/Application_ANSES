#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]

double Cquantile(NumericVector x, double q) 
{
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y[x.size()*(q - 0.000000001)];
}

// [[Rcpp::export]]

// code la fonction de coup

NumericVector costfcpp(DataFrame data,int tstart,int tstop,string init,double q_init,double prec,int Nmax,double lmax)
{
  if((tstart < 0) | (tstart >= data.nrow()))
  {
    stop("Error: tstart must be in data range.");
  }
  if((tstop < 0) | (tstop >= data.nrow()))
  {
    stop("Error: tstop must be in data range.");
  }
  NumericVector y = data[0];
  NumericVector d = data[1];
  y = y[seq(tstart,tstop)];
  d = d[seq(tstart,tstop)];
  int n = y.size();
  NumericVector out(2);
  NumericVector a = y*d;
  a = a[a!=0];
  a = sort_unique(a);
  NumericVector na(a.size());
  for(int i = 0; i < a.size();i++)
  {
    na[i] = sum(y == rep(a[i],n));
  }
  NumericVector ncens;
  int ncensl;
  if(sum(na) == n)
  {
    ncens = rep(0.0,n);
    ncensl = 0;
  }
  else
  {
    ncens = y[d == 0];
    ncensl = ncens.size() ; 
  }
  bool crit = TRUE ;
  int count = 1 ;
  double x;
  if(init == "quantile")
  {
    x = -Cquantile(y,q_init)/log(1-q_init);
  }
  else
  {
    x = 1/mean(y);
  }
  while(crit == TRUE)
  {
    double x1;
    if(na.size() == 0)
    {
      double num = ncensl/x-sum(ncens);
      double den = -ncensl/(x*x);
      x1 = x - num/den;
      // Rcout << num;
    }else if(ncensl == 0)
    {
      double num = sum(na*a*exp(-x*a)/(1-exp(-x*a)));
      double den = -sum(na*a*a*exp(-x*a)/((1-exp(-x*a))*(1-exp(-x*a))));
      x1 = x - num/den;
    }else
    {
      double num = sum(na*a*exp(-x*a)/(1-exp(-x*a)))+ncensl/x-sum(ncens);
      double den = -sum(na*a*a*exp(-x*a)/((1-exp(-x*a))*(1-exp(-x*a))))-ncensl/(x*x);
      x1 = x - num/den;
    }
    if(x1 < 0)
    {
      break;
    }
    if(abs(x1-x) < prec)
    {
      crit = FALSE ;
    }
    x = x1;
    count += 1;
    if(count >= Nmax)
    {
      crit = FALSE ;
    }
    if(x >= lmax)
    {
      crit = FALSE;
    }
  }
  out[0] = x;
  double lv;
  if(na.size() == 0)
  {
    lv = -ncensl*log(x)+x*sum(ncens);
  }
  else if(ncensl == 0)
  {
    lv = sum(-na*log(1-exp(-x*a)));
  }else
  {
    lv = sum(-na*log(1-exp(-x*a)))-ncensl*log(x)+x*sum(ncens);
  }
  out[1] = lv;
  return out;
}

// [[Rcpp::export]]

// code pour la m?thode de recherche optimale

List optcpp(DataFrame sumstat,int Kmax,string init,double q_init,double prec,int Nmax,double lmax)
{
  int n = sumstat.nrow();
  NumericMatrix C(n,n);
  for(int u = 0; u < n-1;u++)
  {
    for(int v = u+1; v < n; v++)
    {
      C(u,v) = costfcpp(sumstat,u,v,init,q_init,prec,Nmax,lmax)[1];
    }
  }
  List L ;
  L.push_back(C);
  if(Kmax != 2)
  {
    for(int i = 0; i < Kmax-2; i++)
    {
      NumericMatrix r = L[i];
      NumericMatrix fs = L[0];
      NumericMatrix C(n,n);
      for(int u = 0; u < n-(i+2);u++)
      {
        for(int v = u+(i+2); v < n; v++)
        {
          C(u,v) = min(as_vector(r(seq(u,u),seq(u+i+1,v-1)))  + as_vector(fs(seq(u+i+2,v),seq(v,v)))); 
        }
      }
      L.push_back(C);
    }
  }
  IntegerVector res = rep(0,Kmax-1);
  List RES = List(res);
  for(int i = 0; i < res.size();i++)
  {
    
    RES[i] = rep(0,i+2);
    int k = i+1;
    as<IntegerVector>(RES[i])[k] = n-1;
    while(k > 0)
    {
      NumericMatrix r = L[k-1];
      NumericMatrix fs = L[0];
      int s =  as<IntegerVector>(RES[i])[k];
      int star = which_min(as_vector(r(seq(0,0),seq(k-1,s-1)))+as_vector(fs(seq(k,s),seq(s,s))));
      as<IntegerVector>(RES[i])[k-1] = seq(k-1,s-1)[star];
      k = k-1;
    }
    
    RES[i] = as<IntegerVector>(RES[i])[Range(0,i)];
  }
  return RES;
}

// [[Rcpp::export]]

// code pour la m?thode de recherche PELT

IntegerVector peltcpp(DataFrame sumstat,int minseg,double pen,string init,double q_init,double prec,int Nmax,double lmax)
{
  int n = sumstat.nrow();
  NumericVector lastchangelike = rep(0.0,n);
  IntegerVector lastchangecpts = rep(0,n);
  // IntegerVector numchangecpts = rep(0,n);
  NumericVector data = sumstat[0];
  lastchangelike[seq(0,minseg-1)] = rep(-pen,minseg);
  NumericVector info = rep(0.0,2);
  for(int i = minseg-1;i < 2*minseg;++i)
  {
    info = costfcpp(sumstat,0,i,init,q_init,prec,Nmax,lmax);
    lastchangelike[i] = info[1];
  }
  int nchecklist = 2;
  IntegerVector checklist(nchecklist);
  checklist[0] = 0;
  checklist[1] = minseg-1;
  for(int tstar = 2*minseg-1 ;tstar < n;++tstar)
  {
    IntegerVector clinter = checklist[checklist <= (tstar-minseg)];
    NumericVector costinter;
    NumericVector restarting;
    for(int i = 0; i < clinter.size();++i)
    {
      if(clinter[i] == 0)
      {
        info = costfcpp(sumstat,clinter[i],tstar,init,q_init,prec,Nmax,lmax);
      }else
      {
        info = costfcpp(sumstat,clinter[i]+1,tstar,init,q_init,prec,Nmax,lmax);
      }
      costinter.push_back(info[1]);
      restarting.push_back(info[0]);
    }
    NumericVector templike = lastchangelike[clinter];
    templike += costinter+pen;
    lastchangelike[tstar] = min(templike);
    lastchangecpts[tstar] = checklist[which_min(templike)];
    // // numchangecpts[tstar] = numchangecpts[lastchangecpts[tstar]]+1;
    LogicalVector coord = rep(true,nchecklist);
    IntegerVector n_c;
    for(int i = 0;i < templike.size();++i)
    {
      if((templike[i]-pen) >= min(templike))
      {
        n_c.push_back(i);
      }
    }
    coord[n_c] = false;
    checklist = checklist[coord];
    checklist.push_back(tstar);
    nchecklist = checklist.size();
  }
  int i = n-1;
  IntegerVector cptsout;
  while(i>1)
  {
    cptsout.push_back(lastchangecpts[i]);
    i = lastchangecpts[i]-1;
  }
  return cptsout;
}





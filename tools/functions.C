#ifndef __RUN1BANA_TOOLS_FUNCTIONS__
#define __RUN1BANA_TOOLS_FUNCTIONS__

//------------------------------------------------------------------------------------------------------------
// Landau core with power-law tails
double landau_crystal_ball(double x, double mean, double a, double b, double alpha1, double alpha2, double n1, double n2) {
  // See: https://github.com/pavel1murat/murat/blob/main/scripts/fit_cb4.C

  // Requirements for normalization are N1 and N2 are > 1:
  if(n1 < 1. || n2 < 1.) return 0.;

  // Evaluate the function
  const double dx = x-mean;
  double val = 0.;
  if (dx < -alpha1) { // Low tail
    const double B1 = -alpha1+n1/(a*b*(1-exp(-b*alpha1)));
    const double A1 = exp(a*(-b*alpha1-exp(-b*alpha1)))*pow(B1+alpha1,n1);
    val  = A1/pow(B1-dx,n1);
  } else if (dx < alpha2) { // Landau core
    val = exp(a*(b*dx-exp(b*dx)));
  } else { // High tail
    const double B2 = -alpha2-n2/(a*b*(1. - exp(b*alpha2)));
    const double A2 = exp(a*(b*alpha2-exp(b*alpha2)))*pow(B2+alpha2,n2);
    val  = A2/pow(B2+dx,n2);
  }

  return val;
}

//------------------------------------------------------------------------------------------------------------
double landau_crystal_ball_func(double* X, double* P) {
  return P[0]*landau_crystal_ball(X[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7]);
}

//------------------------------------------------------------------------------------------------------------
TF1* landau_crystal_ball_tf1() {
  TF1* f = new TF1("landau_crystal_ball", landau_crystal_ball_func, -10., 10., 8);
  f->SetParNames("Norm", "#mu", "a", "b", "#alpha_{1}", "#alpha_{2}", "n_{1}", "n_{2}");
  return f;
}

//------------------------------------------------------------------------------------------------------------
// Gaussian core with power-law tails
double double_crystal_ball(double x, double sigma, double mean, double alpha1, double alpha2, double n1, double n2) {
  //integral from a to b of exp(-x^2/2/s^2) = sqrt(pi/2)*s*erf(x/sqrt(2)/s) eval from a to b
  //bounds are x = -alpha1*sigma and alpha2*sigma
  const double gausVal = sqrt(M_PI/2.)*sigma*(erf(alpha2/sqrt(2.)) - erf(-alpha1/sqrt(2.)));

  //integral (n/a - a - (x-m)/s)^-n = (-s)/(1-n)*(n/a-a-(x-m)/s)^(1-n)
  //from x = 0 to (x-m)/s = -a = (-s)/(1-n)*{(n/a)^(1-n) - (n/a - a + m/s)^(1-n)}
  // double lowTailVal  = abs(pow((n1/alpha1),1.-n1)) - abs(pow((n1/alpha1-alpha1 + mean/sigma),1.-n1));
  double lowTailVal  = pow((n1/alpha1),-n1+1.);//assume --> -inf //- pow((n1/alpha1-alpha1),-n1+1.);
  lowTailVal        *= pow(n1/(alpha1),n1)*exp(-alpha1*alpha1/2.)/(n1-1.)*sigma;

  //integral (n/a - a + (x-m)/s)^-n = (s)/(1-n)*(n/a-a+(x-m)/s)^(1-n)
  //from x = a to (x-m)/s = inf = (s)/(1-n)*{(n/a - a + inf)^(1-n) - (n/a - a + a)^(1-n)} = (-s)/(1-n)*(n/a)^(1-n)
  double highTailVal = -abs(pow((n2/alpha2),1.-n2));
  highTailVal       *= pow(n2/(alpha2),n2)*exp(-alpha2*alpha2/2.)/(1.-n2)*sigma;

  //can be additionally changed to include an (energy dependent)000 efficiency
  const double norm = 1./(lowTailVal + gausVal + highTailVal); //chosen to normalize the pdf
  //  norm = norm / 1.281; // NEED TO FIX, NOT NORMALIZED PROPERLY

  double res = 0.;
  const double arg = (x-mean)/sigma;

  if(arg < -1.*alpha1) { //more than alpha1 sigmas below mean
    res = pow(n1/(alpha1),n1)*exp(-alpha1*alpha1/2.)*pow((n1/alpha1-alpha1-arg),-n1);
  } else if(arg > alpha2) { //more than alpha2 sigmas above mean
    res = pow(n2/alpha2,n2)*exp(-alpha2*alpha2/2.)*pow((n2/alpha2-alpha2+arg),-n2);
  } else { //else, gaussian
    res = exp(-arg*arg/2.);
  }
  res *= norm;

  return res;
}

//------------------------------------------------------------------------------------------------------------
double double_crystal_ball_func(double* X, double* P) {
  return P[0]*double_crystal_ball(X[0], P[1], P[2], P[3], P[4], P[5], P[6]);
}

//------------------------------------------------------------------------------------------------------------
TF1* double_crystal_ball_tf1() {
  TF1* f = new TF1("double_crystal_ball", double_crystal_ball_func, -10., 10., 7);
  f->SetParNames("Norm", "#sigma", "#mu", "#alpha_{1}", "#alpha_{2}", "n_{1}", "n_{2}");
  return f;
}

#endif

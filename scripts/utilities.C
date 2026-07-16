#ifndef __RUN1BANA_SCRIPTS_UTILS__
#define __RUN1BANA_SCRIPTS_UTILS__

double lognormal_pdf(const double x, const double mu, const double SDF) {

  if(x   <= 0.) return 0.;
  if(SDF <= 0.) return 0.;
  if(mu  <= 0.) return 0.;

  const double sigma = std::sqrt(-std::log(SDF));
  const double mu0 = std::log(mu) - 0.5*sigma*sigma;

  const double p = 1./(x*sigma*std::sqrt(2.*M_PI))
    * std::exp(-std::pow(std::log(x) - mu0, 2)
               / (2.*sigma*sigma));

  return p;
}

double lognormal_cdf(const double x, const double mu, const double SDF) {

  if(x   <= 0.) return 0.;
  if(SDF <= 0.) return 0.;
  if(mu  <= 0.) return 0.;

  const double sigma = std::sqrt(-std::log(SDF));
  const double mu0 = std::log(mu) - 0.5*sigma*sigma;

  const double z = (std::log(x) - mu0) / (sigma * std::sqrt(2.));
  return 0.5 * std::erfc(-z);
}

double xlognormal_pdf(const double x, const double mu, const double SDF) {

  if(x   <= 0.) return 0.;
  if(SDF <= 0.) return 0.;
  if(mu  <= 0.) return 0.;

  const double sigma = std::sqrt(-std::log(SDF));
  const double mu0 = std::log(mu) - 0.5*sigma*sigma;

  double shifted_mu = mu0 + sigma * sigma;
  double log_x_minus_mu = std::log(x) - shifted_mu;

  double exponent = -0.5 * (log_x_minus_mu * log_x_minus_mu) / (sigma * sigma);
  double normalization = 1.0 / (x*sigma * std::numbers::sqrt2_v<double> * std::sqrt(M_PI));

  const double p = normalization * std::exp(exponent);
  return p;
}

double xlognormal_cdf(const double x, const double mu, const double SDF) {

  if(x   <= 0.            ) return 0.;
  if(SDF <= 0. || SDF > 1.) return 0.;
  if(mu  <= 0.            ) return 0.;

  const double sigma = std::sqrt(-std::log(SDF));
  const double mu0 = std::log(mu) - 0.5*sigma*sigma;

  // Shifted normal variable for length-biased distribution
  double z = (std::log(x) - (mu0 + sigma * sigma)) / (sigma * std::numbers::sqrt2_v<double>);

  // Standard normal CDF evaluated via complementary error function
  const double p = 0.5 * std::erfc(-z);

  return p;
}

#endif

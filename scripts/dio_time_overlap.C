// Rough evaluation of two DIOs overlapping in time

void dio_time_overlap(long ngen = 1e7) {

  TRandom3 rnd(90);
  constexpr double lifetime = 864.2;
  long noverlap = 0;
  for(long igen = 0; igen < ngen; ++igen) {
    const double t_1 = fmod(rnd.Exp(lifetime) + rnd.Gaus(200., 30.), 1695.);
    const double t_2 = fmod(rnd.Exp(lifetime) + rnd.Gaus(200., 30.), 1695.);
    if(std::abs(t_1 - t_2) < 25.) ++noverlap;
  }

  printf("Rate overlap = %.3g +- %.3g\n", noverlap*1./ngen, std::sqrt(noverlap)/ngen);
}

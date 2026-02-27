#include <cmath>
#include <cstdio>

struct s_rbers {
  float rber[7];
};

double f_erf(double x) {
  double y = 1.0 / (1.0 + 0.3275911 * x);
  return 1 - (((((
        + 1.061405429  * y 
        - 1.453152027) * y 
        + 1.421413741) * y 
        - 0.284496736) * y 
        + 0.254829592) * y) 
        * exp(-x * x);
}

// Returns the probability of x, given the distribution described by mu and
// sigma.
double f_pdf(double x, double mu, double sigma) {
  // Constants
  static const double pi = 3.14159265;
  return exp(-1 * (x - mu) * (x - mu) / (2 * sigma * sigma)) /
         (sigma * sqrt(2 * pi));
}
// Returns the probability of [-inf,x] of a gaussian distribution
double f_cdf(double x, double mu, double sigma) {
  return 0.5 * (1 + f_erf((x - mu) / (sigma * sqrt(2.))));
}

double f_find_sigma(double x, double mu, double target_cdf) {
  int i = 0;
  double sigma = 1.0;
  double cdf = 0.0;
  double error = 1;
  while ((fabs(error) >= 0.00000001) && (i < 10000)) {
    cdf = f_cdf(x, mu, sigma);
    error = target_cdf - cdf;
    sigma *= (1.0 - (1.0 * error));
    i++;
  }
  return sigma;
}

double f_find_x(double mu, double sigma, double target_cdf) {
  int i = 0;
  double x = 0.0;
  double cdf = 0.0;
  double error = 1;
  while ((fabs(error) >= 0.00000001) && (i < 10000)) {
    cdf = f_cdf(x, mu, sigma);
    error = target_cdf - cdf;
    x += (0.1 * error);
    i++;
  }
  return x;
}

s_rbers f_find_error_probabilities(float rber) {
  float vt[7];
  s_rbers rbers;
  float target_cdf;
  float target_rber;
  float mu, sigma;

  rbers.rber[0] = rber;
  vt[0] = 0.0;
  mu = -1.0;

  // for 3 strobes
  target_rber = 0.0005;
  target_cdf = 1.0 - target_rber;
  sigma = f_find_sigma(vt[0], mu, target_cdf);
  target_rber *= 0.10;
  target_cdf = 1.0 - target_rber;
  vt[2] = f_find_x(mu, sigma, target_cdf);
  vt[4] = vt[2];
  vt[6] = vt[4];

  // general for 3, 5, 7 strobes
  vt[1] = -vt[2];
  vt[3] = -vt[4];
  vt[5] = -vt[6];

  target_cdf = 1.0 - rber;
  sigma = f_find_sigma(vt[0], mu, target_cdf);
  rbers.rber[0] = rber;
  rbers.rber[1] = 1.0 - f_cdf(vt[1], mu, sigma);
  rbers.rber[2] = 1.0 - f_cdf(vt[2], mu, sigma);
  rbers.rber[3] = 1.0 - f_cdf(vt[3], mu, sigma);
  rbers.rber[4] = 1.0 - f_cdf(vt[4], mu, sigma);
  rbers.rber[5] = 1.0 - f_cdf(vt[5], mu, sigma);
  rbers.rber[6] = 1.0 - f_cdf(vt[6], mu, sigma);

  // find LLR corresponding to low confidence region
  float P_low_confidence_symbL =
      f_cdf(0, mu, sigma) - f_cdf(-vt[2], mu, sigma); // P(r_low | s==-1)
  float P_low_confidence_symbR =
      f_cdf(vt[2], mu, sigma) - f_cdf(0, mu, sigma); // P(r_low | s==+1)
  float likelihood_ratio_low_confidence =
      P_low_confidence_symbL / P_low_confidence_symbR;
  float llr_low_confidence = log(likelihood_ratio_low_confidence);

  // find LLR corresponding to high confidence region
  float P_high_confidence_symbL = f_cdf(-vt[2], mu, sigma); // P(r_high | s==-1)
  float P_high_confidence_symbR =
      1.0 - f_cdf(vt[2], mu, sigma); // P(r_high | s==+1)
  float likelihood_ratio_high_confidence =
      P_high_confidence_symbL / P_high_confidence_symbR;
  float llr_high_confidence = log(likelihood_ratio_high_confidence);

  printf("\tLLR High Confidence Region: %f\n", llr_high_confidence);
  printf("\tLLR Low Confidence Region: %f\n", llr_low_confidence);
  printf("\tRatio strong-to-weak LLR: %f\n",
         llr_high_confidence / llr_low_confidence);

  return rbers;
}

int main() {
  float rbers[] = {0.016, 0.0165, 0.017, 0.0175, 0.018, 0.0185};
  for (int i = 0; i < sizeof(rbers) / sizeof(float); i++) {
    printf("rber: %f\n", rbers[i]);
    s_rbers error_probabilities = f_find_error_probabilities(rbers[i]);
  }
}
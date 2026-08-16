// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <random>
#include <string>

/**
 * Computes the probability density function (PDF) of the Generalized Kumaraswamy distribution
 * in a numerically stable way using logarithmic calculations.
 *
 * @param x Value at which to evaluate the PDF (must be in (0,1))
 * @param theta Parameter vector [alpha, beta, gamma, delta, lambda]
 * @return The PDF value at x
 */
double gkw_pdf(double x, const arma::vec &theta) {
  if(x <= 0.0 || x >= 1.0) return 0.0;

  const double alpha  = theta(0);
  const double beta   = theta(1);
  const double gamma_ = theta(2);
  const double delta  = theta(3);
  const double lambda = theta(4);

  if(alpha <= 0.0 || beta <= 0.0 || gamma_ <= 0.0 || delta < 0.0 || lambda <= 0.0) return 0.0;

  try {
    const double log_x = std::log(x);
    const double alpha_log_x = alpha * log_x;

    if (!std::isfinite(alpha_log_x)) return 0.0;

    const double x_alpha = std::exp(alpha_log_x);
    const double one_minus_x_alpha = 1.0 - x_alpha;

    if (one_minus_x_alpha <= 0.0 || one_minus_x_alpha >= 1.0) return 0.0;

    double log_pdf = std::log(lambda) + std::log(alpha) + std::log(beta) + (alpha - 1.0) * log_x;
    log_pdf += (beta - 1.0) * std::log(one_minus_x_alpha);

    const double beta_log_one_minus_x_alpha = beta * std::log(one_minus_x_alpha);
    if (!std::isfinite(beta_log_one_minus_x_alpha)) return 0.0;

    const double one_minus_x_alpha_beta = std::exp(beta_log_one_minus_x_alpha);
    const double inner = 1.0 - one_minus_x_alpha_beta;

    if (inner <= 0.0 || inner >= 1.0) return 0.0;

    log_pdf += (gamma_ * lambda - 1.0) * std::log(inner);

    const double lambda_log_inner = lambda * std::log(inner);
    if (!std::isfinite(lambda_log_inner)) return 0.0;

    const double inner_lambda = std::exp(lambda_log_inner);
    const double term3 = 1.0 - inner_lambda;

    if (term3 <= 0.0 || term3 >= 1.0) return 0.0;

    log_pdf += delta * std::log(term3);

    const double log_beta_val = std::lgamma(gamma_) + std::lgamma(delta + 1.0) - std::lgamma(gamma_ + delta + 1.0);
    if (!std::isfinite(log_beta_val)) return 0.0;

    log_pdf -= log_beta_val;

    if (!std::isfinite(log_pdf)) return 0.0;

    const double pdf = std::exp(log_pdf);
    return std::isfinite(pdf) ? pdf : 0.0;
  } catch (...) {
    return 0.0;
  }
}

/**
 * Computes the PDF of the Beta-Kumaraswamy (BKw) distribution
 */
double bkw_pdf(double x, const arma::vec &theta) {
  if(x <= 0.0 || x >= 1.0) return 0.0;

  arma::vec full_theta(5);
  full_theta(0) = theta(0);
  full_theta(1) = theta(1);
  full_theta(2) = theta(2);
  full_theta(3) = theta(3);
  full_theta(4) = 1.0;

  return gkw_pdf(x, full_theta);
}

/**
 * Computes the PDF of the Kumaraswamy-Kumaraswamy (KKw) distribution
 */
double kkw_pdf(double x, const arma::vec &theta) {
  if(x <= 0.0 || x >= 1.0) return 0.0;

  arma::vec full_theta(5);
  full_theta(0) = theta(0);
  full_theta(1) = theta(1);
  full_theta(2) = 1.0;
  full_theta(3) = theta(2);
  full_theta(4) = theta(3);

  return gkw_pdf(x, full_theta);
}

/**
 * Computes the PDF of the Exponentiated Kumaraswamy (EKw) distribution
 */
double ekw_pdf(double x, const arma::vec &theta) {
  if(x <= 0.0 || x >= 1.0) return 0.0;

  arma::vec full_theta(5);
  full_theta(0) = theta(0);
  full_theta(1) = theta(1);
  full_theta(2) = 1.0;
  full_theta(3) = 0.0;  // EKw = GKw(a,b,1,0,l)
  full_theta(4) = theta(2);

  return gkw_pdf(x, full_theta);
}

/**
 * Computes the PDF of the McDonald (MC) distribution
 */
double mc_pdf(double x, const arma::vec &theta) {
  if(x <= 0.0 || x >= 1.0) return 0.0;

  arma::vec full_theta(5);
  full_theta(0) = 1.0;
  full_theta(1) = 1.0;
  full_theta(2) = theta(0);
  full_theta(3) = theta(1);
  full_theta(4) = theta(2);

  return gkw_pdf(x, full_theta);
}

/**
 * Computes the PDF of the Kumaraswamy (Kw) distribution
 */
double kw_pdf(double x, const arma::vec &theta) {
  if(x <= 0.0 || x >= 1.0) return 0.0;

  arma::vec full_theta(5);
  full_theta(0) = theta(0);
  full_theta(1) = theta(1);
  full_theta(2) = 1.0;
  full_theta(3) = 0.0;  // Kw = GKw(a,b,1,0,1)
  full_theta(4) = 1.0;

  return gkw_pdf(x, full_theta);
}

/**
 * Computes the PDF of the Beta distribution
 */
double beta_pdf(double x, const arma::vec &theta) {
  if(x <= 0.0 || x >= 1.0) return 0.0;

  arma::vec full_theta(5);
  full_theta(0) = 1.0;
  full_theta(1) = 1.0;
  full_theta(2) = theta(0);
  full_theta(3) = theta(1);
  full_theta(4) = 1.0;

  return gkw_pdf(x, full_theta);
}

/**
 * Generic PDF dispatcher based on family
 */
double family_pdf(double x, const arma::vec &theta, const std::string &family) {
  if (family == "gkw") {
    return gkw_pdf(x, theta);
  } else if (family == "bkw") {
    return bkw_pdf(x, theta);
  } else if (family == "kkw") {
    return kkw_pdf(x, theta);
  } else if (family == "ekw") {
    return ekw_pdf(x, theta);
  } else if (family == "mc") {
    return mc_pdf(x, theta);
  } else if (family == "kw") {
    return kw_pdf(x, theta);
  } else if (family == "beta") {
    return beta_pdf(x, theta);
  }
  return 0.0;
}

/**
 * Calculate theoretical moment using Simpson's rule with adaptive refinement
 */
double moment_theoretical(int r, const arma::vec &theta, const std::string &family) {
  try {
    const int initial_points = 51;
    const double a = 0.0, b = 1.0;
    const double h = (b - a) / (initial_points - 1);

    double sum = 0.0;

    for (int i = 0; i < initial_points; i++) {
      const double x = a + i * h;
      double weight = 0.0;

      if (i == 0 || i == initial_points - 1)
        weight = 1.0;
      else if (i % 2 == 1)
        weight = 4.0;
      else
        weight = 2.0;

      const double xr = std::pow(x, r);
      const double fx = xr * family_pdf(x, theta, family);

      sum += weight * fx;
    }

    double result = (h/3.0) * sum;

    if (!std::isfinite(result) || std::abs(result) < 1e-10) {
      const int refined_points = 201;
      const double h_refined = (b - a) / (refined_points - 1);

      double sum_refined = 0.0;
      std::vector<double> x_values(refined_points);
      std::vector<double> fx_values(refined_points);

      for (int i = 0; i < refined_points; i++) {
        x_values[i] = a + i * h_refined;
      }

      for (int i = 0; i < refined_points; i++) {
        const double x = x_values[i];
        const double xr = std::pow(x, r);
        fx_values[i] = xr * family_pdf(x, theta, family);
      }

      for (int i = 0; i < refined_points - 1; i++) {
        sum_refined += 0.5 * (fx_values[i] + fx_values[i+1]) * (x_values[i+1] - x_values[i]);
      }

      result = sum_refined;
    }

    if (!std::isfinite(result) || std::abs(result) < 1e-14) {
      double alpha = 1.0, beta = 1.0;

      if (family == "gkw" || family == "bkw" || family == "kkw" || family == "ekw" || family == "kw") {
        alpha = theta(0);
        beta = theta(1);
      }

      result = beta / ((r/alpha) + beta);
    }

    return result;
  } catch (...) {
    return 0.5;
  }
}

/**
 * Objective function
 */
double objective_function(const arma::vec &theta,
                          const arma::vec &sample_moments,
                          const std::string &family) {
  if (arma::any(theta <= 0.0)) {
    return std::numeric_limits<double>::max();
  }

  double error = 0.0;
  bool has_valid_result = false;

  const arma::vec weights = {1.0, 0.8, 0.6, 0.4, 0.2};

  for (int r = 1; r <= 5; r++) {
    const double theor = moment_theoretical(r, theta, family);

    if (!std::isfinite(theor)) {
      continue;
    }

    const double sample_moment = sample_moments(r-1);

    if (std::abs(sample_moment) < 1e-10) {
      error += weights(r-1) * std::pow(theor, 2);
    } else {
      const double rel_err = (theor - sample_moment) / sample_moment;
      error += weights(r-1) * std::pow(rel_err, 2);
    }

    has_valid_result = true;
  }

  if (!has_valid_result || !std::isfinite(error)) {
    return std::numeric_limits<double>::max();
  }

  return error;
}

/**
 * Nelder-Mead optimization
 */
arma::vec optimize_nelder_mead(const arma::vec &initial,
                               const arma::vec &sample_moments,
                               const std::string &family,
                               int max_iter = 1000,
                               double tol = 1e-6) {
  const int n = initial.n_elem;

  arma::mat simplex(n, n+1);
  simplex.col(0) = initial;

  for (int i = 1; i <= n; i++) {
    arma::vec step = initial;
    double step_size = 0.05 * std::max(std::abs(initial(i-1)), 0.1);
    step(i-1) += step_size;

    for (int j = 0; j < n; j++) {
      if (step(j) <= 0.0) step(j) = 0.01;
    }

    simplex.col(i) = step;
  }

  arma::vec function_values(n+1);
  for (int i = 0; i <= n; i++) {
    function_values(i) = objective_function(simplex.col(i), sample_moments, family);
  }

  const double alpha = 1.0;
  const double gamma = 2.0;
  const double rho = 0.5;
  const double sigma = 0.5;

  for (int iter = 0; iter < max_iter; iter++) {
    arma::uvec sorted_indices = arma::sort_index(function_values);
    arma::mat sorted_simplex(n, n+1);
    arma::vec sorted_values(n+1);

    for (int i = 0; i <= n; i++) {
      sorted_simplex.col(i) = simplex.col(sorted_indices(i));
      sorted_values(i) = function_values(sorted_indices(i));
    }

    simplex = sorted_simplex;
    function_values = sorted_values;

    double diameter = 0.0;
    for (int i = 1; i <= n; i++) {
      diameter = std::max(diameter, arma::norm(simplex.col(i) - simplex.col(0)));
    }

    if (diameter < tol) {
      break;
    }

    arma::vec centroid = arma::mean(simplex.cols(0, n-1), 1);

    arma::vec reflected = centroid + alpha * (centroid - simplex.col(n));

    for (int j = 0; j < n; j++) {
      if (reflected(j) <= 0.0) reflected(j) = 0.01;
    }

    double f_reflected = objective_function(reflected, sample_moments, family);

    if (f_reflected < function_values(0)) {
      arma::vec expanded = centroid + gamma * (reflected - centroid);

      for (int j = 0; j < n; j++) {
        if (expanded(j) <= 0.0) expanded(j) = 0.01;
      }

      double f_expanded = objective_function(expanded, sample_moments, family);

      if (f_expanded < f_reflected) {
        simplex.col(n) = expanded;
        function_values(n) = f_expanded;
      } else {
        simplex.col(n) = reflected;
        function_values(n) = f_reflected;
      }
    }
    else if (f_reflected < function_values(n-1)) {
      simplex.col(n) = reflected;
      function_values(n) = f_reflected;
    }
    else {
      arma::vec contracted;
      double f_contracted;

      if (f_reflected < function_values(n)) {
        contracted = centroid + rho * (reflected - centroid);

        for (int j = 0; j < n; j++) {
          if (contracted(j) <= 0.0) contracted(j) = 0.01;
        }

        f_contracted = objective_function(contracted, sample_moments, family);

        if (f_contracted <= f_reflected) {
          simplex.col(n) = contracted;
          function_values(n) = f_contracted;
        } else {
          for (int i = 1; i <= n; i++) {
            simplex.col(i) = simplex.col(0) + sigma * (simplex.col(i) - simplex.col(0));

            for (int j = 0; j < n; j++) {
              if (simplex(j, i) <= 0.0) simplex(j, i) = 0.01;
            }

            function_values(i) = objective_function(simplex.col(i), sample_moments, family);
          }
        }
      } else {
        contracted = centroid - rho * (centroid - simplex.col(n));

        for (int j = 0; j < n; j++) {
          if (contracted(j) <= 0.0) contracted(j) = 0.01;
        }

        f_contracted = objective_function(contracted, sample_moments, family);

        if (f_contracted < function_values(n)) {
          simplex.col(n) = contracted;
          function_values(n) = f_contracted;
        } else {
          for (int i = 1; i <= n; i++) {
            simplex.col(i) = simplex.col(0) + sigma * (simplex.col(i) - simplex.col(0));

            for (int j = 0; j < n; j++) {
              if (simplex(j, i) <= 0.0) simplex(j, i) = 0.01;
            }

            function_values(i) = objective_function(simplex.col(i), sample_moments, family);
          }
        }
      }
    }
  }

  arma::uvec best_idx = arma::sort_index(function_values);
  return simplex.col(best_idx(0));
}

/**
 * Get number of parameters
 */
int get_npar(const std::string &family) {
  if (family == "gkw") return 5;
  if (family == "bkw") return 4;
  if (family == "kkw") return 4;
  if (family == "ekw") return 3;
  if (family == "mc") return 3;
  if (family == "kw") return 2;
  if (family == "beta") return 2;
  return 5;
}

/**
 * Get parameter names
 */
Rcpp::CharacterVector get_param_names(const std::string &family) {
  if (family == "gkw") {
    return Rcpp::CharacterVector::create("alpha", "beta", "gamma", "delta", "lambda");
  } else if (family == "bkw") {
    return Rcpp::CharacterVector::create("alpha", "beta", "gamma", "delta");
  } else if (family == "kkw") {
    return Rcpp::CharacterVector::create("alpha", "beta", "delta", "lambda");
  } else if (family == "ekw") {
    return Rcpp::CharacterVector::create("alpha", "beta", "lambda");
  } else if (family == "mc") {
    return Rcpp::CharacterVector::create("gamma", "delta", "lambda");
  } else if (family == "kw") {
    return Rcpp::CharacterVector::create("alpha", "beta");
  } else if (family == "beta") {
    return Rcpp::CharacterVector::create("gamma", "delta");
  }
  return Rcpp::CharacterVector::create("alpha", "beta", "gamma", "delta", "lambda");
}

/**
 * Generate initial points
 */
std::vector<arma::vec> generate_initial_points(const arma::vec &sample_moments,
                                               const std::string &family,
                                               int n_starts) {
  std::vector<arma::vec> initial_points;

  double m1 = sample_moments(0);
  double m2 = sample_moments(1);
  double var = m2 - m1*m1;

  if (var <= 1e-10) var = 0.01;
  if (m1 <= 0.01) m1 = 0.01;
  if (m1 >= 0.99) m1 = 0.99;

  double alpha_init = std::max(0.1, m1 * (1.0 - m1) / var - m1);
  double beta_init = std::max(0.1, alpha_init * m1 / (1.0 - m1));
  alpha_init = std::min(20.0, std::max(0.1, alpha_init));
  beta_init = std::min(20.0, std::max(0.1, beta_init));

  std::mt19937 gen(42);

  if (family == "gkw") {
    initial_points.push_back(arma::vec({alpha_init, beta_init, 1.0, 0.1, 1.0}));
    initial_points.push_back(arma::vec({2.0, 2.0, 1.0, 0.5, 1.0}));
    initial_points.push_back(arma::vec({1.0, 1.0, 1.0, 0.1, 1.0}));
    initial_points.push_back(arma::vec({4.0, 2.0, 0.8, 0.5, 1.0}));

    std::uniform_real_distribution<double> dist_alpha(0.5, 10.0);
    std::uniform_real_distribution<double> dist_beta(0.5, 10.0);
    std::uniform_real_distribution<double> dist_gamma(0.5, 2.0);
    std::uniform_real_distribution<double> dist_delta(0.1, 1.0);
    std::uniform_real_distribution<double> dist_lambda(0.5, 2.0);

    for (int i = initial_points.size(); i < n_starts; i++) {
      initial_points.push_back(arma::vec({
        dist_alpha(gen), dist_beta(gen), dist_gamma(gen),
        dist_delta(gen), dist_lambda(gen)
      }));
    }
  }
  else if (family == "bkw") {
    initial_points.push_back(arma::vec({alpha_init, beta_init, 1.0, 0.5}));
    initial_points.push_back(arma::vec({2.0, 2.0, 1.0, 0.5}));
    initial_points.push_back(arma::vec({1.0, 1.0, 0.8, 0.3}));
    initial_points.push_back(arma::vec({3.0, 2.0, 1.5, 0.5}));

    std::uniform_real_distribution<double> dist_alpha(0.5, 10.0);
    std::uniform_real_distribution<double> dist_beta(0.5, 10.0);
    std::uniform_real_distribution<double> dist_gamma(0.5, 2.0);
    std::uniform_real_distribution<double> dist_delta(0.1, 1.0);

    for (int i = initial_points.size(); i < n_starts; i++) {
      initial_points.push_back(arma::vec({
        dist_alpha(gen), dist_beta(gen), dist_gamma(gen), dist_delta(gen)
      }));
    }
  }
  else if (family == "kkw") {
    initial_points.push_back(arma::vec({alpha_init, beta_init, 0.5, 1.0}));
    initial_points.push_back(arma::vec({2.0, 2.0, 0.5, 1.0}));
    initial_points.push_back(arma::vec({1.0, 1.0, 0.3, 1.2}));
    initial_points.push_back(arma::vec({3.0, 2.0, 0.7, 1.5}));

    std::uniform_real_distribution<double> dist_alpha(0.5, 10.0);
    std::uniform_real_distribution<double> dist_beta(0.5, 10.0);
    std::uniform_real_distribution<double> dist_delta(0.1, 1.0);
    std::uniform_real_distribution<double> dist_lambda(0.5, 2.0);

    for (int i = initial_points.size(); i < n_starts; i++) {
      initial_points.push_back(arma::vec({
        dist_alpha(gen), dist_beta(gen), dist_delta(gen), dist_lambda(gen)
      }));
    }
  }
  else if (family == "ekw") {
    initial_points.push_back(arma::vec({alpha_init, beta_init, 1.0}));
    initial_points.push_back(arma::vec({2.0, 2.0, 1.0}));
    initial_points.push_back(arma::vec({1.0, 1.0, 1.2}));
    initial_points.push_back(arma::vec({3.0, 2.0, 1.5}));

    std::uniform_real_distribution<double> dist_alpha(0.5, 10.0);
    std::uniform_real_distribution<double> dist_beta(0.5, 10.0);
    std::uniform_real_distribution<double> dist_lambda(0.5, 2.0);

    for (int i = initial_points.size(); i < n_starts; i++) {
      initial_points.push_back(arma::vec({
        dist_alpha(gen), dist_beta(gen), dist_lambda(gen)
      }));
    }
  }
  else if (family == "mc") {
    double gamma_init = m1 > 0.5 ? 2.0 : 1.0;
    double delta_init = m1 < 0.5 ? 2.0 : 1.0;

    initial_points.push_back(arma::vec({gamma_init, delta_init, 1.0}));
    initial_points.push_back(arma::vec({1.0, 1.0, 1.0}));
    initial_points.push_back(arma::vec({2.0, 2.0, 1.2}));
    initial_points.push_back(arma::vec({1.5, 1.5, 1.5}));

    std::uniform_real_distribution<double> dist_gamma(0.5, 5.0);
    std::uniform_real_distribution<double> dist_delta(0.5, 5.0);
    std::uniform_real_distribution<double> dist_lambda(0.5, 2.0);

    for (int i = initial_points.size(); i < n_starts; i++) {
      initial_points.push_back(arma::vec({
        dist_gamma(gen), dist_delta(gen), dist_lambda(gen)
      }));
    }
  }
  else if (family == "kw") {
    initial_points.push_back(arma::vec({alpha_init, beta_init}));
    initial_points.push_back(arma::vec({2.0, 2.0}));
    initial_points.push_back(arma::vec({1.0, 1.0}));
    initial_points.push_back(arma::vec({3.0, 2.0}));

    std::uniform_real_distribution<double> dist_alpha(0.5, 10.0);
    std::uniform_real_distribution<double> dist_beta(0.5, 10.0);

    for (int i = initial_points.size(); i < n_starts; i++) {
      initial_points.push_back(arma::vec({dist_alpha(gen), dist_beta(gen)}));
    }
  }
  else if (family == "beta") {
    double gamma_init = m1 * ((m1 * (1.0 - m1) / var) - 1.0);
    double delta_init = (1.0 - m1) * ((m1 * (1.0 - m1) / var) - 1.0);
    gamma_init = std::min(50.0, std::max(0.1, gamma_init));
    delta_init = std::min(50.0, std::max(0.1, delta_init));

    initial_points.push_back(arma::vec({gamma_init, delta_init}));
    initial_points.push_back(arma::vec({2.0, 2.0}));
    initial_points.push_back(arma::vec({1.0, 1.0}));
    initial_points.push_back(arma::vec({3.0, 2.0}));

    std::uniform_real_distribution<double> dist_gamma(0.5, 10.0);
    std::uniform_real_distribution<double> dist_delta(0.5, 10.0);

    for (int i = initial_points.size(); i < n_starts; i++) {
      initial_points.push_back(arma::vec({dist_gamma(gen), dist_delta(gen)}));
    }
  }

  return initial_points;
}

/**
 * Constrain parameters
 */
arma::vec constrain_parameters(const arma::vec &theta, const std::string &family) {
  arma::vec constrained = theta;

  for (size_t i = 0; i < constrained.n_elem; i++) {
    if (!std::isfinite(constrained(i)) || constrained(i) <= 0.0) {
      constrained(i) = 1.0;
    }
  }

  if (family == "gkw") {
    constrained(0) = std::min(50.0, std::max(0.1, constrained(0)));
    constrained(1) = std::min(50.0, std::max(0.1, constrained(1)));
    constrained(2) = std::min(10.0, std::max(0.1, constrained(2)));
    constrained(3) = std::min(10.0, std::max(0.01, constrained(3)));
    constrained(4) = std::min(20.0, std::max(0.1, constrained(4)));
  }
  else if (family == "bkw") {
    constrained(0) = std::min(50.0, std::max(0.1, constrained(0)));
    constrained(1) = std::min(50.0, std::max(0.1, constrained(1)));
    constrained(2) = std::min(10.0, std::max(0.1, constrained(2)));
    constrained(3) = std::min(10.0, std::max(0.01, constrained(3)));
  }
  else if (family == "kkw") {
    constrained(0) = std::min(50.0, std::max(0.1, constrained(0)));
    constrained(1) = std::min(50.0, std::max(0.1, constrained(1)));
    constrained(2) = std::min(10.0, std::max(0.01, constrained(2)));
    constrained(3) = std::min(20.0, std::max(0.1, constrained(3)));
  }
  else if (family == "ekw") {
    constrained(0) = std::min(50.0, std::max(0.1, constrained(0)));
    constrained(1) = std::min(50.0, std::max(0.1, constrained(1)));
    constrained(2) = std::min(20.0, std::max(0.1, constrained(2)));
  }
  else if (family == "mc") {
    constrained(0) = std::min(20.0, std::max(0.1, constrained(0)));
    constrained(1) = std::min(20.0, std::max(0.1, constrained(1)));
    constrained(2) = std::min(20.0, std::max(0.1, constrained(2)));
  }
  else if (family == "kw") {
    constrained(0) = std::min(50.0, std::max(0.1, constrained(0)));
    constrained(1) = std::min(50.0, std::max(0.1, constrained(1)));
  }
  else if (family == "beta") {
    constrained(0) = std::min(50.0, std::max(0.1, constrained(0)));
    constrained(1) = std::min(50.0, std::max(0.1, constrained(1)));
  }

  return constrained;
}

//' @title Estimate Distribution Parameters Using Method of Moments
//'
//' @description
//' Estimates parameters for various distribution families from the Generalized Kumaraswamy
//' family using the method of moments. The implementation is optimized for numerical
//' stability and computational efficiency through Nelder-Mead optimization and adaptive
//' numerical integration.
//'
//' @param x Numeric vector of observations. All values must be in the open interval (0,1).
//'   Values outside this range will be automatically truncated to avoid numerical issues.
//' @param family Character string specifying the distribution family. Valid options are:
//'   \code{"gkw"} (Generalized Kumaraswamy - 5 parameters),
//'   \code{"bkw"} (Beta-Kumaraswamy - 4 parameters),
//'   \code{"kkw"} (Kumaraswamy-Kumaraswamy - 4 parameters),
//'   \code{"ekw"} (Exponentiated Kumaraswamy - 3 parameters),
//'   \code{"mc"} (McDonald - 3 parameters),
//'   \code{"kw"} (Kumaraswamy - 2 parameters),
//'   \code{"beta"} (Beta - 2 parameters).
//'   The string is case-insensitive. Default is \code{"gkw"}.
//' @param n_starts Integer specifying the number of different initial parameter values
//'   to try during optimization. More starting points increase the probability of finding
//'   the global optimum at the cost of longer computation time. Default is 5.
//'
//' @return Named numeric vector containing the estimated parameters for the specified
//'   distribution family. Parameter names correspond to the distribution specification.
//'   If estimation fails, returns a vector of NA values with appropriate parameter names.
//'
//' @details
//' The function uses the method of moments to estimate distribution parameters by minimizing
//' the weighted sum of squared relative errors between theoretical and sample moments
//' (orders 1 through 5). The optimization employs the Nelder-Mead simplex algorithm,
//' which is derivative-free and particularly robust for this problem.
//'
//' Key implementation features: logarithmic calculations for numerical stability,
//' adaptive numerical integration using Simpson's rule with fallback to trapezoidal rule,
//' multiple random starting points to avoid local minima, decreasing weights for
//' higher-order moments (1.0, 0.8, 0.6, 0.4, 0.2), and automatic parameter constraint
//' enforcement.
//'
//' Parameter Constraints:
//' All parameters are constrained to positive values. Additionally, family-specific
//' constraints are enforced: alpha and beta in (0.1, 50.0), gamma in (0.1, 10.0) for
//' GKw-related families or (0.1, 50.0) for Beta, delta in (0.01, 10.0), and lambda in
//' (0.1, 20.0).
//'
//' The function will issue warnings for empty input vectors, sample sizes less than 10
//' (unreliable estimation), or failure to find valid parameter estimates (returns defaults).
//'
//' @examples
//' \donttest{
//' # Generate sample data from Beta distribution
//' set.seed(123)
//' x <- rbeta(100, shape1 = 2, shape2 = 3)
//'
//' # Estimate Beta parameters
//' params_beta <- gkwgetstartvalues(x, family = "beta")
//' print(params_beta)
//'
//' # Estimate Kumaraswamy parameters
//' params_kw <- gkwgetstartvalues(x, family = "kw")
//' print(params_kw)
//'
//' # Estimate GKw parameters with more starting points
//' params_gkw <- gkwgetstartvalues(x, family = "gkw", n_starts = 10)
//' print(params_gkw)
//' }
//'
//' @references
//' Jones, M. C. (2009). Kumaraswamy's distribution: A beta-type distribution with
//' some tractability advantages. Statistical Methodology, 6(1), 70-81.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector gkwgetstartvalues(const Rcpp::NumericVector &x,
                                     std::string family = "gkw",
                                     int n_starts = 5) {
 std::transform(family.begin(), family.end(), family.begin(), ::tolower);

 if (family != "gkw" && family != "bkw" && family != "kkw" &&
     family != "ekw" && family != "mc" && family != "kw" && family != "beta") {
   Rcpp::stop("Invalid family. Must be one of: 'gkw', 'bkw', 'kkw', 'ekw', 'mc', 'kw', 'beta'");
 }

 if (x.size() == 0) {
   Rcpp::warning("Empty input vector");
   int npar = get_npar(family);
   Rcpp::NumericVector result(npar, Rcpp::NumericVector::get_na());
   result.attr("names") = get_param_names(family);
   return result;
 }

 try {
   // Filter NAs and clamp to (0,1) — std::max(1e-10, NaN) is platform-dependent
   std::vector<double> clean;
   clean.reserve(x.size());
   for (int i = 0; i < x.size(); i++) {
     if (!Rcpp::NumericVector::is_na(x[i]) && std::isfinite(x[i])) {
       clean.push_back(std::max(1e-10, std::min(1.0 - 1e-10, (double)x[i])));
     }
   }
   if (clean.empty()) {
     Rcpp::warning("No valid (non-NA, finite, in (0,1)) observations found.");
     int npar = get_npar(family);
     Rcpp::NumericVector result(npar, Rcpp::NumericVector::get_na());
     result.attr("names") = get_param_names(family);
     return result;
   }
   arma::vec data(clean.data(), clean.size());

   int n = data.n_elem;

   if(n < 10) {
     Rcpp::warning("Insufficient data for robust estimation (n < 10). Results may be unreliable.");
   }

   arma::vec sample_moments(5, arma::fill::zeros);
   for (int r = 1; r <= 5; r++) {
     sample_moments(r-1) = arma::mean(arma::pow(data, r));
   }

   std::vector<arma::vec> initial_points = generate_initial_points(sample_moments, family, n_starts);

   int npar = get_npar(family);
   arma::vec best_theta(npar);
   double best_obj = std::numeric_limits<double>::max();
   bool found_valid_solution = false;

   for (size_t i = 0; i < initial_points.size(); i++) {
     try {
       double obj_val = objective_function(initial_points[i], sample_moments, family);

       if (std::isfinite(obj_val) && obj_val < best_obj) {
         arma::vec optimized = optimize_nelder_mead(initial_points[i], sample_moments, family);

         obj_val = objective_function(optimized, sample_moments, family);

         if (std::isfinite(obj_val) && obj_val < best_obj) {
           best_theta = optimized;
           best_obj = obj_val;
           found_valid_solution = true;
         }
       }
     } catch (...) {
       continue;
     }
   }

   if (!found_valid_solution) {
     Rcpp::warning("Could not find valid parameter estimates. Using defaults.");
     if (family == "gkw") {
       best_theta = arma::vec({1.0, 1.0, 1.0, 0.1, 1.0});
     } else if (family == "bkw") {
       best_theta = arma::vec({1.0, 1.0, 1.0, 0.5});
     } else if (family == "kkw") {
       best_theta = arma::vec({1.0, 1.0, 0.5, 1.0});
     } else if (family == "ekw") {
       best_theta = arma::vec({1.0, 1.0, 1.0});
     } else if (family == "mc") {
       best_theta = arma::vec({1.0, 1.0, 1.0});
     } else if (family == "kw") {
       best_theta = arma::vec({1.0, 1.0});
     } else if (family == "beta") {
       best_theta = arma::vec({1.0, 1.0});
     }
   }

   best_theta = constrain_parameters(best_theta, family);

   Rcpp::NumericVector result(npar);
   for (int i = 0; i < npar; i++) {
     result[i] = best_theta(i);
   }
   result.attr("names") = get_param_names(family);

   return result;
 } catch (std::exception &e) {
   Rcpp::warning("Exception in parameter estimation: %s", e.what());
   int npar = get_npar(family);
   Rcpp::NumericVector result(npar, Rcpp::NumericVector::get_na());
   result.attr("names") = get_param_names(family);
   return result;
 } catch (...) {
   Rcpp::warning("Unknown exception in parameter estimation");
   int npar = get_npar(family);
   Rcpp::NumericVector result(npar, Rcpp::NumericVector::get_na());
   result.attr("names") = get_param_names(family);
   return result;
 }
}

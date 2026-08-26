args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "r_results.csv"

x <- as.numeric(1:100)
pair_y <- 0.5 * x + sin(x / 7)
weights <- 1 + (seq_along(x) %% 7)
quantile_x <- c(8, 1, 5, 3, 9, 2)
quantile_weights <- c(1, 2, 4, 1, 3, 2)
quantile_probs <- c(0, 0.1, 0.25, 0.5, 0.9, 1)
points <- c(0.1, 0.5, 1, 2.5, 10, 100)
log1p_points <- c(-0.9, -1e-12, 0, 1e-12, 1, 1e6)
expm1_points <- c(-100, -1, -1e-12, 1e-12, 1, 10)
log1mexp_points <- c(-100, -10, -1, -1e-12)
softplus_points <- c(-100, -1, 0, 1, 100)
log_weight_points <- c(-1000, -100, -2, 0, 1, 100)
probability_points <- c(1e-12, 0.01, 0.25, 0.5, 0.99, 1 - 1e-12)
log_probability_points <- c(-1, -10, -100, -1000)
order_x <- c(3, 1, 3, 2, 1, 3)
gamma_shapes <- c(0.5, 1, 2.5, 10, 100)
gamma_x <- c(0.1, 1, 3, 8, 110)
beta_x <- c(0.01, 0.1, 0.5, 0.8, 0.99)
beta_shape1 <- c(0.5, 1, 2, 5, 20)
beta_shape2 <- c(2, 4, 3, 1.5, 10)
distribution_x <- c(0.1, 0.5, 1, 2, 5, 10)
student_x <- c(-3, -1, 0, 0.5, 2, 5)
distribution_probs <- c(0.01, 0.1, 0.25, 0.5, 0.9, 0.99)
distribution_df <- c(1, 2, 5, 10, 30, 100)
f_df1 <- c(1, 2, 5, 10, 20, 50)
f_df2 <- c(2, 5, 8, 12, 30, 100)
factorial_points <- c(0, 1, 5, 20, 1000)
choose_n <- c(5, 20, 1000, 1000000)
choose_k <- c(2, 17, 500, 2)
wcheck <- function(value) sum(as.numeric(value) * seq_along(value))
roll_right <- function(value, width) {
  as.numeric(stats::filter(value, rep(1 / width, width), sides = 1))
}
weighted_quantile_ecdf <- function(value, weight, probability) {
  keep <- weight > 0
  value <- value[keep]
  weight <- weight[keep]
  order <- order(value)
  value <- value[order]
  cumulative <- cumsum(weight[order])
  vapply(probability, function(p) value[which(cumulative >= p * sum(weight))[1]], numeric(1))
}
prepare_weighted <- function(value, weight) {
  keep <- weight > 0
  value <- value[keep]
  weight <- weight[keep]
  ord <- order(value)
  list(value = value[ord], weight = weight[ord])
}
weighted_quantile_linear_cdf <- function(value, weight, probability) {
  data <- prepare_weighted(value, weight)
  value <- data$value
  cumulative <- cumsum(data$weight) / sum(data$weight)
  vapply(probability, function(p) {
    if (p <= cumulative[1]) return(value[1])
    high <- which(cumulative > p)[1]
    if (is.na(high)) return(tail(value, 1))
    fraction <- (p - cumulative[high - 1]) / (cumulative[high] - cumulative[high - 1])
    value[high - 1] + fraction * (value[high] - value[high - 1])
  }, numeric(1))
}
weighted_quantile_frequency_type7 <- function(value, weight, probability) {
  data <- prepare_weighted(value, weight)
  value <- data$value
  weight <- data$weight
  cumulative <- cumsum(weight)
  target <- 1 + (sum(weight) - 1) * probability
  vapply(target, function(target_value) {
    i <- which(cumulative >= target_value)[1]
    if (i == 1 || abs(cumulative[i] - target_value) > 1e-12 || i == length(value)) return(value[i])
    fraction <- target_value - floor(target_value)
    (1 - fraction) * value[i] + fraction * value[i + 1]
  }, numeric(1))
}
weighted_quantile_isotone <- function(value, weight, probability) {
  data <- prepare_weighted(value, weight)
  value <- data$value
  weight <- data$weight
  cumulative <- cumsum(weight)
  total <- sum(weight)
  tolerance <- 10 * .Machine$double.eps * max(1, abs(total))
  vapply(probability, function(p) {
    target <- p * total
    i <- which(cumulative >= target - tolerance)[1]
    if (abs(cumulative[i] - target) <= tolerance && i < length(value)) {
      return((weight[i] * value[i] + weight[i + 1] * value[i + 1]) / (weight[i] + weight[i + 1]))
    }
    value[i]
  }, numeric(1))
}
linear_interpolation <- function(knots, value, p) {
  if (p <= knots[1]) return(value[1])
  if (p >= tail(knots, 1)) return(tail(value, 1))
  i <- which(p >= knots[-length(knots)] & p <= knots[-1])[1]
  if (abs(knots[i + 1] - knots[i]) <= .Machine$double.eps * max(1, abs(knots[i]))) return(value[i + 1])
  value[i] + (value[i + 1] - value[i]) * (p - knots[i]) / (knots[i + 1] - knots[i])
}
constant_interpolation <- function(knots, value, p) {
  if (p <= knots[1]) return(value[1])
  i <- which(p < knots[-1])[1]
  if (is.na(i)) tail(value, 1) else value[i]
}
weighted_quantile_survey <- function(value, weight, p, rule) {
  data <- prepare_weighted(value, weight)
  value <- data$value
  weight <- data$weight
  n <- length(value)
  if (n == 1) return(value)
  total <- sum(weight)
  cumulative <- cumsum(weight)
  if (rule %in% c(1, 2, 4, 5, 7)) {
    lower <- tail(which(cumulative <= p * total), 1)
    if (!length(lower)) lower <- 1
    next_value <- min(n, lower + 1)
    lower_gap <- p - cumulative[lower] / total
    upper_gap <- cumulative[next_value] / total - p
    if (rule %in% c(1, 4)) return(if (lower_gap <= 0) value[lower] else value[next_value])
    if (rule %in% c(2, 5)) return(if (lower_gap <= 0) mean(value[c(lower, next_value)]) else value[next_value])
    if (abs(upper_gap + lower_gap) <= .Machine$double.xmin) return(value[lower])
    return(value[lower] + lower_gap * (value[next_value] - value[lower]) / (upper_gap + lower_gap))
  }
  if (rule == 6) {
    unique_value <- unique(value)
    unique_weight <- vapply(unique_value, function(x) sum(weight[value == x]), numeric(1))
    cumulative <- cumsum(unique_weight)
    lower <- tail(which(cumulative <= p * sum(unique_weight)), 1)
    if (!length(lower)) lower <- 1
    next_value <- min(length(unique_value), lower + 1)
    lower_gap <- p - cumulative[lower] / sum(unique_weight)
    return(if (lower_gap <= 0 && lower %% 2 == 0) unique_value[lower] else unique_value[next_value])
  }
  if (rule == 8) return(linear_interpolation((cumulative - 0.5 * weight) / total, value, p))
  if (rule == 9) return(linear_interpolation(cumulative / (total + tail(weight, 1)), value, p))
  if (rule == 3) {
    mean_weight <- total / n
    knots <- (cumsum(weight / mean_weight) + 0.5 - weight / (2 * mean_weight)) / (n + 1)
    return(constant_interpolation(knots, value, p))
  }
  if (rule == 10) return(linear_interpolation(c(0, cumulative[-n] / cumulative[n - 1]), value, p))
  if (rule == 11) {
    knots <- c((2 / 3) * cumulative[1], cumulative[-n] / 3 + 2 * cumulative[-1] / 3)
    return(linear_interpolation(knots / (total + tail(weight, 1) / 3), value, p))
  }
  knots <- c((5 / 8) * cumulative[1], 3 * cumulative[-n] / 8 + 5 * cumulative[-1] / 8)
  linear_interpolation(knots / (total + tail(weight, 1) / 4), value, p)
}

rows <- list()
add <- function(name, fun, reps = 20000L, atol = 1e-11, rtol = 1e-11) {
  tm <- system.time(for (i in seq_len(reps)) value <- as.numeric(fun()))[["elapsed"]]
  rows[[length(rows) + 1L]] <<- data.frame(
    case = name, value = value, seconds = tm, abs_tol = atol, rel_tol = rtol
  )
}

add("digamma_checksum", function() wcheck(digamma(points)))
add("trigamma_checksum", function() wcheck(trigamma(points)))
add("log_beta_checksum", function() wcheck(lbeta(points, rev(points))))
add("rolling_mean_valid", function() wcheck(roll_right(x, 5)[5:100]))
add("rolling_mean_right", function() sum(roll_right(x, 5)[5:100] * (5:100)))
add("log1p_checksum", function() wcheck(log1p(log1p_points)))
add("expm1_checksum", function() wcheck(expm1(expm1_points)))
add("log1mexp_checksum", function() {
  value <- ifelse(log1mexp_points < -log(2),
                  log1p(-exp(log1mexp_points)),
                  log(-expm1(log1mexp_points)))
  wcheck(value)
})
add("log1pexp_checksum", function() {
  value <- ifelse(softplus_points > 0,
                  softplus_points + log1p(exp(-softplus_points)),
                  log1p(exp(softplus_points)))
  wcheck(value)
})
add("log_sum_exp", function() {
  maximum <- max(log_weight_points)
  maximum + log(sum(exp(log_weight_points - maximum)))
})
add("log_mean_exp", function() {
  maximum <- max(log_weight_points)
  maximum + log(mean(exp(log_weight_points - maximum)))
})
add("logistic_checksum", function() wcheck(plogis(softplus_points)))
add("logit_checksum", function() wcheck(qlogis(probability_points)))
add("qnorm_checksum", function() wcheck(qnorm(probability_points)))
add("qnorm_log_probability_checksum", function() wcheck(qnorm(log_probability_points, log.p = TRUE)))
add("student_t_density_checksum", function() wcheck(dt(student_x, distribution_df)), reps = 5000L)
add("student_t_cdf_checksum", function() wcheck(pt(student_x, distribution_df)), reps = 5000L)
add("student_t_log_survival_checksum", function() {
  wcheck(pt(abs(student_x), distribution_df, lower.tail = FALSE, log.p = TRUE))
}, reps = 5000L)
add("student_t_quantile_checksum", function() wcheck(qt(distribution_probs, distribution_df)), reps = 100L)
add("chi_square_density_checksum", function() wcheck(dchisq(distribution_x, distribution_df)), reps = 5000L)
add("chi_square_cdf_checksum", function() wcheck(pchisq(distribution_x, distribution_df)), reps = 5000L)
add("chi_square_log_survival_checksum", function() {
  wcheck(pchisq(distribution_x, distribution_df, lower.tail = FALSE, log.p = TRUE))
}, reps = 5000L)
add("chi_square_quantile_checksum", function() wcheck(qchisq(distribution_probs, distribution_df)), reps = 100L)
add("f_density_checksum", function() wcheck(df(distribution_x, f_df1, f_df2)), reps = 5000L)
add("f_cdf_checksum", function() wcheck(pf(distribution_x, f_df1, f_df2)), reps = 5000L)
add("f_log_survival_checksum", function() {
  wcheck(pf(distribution_x, f_df1, f_df2, lower.tail = FALSE, log.p = TRUE))
}, reps = 5000L)
add("f_quantile_checksum", function() wcheck(qf(distribution_probs, f_df1, f_df2)), reps = 100L)
add("quantile_type7", function() quantile(quantile_x, 0.37, type = 7, names = FALSE))
add("median", function() median(quantile_x))
add("mad", function() mad(quantile_x, constant = 1.4826))
add("regularized_gamma_p_checksum", function() wcheck(pgamma(gamma_x, gamma_shapes)))
add("regularized_gamma_q_checksum", function() wcheck(pgamma(gamma_x, gamma_shapes, lower.tail = FALSE)))
add("regularized_beta_checksum", function() wcheck(pbeta(beta_x, beta_shape1, beta_shape2)))
add("stable_order_checksum", function() wcheck(order(order_x, method = "radix")), reps = 5000L)
add("rolling_variance_valid", function() {
  value <- vapply(seq_len(length(x) - 5 + 1),
                  function(i) var(x[i:(i + 4)]), numeric(1))
  wcheck(value)
}, reps = 1000L)
add("rolling_sd_right", function() {
  value <- vapply(seq_len(length(x) - 5 + 1),
                  function(i) sd(x[i:(i + 4)]), numeric(1))
  sum(value * (5:100))
}, reps = 1000L)
add("rolling_sum_valid", function() {
  value <- vapply(seq_len(length(x) - 5 + 1),
                  function(i) sum(x[i:(i + 4)]), numeric(1))
  wcheck(value)
}, reps = 1000L)
add("rolling_min_right", function() {
  value <- vapply(seq_len(length(x) - 5 + 1),
                  function(i) min(x[i:(i + 4)]), numeric(1))
  sum(value * (5:100))
}, reps = 1000L)
add("rolling_max_valid", function() {
  value <- vapply(seq_len(length(x) - 5 + 1),
                  function(i) max(x[i:(i + 4)]), numeric(1))
  wcheck(value)
}, reps = 1000L)
add("covariance", function() cov(x, pair_y))
add("correlation", function() cor(x, pair_y))
add("weighted_mean", function() weighted.mean(x, weights), reps = 5000L)
add("weighted_variance_ml", function() {
  cov.wt(cbind(x), wt = weights, method = "ML")$cov[1, 1]
}, reps = 5000L)
add("weighted_covariance_unbiased", function() {
  cov.wt(cbind(x, pair_y), wt = weights, method = "unbiased")$cov[1, 2]
}, reps = 5000L)
add("weighted_correlation", function() {
  value <- cov.wt(cbind(x, pair_y), wt = weights, method = "ML")$cov
  value[1, 2] / sqrt(value[1, 1] * value[2, 2])
}, reps = 5000L)
add("weighted_quantile_ecdf", function() {
  wcheck(weighted_quantile_ecdf(quantile_x, quantile_weights, quantile_probs))
}, reps = 5000L)
add("weighted_quantile_linear_cdf", function() {
  wcheck(weighted_quantile_linear_cdf(quantile_x, quantile_weights, quantile_probs))
}, reps = 5000L)
add("weighted_quantile_frequency_type7", function() {
  wcheck(weighted_quantile_frequency_type7(quantile_x, quantile_weights, quantile_probs))
}, reps = 5000L)
add("weighted_quantile_isotone", function() {
  wcheck(weighted_quantile_isotone(quantile_x, quantile_weights, quantile_probs))
}, reps = 5000L)
add("weighted_quantile_survey_rules", function() {
  sum(seq_len(12) * vapply(seq_len(12), function(rule) {
    weighted_quantile_survey(quantile_x, quantile_weights, 0.37, rule)
  }, numeric(1)))
}, reps = 5000L)
add("rolling_covariance_valid", function() {
  value <- vapply(seq_len(length(x) - 5 + 1),
                  function(i) cov(x[i:(i + 4)], pair_y[i:(i + 4)]), numeric(1))
  wcheck(value)
}, reps = 1000L)
add("rolling_correlation_right", function() {
  value <- vapply(seq_len(length(x) - 5 + 1),
                  function(i) cor(x[i:(i + 4)], pair_y[i:(i + 4)]), numeric(1))
  sum(value * (5:100))
}, reps = 1000L)
add("log_factorial_checksum", function() wcheck(lfactorial(factorial_points)))
add("log_choose_checksum", function() wcheck(lchoose(choose_n, choose_k)))

write.csv(do.call(rbind, rows), out, row.names = FALSE, quote = FALSE)

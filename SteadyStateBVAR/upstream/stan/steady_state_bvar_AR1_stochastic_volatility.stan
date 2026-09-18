functions {
  matrix kron(matrix A, matrix B) { //kronecker product, stan  does not have built in
  matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
  int m;
  int n;
  int p;
  int q;
  m = rows(A);
  n = cols(A);
  p = rows(B);
  q = cols(B);
  for (i in 1:m) {
    for (j in 1:n) {
      int row_start;
      int row_end;
      int col_start;
      int col_end;
      row_start = (i - 1) * p + 1;
      row_end = (i - 1) * p + p;
      col_start = (j - 1) * q + 1;
      col_end = (j - 1) * q + q;
      C[row_start:row_end, col_start:col_end] = A[i, j] * B;
    }
  }
  return C;
  }
}

data {
  int<lower=1> p; //lag order
  int<lower=p+1> N; //number of observations
  int<lower=2> k; //number of endogenous variables
  int<lower=1> q; //number of deterministic variables
  matrix[N, k] Y; //endogenous variables (y_t)
  matrix[N, q] D; //deterministic variables (d_t)
  matrix[N, k*p] W; //lagged endogenous variables
  matrix[N, q*p] Q; //lagged deterministic variables
  vector[k*p*k] theta_beta; //vec_beta prior mean
  matrix[k*p*k, k*p*k] Omega_beta; //vec_beta prior covariance matrix
  vector[k*q] theta_Psi; //vec_Psi prior mean
  matrix[k*q, k*q] Omega_Psi; //vec_Psi prior covariance matrix
  vector[k*(k-1)/2] theta_A; //a prior mean
  matrix[k*(k-1)/2, k*(k-1)/2] Omega_A; //a prior covariance matrix
  vector[k] theta_gamma_0; //gamma_0 prior mean
  matrix[k, k] Omega_gamma_0; //gamma_0 prior covariance matrix
  vector[k] theta_gamma_1; //gamma_1 prior mean
  matrix[k, k] Omega_gamma_1; //gamma_1 prior covariance matrix
  vector[k] theta_log_lambda_1; //log lambda initial condition prior mean
  matrix[k, k] Omega_log_lambda_1; //log lambda initial condition prior covariance matrix
  int<lower=k> m_Phi; //Phi prior degrees of freedom
  matrix[k,k] V_Phi; //Phi prior scale matrix
  int<lower=1> H; // Forecast horizon
  matrix[H, q] d_pred; //future deterministic variables
}

transformed data {
    matrix[p, p] I_p = diag_matrix(rep_vector(1, p)); // Identity matrix
    vector[k] Omega_log_lambda_1_sqrt_diag = sqrt(diagonal(Omega_log_lambda_1));
    vector[k*p*k] Omega_beta_sqrt_diag = sqrt(diagonal(Omega_beta));
    vector[k*q] Omega_psi_sqrt_diag = sqrt(diagonal(Omega_Psi));
    vector[k*(k-1)/2] Omega_A_sqrt_diag = sqrt(diagonal(Omega_A));
    vector[k] Omega_gamma_0_sqrt_diag = sqrt(diagonal(Omega_gamma_0));
    vector[k] Omega_gamma_1_sqrt_diag = sqrt(diagonal(Omega_gamma_1));
}

parameters {
  matrix[k*p, k] beta; //beta' = (Pi_1,...,Pi_p)
  matrix[k, q] Psi; //Psi * d_t = steady state
  matrix[N, k] z; //standard normals
  vector[k] gamma_0; //log volatility intercept
  vector<lower=-1, upper=1>[k] gamma_1; //log volatility slope
  cov_matrix[k] Phi; //log volatility innovation covariance matrix
  vector[k*(k-1)/2] a; //free parameters in A
}

transformed parameters {
  matrix[k,k] A;
  matrix[k,k] Ainv;
  array[N] matrix[k,k] Sigma_u; //time varying covariance matrix of reduced form errors u_t
  array[N] matrix[k,k] L_Sigma_u;
  matrix[N, k] log_lambda; //log volatilities
  matrix[k,k] L_Phi;

  // construct A (1's on diagonal, and then free parameters on lower triangular)
  A = diag_matrix(rep_vector(1,k));
  {
    int idx = 1;
    for (i in 2:k) {
      for (j in 1:(i-1)) {
        A[i,j] = a[idx];
        idx += 1;
      }
    }
  }
  Ainv = inverse(A);
  L_Phi = cholesky_decompose(Phi);
  
  log_lambda[1] = (theta_log_lambda_1 + Omega_log_lambda_1_sqrt_diag .* z[1]')';
  for (t in 2:N) {
    log_lambda[t] = (gamma_0 + gamma_1 .* log_lambda[t-1]' + L_Phi * z[t]')';
  }
  
  
  for (t in 1:N) {
    matrix[k,k] sqrt_Lambda_t = diag_matrix(sqrt(exp(log_lambda[t]')));
    L_Sigma_u[t] = Ainv * sqrt_Lambda_t;
    Sigma_u[t] = tcrossprod(L_Sigma_u[t]);
  }
}

model {
  to_vector(z) ~ std_normal();

  matrix[p*q, p*k] Kron_I_p_Psi = kron(I_p, Psi');
  for(t in 1:N){
      vector[k] u_t = (Y[t] - (D[t]*Psi' + (W[t]-Q[t]*Kron_I_p_Psi)*beta))';
      u_t ~ multi_normal_cholesky(rep_vector(0,k), L_Sigma_u[t]);
  }

  to_vector(beta) ~ normal(theta_beta, Omega_beta_sqrt_diag);
  to_vector(Psi)  ~ normal(theta_Psi, Omega_psi_sqrt_diag);
  a               ~ normal(theta_A, Omega_A_sqrt_diag);
  gamma_0         ~ normal(theta_gamma_0, Omega_gamma_0_sqrt_diag);
  gamma_1         ~ normal(theta_gamma_1, Omega_gamma_1_sqrt_diag);
  Phi             ~ inv_wishart(m_Phi, V_Phi);
}

generated quantities {
  
  matrix[N, k] mu = D * Psi';
  matrix[H, k] mu_pred = d_pred * Psi';

  array[p] matrix[k, k] Pi;
  for (i in 1:p) {
    Pi[i] = (beta[((i - 1) * k + 1):(i * k), :])'; //extract Pi_1, ..., Pi_p
  }
  
  matrix[H,k] y_pred;
  matrix[H,k] log_lambda_pred;
  array[H] matrix[k,k] Sigma_u_pred;
  vector[k] nu;
  
  nu = multi_normal_cholesky_rng(rep_vector(0,k), L_Phi);
  
  for (i in 1:k) {
    log_lambda_pred[1,i] = (gamma_0[i] + gamma_1[i] * log_lambda[N, i] + nu[i]);
    }
  
  for (h in 2:H) {
    nu = multi_normal_cholesky_rng(rep_vector(0,k), L_Phi);
    for (i in 1:k) {
      log_lambda_pred[h,i] = (gamma_0[i] + gamma_1[i] * log_lambda_pred[h-1, i] + nu[i]);
      }
    }

  
  for (h in 1:H) {
    
    matrix[k,k] sqrt_Lambda;
    sqrt_Lambda = diag_matrix(sqrt(exp(log_lambda_pred[h]')));
    matrix[k,k] L_Sigma_u_pred_h = Ainv * sqrt_Lambda;
    
    vector[k] epsilon;
    for (i in 1:k) {
      epsilon[i] = std_normal_rng();
    }
    vector[k] u_t = L_Sigma_u_pred_h * epsilon;
    
    Sigma_u_pred[h] = tcrossprod(L_Sigma_u_pred_h);
    
    vector[k] yhat_t = (d_pred[h]*Psi')';

    if (h > 1) {
      for (i in 1:min(h-1, p)) {
        yhat_t += to_vector((y_pred[h-i] - d_pred[h-i]*Psi') * Pi[i]');
      }
    }

    if (h <= p) {
      for (i in h:p) {
        yhat_t += to_vector((Y[N + h - i] - D[N + h - i]*Psi') * Pi[i]');
      }
    }
    y_pred[h] = (yhat_t + u_t)';
  }
}

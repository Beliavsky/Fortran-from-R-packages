#include <RcppArmadillo.h>
#include "m_exp.h"
#include "EM_LL_UNI.h"
// [[Rcpp::depends(RcppArmadillo)]]


//' EM step for the mPH class with right-censoring, for different marginal 
//'  sub-intensity matrices
//'
//' @param alpha Common initial distribution vector.
//' @param S_list List of marginal sub-intensity matrices.
//' @param y Matrix of marginal observations.
//' @param delta Matrix with right-censoring indications (1 uncensored, 0 right-censored).
//' @param h Tolerance of uniformization.
//'
// [[Rcpp::export]]
void EM_step_mPH_rc(arma::vec & alpha, Rcpp::List & S_list, const arma::mat y, const arma::mat delta, double h) {
  unsigned p{alpha.size()}; // dimension of the distribution
  unsigned n{y.n_rows}; // number of observaions
  
  unsigned d{y.n_cols}; // number of uncensored marginals
  
  arma::vec e(p); 
  e.ones(); // vector of ones
  
  arma::field<arma::cube> matrix_integrals(d,p);
  
  arma::field<arma::cube> list_exp(d);
  arma::cube matrix_exp(p,p,n);
  
  arma::mat J(2 * p, 2 * p); // Van Loan matrix
  arma::mat aux_exp(p,p); //for extracting exp(Ti*xi) from uniformization
  arma::mat aux_int(p,p); //for extracting the matrix integrals from uniformization
  std::vector<arma::mat> aux_vect;
  
  // Compute all matrix exponentials necessary for the E-step
  for (int i{0}; i < d; ++i) {
    arma::mat S = S_list[i];
    arma::vec s = -S * e;
    
    matrix_exp.zeros();
    aux_exp.zeros();
    aux_int.zeros();
    
    for (int j{0}; j < p; ++j) {
      arma::rowvec e_k(p); 
      e_k.zeros();
      e_k(j) = 1;
      
      arma::vec obs = y.col(i);
      arma::cube m_exp(p, p, n); 
      m_exp.zeros();
      
      for(int m{0}; m < n; ++m){
        int rc{static_cast<int>(delta(m,i))};
        if (rc == 1) {
          J = matrix_vanloan(S, S, s * e_k); //Van Loan for uncensored case
        } else {
          J = matrix_vanloan(S, S, e * e_k); //Van Loan for right-censored case
        }
        aux_vect.clear();
        
        double x{obs(m)};
        
        //Van Loan approach with Uniformization
        double a = max_diagonal(J * (-1));
        
        int M{find_n(h, 1)};
        
        vector_of_matrices(aux_vect, J, a, M);
        
        if (x * a <= 1.0) {
          J = m_exp_sum(x, M, aux_vect, a);
        }
        else {
          int n{};
          n = std::log(a * x) / std::log(2.0);
          ++n;
          
          J = m_exp_sum(x / pow(2.0, n), M, aux_vect, a);
          
          pow2_matrix(n, J);
        }
        
        //Extract uniformization results
        for (int k{0}; k < p; ++k) {
          for (int l{0}; l < p; ++l) {
            aux_exp(k,l) = J(k,l); //exp(Si*x)
            aux_int(k,l) = J(k,l + p);
          }
        }
        m_exp.slice(m) = aux_int;
        matrix_exp.slice(m) = aux_exp;
      }
      
      matrix_integrals(i,j) = m_exp;
    }
    list_exp(i) = matrix_exp;
  }
  
  //////////////
  std::vector<std::vector<std::vector<std::vector<double>>>> a_kij(n,std::vector<std::vector<std::vector<double>>>(p, std::vector<std::vector<double>>(d, std::vector<double>(p))));
  //double a_kij [n][p][d][p];
  for (int k{0}; k < p; ++k) {
    for (int j{0}; j < p; ++j) {
      for (int i{0}; i < d; ++i) {
        for (int m{0}; m < n; ++m) {
          arma::cube step_cube = list_exp(i);
          
          a_kij[m][k][i][j] = step_cube(k,j,m);
          // a_kij[m][k][i][j]=matrix_integrals[[i]][[1]][k,j,m] // check if it works
        }
      }
    }
  }
  
  //////////////
  arma::cube a_ki(n,p,d);
  for (int k{0}; k < p; ++k) {
    for (int i{0}; i < d; ++i) {
      arma::mat S = S_list[i];
      arma::vec ti = -S * e;
      
      for (int m{0}; m < n; ++m) {
        //a_kij[m][k][i][]
        arma::vec step_vec(p); 
        step_vec.zeros();
        for (int j{0}; j < p; ++j) {
          step_vec(j) = a_kij[m][k][i][j];
        }
        int rc{static_cast<int>(delta(m,i))};
        
        if (rc == 1) {
          a_ki(m,k,i) = arma::sum(ti % step_vec);
        } else{
          a_ki(m,k,i) = arma::sum(step_vec);
        }
        // a_ki(m,k,i)= arma::sum(s % a_kij[m][k][i][]); //a_kij[m][k][i][] is a vector
      }
    }
  }
  
  //////////////
  arma::cube a_k_minus_i(n, p, d);
  for (int k{0}; k < p; ++k) {
    for (int i{0}; i < d; ++i) {
      for (int m{0}; m < n; ++m) {
        arma::cube copy = a_ki; 
        copy.shed_slice(i);
        arma::mat step_mat = copy.row(m);
        
        if (d > 2) {
          arma::rowvec step_vec = step_mat.row(k);
          a_k_minus_i(m,k,i) = arma::prod(step_vec);
        } else {
          double step_vec = step_mat[k];
          a_k_minus_i(m,k,i) = step_vec;
        }
        // a_k_minus_i(m,k,i)=arma::prod(step_vec);
      }
    }
  }
  
  //////////////
  arma::mat a_k(n,p);
  for (int k{0}; k < p; ++k) {
    for (int m{0}; m < n; ++m) {
      arma::mat step_mat = a_ki.row(m);
      arma::rowvec step_vec = step_mat.row(k);
      
      a_k(m,k) = arma::prod(step_vec);
      //a_k[m, k] <- prod(a_ki[m, k, ])
    }
  }
  
  //////////////
  arma::cube a_tilde_ki(n,p,d);
  for (int k{0}; k < p; ++k) {
    for (int i{0}; i < d; ++i) {
      for (int m{0}; m < n; ++m) {
        arma::mat step_mat = a_k_minus_i.slice(i);
        arma::rowvec step_rowvec = step_mat.row(m);
        
        arma::vec step_vec(p); 
        step_vec.zeros();
        for (int j{0}; j < p; ++j){
          step_vec(j) = a_kij[m][j][i][k];
        }
        int rc{static_cast<int>(delta(m,i))};
        if (rc == 1) {
          a_tilde_ki(m,k,i) = arma::sum(alpha % step_vec % step_rowvec.t());
        } else {
          a_tilde_ki(m,k,i) = 0;
        }
        
        // a_tilde_ki(m,k,i)=arma::sum(alpha % a_kij[m][][i][k] % a_k_minus_i(m, , i));
      }
    }
  }
  
  ////////////
  arma::vec a(n);
  for (int m{0}; m < n; ++m) {
    a(m) = arma::sum(alpha.t() % a_k.row(m));
  }
  
  //////////////
  std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> b_skij(n,std::vector<std::vector<std::vector<std::vector<double>>>>(p, std::vector<std::vector<std::vector<double>>>(p, std::vector<std::vector<double>>(d, std::vector<double>(p)))));
  //double b_skij[n][p][p][d][p];
  for (int s{0}; s < p; ++s) {
    for (int k{0}; k < p; ++k) {
      for (int j{0}; j < p; ++j) {
        for (int i{0}; i < d; ++i) {
          for (int m{0}; m < n; ++m) {
            arma::cube step_cube = matrix_integrals(i,j);
            
            b_skij[m][s][k][i][j] = step_cube(s,k,m);
            
            // b_skij[m][s][k][i][j]=matrix_integrals[[i]][[j]][s,k+p,m];
          }
        }
      }
    }
  }
  
  //////////////
  std::vector<std::vector<std::vector<std::vector<double>>>> b_ski(n,std::vector<std::vector<std::vector<double>>>(p, std::vector<std::vector<double>>(p, std::vector<double>(d))));
  //double b_ski[n][p][p][d];
  for (int s{0}; s < p; ++s) {
    for (int k{0}; k < p; ++k) {
      for (int i{0}; i < d; ++i) {
        for (int m{0}; m < n; ++m) {
          //a_k_minus_i(m,,i)
          arma::mat step_mat = a_k_minus_i.slice(i);
          arma::rowvec step_rowvec = step_mat.row(m);
          
          //b_skij[m][s][k][i][]
          arma::vec step_vec(p); 
          step_vec.zeros();
          for (int j{0}; j < p; ++j) {
            step_vec(j) = b_skij[m][s][k][i][j];
          }
          
          b_ski[m][s][k][i] = arma::sum(alpha % step_rowvec.t() % step_vec);
          
          // b_ski[m][s][k][i][j]=arma::sum(alpha % a_k_minus_i(m,,i) % b_skij[m][s][k][i][]);
          
        }
      }
    }
  }
  
  //////////////////////////
  // E step
  arma::vec EB_k(p);
  
  for (int k{0}; k < p; ++k) {
    EB_k(k) = alpha(k) * arma::sum(a_k.col(k) / a);
  }
  
  arma::mat EZ_ki(p,d);
  for (int k{0}; k < p; ++k) {
    for (int i{0}; i < d; ++i) {
      
      //b_ski[][k][k][i]
      arma::vec step_vec(n); step_vec.zeros();
      for (int m{0}; m < n; ++m) {
        step_vec(m) = b_ski[m][k][k][i];
      }
      
      EZ_ki(k,i) = arma::sum(step_vec / a);
      
      // EZ_ki(k,i)=arma::sum(b_ski[][k][k][i]/a);
    }
  }
  
  
  arma::cube EN_ksi(p,p,d);
  for (int i{0}; i < d; ++i) {
    for (int s{0}; s < p; ++s) {
      for (int k{0}; k < p; ++k) {
        arma::mat step_S = S_list[i];
        
        double t_ksi = step_S(k,s);
        
        //b_ski[][s][k][i]
        arma::vec step_vec(n); 
        step_vec.zeros();
        for (int m{0}; m < n; ++m) {
          step_vec(m) = b_ski[m][s][k][i];
        }
        
        EN_ksi(k,s,i) = t_ksi * arma::sum(step_vec / a);
        // EN_ksi(k,s,i)=t_ksi*arma::sum(b_ski[][s][k][i]/a);
      }
    }
  }
  
  arma::mat EN_ki(p,d);
  for (int i{0}; i < d; ++i) {
    for (int k{0}; k < p; ++k) {
      arma::mat S = S_list[i];
      arma::vec t_ki = -S * e;
      
      arma::mat step_mat = a_tilde_ki.slice(i);
      arma::vec step_vec = step_mat.col(k);
      
      EN_ki(k,i) = t_ki(k) * arma::sum(step_vec / a);
      // EN_ki(k,i)=t_ki(k)*arma::sum(a_tilde_ki(,k,i)/a);
    }
  }
  
  //////////////////////////
  // M step
  for (int k{0}; k < p; ++k) {
    alpha(k) = EB_k(k) / n;
  }
  
  arma::cube T(p,p,d);
  for (int i{0}; i < d; ++i) {
    for (int s{0}; s < p; ++s) {
      for (int k{0}; k < p; ++k) {
        T(k,s,i) = EN_ksi(k, s, i) / EZ_ki(k,i);
      }
    }
  }
  
  arma::mat t(p,d);
  for (int i{0}; i < d; ++i) {
    for (int k{0}; k < p; ++k) {
      t(k,i) = EN_ki(k,i) / EZ_ki(k,i);
    }
  }
  
  for (int k{0}; k < p; ++k) {
    for (int i{0}; i < d; ++i) {
      arma::mat step_mat = T.slice(i);
      arma::mat copy = step_mat; 
      copy.shed_col(k);
      arma::rowvec step_vec = copy.row(k);
      
      T(k,k,i) = -arma::sum(step_vec) - t(k,i);
    }
  }
  
  for (int i{0}; i < d; ++i) {
    S_list[i] = T.slice(i);
  }
}

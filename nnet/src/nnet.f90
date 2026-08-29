! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module nnet
use r_compat, only: dp
use nnet_types, only: nnet_model_t, multinom_model_t
use nnet_core, only: build_network, nnet_weight_count, nnet_objective, nnet_gradient, nnet_objective_gradient, &
   nnet_hessian_exact, nnet_predict_raw
use nnet_fit_mod, only: nnet_fit, nnet_refit, nnet_predict, nnet_predict_class
use nnet_multinom, only: multinom_fit_labels, multinom_fit_counts, multinom_predict_proba, multinom_predict_class, &
   multinom_information, multinom_covariance, multinom_loglik
use nnet_utils, only: class_ind, which_is_max, which_is_max_rows, summarize_rows
implicit none
public
end module nnet

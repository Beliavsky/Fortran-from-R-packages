module polca
   use polca_kinds, only : dp
   use polca_types, only : polca_model, polca_simulation
   use polca_core, only : polca_fit, polca_item_likelihood, polca_postclass, polca_probhat, &
        polca_update_prior, polca_beta_derivatives, polca_standard_errors
   use polca_api, only : polca_posterior, polca_predcell, polca_entropy, polca_reorder_probs, &
        polca_table_1d, polca_table_2d, polca_simdata, polca_rmulti, polca_coef, polca_vcov
   implicit none
   public
end module polca

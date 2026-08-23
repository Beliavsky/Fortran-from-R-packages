module mixspe
  use mvtnorm_kinds, only : dp
  use mixspe_distributions, only : dpe, log_dpe, rpe, dspe, log_dspe, rspe, cov_pe
  use mixspe_mixture, only : spe_model, em_fit, emgr_fit, model_num_parameters, map_labels
  implicit none
  public
end module

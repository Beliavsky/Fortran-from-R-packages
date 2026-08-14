module kriginv
  use kriginv_kinds, only : dp
  use kriginv_model, only : krig_model, krig_prediction, km_control, scaling_axis, init_krig_model, &
                            init_dice_krig_model, fit_krig_model, update_krig_model, predict_nobias_km, &
                            covariance_matrix, covariance_cross, posterior_covariance
  use kriginv_criteria, only : excursion_probability, vorob_threshold, ranjan_optim, bichon_optim, tmse_optim, &
                               tsee_optim, predict_update_km_parallel, sur_optim_parallel, jn_optim_parallel, &
                               timse_optim_parallel, vorob_optim_parallel, vorobvol_optim_parallel, &
                               sur_optim_parallel2, jn_optim_parallel2, timse_optim_parallel2, &
                               vorob_optim_parallel2, vorobvol_optim_parallel2, compute_real_volume_constant
  use kriginv_update, only : precomputed_update_data, precompute_update_data, compute_quick_krigcov
  use kriginv_integration, only : integration_control, integration_result, integration_design
  use kriginv_optimize, only : optimizer_control
  use kriginv_maximize, only : optimization_result, max_infill_criterion, max_sur_parallel, max_timse_parallel, &
                               max_vorob_parallel, max_futurevol_parallel
  use kriginv_egi, only : egi_result, egi, egi_parallel
  implicit none
  public
end module kriginv

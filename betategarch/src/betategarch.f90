! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2011-2025 Genaro Sucarrat
! Copyright (C) 2026 contributors to the Modern Fortran translation
!
! This file is part of betategarch-modern-fortran.
! It is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License version 2 only.

module betategarch
  use betategarch_kinds, only : dp
  use betategarch_rng, only : set_random_seed
  use skew_t_mod, only : skew_t_random, skew_t_pdf, skew_t_logpdf, skew_t_mean, &
    skew_t_variance, skew_t_skewness, skew_t_kurtosis, skew_t_raw_moment
  use tegarch_mod, only : tegarch_parameters, tegarch_filter_result, tegarch_fit_result, &
    tegarch_filter, tegarch_simulate, tegarch_loglik, tegarch_fit, tegarch_forecast, &
    tegarch_default_parameters, tegarch_default_bounds, tegarch_free_parameter_count, &
    tegarch_params_to_free, tegarch_free_to_params, tegarch_bic_per_observation, &
    tegarch_standard_errors
  implicit none
  public
end module betategarch

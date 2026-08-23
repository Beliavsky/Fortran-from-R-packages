! SPDX-License-Identifier: GPL-3.0-or-later
module sgt
  use sgt_kinds, only : dp
  use sgt_distribution, only : dsgt, psgt, qsgt, rsgt, sgt_logpdf, sgt_pdf, &
    sgt_cdf, sgt_quantile, sgt_mean_shift, sgt_variance_scale
  use sgt_mle_mod, only : sgt_params, sgt_mle_result, sgt_observation_model, &
    sgt_mle_model, sgt_mle_constant
  implicit none
  public
end module sgt

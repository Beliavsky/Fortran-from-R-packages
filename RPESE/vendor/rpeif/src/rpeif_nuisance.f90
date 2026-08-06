! SPDX-License-Identifier: GPL-3.0-or-later
module rpeif_nuisance
  use rpeif_kinds, only : dp
  use rpeif_types, only : nuisance_parameters, rpeif_success, rpeif_invalid_argument, &
    rpeif_numerical_failure
  use rpeif_stats, only : normal_pdf, normal_cdf, normal_quantile
  implicit none
  private
  public :: nuisance_parameters_fn
contains
  subroutine nuisance_parameters_fn(pars, mu, sd, threshold, alpha, beta, status)
    type(nuisance_parameters), intent(out) :: pars
    real(dp), intent(in), optional :: mu, sd, threshold, alpha, beta
    integer, intent(out), optional :: status
    real(dp) :: d, z_alpha, z_beta, denominator
    integer :: stat

    stat = rpeif_success
    if (present(mu)) pars%mu = mu
    if (present(sd)) pars%sd = sd
    if (present(threshold)) pars%c = threshold
    if (present(alpha)) pars%alpha = alpha
    if (present(beta)) pars%beta = beta

    if (pars%sd <= 0.0_dp .or. pars%alpha <= 0.0_dp .or. pars%alpha >= 1.0_dp .or. &
        pars%beta <= 0.0_dp .or. pars%beta >= 1.0_dp) then
      stat = rpeif_invalid_argument
      if (present(status)) status = stat
      return
    end if

    z_alpha = normal_quantile(pars%alpha)
    pars%q_alpha = pars%mu + z_alpha * pars%sd
    pars%es_alpha = -pars%mu + normal_pdf(z_alpha) * pars%sd / pars%alpha

    d = (pars%c - pars%mu) / pars%sd
    pars%lpm1 = (d * normal_cdf(d) + normal_pdf(d)) * pars%sd
    pars%lpm2 = ((d * d + 1.0_dp) * normal_cdf(d) + d * normal_pdf(d)) * pars%sd ** 2

    pars%semisd = pars%sd / sqrt(2.0_dp)
    pars%semimean = -normal_pdf(0.0_dp) * pars%sd
    pars%fq_alpha = normal_pdf(pars%q_alpha, pars%mu, pars%sd)

    denominator = max(abs(pars%semisd), tiny(1.0_dp))
    pars%dsr = pars%mu / denominator
    denominator = max(abs(pars%es_alpha), tiny(1.0_dp))
    pars%es_ratio = pars%mu / denominator
    pars%upm1 = pars%lpm1 + pars%mu - pars%c
    denominator = max(abs(pars%lpm1), tiny(1.0_dp))
    pars%omega = pars%upm1 / denominator

    z_beta = normal_quantile(1.0_dp - pars%beta)
    pars%q_beta = pars%mu + z_beta * pars%sd
    pars%eg_beta = pars%mu + normal_pdf(z_beta) * pars%sd / pars%beta
    denominator = max(abs(pars%es_alpha), tiny(1.0_dp))
    pars%rachev_ratio = pars%eg_beta / denominator

    denominator = max(sqrt(max(pars%lpm2, 0.0_dp)), tiny(1.0_dp))
    pars%sor_c = pars%mu / denominator
    denominator = max(abs(pars%semisd), tiny(1.0_dp))
    pars%sor_mu = pars%mu / denominator
    pars%sr = pars%mu / pars%sd
    if (abs(pars%q_alpha) <= tiny(1.0_dp)) then
      pars%var_ratio = 0.0_dp
      stat = rpeif_numerical_failure
    else
      pars%var_ratio = -pars%mu / pars%q_alpha
    end if

    if (present(status)) status = stat
  end subroutine nuisance_parameters_fn
end module rpeif_nuisance

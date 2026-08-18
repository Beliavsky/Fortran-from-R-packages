! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of RMKdiscrete 0.1 by Robert M. Kirkpatrick.
! See LICENSE, COPYING, and TRANSLATION_NOTES.md for provenance.
module rmkdiscrete
  use rmkdiscrete_kinds, only : dp
  use rmkdiscrete_lgp, only : lgp_summary, lgp_from_mu_sigma2, lgp_from_theta_lambda, &
    lgp_from_mu_theta, lgp_from_sigma2_lambda, lgp_from_sigma2_theta, lgp_from_mu_lambda, &
    lgp_findmax, lgp_get_nc, dlgp, plgp, qlgp, rlgp, rlgp_sample, slgp
  use rmkdiscrete_negbin, only : dnegbin, negbin_from_nu_p, negbin_from_mu_sigma2, &
    negbin_from_mu_nu, negbin_from_mu_p, negbin_from_sigma2_p, rnegbin, rnegbin_sample
  use rmkdiscrete_bivariate, only : log_moments2, dbilgp, rbilgp, rbilgp_sample, bilgp_logmv, &
    dbinegbin, rbinegbin, rbinegbin_sample, binegbin_logmv
  use rmkdiscrete_manaclash, only : dmanaclash_dmg, dmanaclash_xyn, dmanaclash_net, rmanaclash, rmanaclash_sample
  implicit none
  public
end module rmkdiscrete

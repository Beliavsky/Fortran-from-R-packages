! Umbrella module for the gamlss.dist computational translation.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_dist
   use gamlss_kinds
   use gamlss_links
   use gamlss_base, dzipf_base => dzipf, pzipf_base => pzipf, &
      qzipf_base => qzipf, rzipf_base => rzipf, &
      dpareto1_base => dpareto1, ppareto1_base => ppareto1, qpareto1_base => qpareto1, &
      dpareto2_base => dpareto2, ppareto2_base => ppareto2, qpareto2_base => qpareto2
   use gamlss_student_t
   use gamlss_continuous
   use gamlss_discrete
   use gamlss_boxcox
   use gamlss_continuous_v02
   use gamlss_discrete_v02
   use gamlss_continuous_v03
   use gamlss_discrete_v03
   use gamlss_flexible_v03
   use gamlss_fit
   use gamlss_fit_v03
   implicit none
   public
end module gamlss_dist

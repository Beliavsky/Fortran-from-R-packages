! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from Risk 1.0 by Saralees Nadarajah and Stephen Chan.
! Copyright (c) 2017 Saralees Nadarajah and Stephen Chan.
module risk
   use risk_kinds, only : dp
   use risk_distributions, only : continuous_distribution, callback_distribution, &
      normal_distribution, lognormal_distribution, uniform_distribution, &
      exponential_distribution, logistic_distribution, student_t_distribution
   use risk_measures
   implicit none
   public
end module risk

! SPDX-License-Identifier: LGPL-2.0-or-later
module survival_types
  use survival_kinds, only : dp
  implicit none
  private
  public :: survfit_result, coxph_result, aft_result, concordance_result, survdiff_result

  type :: survfit_result
     real(dp), allocatable :: time(:), n_risk(:), n_event(:), n_censor(:)
     real(dp), allocatable :: survival(:), cumhaz(:), std_err(:), std_chaz(:)
  end type

  type :: coxph_result
     real(dp), allocatable :: coef(:), var(:,:), score(:), means(:)
     real(dp) :: loglik_initial = 0.0_dp, loglik = 0.0_dp, score_test = 0.0_dp
     integer :: iterations = 0, rank = 0
     logical :: converged = .false.
  end type

  type :: aft_result
     real(dp), allocatable :: coef(:), var(:,:)
     real(dp) :: scale = 1.0_dp, loglik = 0.0_dp
     integer :: iterations = 0
     logical :: converged = .false.
  end type

  type :: concordance_result
     real(dp) :: concordant=0.0_dp, discordant=0.0_dp, tied_risk=0.0_dp, tied_time=0.0_dp
     real(dp) :: cindex=0.0_dp
  end type

  type :: survdiff_result
     real(dp), allocatable :: observed(:), expected(:), variance(:,:)
     real(dp) :: chisq=0.0_dp
  end type
end module survival_types

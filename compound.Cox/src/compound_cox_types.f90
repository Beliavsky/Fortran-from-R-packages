! SPDX-License-Identifier: GPL-2.0-only
module compound_cox_types
  use compound_cox_kinds, only : dp
  implicit none
  type :: cg_result
    real(dp) :: tau=0.0_dp, median=0.0_dp
    real(dp), allocatable :: time(:), n_risk(:), surv(:)
  end type
  type :: univariate_result
    real(dp), allocatable :: beta(:), se(:), z(:), p(:)
  end type
  type :: depend_cox_result
    real(dp) :: beta=0.0_dp, se=0.0_dp, z=0.0_dp, p=1.0_dp
    real(dp) :: beta_censor=0.0_dp, se_censor=0.0_dp, z_censor=0.0_dp, p_censor=1.0_dp
    real(dp), allocatable :: baseline(:), censor_baseline(:)
    logical :: converged=.false.
  end type
  type :: depend_cv_result
    real(dp), allocatable :: beta(:), se(:), z(:), p(:)
    real(dp) :: alpha=0.0_dp, c_index=0.0_dp
  end type
  type :: compound_result
    real(dp) :: a=0.0_dp
    real(dp), allocatable :: beta(:), se(:), lower95(:), upper95(:)
    real(dp), allocatable :: sigma(:,:), v(:,:), h_dot(:)
    real(dp) :: hessian_cv=0.0_dp
    logical :: converged=.false.
  end type
  type :: selection_result
    integer :: n_genes=0, n_selected=0
    real(dp), allocatable :: beta(:), z(:), p(:)
    integer, allocatable :: selected(:)
    real(dp) :: cvl=0.0_dp, rcvl1=0.0_dp, rcvl2=0.0_dp
    real(dp) :: c_index_no_cv=0.0_dp, c_index_incomplete=0.0_dp, c_index_full=0.0_dp
    real(dp) :: false_selected=-1.0_dp, fdr_formula=0.0_dp, fdr_permutation=-1.0_dp
  end type
  type :: cg_test_result
    real(dp) :: survival_diff=0.0_dp, rmstd=0.0_dp, p_value=1.0_dp
    real(dp) :: l1_distance=0.0_dp, integrated_l1=0.0_dp, l1_p_value=1.0_dp
    real(dp) :: tau=0.0_dp
    integer :: n_good=0, n_poor=0, events_good=0, events_poor=0
    real(dp) :: rmst_good=0.0_dp, rmst_poor=0.0_dp, mean_pi_good=0.0_dp, mean_pi_poor=0.0_dp
  end type
  type :: factorial_result
    real(dp), allocatable :: estimate(:), se(:), lower(:), upper(:), p(:), variance(:,:)
    real(dp) :: f_stat=0.0_dp, p_simu=1.0_dp, p_anal=1.0_dp, df=0.0_dp
    real(dp) :: c_simu(3)=0.0_dp, c_anal(3)=0.0_dp
  end type
end module compound_cox_types

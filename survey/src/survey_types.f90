! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_types
  use r_quantiles, only : r_qrule_math, r_qrule_school, r_qrule_shahvaish
  use r_quantiles, only : r_qrule_hf1, r_qrule_hf2, r_qrule_hf3, r_qrule_hf4
  use r_quantiles, only : r_qrule_hf5, r_qrule_hf6, r_qrule_hf7, r_qrule_hf8, r_qrule_hf9
  use survey_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: LONELY_FAIL=1, LONELY_REMOVE=2, LONELY_CERTAINTY=3, &
                                LONELY_ADJUST=4, LONELY_AVERAGE=5
  integer, parameter, public :: FAMILY_GAUSSIAN=1, FAMILY_BINOMIAL=2, FAMILY_POISSON=3
  integer, parameter, public :: LINK_IDENTITY=1, LINK_LOGIT=2, LINK_LOG=3
  integer, parameter, public :: OLR_LOGISTIC=1, OLR_PROBIT=2, OLR_CLOGLOG=3, OLR_CAUCHIT=4
  integer, parameter, public :: FACT_N_NONE=0, FACT_N_SAMPLE=1, FACT_N_DEGF=2, &
                                FACT_N_EFFECTIVE=3, FACT_N_MIN_EFFECTIVE=4
  integer, parameter, public :: QRULE_MATH = r_qrule_math
  integer, parameter, public :: QRULE_SCHOOL = r_qrule_school
  integer, parameter, public :: QRULE_SHAHVAISH = r_qrule_shahvaish
  integer, parameter, public :: QRULE_HF1 = r_qrule_hf1
  integer, parameter, public :: QRULE_HF2 = r_qrule_hf2
  integer, parameter, public :: QRULE_HF3 = r_qrule_hf3
  integer, parameter, public :: QRULE_HF4 = r_qrule_hf4
  integer, parameter, public :: QRULE_HF5 = r_qrule_hf5
  integer, parameter, public :: QRULE_HF6 = r_qrule_hf6
  integer, parameter, public :: QRULE_HF7 = r_qrule_hf7
  integer, parameter, public :: QRULE_HF8 = r_qrule_hf8
  integer, parameter, public :: QRULE_HF9 = r_qrule_hf9

  type, public :: survey_design_t
    integer :: n = 0
    integer :: stages = 0
    logical :: has_strata = .false.
    logical :: ultimate_cluster = .false.
    logical :: adjust_domain_lonely = .false.
    integer :: lonely_psu = LONELY_REMOVE
    real(dp), allocatable :: weight(:)
    integer, allocatable :: cluster(:,:)
    integer, allocatable :: strata(:,:)
    real(dp), allocatable :: samp_size(:,:)
    real(dp), allocatable :: pop_size(:,:)
  end type survey_design_t

  type, public :: rep_design_t
    integer :: n = 0
    integer :: r = 0
    real(dp) :: scale = 1.0_dp
    logical :: mse = .false.
    real(dp), allocatable :: weight(:)
    real(dp), allocatable :: repweights(:,:)
    real(dp), allocatable :: rscales(:)
  end type rep_design_t

  type, public :: svystat_t
    real(dp), allocatable :: estimate(:)
    real(dp), allocatable :: variance(:,:)
    real(dp), allocatable :: influence(:,:)
  end type svystat_t

  type, public :: ratio_result_t
    real(dp), allocatable :: ratio(:,:)
    real(dp), allocatable :: variance(:,:)
    real(dp), allocatable :: influence(:,:)
  end type ratio_result_t

  type, public :: glm_result_t
    real(dp), allocatable :: coef(:)
    real(dp), allocatable :: vcov(:,:)
    real(dp), allocatable :: naive_vcov(:,:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: residual(:)
    real(dp) :: deviance = 0.0_dp
    integer :: iterations = 0
    integer :: rank = 0
    integer :: df_residual = 0
    logical :: converged = .false.
  end type glm_result_t

  type, public :: quantile_result_t
    real(dp), allocatable :: quantile(:)
    real(dp), allocatable :: se(:)
    real(dp), allocatable :: lower(:)
    real(dp), allocatable :: upper(:)
  end type quantile_result_t

  type, public :: survival_curve_t
    real(dp), allocatable :: time(:)
    real(dp), allocatable :: survival(:)
    real(dp), allocatable :: hazard(:)
    real(dp), allocatable :: variance(:)
  end type survival_curve_t


  type, public :: cox_result_t
    real(dp), allocatable :: coef(:)
    real(dp), allocatable :: vcov(:,:)
    real(dp), allocatable :: naive_vcov(:,:)
    real(dp) :: loglik = 0.0_dp
    integer :: iterations = 0
    integer :: rank = 0
    logical :: converged = .false.
  end type cox_result_t

  type, public :: chisq_result_t
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    real(dp) :: numerator_df = 0.0_dp
    real(dp) :: denominator_df = 0.0_dp
    integer :: rank = 0
  end type chisq_result_t

  type, public :: ivreg_result_t
    real(dp), allocatable :: coef(:)
    real(dp), allocatable :: vcov(:,:)
    real(dp), allocatable :: naive_vcov(:,:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: residual(:)
    integer :: rank = 0
    integer :: df_residual = 0
  end type ivreg_result_t

  type, public :: aft_result_t
    real(dp), allocatable :: coef(:)
    real(dp), allocatable :: vcov(:,:)
    real(dp), allocatable :: naive_vcov(:,:)
    real(dp) :: scale = 1.0_dp
    real(dp) :: loglik = 0.0_dp
    integer :: iterations = 0
    logical :: converged = .false.
  end type aft_result_t

  type, public :: mle_result_t
    real(dp), allocatable :: par(:)
    real(dp), allocatable :: vcov(:,:)
    real(dp), allocatable :: model_vcov(:,:)
    real(dp), allocatable :: influence(:,:)
    real(dp) :: loglik = 0.0_dp
    integer :: evaluations = 0
    integer :: status = 0
    logical :: converged = .false.
  end type mle_result_t

  type, public :: nls_result_t
    real(dp), allocatable :: coef(:)
    real(dp), allocatable :: vcov(:,:)
    real(dp), allocatable :: naive_vcov(:,:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: residual(:)
    real(dp) :: rss = 0.0_dp
    integer :: iterations = 0
    logical :: converged = .false.
  end type nls_result_t

  type, public :: pca_result_t
    real(dp), allocatable :: center(:)
    real(dp), allocatable :: scale(:)
    real(dp), allocatable :: rotation(:,:)
    real(dp), allocatable :: sdev(:)
  end type pca_result_t

  type, public :: olr_result_t
    real(dp), allocatable :: coef(:)
    real(dp), allocatable :: zeta(:)
    real(dp), allocatable :: vcov(:,:)
    real(dp), allocatable :: naive_vcov(:,:)
    real(dp), allocatable :: fitted(:,:)
    real(dp) :: deviance = 0.0_dp
    integer :: iterations = 0
    integer :: method = OLR_LOGISTIC
    logical :: converged = .false.
  end type olr_result_t

  type, public :: loglin_result_t
    real(dp) :: intercept = 0.0_dp
    real(dp), allocatable :: coef(:)
    real(dp), allocatable :: vcov(:,:)
    real(dp), allocatable :: fitted_prob(:)
    real(dp), allocatable :: cell_prob(:)
    real(dp), allocatable :: cell_vcov(:,:)
    real(dp) :: deviance = 0.0_dp
    integer :: df_residual = 0
    integer :: df_null = 0
    logical :: converged = .false.
  end type loglin_result_t

  type, public :: loglin_test_t
    real(dp) :: deviance = 0.0_dp
    real(dp) :: score = 0.0_dp
    real(dp) :: p_deviance_satterthwaite = 1.0_dp
    real(dp) :: p_deviance_saddle = 1.0_dp
    real(dp) :: p_score_satterthwaite = 1.0_dp
    real(dp) :: p_score_saddle = 1.0_dp
    real(dp), allocatable :: lambda(:)
    integer :: df = 0
  end type loglin_test_t

  type, public :: factor_result_t
    real(dp), allocatable :: loadings(:,:)
    real(dp), allocatable :: uniqueness(:)
    real(dp), allocatable :: communalities(:)
    real(dp) :: criterion = 0.0_dp
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    real(dp) :: effective_n = 0.0_dp
    integer :: factors = 0
    integer :: dof = 0
    integer :: iterations = 0
    logical :: converged = .false.
  end type factor_result_t


  type, public :: model_test_t
    real(dp) :: statistic = 0.0_dp
    real(dp) :: p_value = 1.0_dp
    integer :: df = 0
    integer :: ddf = 0
    real(dp), allocatable :: lambda(:)
  end type model_test_t

  type, public :: phase_variance_t
    real(dp), allocatable :: variance(:,:)
    real(dp), allocatable :: phase(:,:,:)
  end type phase_variance_t

  type, public :: multiframe_design_t
    integer :: n1 = 0
    integer :: n2 = 0
    real(dp), allocatable :: frame_weight1(:)
    real(dp), allocatable :: frame_weight2(:)
    real(dp), allocatable :: design_weight1(:)
    real(dp), allocatable :: design_weight2(:)
  end type multiframe_design_t

end module survey_types

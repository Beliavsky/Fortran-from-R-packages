module bayesm_types
  use bayesm_kinds, only: dp
  implicit none
  private
  public :: normal_component, normal_mixture, reg_data, mnl_data
  public :: unireg_result, multireg_draw, mixture_step_result, nmix_result
  public :: probit_result, matrix_mcmc_result, hier_linear_result, hier_mnl_result
  public :: negbin_result, iv_result, dp_mixture_result, sur_result, ordprobit_result
  public :: scale_usage_result, blp_result, mnl_metrop_result, hier_negbin_result, wishart_result


  type :: wishart_result
    real(dp), allocatable :: w(:,:)
    real(dp), allocatable :: iw(:,:)
    real(dp), allocatable :: c(:,:)
    real(dp), allocatable :: ci(:,:)
  end type wishart_result

  type :: normal_component
    real(dp), allocatable :: mu(:)
    real(dp), allocatable :: rooti(:,:)
  end type normal_component

  type :: normal_mixture
    real(dp), allocatable :: p(:)
    type(normal_component), allocatable :: comp(:)
  end type normal_mixture

  type :: reg_data
    real(dp), allocatable :: y(:)
    real(dp), allocatable :: x(:,:)
  end type reg_data

  type :: mnl_data
    integer, allocatable :: y(:)
    real(dp), allocatable :: x(:,:)
  end type mnl_data

  type :: unireg_result
    real(dp), allocatable :: betadraw(:,:)
    real(dp), allocatable :: sigmasqdraw(:)
  end type unireg_result

  type :: multireg_draw
    real(dp), allocatable :: b(:,:)
    real(dp), allocatable :: sigma(:,:)
  end type multireg_draw

  type :: mixture_step_result
    real(dp), allocatable :: p(:)
    integer, allocatable :: z(:)
    type(normal_component), allocatable :: comp(:)
  end type mixture_step_result

  type :: nmix_result
    real(dp), allocatable :: probdraw(:,:)
    integer, allocatable :: zdraw(:,:)
    type(normal_mixture), allocatable :: mixdraw(:)
  end type nmix_result

  type :: probit_result
    real(dp), allocatable :: betadraw(:,:)
    real(dp), allocatable :: sigmadraw(:,:,:)
    real(dp), allocatable :: ddraw(:,:)
    real(dp), allocatable :: llike(:)
    real(dp) :: accept = 0.0_dp
  end type probit_result

  type :: matrix_mcmc_result
    real(dp), allocatable :: draw1(:,:)
    real(dp), allocatable :: draw2(:,:)
    real(dp), allocatable :: draw3(:,:)
    real(dp), allocatable :: stat(:)
  end type matrix_mcmc_result

  type :: hier_linear_result
    real(dp), allocatable :: betadraw(:,:,:)
    real(dp), allocatable :: deltadraw(:,:,:)
    real(dp), allocatable :: vbetadraw(:,:,:)
    real(dp), allocatable :: taudraw(:,:)
    real(dp), allocatable :: probdraw(:,:)
    integer, allocatable :: zdraw(:,:)
    type(normal_mixture), allocatable :: mixdraw(:)
  end type hier_linear_result

  type :: hier_mnl_result
    real(dp), allocatable :: betadraw(:,:,:)
    real(dp), allocatable :: deltadraw(:,:,:)
    real(dp), allocatable :: vbetadraw(:,:,:)
    real(dp), allocatable :: probdraw(:,:)
    integer, allocatable :: zdraw(:,:)
    real(dp), allocatable :: llike(:)
    real(dp), allocatable :: reject(:)
  end type hier_mnl_result


  type :: hier_negbin_result
    real(dp), allocatable :: betadraw(:,:,:)
    real(dp), allocatable :: alphadraw(:)
    real(dp), allocatable :: vbetadraw(:,:,:)
    real(dp), allocatable :: deltadraw(:,:,:)
    real(dp), allocatable :: llike(:)
    real(dp) :: beta_accept = 0.0_dp
    real(dp) :: alpha_accept = 0.0_dp
  end type hier_negbin_result

  type :: negbin_result
    real(dp), allocatable :: betadraw(:,:)
    real(dp), allocatable :: alphadraw(:)
    real(dp), allocatable :: llike(:)
    real(dp) :: beta_accept = 0.0_dp
    real(dp) :: alpha_accept = 0.0_dp
  end type negbin_result

  type :: iv_result
    real(dp), allocatable :: betadraw(:,:)
    real(dp), allocatable :: deltadraw(:,:)
    real(dp), allocatable :: sigmadraw(:,:,:)
    real(dp), allocatable :: alphadraw(:)
  end type iv_result

  type :: dp_mixture_result
    real(dp), allocatable :: alphadraw(:)
    integer, allocatable :: zdraw(:,:)
    integer, allocatable :: ncompdraw(:)
    type(normal_mixture), allocatable :: mixdraw(:)
  end type dp_mixture_result

  type :: sur_result
    real(dp), allocatable :: betadraw(:,:)
    real(dp), allocatable :: sigmadraw(:,:,:)
  end type sur_result

  type :: ordprobit_result
    real(dp), allocatable :: betadraw(:,:)
    real(dp), allocatable :: ddraw(:,:)
    real(dp), allocatable :: llike(:)
    real(dp) :: accept = 0.0_dp
  end type ordprobit_result


  type :: mnl_metrop_result
    real(dp), allocatable :: betadraw(:,:)
    real(dp), allocatable :: loglike(:)
    integer :: naccept = 0
  end type mnl_metrop_result

  type :: scale_usage_result
    real(dp), allocatable :: mudraw(:,:)
    real(dp), allocatable :: sigmadraw(:,:,:)
    real(dp), allocatable :: taudraw(:,:)
    real(dp), allocatable :: sdraw(:,:)
    real(dp), allocatable :: lambdadraw(:,:,:)
    real(dp), allocatable :: edraw(:)
  end type scale_usage_result

  type :: blp_result
    real(dp), allocatable :: thetadraw(:,:)
    real(dp), allocatable :: rdraw(:,:)
    real(dp), allocatable :: tausqdraw(:)
    real(dp), allocatable :: omegadraw(:,:,:)
    real(dp), allocatable :: deltadraw(:,:)
    real(dp), allocatable :: llike(:)
    real(dp) :: accept = 0.0_dp
  end type blp_result
end module bayesm_types

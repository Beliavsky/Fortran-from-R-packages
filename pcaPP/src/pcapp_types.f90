module pcapp_types
  use pcapp_kinds, only: dp
  implicit none
  private

  type, public :: scale_result
    real(dp), allocatable :: x(:,:)
    real(dp), allocatable :: center(:)
    real(dp), allocatable :: scale(:)
  end type scale_result

  type, public :: median_result
    real(dp), allocatable :: par(:)
    real(dp) :: value = 0.0_dp
    integer :: code = 0
    integer :: iterations = 0
  end type median_result

  type, public :: pca_result
    real(dp), allocatable :: loadings(:,:)
    real(dp), allocatable :: sdev(:)
    real(dp), allocatable :: scores(:,:)
    real(dp), allocatable :: center(:)
    real(dp), allocatable :: scale(:)
    real(dp), allocatable :: objective(:)
    real(dp), allocatable :: lambda(:)
    integer :: k = 0
    integer :: n_obs = 0
  end type pca_result

  type, public :: covariance_result
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: center(:)
  end type covariance_result

  type, public :: tuning_result
    type(pca_result) :: pc
    real(dp), allocatable :: lambda_grid(:)
    real(dp), allocatable :: criterion(:)
    integer :: best_index = 0
  end type tuning_result
end module pcapp_types

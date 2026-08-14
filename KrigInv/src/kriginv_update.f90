module kriginv_update
  use kriginv_kinds, only : dp
  use kriginv_model, only : krig_model, posterior_covariance
  implicit none
  private
  type, public :: precomputed_update_data
    real(dp), allocatable :: integration_points(:,:)
  end type precomputed_update_data
  public :: precompute_update_data, compute_quick_krigcov
contains
  function precompute_update_data(model,integration_points) result(data)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: integration_points(:,:)
    type(precomputed_update_data) :: data
    if(size(integration_points,2)/=model%d) error stop 'kriginv: integration point dimension mismatch'
    data%integration_points=integration_points
  end function precompute_update_data

  function compute_quick_krigcov(model,integration_points,xnew,precalc) result(kn)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: integration_points(:,:),xnew(:,:)
    type(precomputed_update_data), intent(in), optional :: precalc
    real(dp), allocatable :: kn(:,:)
    if(present(precalc)) then
      if(size(precalc%integration_points,1)/=size(integration_points,1)) &
        error stop 'kriginv: precomputed integration data mismatch'
    end if
    kn=posterior_covariance(model,integration_points,xnew)
  end function compute_quick_krigcov
end module kriginv_update

module mixsqp
  use mixsqp_kinds, only : dp
  use mixsqp_types, only : mixsqp_control, mixsqp_result
  use mixsqp_utils, only : mixobjective_raw => mixobjective, normalize_likelihoods, &
    normalize_loglikelihoods, set_seed
  use mixsqp_em, only : mixem_update
  use mixsqp_solver, only : active_set_qp
  use mixsqp_highlevel, only : fit_mixsqp, mixsqp_default_control
  use mixsqp_simulate, only : simulate_mix_data
  implicit none
  private
  public :: dp, mixsqp_control, mixsqp_result
  public :: fit_mixsqp, mixsqp_default_control, mixobjective
  public :: normalize_likelihoods, normalize_loglikelihoods, set_seed
  public :: mixem_update, active_set_qp, simulate_mix_data
contains
  function mixobjective(L,x,w) result(f)
    real(dp), intent(in) :: L(:,:),x(:)
    real(dp), intent(in), optional :: w(:)
    real(dp) :: f
    real(dp), allocatable :: xn(:),wn(:)
    if (size(x)/=size(L,2) .or. any(x<0.0_dp) .or. sum(x)<=0.0_dp) then
      f=huge(1.0_dp); return
    end if
    allocate(xn(size(x)),wn(size(L,1)))
    xn=x/sum(x)
    if (present(w)) then
      if (size(w)/=size(L,1) .or. any(w<0.0_dp) .or. sum(w)<=0.0_dp) then
        f=huge(1.0_dp); return
      end if
      wn=w/sum(w)
    else
      wn=1.0_dp/real(size(L,1),dp)
    end if
    f=mixobjective_raw(L,xn,wn)
  end function mixobjective
end module mixsqp

! SPDX-License-Identifier: GPL-3.0-only
module svdnf_grids
  use svdnf_kinds, only : dp
  use svdnf_types, only : svm_dynamics, grid_type, model_duffie_pan_singleton, &
    model_bates, model_heston, model_pitt_malik_doucet, model_taylor, &
    model_taylor_leverage, model_capm_sv
  implicit none
  private
  public :: grid_maker, validate_grid

contains

  function grid_maker(dynamics, n, k, r) result(grids)
    type(svm_dynamics), intent(in) :: dynamics
    integer, intent(in), optional :: n, k, r
    type(grid_type) :: grids
    integer :: nn, kk, rr, i, nleft
    real(dp) :: lower, upper, stationary_sd, amplitude
    real(dp), allocatable :: half_grid(:), combined(:)

    nn = 50
    kk = 20
    rr = 1
    if (present(n)) nn = n
    if (present(k)) kk = k
    if (present(r)) rr = r
    nn = max(nn,2)
    kk = max(kk,1)
    rr = max(rr,0)

    select case (dynamics%model_id)
    case (model_duffie_pan_singleton, model_bates, model_heston)
      stationary_sd = sqrt(max(0.5_dp*dynamics%theta*dynamics%sigma**2/dynamics%kappa,0.0_dp))
      lower = sqrt(max(dynamics%theta-(3.0_dp+log(real(nn,dp)))*stationary_sd,1.0e-7_dp))
      upper = max(sqrt(max(dynamics%theta+(3.0_dp+log(real(nn,dp)))*stationary_sd,0.0_dp)),sqrt(0.05_dp))
      allocate(grids%var_mid_points(nn))
      do i = 1, nn
        grids%var_mid_points(i) = (lower + real(i-1,dp)*(upper-lower)/real(nn-1,dp))**2
      end do
      if (dynamics%model_id == model_heston) then
        grids%jump_counts = [0]
        grids%jump_mid_points = [0.0_dp]
      else
        allocate(grids%jump_counts(rr+1))
        grids%jump_counts = [(i-1,i=1,rr+1)]
        if (dynamics%model_id == model_duffie_pan_singleton) then
          allocate(grids%jump_mid_points(kk))
          upper = max((3.0_dp+log(real(kk,dp)))*sqrt(real(max(rr,1),dp))*dynamics%nu,1.0e-6_dp)
          do i = 1, kk
            grids%jump_mid_points(i) = 1.0e-6_dp + real(i-1,dp)*(upper-1.0e-6_dp)/real(max(kk-1,1),dp)
          end do
        else
          grids%jump_mid_points = [0.0_dp]
        end if
      end if
    case (model_taylor, model_taylor_leverage, model_pitt_malik_doucet, model_capm_sv)
      stationary_sd = dynamics%sigma/sqrt(max(1.0_dp-dynamics%phi**2,epsilon(1.0_dp)))
      nleft = nn/2
      allocate(half_grid(nleft+1))
      amplitude = sqrt(max((3.0_dp+log(real(nn,dp)))*stationary_sd,0.0_dp))
      do i = 1, nleft+1
        half_grid(i) = (real(i-1,dp)*amplitude/real(max(nleft,1),dp))**2
      end do
      allocate(combined(2*nleft+1))
      do i = 1, nleft
        combined(i) = -half_grid(nleft-i+2)
      end do
      combined(nleft+1:) = half_grid
      combined = combined + dynamics%theta
      allocate(grids%var_mid_points(nn))
      grids%var_mid_points = combined(1:nn)
      call sort_ascending(grids%var_mid_points)
      if (dynamics%model_id == model_pitt_malik_doucet) then
        grids%jump_counts = [0,1]
      else
        grids%jump_counts = [0]
      end if
      grids%jump_mid_points = [0.0_dp]
    case default
      allocate(grids%var_mid_points(nn))
      do i = 1, nn
        grids%var_mid_points(i) = -2.0_dp + 4.0_dp*real(i-1,dp)/real(nn-1,dp)
      end do
      grids%jump_counts = [0]
      grids%jump_mid_points = [0.0_dp]
    end select
  end function grid_maker

  subroutine validate_grid(grids, ok, message)
    type(grid_type), intent(in) :: grids
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    integer :: i
    ok = .false.
    message = ''
    if (.not. allocated(grids%var_mid_points) .or. size(grids%var_mid_points) < 2) then
      message = 'The volatility grid must have at least two points.'
      return
    end if
    do i = 2, size(grids%var_mid_points)
      if (grids%var_mid_points(i) <= grids%var_mid_points(i-1)) then
        message = 'The volatility grid must be strictly increasing.'
        return
      end if
    end do
    if (.not. allocated(grids%jump_counts) .or. size(grids%jump_counts) < 1) then
      message = 'The jump-count grid is empty.'
      return
    end if
    if (.not. allocated(grids%jump_mid_points) .or. size(grids%jump_mid_points) < 1) then
      message = 'The jump-size grid is empty.'
      return
    end if
    ok = .true.
  end subroutine validate_grid

  pure subroutine sort_ascending(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j
    real(dp) :: key
    do i = 2, size(x)
      key = x(i)
      j = i-1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j+1) = x(j)
        j = j-1
      end do
      x(j+1) = key
    end do
  end subroutine sort_ascending

end module svdnf_grids

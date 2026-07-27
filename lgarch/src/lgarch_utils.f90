! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2014-2025 Genaro Sucarrat (original R/C package)
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is distributed under the GNU General Public License version 2 only.
module lgarch_utils
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
  use lgarch_kinds, only : dp
  implicit none
  private
  public :: lag_vector, lag_matrix, diff_vector, diff_matrix
  public :: mean_value, sample_variance, all_finite, safe_log_mean_exp

  interface lag_values
    module procedure lag_vector
    module procedure lag_matrix
  end interface
  interface diff_values
    module procedure diff_vector
    module procedure diff_matrix
  end interface
  public :: lag_values, diff_values
contains
  function lag_vector(x, k, pad, pad_value) result(y)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: k
    logical, intent(in), optional :: pad
    real(dp), intent(in), optional :: pad_value
    real(dp), allocatable :: y(:)
    logical :: do_pad
    real(dp) :: pv
    if (k < 1) error stop "lag_vector: k must be positive"
    if (k >= size(x)) then
      if (present(pad)) then
        if (.not. pad) then
          allocate(y(0)); return
        end if
      end if
    end if
    do_pad = .true.; if (present(pad)) do_pad = pad
    pv = ieee_value(0.0_dp, ieee_quiet_nan); if (present(pad_value)) pv = pad_value
    if (do_pad) then
      allocate(y(size(x))); y = pv
      if (k < size(x)) y(k+1:) = x(:size(x)-k)
    else
      allocate(y(max(0,size(x)-k)))
      if (size(y) > 0) y = x(:size(x)-k)
    end if
  end function lag_vector

  function lag_matrix(x, k, pad, pad_value) result(y)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: k
    logical, intent(in), optional :: pad
    real(dp), intent(in), optional :: pad_value
    real(dp), allocatable :: y(:,:)
    logical :: do_pad
    real(dp) :: pv
    if (k < 1) error stop "lag_matrix: k must be positive"
    do_pad = .true.; if (present(pad)) do_pad = pad
    pv = ieee_value(0.0_dp, ieee_quiet_nan); if (present(pad_value)) pv = pad_value
    if (do_pad) then
      allocate(y(size(x,1),size(x,2))); y = pv
      if (k < size(x,1)) y(k+1:,:) = x(:size(x,1)-k,:)
    else
      allocate(y(max(0,size(x,1)-k),size(x,2)))
      if (size(y,1) > 0) y = x(:size(x,1)-k,:)
    end if
  end function lag_matrix

  function diff_vector(x, lag, pad, pad_value) result(y)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: lag
    logical, intent(in), optional :: pad
    real(dp), intent(in), optional :: pad_value
    real(dp), allocatable :: y(:)
    logical :: do_pad
    real(dp) :: pv
    if (lag < 1) error stop "diff_vector: lag must be positive"
    do_pad = .true.; if (present(pad)) do_pad = pad
    pv = ieee_value(0.0_dp, ieee_quiet_nan); if (present(pad_value)) pv = pad_value
    if (do_pad) then
      allocate(y(size(x))); y = pv
      if (lag < size(x)) y(lag+1:) = x(lag+1:) - x(:size(x)-lag)
    else
      allocate(y(max(0,size(x)-lag)))
      if (size(y) > 0) y = x(lag+1:) - x(:size(x)-lag)
    end if
  end function diff_vector

  function diff_matrix(x, lag, pad, pad_value) result(y)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: lag
    logical, intent(in), optional :: pad
    real(dp), intent(in), optional :: pad_value
    real(dp), allocatable :: y(:,:)
    logical :: do_pad
    real(dp) :: pv
    if (lag < 1) error stop "diff_matrix: lag must be positive"
    do_pad = .true.; if (present(pad)) do_pad = pad
    pv = ieee_value(0.0_dp, ieee_quiet_nan); if (present(pad_value)) pv = pad_value
    if (do_pad) then
      allocate(y(size(x,1),size(x,2))); y = pv
      if (lag < size(x,1)) y(lag+1:,:) = x(lag+1:,:) - x(:size(x,1)-lag,:)
    else
      allocate(y(max(0,size(x,1)-lag),size(x,2)))
      if (size(y,1) > 0) y = x(lag+1:,:) - x(:size(x,1)-lag,:)
    end if
  end function diff_matrix

  pure real(dp) function mean_value(x) result(ans)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      ans = 0.0_dp
    else
      ans = sum(x)/real(size(x),dp)
    end if
  end function mean_value

  pure real(dp) function sample_variance(x) result(ans)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if (size(x) < 2) then
      ans = 0.0_dp
    else
      m = sum(x)/real(size(x),dp)
      ans = sum((x-m)**2)/real(size(x)-1,dp)
    end if
  end function sample_variance

  pure logical function all_finite(x) result(ok)
    real(dp), intent(in) :: x(:)
    ok = all(ieee_is_finite(x))
  end function all_finite

  pure real(dp) function safe_log_mean_exp(x) result(ans)
    real(dp), intent(in) :: x(:)
    real(dp) :: xmax
    if (size(x) == 0) then
      ans = -huge(1.0_dp)
      return
    end if
    xmax = maxval(x)
    ans = xmax + log(sum(exp(x-xmax))/real(size(x),dp))
  end function safe_log_mean_exp
end module lgarch_utils

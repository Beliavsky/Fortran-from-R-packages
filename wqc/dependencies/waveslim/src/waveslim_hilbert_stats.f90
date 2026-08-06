! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
module waveslim_hilbert_stats
  use waveslim_kinds, only : dp
  use waveslim_types, only : complex_wavelet_transform, real_vector
  implicit none
  private

  type, public :: seasonal_coherence_result
    type(real_vector), allocatable :: coherence(:)
    type(real_vector), allocatable :: variance(:)
  end type seasonal_coherence_result

  public :: modhwt_coherence, modhwt_phase
  public :: modhwt_coherence_seasonal, modhwt_phase_seasonal

contains

  function modhwt_coherence(x, y, filter_length) result(coh)
    type(complex_wavelet_transform), intent(in) :: x, y
    integer, intent(in), optional :: filter_length
    type(real_vector), allocatable :: coh(:)
    complex(dp), allocatable :: cross(:)
    real(dp), allocatable :: co(:), quad(:), xs(:), ys(:)
    integer :: j, levels, width
    levels = min(x%levels(), y%levels())
    width = 0
    if (present(filter_length)) width = max(0, filter_length)
    allocate(coh(levels))
    do j = 1, levels
      cross = x%detail(j)%values*conjg(y%detail(j)%values)
      co = moving_average(real(cross,dp), width+1)
      quad = moving_average(-aimag(cross), width+1)
      xs = moving_average(abs(x%detail(j)%values)**2, width+1)
      ys = moving_average(abs(y%detail(j)%values)**2, width+1)
      coh(j)%values = (co*co+quad*quad)/max(xs*ys,tiny(1.0_dp))
      coh(j)%values = max(0.0_dp,min(1.0_dp,coh(j)%values))
    end do
  end function modhwt_coherence

  function modhwt_phase(x, y, filter_length) result(phase)
    type(complex_wavelet_transform), intent(in) :: x, y
    integer, intent(in), optional :: filter_length
    type(real_vector), allocatable :: phase(:)
    complex(dp), allocatable :: cross(:)
    real(dp), allocatable :: co(:), quad(:)
    integer :: j, levels, width
    levels = min(x%levels(), y%levels())
    width = 0
    if (present(filter_length)) width = max(0, filter_length)
    allocate(phase(levels))
    do j = 1, levels
      cross = x%detail(j)%values*conjg(y%detail(j)%values)
      co = moving_average(real(cross,dp), width+1)
      quad = moving_average(-aimag(cross), width+1)
      phase(j)%values = atan2(-quad,co)
    end do
  end function modhwt_phase

  function modhwt_coherence_seasonal(x, y, season) result(ans)
    type(complex_wavelet_transform), intent(in) :: x, y
    integer, intent(in) :: season
    type(seasonal_coherence_result) :: ans
    complex(dp), allocatable :: cross(:)
    real(dp), allocatable :: co(:), quad(:), xs(:), ys(:)
    integer :: j, levels
    levels = min(x%levels(),y%levels())
    allocate(ans%coherence(levels),ans%variance(levels))
    do j = 1, levels
      cross = x%detail(j)%values*conjg(y%detail(j)%values)
      co = seasonal_mean(real(cross,dp),season)
      quad = seasonal_mean(-aimag(cross),season)
      xs = seasonal_mean(abs(x%detail(j)%values)**2,season)
      ys = seasonal_mean(abs(y%detail(j)%values)**2,season)
      ans%coherence(j)%values = (co*co+quad*quad)/max(xs*ys,tiny(1.0_dp))
      ans%coherence(j)%values = max(0.0_dp,min(1.0_dp,ans%coherence(j)%values))
      ! A conservative delta-method proxy; the exact R code uses
      ! seasonal HAC terms.
      ans%variance(j)%values = 4.0_dp*ans%coherence(j)%values* &
        max(0.0_dp,1.0_dp-ans%coherence(j)%values)/ &
        real(max(1,size(cross)/max(1,season)),dp)
    end do
  end function modhwt_coherence_seasonal

  function modhwt_phase_seasonal(x, y, season) result(phase)
    type(complex_wavelet_transform), intent(in) :: x, y
    integer, intent(in) :: season
    type(real_vector), allocatable :: phase(:)
    complex(dp), allocatable :: cross(:)
    real(dp), allocatable :: co(:), quad(:)
    integer :: j, levels
    levels = min(x%levels(),y%levels())
    allocate(phase(levels))
    do j = 1, levels
      cross = x%detail(j)%values*conjg(y%detail(j)%values)
      co = seasonal_mean(real(cross,dp),season)
      quad = seasonal_mean(-aimag(cross),season)
      phase(j)%values = atan2(-quad,co)
    end do
  end function modhwt_phase_seasonal

  function moving_average(x, width) result(y)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: width
    real(dp), allocatable :: y(:)
    integer :: i, lo, hi, half
    allocate(y(size(x)))
    half = max(0,width-1)/2
    do i = 1, size(x)
      lo = max(1,i-half)
      hi = min(size(x),lo+max(1,width)-1)
      lo = max(1,hi-max(1,width)+1)
      y(i) = sum(x(lo:hi))/real(hi-lo+1,dp)
    end do
  end function moving_average

  function seasonal_mean(x, season) result(y)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: season
    real(dp), allocatable :: y(:)
    integer :: i, k, count, s
    s = max(1,season)
    allocate(y(s))
    y = 0.0_dp
    do k = 1, s
      count = 0
      do i = k, size(x), s
        y(k) = y(k)+x(i)
        count = count+1
      end do
      if (count > 0) y(k) = y(k)/real(count,dp)
    end do
  end function seasonal_mean

end module waveslim_hilbert_stats

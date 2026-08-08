! SPDX-License-Identifier: MIT
module gradient_stats
  use gradient_kinds, only : dp
  implicit none
  private
  public :: vector_median, vector_sd
contains
  function vector_median(x) result(m)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    real(dp), allocatable :: y(:)
    real(dp) :: tmp
    integer :: i, j, n
    n = size(x)
    if (n == 0) then
      m = 0.0_dp
      return
    end if
    allocate(y(n))
    y = x
    do i = 2, n
      tmp = y(i)
      j = i - 1
      do while (j >= 1)
        if (y(j) <= tmp) exit
        y(j+1) = y(j)
        j = j - 1
      end do
      y(j+1) = tmp
    end do
    if (mod(n,2) == 1) then
      m = y((n+1)/2)
    else
      m = 0.5_dp*(y(n/2)+y(n/2+1))
    end if
  end function vector_median

  function vector_sd(x) result(s)
    real(dp), intent(in) :: x(:)
    real(dp) :: s, mu
    integer :: n
    n = size(x)
    if (n <= 1) then
      s = 0.0_dp
      return
    end if
    mu = sum(x)/real(n,dp)
    s = sqrt(sum((x-mu)**2)/real(n-1,dp))
  end function vector_sd
end module gradient_stats

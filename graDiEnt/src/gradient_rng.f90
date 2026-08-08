! SPDX-License-Identifier: MIT
module gradient_rng
  use gradient_kinds, only : dp
  implicit none
  private
  public :: seed_rng, random_normal, random_normal_vector
  public :: random_uniform_range, sample_without_replacement
contains
  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)

    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(seed + 104729*i + 37*i*i, huge(1)-1)
      if (put(i) <= 0) put(i) = i + 17
    end do
    call random_seed(put=put)
  end subroutine seed_rng

  function random_normal() result(z)
    real(dp) :: z
    real(dp) :: u1, u2
    real(dp), parameter :: twopi = 6.2831853071795864769252867665590058_dp

    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    z = sqrt(-2.0_dp*log(u1))*cos(twopi*u2)
  end function random_normal

  subroutine random_normal_vector(x)
    real(dp), intent(out) :: x(:)
    integer :: i
    do i = 1, size(x)
      x(i) = random_normal()
    end do
  end subroutine random_normal_vector

  function random_uniform_range(lo, hi) result(x)
    real(dp), intent(in) :: lo, hi
    real(dp) :: x, u
    call random_number(u)
    x = lo + (hi-lo)*u
  end function random_uniform_range

  subroutine sample_without_replacement(pool, k, sample, status)
    integer, intent(in) :: pool(:)
    integer, intent(in) :: k
    integer, intent(out) :: sample(:)
    integer, intent(out) :: status
    integer, allocatable :: work(:)
    integer :: i, j, tmp, n
    real(dp) :: u

    n = size(pool)
    status = 0
    if (k < 0 .or. k > n .or. size(sample) /= k) then
      status = 1
      return
    end if
    allocate(work(n))
    work = pool
    do i = 1, k
      call random_number(u)
      j = i + int(u*real(n-i+1,dp))
      if (j > n) j = n
      tmp = work(i)
      work(i) = work(j)
      work(j) = tmp
      sample(i) = work(i)
    end do
  end subroutine sample_without_replacement
end module gradient_rng

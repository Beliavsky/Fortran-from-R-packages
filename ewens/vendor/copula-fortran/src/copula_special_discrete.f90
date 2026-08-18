! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan
module copula_special_discrete
  use copula_kinds, only : dp, i8
  use copula_random, only : seed_random, random_uniform
  implicit none
  private
  public :: stirling_first, stirling_second, eulerian_number
  public :: sibuya_pmf, random_sibuya, logseries_pmf, random_logseries
contains
  integer(i8) function stirling_first(n, k) result(value)
    integer, intent(in) :: n, k
    integer(i8), allocatable :: table(:,:)
    integer :: i, j
    if (n < 0 .or. k < 0 .or. k > n) then
      value = 0_i8
      return
    end if
    allocate(table(0:n,0:n))
    table = 0_i8
    table(0,0) = 1_i8
    do i = 1, n
      do j = 1, i
        table(i,j) = table(i-1,j-1)+int(i-1,i8)*table(i-1,j)
      end do
    end do
    value = table(n,k)
  end function stirling_first

  integer(i8) function stirling_second(n, k) result(value)
    integer, intent(in) :: n, k
    integer(i8), allocatable :: table(:,:)
    integer :: i, j
    if (n < 0 .or. k < 0 .or. k > n) then
      value = 0_i8
      return
    end if
    allocate(table(0:n,0:n))
    table = 0_i8
    table(0,0) = 1_i8
    do i = 1, n
      do j = 1, i
        table(i,j) = table(i-1,j-1)+int(j,i8)*table(i-1,j)
      end do
    end do
    value = table(n,k)
  end function stirling_second

  integer(i8) function eulerian_number(n, k) result(value)
    integer, intent(in) :: n, k
    integer(i8), allocatable :: table(:,:)
    integer :: i, j
    if (n < 0 .or. k < 0 .or. k >= max(1,n)) then
      if (n == 0 .and. k == 0) then
        value = 1_i8
      else
        value = 0_i8
      end if
      return
    end if
    allocate(table(0:n,0:n))
    table = 0_i8
    table(0,0) = 1_i8
    do i = 1, n
      do j = 0, i-1
        if (j > 0) table(i,j) = table(i,j)+int(i-j,i8)*table(i-1,j-1)
        table(i,j) = table(i,j)+int(j+1,i8)*table(i-1,j)
      end do
    end do
    value = table(n,k)
  end function eulerian_number

  real(dp) function sibuya_pmf(k, alpha) result(value)
    integer, intent(in) :: k
    real(dp), intent(in) :: alpha
    integer :: j
    if (k < 1 .or. alpha <= 0.0_dp .or. alpha > 1.0_dp) then
      value = 0.0_dp
      return
    end if
    value = alpha
    do j = 2, k
      value = value*(real(j-1,dp)-alpha)/real(j,dp)
    end do
  end function sibuya_pmf

  integer function random_sibuya(alpha, seed) result(k)
    real(dp), intent(in) :: alpha
    integer(i8), intent(in), optional :: seed
    real(dp) :: target, cumulative, mass
    if (present(seed)) call seed_random(seed)
    target = random_uniform()
    k = 1
    mass = alpha
    cumulative = mass
    do while (target > cumulative .and. k < 1000000)
      k = k+1
      mass = mass*(real(k-1,dp)-alpha)/real(k,dp)
      cumulative = cumulative+mass
    end do
  end function random_sibuya

  real(dp) function logseries_pmf(k, probability) result(value)
    integer, intent(in) :: k
    real(dp), intent(in) :: probability
    if (k < 1 .or. probability <= 0.0_dp .or. probability >= 1.0_dp) then
      value = 0.0_dp
    else
      value = -probability**k/(real(k,dp)*log(1.0_dp-probability))
    end if
  end function logseries_pmf

  integer function random_logseries(probability, seed) result(k)
    real(dp), intent(in) :: probability
    integer(i8), intent(in), optional :: seed
    real(dp) :: target, cumulative, mass
    if (present(seed)) call seed_random(seed)
    target = random_uniform()
    k = 1
    mass = logseries_pmf(1,probability)
    cumulative = mass
    do while (target > cumulative .and. k < 1000000)
      k = k+1
      mass = mass*probability*real(k-1,dp)/real(k,dp)
      cumulative = cumulative+mass
    end do
  end function random_logseries
end module copula_special_discrete

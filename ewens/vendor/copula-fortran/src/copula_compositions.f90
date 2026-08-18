! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan
module copula_compositions
  use copula_kinds, only : dp, i8
  use copula_types, only : copula_model
  use copula_families, only : copula_cdf, copula_density
  use copula_simulation, only : random_copula
  use copula_random, only : seed_random, random_uniform
  implicit none
  private
  public :: mixture_cdf, mixture_density, random_mixture
  public :: khoudraji_cdf, khoudraji_density, nested_clayton_cdf
contains
  real(dp) function mixture_cdf(u, models, weights) result(value)
    real(dp), intent(in) :: u(:)
    type(copula_model), intent(in) :: models(:)
    real(dp), intent(in) :: weights(:)
    integer :: i
    value = 0.0_dp
    if (size(models) /= size(weights) .or. sum(weights) <= 0.0_dp) return
    do i = 1, size(models)
      value = value+weights(i)*copula_cdf(u,models(i))
    end do
    value = value/sum(weights)
  end function mixture_cdf

  real(dp) function mixture_density(u, models, weights) result(value)
    real(dp), intent(in) :: u(:)
    type(copula_model), intent(in) :: models(:)
    real(dp), intent(in) :: weights(:)
    integer :: i
    value = 0.0_dp
    if (size(models) /= size(weights) .or. sum(weights) <= 0.0_dp) return
    do i = 1, size(models)
      value = value+weights(i)*copula_density(u,models(i))
    end do
    value = value/sum(weights)
  end function mixture_density

  subroutine random_mixture(n, models, weights, samples, ok, seed)
    integer, intent(in) :: n
    type(copula_model), intent(in) :: models(:)
    real(dp), intent(in) :: weights(:)
    real(dp), allocatable, intent(out) :: samples(:,:)
    logical, intent(out) :: ok
    integer(i8), intent(in), optional :: seed
    real(dp), allocatable :: one(:,:), cumulative(:)
    real(dp) :: total, draw
    integer :: i, j, selected, d
    ok = n >= 1 .and. size(models) == size(weights) .and. size(models) >= 1 .and. all(weights >= 0.0_dp)
    if (.not. ok) then
      allocate(samples(0,0))
      return
    end if
    d = models(1)%dimension
    if (any([(models(i)%dimension /= d,i=1,size(models))])) then
      ok = .false.
      allocate(samples(0,0))
      return
    end if
    if (present(seed)) call seed_random(seed)
    total = sum(weights)
    if (total <= 0.0_dp) then
      ok = .false.
      allocate(samples(0,0))
      return
    end if
    allocate(samples(n,d),cumulative(size(weights)))
    cumulative(1) = weights(1)/total
    do i = 2, size(weights)
      cumulative(i) = cumulative(i-1)+weights(i)/total
    end do
    do i = 1, n
      draw = random_uniform()
      selected = size(models)
      do j = 1, size(models)
        if (draw <= cumulative(j)) then
          selected = j
          exit
        end if
      end do
      call random_copula(1,models(selected),one,ok)
      if (.not. ok) return
      samples(i,:) = one(1,:)
    end do
  end subroutine random_mixture

  real(dp) function khoudraji_cdf(u, base, shape) result(value)
    real(dp), intent(in) :: u(2), shape(2)
    type(copula_model), intent(in) :: base
    real(dp) :: transformed(2)
    transformed = u**shape
    value = u(1)**(1.0_dp-shape(1))*u(2)**(1.0_dp-shape(2))*copula_cdf(transformed,base)
  end function khoudraji_cdf

  real(dp) function khoudraji_density(u, base, shape) result(value)
    real(dp), intent(in) :: u(2), shape(2)
    type(copula_model), intent(in) :: base
    real(dp) :: h1, h2
    h1 = min(2.0e-5_dp,0.25_dp*min(u(1),1.0_dp-u(1)))
    h2 = min(2.0e-5_dp,0.25_dp*min(u(2),1.0_dp-u(2)))
    value = (khoudraji_cdf([u(1)+h1,u(2)+h2],base,shape) - &
      khoudraji_cdf([u(1)+h1,u(2)-h2],base,shape) - &
      khoudraji_cdf([u(1)-h1,u(2)+h2],base,shape) + &
      khoudraji_cdf([u(1)-h1,u(2)-h2],base,shape))/(4.0_dp*h1*h2)
  end function khoudraji_density

  real(dp) function nested_clayton_cdf(u, child_dimension, root_theta, child_theta) result(value)
    real(dp), intent(in) :: u(:)
    integer, intent(in) :: child_dimension
    real(dp), intent(in) :: root_theta, child_theta
    real(dp) :: child_value, total
    integer :: d
    d = size(u)
    if (child_dimension < 2 .or. child_dimension >= d .or. root_theta <= 0.0_dp .or. &
        child_theta < root_theta .or. any(u <= 0.0_dp) .or. any(u > 1.0_dp)) then
      value = 0.0_dp
      return
    end if
    child_value = (sum(u(1:child_dimension)**(-child_theta))-real(child_dimension,dp)+1.0_dp)** &
      (-1.0_dp/child_theta)
    total = child_value**(-root_theta)+sum(u(child_dimension+1:d)**(-root_theta)) - real(d-child_dimension,dp)
    value = total**(-1.0_dp/root_theta)
  end function nested_clayton_cdf
end module copula_compositions

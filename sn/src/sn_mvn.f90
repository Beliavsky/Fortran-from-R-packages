! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
module sn_mvn
  use sn_kinds, only : dp, pi, tiny_dp
  use sn_math, only : normal_cdf, normal_quantile, adaptive_simpson, &
                      chi_square_quantile, clamp_probability
  use sn_linalg, only : cholesky_lower, covariance_to_correlation, rmvn
  use sn_rng, only : sn_rng_state
  use sn_status, only : sn_ok, sn_dimension_mismatch, sn_not_positive_definite, &
                        sn_no_convergence
  implicit none
  private

  public :: bivariate_normal_cdf, mvn_cdf, mvt_cdf
  public :: rtmvn_lower

contains

  real(dp) function bivariate_normal_cdf(h, k, rho) result(value)
    real(dp), intent(in) :: h, k, rho
    real(dp) :: r, base, integ

    if (h < -38.0_dp .or. k < -38.0_dp) then
      value = 0.0_dp
      return
    end if
    if (h > 38.0_dp) then
      value = normal_cdf(k)
      return
    end if
    if (k > 38.0_dp) then
      value = normal_cdf(h)
      return
    end if
    r = max(-1.0_dp+1.0e-12_dp,min(1.0_dp-1.0e-12_dp,rho))
    if (abs(r) <= 10.0_dp*epsilon(1.0_dp)) then
      value = normal_cdf(h)*normal_cdf(k)
      return
    end if
    base = normal_cdf(h)*normal_cdf(k)
    integ = adaptive_simpson(integrand,0.0_dp,r,1.0e-11_dp,24)
    value = clamp_probability(base+integ/(2.0_dp*pi))
  contains
    real(dp) function integrand(t) result(y)
      real(dp), intent(in) :: t
      real(dp) :: den
      den = max(1.0_dp-t*t,tiny_dp)
      y = exp(-(h*h-2.0_dp*h*k*t+k*k)/(2.0_dp*den))/sqrt(den)
    end function integrand
  end function bivariate_normal_cdf

  real(dp) function mvn_cdf(upper, mean, cov, info, samples) result(value)
    real(dp), intent(in) :: upper(:), mean(:), cov(:,:)
    integer, intent(out), optional :: info
    integer, intent(in), optional :: samples
    real(dp), allocatable :: l(:,:), sd(:), cor(:,:), z(:), shift(:)
    real(dp) :: prod, bound, u, total
    integer :: d, ierr, ns, i, j, s

    d = size(upper)
    if (size(mean) /= d .or. size(cov,1) /= d .or. size(cov,2) /= d) then
      value = 0.0_dp
      if (present(info)) info = sn_dimension_mismatch
      return
    end if
    if (d == 0) then
      value = 1.0_dp
      if (present(info)) info = sn_ok
      return
    else if (d == 1) then
      if (cov(1,1) <= 0.0_dp) then
        value = 0.0_dp
        if (present(info)) info = sn_not_positive_definite
      else
        value = normal_cdf((upper(1)-mean(1))/sqrt(cov(1,1)))
        if (present(info)) info = sn_ok
      end if
      return
    else if (d == 2) then
      call covariance_to_correlation(cov,cor,sd,ierr)
      if (ierr /= sn_ok) then
        value = 0.0_dp
        if (present(info)) info = ierr
        return
      end if
      value = bivariate_normal_cdf((upper(1)-mean(1))/sd(1), &
                                   (upper(2)-mean(2))/sd(2),cor(1,2))
      if (present(info)) info = sn_ok
      return
    end if

    call cholesky_lower(cov,l,ierr)
    if (ierr /= sn_ok) then
      value = 0.0_dp
      if (present(info)) info = ierr
      return
    end if
    ns = 32768
    if (present(samples)) ns = max(256,samples)
    allocate(z(d),shift(d-1))
    do j=1,d-1
      shift(j) = modulo(0.6180339887498948482_dp*real(j*j+3*j+1,dp),1.0_dp)
    end do
    total = 0.0_dp
    do s=1,ns
      z = 0.0_dp
      prod = 1.0_dp
      do i=1,d
        bound = (upper(i)-mean(i)-dot_product(l(i,1:i-1),z(1:i-1)))/l(i,i)
        bound = normal_cdf(bound)
        prod = prod*bound
        if (prod <= tiny_dp) exit
        if (i < d) then
          u = modulo(halton(s,prime(i))+shift(i),1.0_dp)
          u = max(tiny_dp,min(1.0_dp-epsilon(1.0_dp),u*bound))
          z(i) = normal_quantile(u)
        end if
      end do
      total = total+prod
      ! Antithetic point.
      z = 0.0_dp
      prod = 1.0_dp
      do i=1,d
        bound = (upper(i)-mean(i)-dot_product(l(i,1:i-1),z(1:i-1)))/l(i,i)
        bound = normal_cdf(bound)
        prod = prod*bound
        if (prod <= tiny_dp) exit
        if (i < d) then
          u = 1.0_dp-modulo(halton(s,prime(i))+shift(i),1.0_dp)
          u = max(tiny_dp,min(1.0_dp-epsilon(1.0_dp),u*bound))
          z(i) = normal_quantile(u)
        end if
      end do
      total = total+prod
    end do
    value = clamp_probability(total/(2.0_dp*real(ns,dp)))
    if (present(info)) info = sn_ok
  end function mvn_cdf

  real(dp) function mvt_cdf(upper, mean, scale, nu, info, samples) result(value)
    real(dp), intent(in) :: upper(:), mean(:), scale(:,:), nu
    integer, intent(out), optional :: info
    integer, intent(in), optional :: samples
    real(dp), allocatable :: l(:,:), z(:), x(:), shift(:)
    real(dp) :: v, u
    integer :: d, ns, ierr, j, s, count

    d = size(upper)
    if (size(mean) /= d .or. size(scale,1) /= d .or. size(scale,2) /= d .or. nu <= 0.0_dp) then
      value = 0.0_dp
      if (present(info)) info = sn_dimension_mismatch
      return
    end if
    if (d == 1) then
      use_univariate: block
        use sn_math, only : student_t_cdf
        value = student_t_cdf((upper(1)-mean(1))/sqrt(scale(1,1)),nu)
      end block use_univariate
      if (present(info)) info = sn_ok
      return
    end if
    call cholesky_lower(scale,l,ierr)
    if (ierr /= sn_ok) then
      value = 0.0_dp
      if (present(info)) info = ierr
      return
    end if
    ns = 65536
    if (present(samples)) ns = max(1024,samples)
    allocate(z(d),x(d),shift(d+1))
    do j=1,d+1
      shift(j) = modulo(0.7548776662466927_dp*real(j*j+7*j+3,dp),1.0_dp)
    end do
    count = 0
    do s=1,ns
      do j=1,d
        u = modulo(halton(s,prime(j))+shift(j),1.0_dp)
        z(j) = normal_quantile(max(tiny_dp,min(1.0_dp-epsilon(1.0_dp),u)))
      end do
      u = modulo(halton(s,prime(d+1))+shift(d+1),1.0_dp)
      v = chi_square_quantile(max(tiny_dp,min(1.0_dp-epsilon(1.0_dp),u)),nu)/nu
      x = mean+matmul(l,z)/sqrt(v)
      if (all(x <= upper)) count = count+1
    end do
    value = real(count,dp)/real(ns,dp)
    if (present(info)) info = sn_ok
  end function mvt_cdf

  subroutine rtmvn_lower(rng, n, mean, cov, lower, x, info, max_attempts)
    type(sn_rng_state), intent(inout) :: rng
    integer, intent(in) :: n
    real(dp), intent(in) :: mean(:), cov(:,:), lower(:)
    real(dp), allocatable, intent(out) :: x(:,:)
    integer, intent(out) :: info
    integer, intent(in), optional :: max_attempts
    real(dp), allocatable :: candidate(:,:)
    integer :: d, accepted, attempts, limit, ierr

    d = size(mean)
    if (size(lower) /= d .or. size(cov,1) /= d .or. size(cov,2) /= d .or. n < 0) then
      allocate(x(0,0))
      info = sn_dimension_mismatch
      return
    end if
    limit = max(10000,1000*n)
    if (present(max_attempts)) limit = max_attempts
    allocate(x(n,d))
    accepted = 0
    attempts = 0
    do while (accepted < n .and. attempts < limit)
      call rmvn(rng,1,mean,cov,candidate,ierr)
      if (ierr /= sn_ok) then
        info = ierr
        return
      end if
      attempts = attempts+1
      if (all(candidate(1,:) >= lower)) then
        accepted = accepted+1
        x(accepted,:) = candidate(1,:)
      end if
    end do
    if (accepted < n) then
      if (accepted > 0) x(accepted+1:n,:) = 0.0_dp
      info = sn_no_convergence
    else
      info = sn_ok
    end if
  end subroutine rtmvn_lower

  pure real(dp) function halton(index, base) result(value)
    integer, intent(in) :: index, base
    integer :: i
    real(dp) :: f
    i = index
    f = 1.0_dp
    value = 0.0_dp
    do while (i > 0)
      f = f/real(base,dp)
      value = value+f*real(mod(i,base),dp)
      i = i/base
    end do
  end function halton

  pure integer function prime(k) result(p)
    integer, intent(in) :: k
    integer, parameter :: primes(64) = [ &
       2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89, &
       97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181, &
       191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281, &
       283,293,307,311]
    if (k <= size(primes)) then
      p = primes(max(1,k))
    else
      p = primes(size(primes))
    end if
  end function prime

end module sn_mvn

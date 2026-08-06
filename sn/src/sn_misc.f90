! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
module sn_misc
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use sn_kinds, only : dp, pi, tiny_dp
  use sn_status, only : sn_ok, sn_invalid_argument, sn_no_convergence
  use sn_math, only : normal_pdf, normal_cdf, student_t_pdf, student_t_cdf, &
                      adaptive_simpson, clamp_probability
  use sn_rng, only : sn_rng_state
  implicit none
  private

  abstract interface
    pure real(dp) function scalar_density(x)
      import dp
      real(dp), intent(in) :: x
    end function scalar_density

    pure real(dp) function perturbation_function(x)
      import dp
      real(dp), intent(in) :: x
    end function perturbation_function

    subroutine scalar_random(rng,x)
      import dp, sn_rng_state
      type(sn_rng_state), intent(inout) :: rng
      real(dp), intent(out) :: x
    end subroutine scalar_random
  end interface

  public :: d_symm_modulated, r_symm_modulated
  public :: pprodn2, pprodt2, qprodt2, q_penalty
  public :: galton_moors_to_alpha_nu

contains

  real(dp) function d_symm_modulated(x, base_density, perturbation, log_pdf) result(value)
    real(dp), intent(in) :: x
    procedure(scalar_density) :: base_density
    procedure(perturbation_function) :: perturbation
    logical, intent(in), optional :: log_pdf
    real(dp) :: d, g
    logical :: give_log
    give_log = .false.
    if (present(log_pdf)) give_log = log_pdf
    d = base_density(x)
    g = min(1.0_dp,max(0.0_dp,perturbation(x)))
    if (give_log) then
      if (d <= 0.0_dp .or. g <= 0.0_dp) then
        value = -huge(1.0_dp)
      else
        value = log(2.0_dp)+log(d)+log(g)
      end if
    else
      value = 2.0_dp*d*g
    end if
  end function d_symm_modulated

  subroutine r_symm_modulated(rng, n, base_random, perturbation, x, info, max_trials)
    type(sn_rng_state), intent(inout) :: rng
    integer, intent(in) :: n
    procedure(scalar_random) :: base_random
    procedure(perturbation_function) :: perturbation
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out), optional :: info
    integer, intent(in), optional :: max_trials
    integer :: i, trials, limit
    real(dp) :: candidate, g
    if (n < 0) then
      allocate(x(0))
      if (present(info)) info = sn_invalid_argument
      return
    end if
    allocate(x(n))
    limit = 100000
    if (present(max_trials)) limit = max_trials
    do i=1,n
      trials = 0
      do
        call base_random(rng,candidate)
        g = min(1.0_dp,max(0.0_dp,perturbation(candidate)))
        trials = trials+1
        if (rng%uniform() <= g) exit
        if (trials >= limit) then
          x(i:) = ieee_value(0.0_dp,ieee_quiet_nan)
          if (present(info)) info = sn_no_convergence
          return
        end if
      end do
      x(i) = candidate
    end do
    if (present(info)) info = sn_ok
  end subroutine r_symm_modulated

  real(dp) function pprodn2(q, rho, tol) result(value)
    real(dp), intent(in) :: q, rho
    real(dp), intent(in), optional :: tol
    real(dp) :: eps, rr, tail
    eps = 1.0e-9_dp
    if (present(tol)) eps = tol
    if (abs(rho) >= 1.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    rr = sqrt(max(tiny_dp,1.0_dp-rho*rho))
    tail = max(10.0_dp,sqrt(max(0.0_dp,2.0_dp*log(1.0_dp/max(eps,tiny_dp)))))
    value = adaptive_simpson(integrand,-tail,-sqrt(tiny_dp),eps,24)+ &
            adaptive_simpson(integrand,sqrt(tiny_dp),tail,eps,24)
    value = clamp_probability(value)
  contains
    real(dp) function integrand(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: z
      if (abs(x) <= tiny_dp) then
        y = 0.0_dp
      else
        z = (q/x-rho*x)/rr
        if (x > 0.0_dp) then
          y = normal_pdf(x)*normal_cdf(z)
        else
          y = normal_pdf(x)*(1.0_dp-normal_cdf(z))
        end if
      end if
    end function integrand
  end function pprodn2

  real(dp) function pprodt2(q, rho, nu, tol) result(value)
    real(dp), intent(in) :: q, rho, nu
    real(dp), intent(in), optional :: tol
    real(dp) :: eps, rr, tail
    eps = 2.0e-8_dp
    if (present(tol)) eps = tol
    if (abs(rho) >= 1.0_dp .or. nu <= 0.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    rr = sqrt(max(tiny_dp,1.0_dp-rho*rho))
    tail = max(50.0_dp,20.0_dp*sqrt(max(1.0_dp,nu/max(nu-2.0_dp,0.25_dp))))
    value = adaptive_simpson(integrand,-tail,-sqrt(tiny_dp),eps,28)+ &
            adaptive_simpson(integrand,sqrt(tiny_dp),tail,eps,28)
    value = clamp_probability(value)
  contains
    real(dp) function integrand(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: scale, z
      if (abs(x) <= tiny_dp) then
        y = 0.0_dp
      else
        scale = sqrt((nu+x*x)/(nu+1.0_dp))*rr
        z = (q/x-rho*x)/max(scale,tiny_dp)
        if (x > 0.0_dp) then
          y = student_t_pdf(x,nu)*student_t_cdf(z,nu+1.0_dp)
        else
          y = student_t_pdf(x,nu)*(1.0_dp-student_t_cdf(z,nu+1.0_dp))
        end if
      end if
    end function integrand
  end function pprodt2

  real(dp) function qprodt2(p, rho, nu, tol, info) result(value)
    real(dp), intent(in) :: p, rho, nu
    real(dp), intent(in), optional :: tol
    integer, intent(out), optional :: info
    real(dp) :: eps, lo, hi, mid, pmid
    integer :: iter
    eps = 1.0e-7_dp
    if (present(tol)) eps = tol
    if (p < 0.0_dp .or. p > 1.0_dp .or. abs(rho)>=1.0_dp .or. nu<=0.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      if (present(info)) info = sn_invalid_argument
      return
    end if
    lo = -1.0_dp
    hi = 1.0_dp
    do while (pprodt2(lo,rho,nu,0.2_dp*eps) > p)
      lo = 2.0_dp*lo
      if (lo < -1.0e12_dp) exit
    end do
    do while (pprodt2(hi,rho,nu,0.2_dp*eps) < p)
      hi = 2.0_dp*hi
      if (hi > 1.0e12_dp) exit
    end do
    do iter=1,120
      mid = 0.5_dp*(lo+hi)
      pmid = pprodt2(mid,rho,nu,0.2_dp*eps)
      if (abs(pmid-p) <= eps) then
        lo = mid
        hi = mid
        exit
      else if (pmid < p) then
        lo = mid
      else
        hi = mid
      end if
    end do
    value = 0.5_dp*(lo+hi)
    if (present(info)) info = merge(sn_ok,sn_no_convergence,iter<=120)
  end function qprodt2

  pure real(dp) function q_penalty(alpha, nu) result(value)
    real(dp), intent(in) :: alpha
    real(dp), intent(in), optional :: nu
    real(dp) :: e1, e2, c1, c2, df, euler_gamma
    e1 = 1.0_dp/3.0_dp
    e2 = 0.2854166_dp
    euler_gamma = 0.57721_dp
    if (present(nu)) then
      df = nu
      if (df < huge(1.0_dp)/10.0_dp) then
        e1 = e1*(df+2.0_dp)*(df+3.0_dp)/(df+1.0_dp)**2
        e2 = e2*(1.0_dp+4.0_dp/(df+euler_gamma))
      end if
    end if
    c1 = 1.0_dp/(4.0_dp*e2)
    c2 = e2/e1
    value = c1*log(1.0_dp+c2*alpha*alpha)
  end function q_penalty

  subroutine galton_moors_to_alpha_nu(galton, moors, alpha, nu, info)
    real(dp), intent(in) :: galton, moors
    real(dp), intent(out) :: alpha, nu
    integer, intent(out), optional :: info
    real(dp) :: a, df, best, score, sg, sm
    integer :: ia, idf
    ! A deterministic grid inversion of the quantile summaries. This mirrors
    ! the role of the package's spline-based preliminary estimator without
    ! depending on serialized R interpolation objects.
    best = huge(1.0_dp)
    alpha = 0.0_dp
    nu = 20.0_dp
    do idf=0,78
      df = 1.0_dp+0.5_dp*real(idf,dp)
      do ia=-80,80
        a = 0.1_dp*real(ia,dp)
        call approximate_summaries(a,df,sg,sm)
        score = (sg-galton)**2+0.25_dp*(sm-moors)**2
        if (score < best) then
          best = score
          alpha = a
          nu = df
        end if
      end do
    end do
    if (present(info)) info = sn_ok
  contains
    subroutine approximate_summaries(a0,df0,g,m)
      use sn_univariate, only : qst
      real(dp), intent(in) :: a0,df0
      real(dp), intent(out) :: g,m
      real(dp) :: q1,q2,q3,e1,e2,e3,e4,e5,e6,e7
      q1=qst(0.25_dp,alpha=a0,nu=df0)
      q2=qst(0.50_dp,alpha=a0,nu=df0)
      q3=qst(0.75_dp,alpha=a0,nu=df0)
      e1=qst(0.125_dp,alpha=a0,nu=df0)
      e2=q1
      e3=qst(0.375_dp,alpha=a0,nu=df0)
      e4=q2
      e5=qst(0.625_dp,alpha=a0,nu=df0)
      e6=q3
      e7=qst(0.875_dp,alpha=a0,nu=df0)
      g=(q3+q1-2.0_dp*q2)/max(q3-q1,tiny_dp)
      m=((e7-e5)+(e3-e1))/max(e6-e2,tiny_dp)
    end subroutine approximate_summaries
  end subroutine galton_moors_to_alpha_nu

end module sn_misc

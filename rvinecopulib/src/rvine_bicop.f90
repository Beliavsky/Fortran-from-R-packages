! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
module rvine_bicop
  use rvine_kinds, only : dp, pi, eps_prob, clamp_prob
  use rvine_math, only : normal_cdf, normal_pdf, normal_quantile, student_cdf, &
                         student_pdf, student_quantile, gauss_legendre_rule, random_uniform
  implicit none
  private

  integer, parameter, public :: bicop_indep = 0
  integer, parameter, public :: bicop_gaussian = 1
  integer, parameter, public :: bicop_student = 2
  integer, parameter, public :: bicop_clayton = 3
  integer, parameter, public :: bicop_gumbel = 4
  integer, parameter, public :: bicop_frank = 5
  integer, parameter, public :: bicop_joe = 6
  integer, parameter, public :: bicop_bb1 = 7
  integer, parameter, public :: bicop_bb6 = 8
  integer, parameter, public :: bicop_bb7 = 9
  integer, parameter, public :: bicop_bb8 = 10
  integer, parameter, public :: bicop_tawn = 11

  type, public :: bicop_model
    integer :: family = bicop_indep
    integer :: rotation = 0
    integer :: npar = 0
    real(dp) :: parameters(3) = 0.0_dp
    real(dp) :: loglik = 0.0_dp
    integer :: nobs = 0
  contains
    procedure :: pdf => bicop_pdf
    procedure :: cdf => bicop_cdf
    procedure :: hfunc1 => bicop_hfunc1
    procedure :: hfunc2 => bicop_hfunc2
    procedure :: hinv1 => bicop_hinv1
    procedure :: hinv2 => bicop_hinv2
    procedure :: tau => bicop_tau
    procedure :: simulate => bicop_simulate
    procedure :: aic => bicop_aic
    procedure :: bic => bicop_bic
    procedure :: name => bicop_name
  end type bicop_model

  public :: make_bicop, family_parameter_count, family_bounds

contains

  function make_bicop(family, rotation, parameters) result(model)
    integer, intent(in) :: family
    integer, intent(in), optional :: rotation
    real(dp), intent(in), optional :: parameters(:)
    type(bicop_model) :: model
    model%family = family
    model%npar = family_parameter_count(family)
    if (present(rotation)) model%rotation = rotation
    if (present(parameters)) then
      model%parameters(1:min(3,size(parameters))) = parameters(1:min(3,size(parameters)))
    else
      call default_parameters(model)
    end if
  end function make_bicop

  subroutine default_parameters(model)
    type(bicop_model), intent(inout) :: model
    select case (model%family)
    case (bicop_gaussian)
      model%parameters(1) = 0.0_dp
    case (bicop_student)
      model%parameters(1:2) = [0.0_dp, 5.0_dp]
    case (bicop_clayton)
      model%parameters(1) = 1.0_dp
    case (bicop_gumbel)
      model%parameters(1) = 1.5_dp
    case (bicop_frank)
      model%parameters(1) = 2.0_dp
    case (bicop_joe)
      model%parameters(1) = 1.5_dp
    case (bicop_bb1)
      model%parameters(1:2) = [1.0_dp, 1.5_dp]
    case (bicop_bb6)
      model%parameters(1:2) = [1.5_dp, 1.5_dp]
    case (bicop_bb7)
      model%parameters(1:2) = [1.5_dp, 1.0_dp]
    case (bicop_bb8)
      model%parameters(1:2) = [1.5_dp, 0.8_dp]
    case (bicop_tawn)
      model%parameters(1:3) = [1.0_dp, 1.0_dp, 1.5_dp]
    end select
  end subroutine default_parameters

  pure integer function family_parameter_count(family) result(n)
    integer, intent(in) :: family
    select case (family)
    case (bicop_indep)
      n = 0
    case (bicop_gaussian, bicop_clayton, bicop_gumbel, bicop_frank, bicop_joe)
      n = 1
    case (bicop_student, bicop_bb1, bicop_bb6, bicop_bb7, bicop_bb8)
      n = 2
    case (bicop_tawn)
      n = 3
    case default
      n = 0
    end select
  end function family_parameter_count

  subroutine family_bounds(family, lower, upper)
    integer, intent(in) :: family
    real(dp), intent(out) :: lower(3), upper(3)
    lower = 0.0_dp
    upper = 0.0_dp
    select case (family)
    case (bicop_gaussian)
      lower(1) = -0.995_dp; upper(1) = 0.995_dp
    case (bicop_student)
      lower(1:2) = [-0.995_dp, 2.01_dp]
      upper(1:2) = [0.995_dp, 50.0_dp]
    case (bicop_clayton)
      lower(1) = 1.0e-4_dp; upper(1) = 30.0_dp
    case (bicop_gumbel)
      lower(1) = 1.0_dp; upper(1) = 30.0_dp
    case (bicop_frank)
      lower(1) = -35.0_dp; upper(1) = 35.0_dp
    case (bicop_joe)
      lower(1) = 1.0_dp; upper(1) = 30.0_dp
    case (bicop_bb1)
      lower(1:2) = [1.0e-4_dp, 1.0_dp]
      upper(1:2) = [7.0_dp, 7.0_dp]
    case (bicop_bb6)
      lower(1:2) = [1.0_dp, 1.0_dp]
      upper(1:2) = [6.0_dp, 8.0_dp]
    case (bicop_bb7)
      lower(1:2) = [1.0_dp, 0.01_dp]
      upper(1:2) = [6.0_dp, 25.0_dp]
    case (bicop_bb8)
      lower(1:2) = [1.0_dp, 1.0e-4_dp]
      upper(1:2) = [8.0_dp, 1.0_dp]
    case (bicop_tawn)
      lower(1:3) = [0.0_dp, 0.0_dp, 1.0_dp]
      upper(1:3) = [1.0_dp, 1.0_dp, 30.0_dp]
    end select
  end subroutine family_bounds

  pure function bicop_name(self) result(name)
    class(bicop_model), intent(in) :: self
    character(len=16) :: name
    select case (self%family)
    case (bicop_indep); name = 'independence'
    case (bicop_gaussian); name = 'gaussian'
    case (bicop_student); name = 'student'
    case (bicop_clayton); name = 'clayton'
    case (bicop_gumbel); name = 'gumbel'
    case (bicop_frank); name = 'frank'
    case (bicop_joe); name = 'joe'
    case (bicop_bb1); name = 'bb1'
    case (bicop_bb6); name = 'bb6'
    case (bicop_bb7); name = 'bb7'
    case (bicop_bb8); name = 'bb8'
    case (bicop_tawn); name = 'tawn'
    case default; name = 'unknown'
    end select
  end function bicop_name

  pure real(dp) function bicop_aic(self) result(value)
    class(bicop_model), intent(in) :: self
    value = -2.0_dp*self%loglik + 2.0_dp*real(self%npar,dp)
  end function bicop_aic

  pure real(dp) function bicop_bic(self) result(value)
    class(bicop_model), intent(in) :: self
    if (self%nobs > 0) then
      value = -2.0_dp*self%loglik + log(real(self%nobs,dp))*real(self%npar,dp)
    else
      value = huge(1.0_dp)
    end if
  end function bicop_bic

  real(dp) function bicop_cdf(self, u, v) result(value)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: u, v
    real(dp) :: x, y
    x = clamp_prob(u)
    y = clamp_prob(v)
    select case (self%rotation)
    case (0)
      value = base_cdf(self, x, y)
    case (90)
      value = y - base_cdf(self, 1.0_dp-x, y)
    case (180)
      value = x + y - 1.0_dp + base_cdf(self, 1.0_dp-x, 1.0_dp-y)
    case (270)
      value = x - base_cdf(self, x, 1.0_dp-y)
    case default
      value = base_cdf(self, x, y)
    end select
    value = min(min(x,y), max(max(0.0_dp,x+y-1.0_dp), value))
  end function bicop_cdf

  real(dp) function bicop_pdf(self, u, v) result(value)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: u, v
    real(dp) :: x, y
    x = clamp_prob(u)
    y = clamp_prob(v)
    select case (self%rotation)
    case (0)
      value = base_pdf(self, x, y)
    case (90)
      value = base_pdf(self, 1.0_dp-x, y)
    case (180)
      value = base_pdf(self, 1.0_dp-x, 1.0_dp-y)
    case (270)
      value = base_pdf(self, x, 1.0_dp-y)
    case default
      value = base_pdf(self, x, y)
    end select
    value = max(tiny(1.0_dp), value)
  end function bicop_pdf

  real(dp) function bicop_hfunc1(self, u, v) result(value)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: u, v
    real(dp) :: x, y
    x = clamp_prob(u)
    y = clamp_prob(v)
    select case (self%rotation)
    case (0)
      value = base_h1(self, x, y)
    case (90)
      value = base_h1(self, 1.0_dp-x, y)
    case (180)
      value = 1.0_dp - base_h1(self, 1.0_dp-x, 1.0_dp-y)
    case (270)
      value = 1.0_dp - base_h1(self, x, 1.0_dp-y)
    case default
      value = base_h1(self, x, y)
    end select
    value = clamp_prob(value)
  end function bicop_hfunc1

  real(dp) function bicop_hfunc2(self, u, v) result(value)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: u, v
    real(dp) :: x, y
    x = clamp_prob(u)
    y = clamp_prob(v)
    select case (self%rotation)
    case (0)
      value = base_h2(self, x, y)
    case (90)
      value = 1.0_dp - base_h2(self, 1.0_dp-x, y)
    case (180)
      value = 1.0_dp - base_h2(self, 1.0_dp-x, 1.0_dp-y)
    case (270)
      value = base_h2(self, x, 1.0_dp-y)
    case default
      value = base_h2(self, x, y)
    end select
    value = clamp_prob(value)
  end function bicop_hfunc2

  real(dp) function bicop_hinv1(self, u, p) result(v)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: u, p
    real(dp) :: lo, hi, mid, target
    integer :: iter
    lo = eps_prob
    hi = 1.0_dp - eps_prob
    target = clamp_prob(p)
    do iter = 1, 70
      mid = 0.5_dp*(lo+hi)
      if (self%hfunc1(u,mid) < target) then
        lo = mid
      else
        hi = mid
      end if
    end do
    v = 0.5_dp*(lo+hi)
  end function bicop_hinv1

  real(dp) function bicop_hinv2(self, p, v) result(u)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: p, v
    real(dp) :: lo, hi, mid, target
    integer :: iter
    lo = eps_prob
    hi = 1.0_dp - eps_prob
    target = clamp_prob(p)
    do iter = 1, 70
      mid = 0.5_dp*(lo+hi)
      if (self%hfunc2(mid,v) < target) then
        lo = mid
      else
        hi = mid
      end if
    end do
    u = 0.5_dp*(lo+hi)
  end function bicop_hinv2

  subroutine bicop_simulate(self, n, sample)
    class(bicop_model), intent(in) :: self
    integer, intent(in) :: n
    real(dp), intent(out) :: sample(2,n)
    integer :: i
    real(dp) :: u, w
    do i = 1, n
      u = random_uniform()
      w = random_uniform()
      sample(1,i) = u
      sample(2,i) = self%hinv1(u,w)
    end do
  end subroutine bicop_simulate

  real(dp) function bicop_tau(self) result(tau)
    class(bicop_model), intent(in) :: self
    real(dp) :: sgn, a, b, xm, xl, xx
    real(dp), allocatable :: nodes(:), weights(:)
    integer :: i
    sgn = 1.0_dp
    if (self%rotation == 90 .or. self%rotation == 270) sgn = -1.0_dp
    select case (self%family)
    case (bicop_indep)
      tau = 0.0_dp
    case (bicop_gaussian, bicop_student)
      tau = 2.0_dp/pi * asin(max(-1.0_dp,min(1.0_dp,self%parameters(1))))
    case (bicop_clayton)
      tau = self%parameters(1)/(self%parameters(1)+2.0_dp)
    case (bicop_gumbel)
      tau = 1.0_dp - 1.0_dp/self%parameters(1)
    case (bicop_bb1)
      tau = 1.0_dp - 2.0_dp / &
            (self%parameters(2)*(self%parameters(1)+2.0_dp))
    case (bicop_tawn)
      call gauss_legendre_rule(80,nodes,weights)
      tau = 0.0_dp
      do i=1,size(nodes)
        xx=0.5_dp*(nodes(i)+1.0_dp)
        tau=tau+0.5_dp*weights(i)*xx*(1.0_dp-xx)* &
            tawn_pickands_d2(self,xx)/max(tawn_pickands(self,xx),tiny(1.0_dp))
      end do
    case default
      call gauss_legendre_rule(80,nodes,weights)
      a=1.0e-8_dp; b=1.0_dp-1.0e-8_dp
      xm=0.5_dp*(a+b); xl=0.5_dp*(b-a); tau=0.0_dp
      do i=1,size(nodes)
        xx=xm+xl*nodes(i)
        tau=tau+weights(i)*generator(self,xx)/generator_d1(self,xx)
      end do
      tau=1.0_dp+4.0_dp*xl*tau
    end select
    tau = sgn*tau
  end function bicop_tau

  real(dp) function base_cdf(self, u, v) result(c)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: u, v
    real(dp) :: s, a, b, xm, xl, xx
    real(dp), allocatable :: nodes(:), weights(:)
    integer :: i
    select case (self%family)
    case (bicop_indep)
      c = u*v
    case (bicop_gaussian, bicop_student)
      call gauss_legendre_rule(72,nodes,weights)
      a=eps_prob; b=u; xm=0.5_dp*(a+b); xl=0.5_dp*(b-a); c=0.0_dp
      do i=1,size(nodes)
        xx=xm+xl*nodes(i)
        c=c+weights(i)*base_h1(self,xx,v)
      end do
      c=xl*c
    case (bicop_tawn)
      c = tawn_cdf(self,u,v)
    case default
      s = generator(self,u) + generator(self,v)
      c = generator_inv(self,s)
    end select
  end function base_cdf

  real(dp) function base_pdf(self, u, v) result(c)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: u, v
    real(dp) :: x, y, rho, nu, q, logjoint, logmarg, cc, p1, p2, p3
    select case (self%family)
    case (bicop_indep)
      c = 1.0_dp
    case (bicop_gaussian)
      rho = self%parameters(1)
      x = normal_quantile(u)
      y = normal_quantile(v)
      c = exp((2.0_dp*rho*x*y-rho*rho*(x*x+y*y)) / &
              (2.0_dp*(1.0_dp-rho*rho))) / sqrt(1.0_dp-rho*rho)
    case (bicop_student)
      rho = self%parameters(1)
      nu = self%parameters(2)
      x = student_quantile(u,nu)
      y = student_quantile(v,nu)
      q = (x*x - 2.0_dp*rho*x*y + y*y)/(1.0_dp-rho*rho)
      logjoint = log_gamma(0.5_dp*(nu+2.0_dp)) - log_gamma(0.5_dp*nu) - &
                 log(nu*pi) - 0.5_dp*log(1.0_dp-rho*rho) - &
                 0.5_dp*(nu+2.0_dp)*log(1.0_dp + q/nu)
      logmarg = log(student_pdf(x,nu)) + log(student_pdf(y,nu))
      c = exp(logjoint-logmarg)
    case (bicop_tawn)
      c = tawn_pdf(self,u,v)
    case default
      cc = generator_inv(self,generator(self,u)+generator(self,v))
      p1 = generator_d1(self,u)
      p2 = generator_d1(self,v)
      p3 = generator_d1(self,cc)
      c = -p1*p2*generator_d2(self,cc)/(p3*p3*p3)
    end select
    if (.not. (c > 0.0_dp)) c = tiny(1.0_dp)
  end function base_pdf

  real(dp) function base_h1(self, u, v) result(h)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: u, v
    real(dp) :: x, y, rho, nu, cc
    select case (self%family)
    case (bicop_indep)
      h = v
    case (bicop_gaussian)
      rho = self%parameters(1)
      x = normal_quantile(u)
      y = normal_quantile(v)
      h = normal_cdf((y-rho*x)/sqrt(1.0_dp-rho*rho))
    case (bicop_student)
      rho = self%parameters(1)
      nu = self%parameters(2)
      x = student_quantile(u,nu)
      y = student_quantile(v,nu)
      h = student_cdf((y-rho*x)*sqrt((nu+1.0_dp)/ &
          ((nu+x*x)*(1.0_dp-rho*rho))),nu+1.0_dp)
    case (bicop_tawn)
      h = tawn_h1(self,u,v)
    case default
      cc = generator_inv(self,generator(self,u)+generator(self,v))
      h = generator_d1(self,u)/generator_d1(self,cc)
    end select
  end function base_h1

  real(dp) function base_h2(self, u, v) result(h)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: u, v
    real(dp) :: x, y, rho, nu, cc
    select case (self%family)
    case (bicop_indep)
      h = u
    case (bicop_gaussian)
      rho = self%parameters(1)
      x = normal_quantile(u)
      y = normal_quantile(v)
      h = normal_cdf((x-rho*y)/sqrt(1.0_dp-rho*rho))
    case (bicop_student)
      rho = self%parameters(1)
      nu = self%parameters(2)
      x = student_quantile(u,nu)
      y = student_quantile(v,nu)
      h = student_cdf((x-rho*y)*sqrt((nu+1.0_dp)/ &
          ((nu+y*y)*(1.0_dp-rho*rho))),nu+1.0_dp)
    case (bicop_tawn)
      h = tawn_h2(self,u,v)
    case default
      cc = generator_inv(self,generator(self,u)+generator(self,v))
      h = generator_d1(self,v)/generator_d1(self,cc)
    end select
  end function base_h2

  real(dp) function generator(self, u) result(phi)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: u
    real(dp) :: th, de, a
    th = self%parameters(1)
    de = self%parameters(2)
    select case (self%family)
    case (bicop_clayton)
      phi = u**(-th)-1.0_dp
    case (bicop_gumbel)
      phi = (-log(u))**th
    case (bicop_frank)
      if (abs(th) < 1.0e-7_dp) then
        phi = -log(u)
      else
        phi = -log((exp(-th*u)-1.0_dp)/(exp(-th)-1.0_dp))
      end if
    case (bicop_joe)
      phi = -log(1.0_dp-(1.0_dp-u)**th)
    case (bicop_bb1)
      phi = (u**(-th)-1.0_dp)**de
    case (bicop_bb6)
      phi = (-log(1.0_dp-(1.0_dp-u)**th))**de
    case (bicop_bb7)
      phi = (1.0_dp-(1.0_dp-u)**th)**(-de)-1.0_dp
    case (bicop_bb8)
      a = 1.0_dp-(1.0_dp-de*u)**th
      phi = -log(a/(1.0_dp-(1.0_dp-de)**th))
    case default
      phi = -log(u)
    end select
  end function generator

  real(dp) function generator_inv(self, s) result(u)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: s
    real(dp) :: th, de, a
    th = self%parameters(1)
    de = self%parameters(2)
    select case (self%family)
    case (bicop_clayton)
      u = (1.0_dp+s)**(-1.0_dp/th)
    case (bicop_gumbel)
      u = exp(-s**(1.0_dp/th))
    case (bicop_frank)
      if (abs(th) < 1.0e-7_dp) then
        u = exp(-s)
      else
        u = -log(1.0_dp + exp(-s)*(exp(-th)-1.0_dp))/th
      end if
    case (bicop_joe)
      u = 1.0_dp-(1.0_dp-exp(-s))**(1.0_dp/th)
    case (bicop_bb1)
      u = (s**(1.0_dp/de)+1.0_dp)**(-1.0_dp/th)
    case (bicop_bb6)
      u = 1.0_dp-(1.0_dp-exp(-s**(1.0_dp/de)))**(1.0_dp/th)
    case (bicop_bb7)
      u = 1.0_dp-(1.0_dp-(1.0_dp+s)**(-1.0_dp/de))**(1.0_dp/th)
    case (bicop_bb8)
      a = exp(-s)*((1.0_dp-de)**th-1.0_dp)
      u = (1.0_dp-(1.0_dp+a)**(1.0_dp/th))/de
    case default
      u = exp(-s)
    end select
    u = clamp_prob(u)
  end function generator_inv

  real(dp) function generator_d1(self, u) result(d1)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: u
    real(dp) :: th, de, a, b
    th = self%parameters(1)
    de = self%parameters(2)
    select case (self%family)
    case (bicop_clayton)
      d1 = -th*u**(-th-1.0_dp)
    case (bicop_gumbel)
      d1 = -th*(-log(u))**(th-1.0_dp)/u
    case (bicop_frank)
      if (abs(th) < 1.0e-7_dp) then
        d1 = -1.0_dp/u
      else
        d1 = th/(1.0_dp-exp(th*u))
      end if
    case (bicop_joe)
      a = 1.0_dp-u
      b = 1.0_dp-a**th
      d1 = -th*a**(th-1.0_dp)/b
    case (bicop_bb1)
      d1 = -de*th*u**(-1.0_dp-th)*(u**(-th)-1.0_dp)**(de-1.0_dp)
    case (bicop_bb6)
      a = log(1.0_dp-(1.0_dp-u)**th)
      d1 = de*th*(-a)**(de-1.0_dp)*(1.0_dp-u)**(th-1.0_dp) / &
           ((1.0_dp-u)**th-1.0_dp)
    case (bicop_bb7)
      a = 1.0_dp-(1.0_dp-u)**th
      d1 = -de*th*a**(-1.0_dp-de)*(1.0_dp-u)**(th-1.0_dp)
    case (bicop_bb8)
      a = 1.0_dp-de*u
      d1 = -de*th*a**(th-1.0_dp)/(1.0_dp-a**th)
    case default
      d1 = -1.0_dp/u
    end select
  end function generator_d1

  real(dp) function generator_d2(self, u) result(d2)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: u
    real(dp) :: th, de, a, b, tmp, res
    th = self%parameters(1)
    de = self%parameters(2)
    select case (self%family)
    case (bicop_clayton)
      d2 = th*(th+1.0_dp)*u**(-th-2.0_dp)
    case (bicop_gumbel)
      a = -log(u)
      d2 = th*a**(th-2.0_dp)*(a+th-1.0_dp)/(u*u)
    case (bicop_frank)
      if (abs(th) < 1.0e-7_dp) then
        d2 = 1.0_dp/(u*u)
      else
        d2 = th*th*exp(th*u)/(1.0_dp-exp(th*u))**2
      end if
    case (bicop_joe)
      a = 1.0_dp-u
      b = 1.0_dp-a**th
      d2 = th*(th-1.0_dp)*a**(th-2.0_dp)/b + &
           th*th*a**(2.0_dp*th-2.0_dp)/(b*b)
    case (bicop_bb1)
      res = de*th*(u**(-th)-1.0_dp)**de
      d2 = res*(1.0_dp+de*th-(1.0_dp+th)*u**th) / &
           ((u**th-1.0_dp)**2*u*u)
    case (bicop_bb6)
      tmp = (1.0_dp-u)**th
      a = log(1.0_dp-tmp)
      res = (-a)**(de-2.0_dp)
      res = res*((de-1.0_dp)*th*tmp-(tmp+th-1.0_dp)*a)
      d2 = res*de*th*(1.0_dp-u)**(th-2.0_dp)/(tmp-1.0_dp)**2
    case (bicop_bb7)
      tmp = (1.0_dp-u)**th
      res = de*th*(1.0_dp-tmp)**(-2.0_dp-de)*(1.0_dp-u)**(th-2.0_dp)
      d2 = res*(th-1.0_dp+(1.0_dp+de*th)*tmp)
    case (bicop_bb8)
      tmp = (1.0_dp-de*u)**th
      res = de*de*th*(1.0_dp-de*u)**(th-2.0_dp)
      d2 = res*(th-1.0_dp+tmp)/(tmp-1.0_dp)**2
    case default
      d2 = 1.0_dp/(u*u)
    end select
  end function generator_d2

  pure real(dp) function tawn_pickands(self,t) result(a)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: t
    real(dp) :: p1,p2,th,temp
    p1=self%parameters(1); p2=self%parameters(2); th=self%parameters(3)
    temp=(p2*t)**th+(p1*(1.0_dp-t))**th
    a=(1.0_dp-p1)*(1.0_dp-t)+(1.0_dp-p2)*t+temp**(1.0_dp/th)
  end function tawn_pickands

  pure real(dp) function tawn_pickands_d1(self,t) result(a1)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: t
    real(dp) :: p1,p2,th,temp,temp2
    p1=self%parameters(1); p2=self%parameters(2); th=self%parameters(3)
    temp=(p2*t)**th+(p1*(1.0_dp-t))**th
    temp2=p2*(p2*t)**(th-1.0_dp)-p1*(p1*(1.0_dp-t))**(th-1.0_dp)
    a1=p1-p2+temp**(1.0_dp/th-1.0_dp)*temp2
  end function tawn_pickands_d1

  pure real(dp) function tawn_pickands_d2(self,t) result(a2)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: t
    real(dp) :: p1,p2,th,temp,temp2,temp3
    p1=self%parameters(1); p2=self%parameters(2); th=self%parameters(3)
    temp=(p2*t)**th+(p1*(1.0_dp-t))**th
    temp2=p2*(p2*t)**(th-1.0_dp)-p1*(p1*(1.0_dp-t))**(th-1.0_dp)
    temp3=p2*p2*(p2*t)**(th-2.0_dp)+p1*p1*(p1*(1.0_dp-t))**(th-2.0_dp)
    a2=(1.0_dp-th)*temp**(1.0_dp/th-2.0_dp)*temp2*temp2 + &
       temp**(1.0_dp/th-1.0_dp)*(th-1.0_dp)*temp3
  end function tawn_pickands_d2

  real(dp) function tawn_cdf(self,u,v) result(c)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: u,v
    real(dp) :: t,s
    s=log(u)+log(v)
    t=log(v)/s
    c=exp(s*tawn_pickands(self,t))
  end function tawn_cdf

  real(dp) function tawn_h1(self,u,v) result(h)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: u,v
    real(dp) :: t,a,a1,c
    t=log(v)/log(u*v)
    a=tawn_pickands(self,t)
    a1=tawn_pickands_d1(self,t)
    c=exp(log(u*v)*a)
    h=c*(a-t*a1)/u
  end function tawn_h1

  real(dp) function tawn_h2(self,u,v) result(h)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: u,v
    real(dp) :: t,a,a1,c
    t=log(v)/log(u*v)
    a=tawn_pickands(self,t)
    a1=tawn_pickands_d1(self,t)
    c=exp(log(u*v)*a)
    h=c*(a+(1.0_dp-t)*a1)/v
  end function tawn_h2

  real(dp) function tawn_pdf(self,u,v) result(c)
    class(bicop_model), intent(in) :: self
    real(dp), intent(in) :: u,v
    real(dp) :: t,a,a1,a2,s,term
    s=log(u*v)
    t=log(v)/s
    a=tawn_pickands(self,t)
    a1=tawn_pickands_d1(self,t)
    a2=tawn_pickands_d2(self,t)
    term=a*a+(1.0_dp-2.0_dp*t)*a1*a - &
         (1.0_dp-t)*t*(a1*a1+a2/s)
    c=exp(s*a)*term/(u*v)
  end function tawn_pdf

end module rvine_bicop

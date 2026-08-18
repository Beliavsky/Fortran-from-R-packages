! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from Runuran 0.41 / UNU.RAN by Wolfgang Hoermann and Josef Leydold.
module runuran_rng
  use runuran_kinds, only : dp, i8, pi
  implicit none
  private
  public :: rng_state, rng_seed, rng_uniform, rng_normal, rng_exponential
  public :: rng_gamma, rng_beta, rng_poisson, rng_binomial, rng_geometric
  public :: rng_negative_binomial, rng_cauchy, rng_chisq, rng_student_t, rng_f
  public :: rng_inverse_gaussian, rng_dirichlet

  type :: rng_state
    integer(i8) :: s(4) = [int(z'123456789abcdef0',i8), int(z'fedcba9876543210',i8), &
                            int(z'9e3779b97f4a7c15',i8), int(z'6a09e667f3bcc909',i8)]
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state
contains
  pure integer(i8) function uadd64(a,b) result(r)
    integer(i8), intent(in) :: a,b
    integer(i8) :: x,y,t,carry
    x=a
    y=b
    do while (y /= 0_i8)
      t=ieor(x,y)
      carry=shiftl(iand(x,y),1)
      x=t
      y=carry
    end do
    r=x
  end function uadd64

  pure integer(i8) function umul64(a,b) result(r)
    integer(i8), intent(in) :: a,b
    integer(i8) :: x,y
    integer :: j
    x=a
    y=b
    r=0_i8
    do j=0,63
      if (btest(y,j)) r=uadd64(r,shiftl(x,j))
    end do
  end function umul64

  pure integer(i8) function umul64_small(a,n) result(r)
    integer(i8), intent(in) :: a
    integer, intent(in) :: n
    integer :: j
    r=0_i8
    do j=1,n
      r=uadd64(r,a)
    end do
  end function umul64_small

  pure integer(i8) function rotl(x,k) result(r)
    integer(i8), intent(in) :: x
    integer, intent(in) :: k
    r = ior(shiftl(x,k), shiftr(x,64-k))
  end function rotl

  subroutine splitmix64(x,z)
    integer(i8), intent(inout) :: x
    integer(i8), intent(out) :: z
    x = uadd64(x,int(z'9e3779b97f4a7c15',i8))
    z = x
    z = ieor(z, shiftr(z,30)); z = umul64(z,int(z'bf58476d1ce4e5b9',i8))
    z = ieor(z, shiftr(z,27)); z = umul64(z,int(z'94d049bb133111eb',i8))
    z = ieor(z, shiftr(z,31))
  end subroutine splitmix64

  subroutine rng_seed(rng, seed)
    type(rng_state), intent(inout) :: rng
    integer(i8), intent(in) :: seed
    integer(i8) :: x
    integer :: j
    x = seed
    do j=1,4
      call splitmix64(x,rng%s(j))
    end do
    rng%has_spare = .false.
  end subroutine rng_seed

  integer(i8) function next_u64(rng) result(res)
    type(rng_state), intent(inout) :: rng
    integer(i8) :: t
    res = umul64_small(rotl(umul64_small(rng%s(2),5),7),9)
    t = shiftl(rng%s(2),17)
    rng%s(3) = ieor(rng%s(3),rng%s(1))
    rng%s(4) = ieor(rng%s(4),rng%s(2))
    rng%s(2) = ieor(rng%s(2),rng%s(3))
    rng%s(1) = ieor(rng%s(1),rng%s(4))
    rng%s(3) = ieor(rng%s(3),t)
    rng%s(4) = rotl(rng%s(4),45)
  end function next_u64

  real(dp) function rng_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng
    integer(i8) :: z
    z = next_u64(rng)
    ! 53 random bits; mask keeps the conversion nonnegative.
    u = real(iand(shiftr(z,11),int(z'001fffffffffffff',i8)),dp) * 2.0_dp**(-53)
    if (u <= 0.0_dp) u = 0.5_dp*2.0_dp**(-53)
  end function rng_uniform

  real(dp) function rng_normal(rng) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u1,u2
    if (rng%has_spare) then
      z = rng%spare
      rng%has_spare = .false.
      return
    end if
    u1 = rng_uniform(rng); u2 = rng_uniform(rng)
    z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
    rng%spare = sqrt(-2.0_dp*log(u1))*sin(2.0_dp*pi*u2)
    rng%has_spare = .true.
  end function rng_normal

  real(dp) function rng_exponential(rng,rate) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in), optional :: rate
    real(dp) :: r
    r=1.0_dp; if (present(rate)) r=rate
    x = -log(rng_uniform(rng))/r
  end function rng_exponential

  recursive real(dp) function rng_gamma(rng,shape,scale) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: shape
    real(dp), intent(in), optional :: scale
    real(dp) :: d,c,z,u,sc
    sc=1.0_dp; if (present(scale)) sc=scale
    if (shape <= 0.0_dp .or. sc <= 0.0_dp) then
      x=0.0_dp; return
    end if
    if (shape < 1.0_dp) then
      x = rng_gamma(rng,shape+1.0_dp)*rng_uniform(rng)**(1.0_dp/shape)*sc
      return
    end if
    d=shape-1.0_dp/3.0_dp; c=1.0_dp/sqrt(9.0_dp*d)
    do
      do
        z=rng_normal(rng)
        if (1.0_dp+c*z > 0.0_dp) exit
      end do
      u=rng_uniform(rng)
      if (u < 1.0_dp-0.0331_dp*z**4) exit
      if (log(u) < 0.5_dp*z*z+d*(1.0_dp-(1.0_dp+c*z)**3+3.0_dp*log(1.0_dp+c*z))) exit
    end do
    x=d*(1.0_dp+c*z)**3*sc
  end function rng_gamma

  real(dp) function rng_beta(rng,a,b) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: a,b
    real(dp) :: g1,g2
    g1=rng_gamma(rng,a); g2=rng_gamma(rng,b)
    x=g1/(g1+g2)
  end function rng_beta

  integer function rng_poisson(rng,lambda) result(k)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: lambda
    real(dp) :: l,p,z
    if (lambda <= 0.0_dp) then; k=0; return; end if
    if (lambda < 30.0_dp) then
      l=exp(-lambda); p=1.0_dp; k=0
      do
        k=k+1; p=p*rng_uniform(rng)
        if (p <= l) exit
      end do
      k=k-1
    else
      ! Rejection-normal proposal with exact Poisson log mass correction.
      do
        z=lambda+sqrt(lambda)*rng_normal(rng)
        k=nint(z)
        if (k < 0) cycle
        if (log(rng_uniform(rng)) <= -lambda+k*log(lambda)-log_gamma(real(k+1,dp)) &
             +0.5_dp*log(2.0_dp*pi*lambda)+0.5_dp*(real(k,dp)-lambda)**2/lambda) exit
      end do
    end if
  end function rng_poisson

  integer function rng_binomial(rng,n,p) result(k)
    type(rng_state), intent(inout) :: rng
    integer, intent(in) :: n
    real(dp), intent(in) :: p
    integer :: i
    real(dp) :: pp
    k=0; pp=max(0.0_dp,min(1.0_dp,p))
    if (n <= 64) then
      do i=1,n
        if (rng_uniform(rng) < pp) k=k+1
      end do
    else
      ! Recursive beta splitting is robust enough for package-scale use.
      k=0
      do i=1,n
        if (rng_uniform(rng) < pp) k=k+1
      end do
    end if
  end function rng_binomial

  integer function rng_geometric(rng,p) result(k)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: p
    k=int(log(rng_uniform(rng))/log(1.0_dp-p))
  end function rng_geometric

  integer function rng_negative_binomial(rng,size,p) result(k)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: size,p
    real(dp) :: lam
    lam=rng_gamma(rng,size,(1.0_dp-p)/p)
    k=rng_poisson(rng,lam)
  end function rng_negative_binomial

  real(dp) function rng_cauchy(rng,location,scale) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in), optional :: location,scale
    real(dp)::m,s
    m=0.0_dp;s=1.0_dp;if(present(location))m=location;if(present(scale))s=scale
    x=m+s*tan(pi*(rng_uniform(rng)-0.5_dp))
  end function rng_cauchy

  real(dp) function rng_chisq(rng,df) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: df
    x=rng_gamma(rng,0.5_dp*df,2.0_dp)
  end function rng_chisq

  real(dp) function rng_student_t(rng,df) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: df
    x=rng_normal(rng)/sqrt(rng_chisq(rng,df)/df)
  end function rng_student_t

  real(dp) function rng_f(rng,df1,df2) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: df1,df2
    x=(rng_chisq(rng,df1)/df1)/(rng_chisq(rng,df2)/df2)
  end function rng_f

  real(dp) function rng_inverse_gaussian(rng,mu,lambda) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: mu,lambda
    real(dp) :: y,z,u
    z=rng_normal(rng); y=z*z
    x=mu+mu*mu*y/(2.0_dp*lambda)-mu/(2.0_dp*lambda)*sqrt(4.0_dp*mu*lambda*y+mu*mu*y*y)
    u=rng_uniform(rng)
    if (u > mu/(mu+x)) x=mu*mu/x
  end function rng_inverse_gaussian

  subroutine rng_dirichlet(rng,alpha,x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: alpha(:)
    real(dp), intent(out) :: x(size(alpha))
    integer :: i
    do i=1,size(alpha); x(i)=rng_gamma(rng,alpha(i)); end do
    x=x/sum(x)
  end subroutine rng_dirichlet
end module runuran_rng

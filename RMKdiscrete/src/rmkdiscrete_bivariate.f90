! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of RMKdiscrete 0.1 by Robert M. Kirkpatrick.
! See LICENSE, COPYING, and TRANSLATION_NOTES.md for provenance.
module rmkdiscrete_bivariate
  use rmkdiscrete_kinds, only : dp
  use rmkdiscrete_math, only : qnan, ninf, real_equal, kahan_add, dnbinom_prob, rnbinom_prob
  use rmkdiscrete_lgp, only : dlgp, lgp_get_nc, lgp_findmax, rlgp
  implicit none
  private

  type, public :: log_moments2
    real(dp) :: ey1=0.0_dp
    real(dp) :: ey2=0.0_dp
    real(dp) :: vy1=0.0_dp
    real(dp) :: vy2=0.0_dp
    real(dp) :: cov=0.0_dp
  end type log_moments2

  public :: dbilgp, rbilgp, rbilgp_sample, bilgp_logmv
  public :: dbinegbin, rbinegbin, rbinegbin_sample, binegbin_logmv

contains

  real(dp) function dbilgp(x,y,theta,lambda,nc,give_log) result(v)
    integer, intent(in) :: x,y
    real(dp), intent(in) :: theta(3),lambda(3)
    real(dp), intent(in), optional :: nc(3)
    logical, intent(in), optional :: give_log
    real(dp) :: z(3), s,c,term
    integer :: u,umax
    logical :: gl
    gl=.false.
    if(present(give_log)) gl=give_log
    if (x < 0 .or. y < 0 .or. any(theta<0.0_dp) .or. any(abs(lambda)>1.0_dp)) then
      if (x < 0 .or. y < 0) then
        v=merge(ninf(),0.0_dp,gl)
      else
        v=qnan()
      end if
      return
    end if
    if(present(nc)) then
      z=nc
    else
      z(1)=lgp_get_nc(theta(1),lambda(1))
      z(2)=lgp_get_nc(theta(2),lambda(2))
      z(3)=lgp_get_nc(theta(3),lambda(3))
    end if
    if(real_equal(theta(1),0.0_dp)) then
      if(gl) then
        v=dlgp(x,theta(2),lambda(2),z(2),.true.)+dlgp(y,theta(3),lambda(3),z(3),.true.)
      else
        v=dlgp(x,theta(2),lambda(2),z(2))*dlgp(y,theta(3),lambda(3),z(3))
      end if
      return
    end if
    umax=min(x,y)
    if(lambda(1)<0.0_dp) umax=min(umax,int(lgp_findmax(theta(1),lambda(1))))
    s=0.0_dp
    c=0.0_dp
    do u=0,umax
      term=exp(dlgp(x-u,theta(2),lambda(2),z(2),.true.) + &
               dlgp(y-u,theta(3),lambda(3),z(3),.true.) + &
               dlgp(u,theta(1),lambda(1),z(1),.true.))
      call kahan_add(s,c,term)
    end do
    if(gl) then
      v=merge(log(s),ninf(),s>0.0_dp)
    else
      v=s
    end if
  end function dbilgp

  subroutine rbilgp(theta,lambda,out)
    real(dp), intent(in) :: theta(3),lambda(3)
    integer, intent(out) :: out(2)
    integer :: u
    u=rlgp(theta(1),lambda(1))
    out(1)=u+rlgp(theta(2),lambda(2))
    out(2)=u+rlgp(theta(3),lambda(3))
  end subroutine rbilgp

  function rbilgp_sample(n,theta,lambda) result(out)
    integer, intent(in) :: n
    real(dp), intent(in) :: theta(3),lambda(3)
    integer, allocatable :: out(:,:)
    integer :: i
    allocate(out(max(0,n),2))
    do i=1,n
      call rbilgp(theta,lambda,out(i,:))
    end do
  end function rbilgp_sample

  real(dp) function lgp_conv_pmf(x,t0,t1,l0,l1,n0,n1) result(v)
    integer,intent(in)::x
    real(dp),intent(in)::t0,t1,l0,l1,n0,n1
    integer::u,umax
    real(dp)::s,c,term
    if(real_equal(l0,l1) .and. l0>=0.0_dp) then
      v=dlgp(x,t0+t1,l0,1.0_dp)
      return
    end if
    umax=x
    if(l0<0.0_dp) umax=min(umax,int(lgp_findmax(t0,l0)))
    s=0.0_dp
    c=0.0_dp
    do u=0,umax
      term=exp(dlgp(x-u,t1,l1,n1,.true.)+dlgp(u,t0,l0,n0,.true.))
      call kahan_add(s,c,term)
    end do
    v=s
  end function lgp_conv_pmf

  function bilgp_logmv(theta,lambda,const_add,tol,nc) result(out)
    real(dp),intent(in)::theta(3),lambda(3)
    real(dp),intent(in),optional::const_add,tol,nc(3)
    type(log_moments2)::out
    real(dp)::ca,eps,z(3),p,old,ex,ey,ex2,ey2,exy
    integer::i,j,x,y,imax,jmax
    logical::past
    ca=1.0_dp
    if(present(const_add)) ca=const_add
    eps=1.0e-14_dp
    if(present(tol)) eps=tol
    if(ca<=0.0_dp .or. eps<=0.0_dp .or. any(theta<0.0_dp) .or. any(abs(lambda)>1.0_dp)) then
      out%ey1=qnan()
      out%ey2=qnan()
      out%vy1=qnan()
      out%vy2=qnan()
      out%cov=qnan()
      return
    end if
    if(present(nc)) then
      z=nc
    else
      do i=1,3
      z(i)=lgp_get_nc(theta(i),lambda(i))
      end do
    end if
    ex=0.0_dp
    ex2=0.0_dp
    old=0.0_dp
    past=.false.
    imax=0
    do i=0,1000000
      p=lgp_conv_pmf(i,theta(1),theta(2),lambda(1),lambda(2),z(1),z(2))
      if(p<old) past=.true.
      ex=ex+p*log(real(i,dp)+ca)
      ex2=ex2+p*log(real(i,dp)+ca)**2
      imax=i
      if(past .and. p*log(real(i,dp)+ca)**2<eps) exit
      old=p
    end do
    if(real_equal(theta(2),theta(3)) .and. real_equal(lambda(2),lambda(3))) then
      ey=ex
      ey2=ex2
      jmax=imax
    else
      ey=0.0_dp
      ey2=0.0_dp
      old=0.0_dp
      past=.false.
      jmax=0
      do j=0,1000000
        p=lgp_conv_pmf(j,theta(1),theta(3),lambda(1),lambda(3),z(1),z(3))
        if(p<old) past=.true.
        ey=ey+p*log(real(j,dp)+ca)
        ey2=ey2+p*log(real(j,dp)+ca)**2
        jmax=j
        if(past .and. p*log(real(j,dp)+ca)**2<eps) exit
        old=p
      end do
    end if
    exy=0.0_dp
    do x=0,imax
      do y=0,jmax
        exy=exy+dbilgp(x,y,theta,lambda,z)*log(real(x,dp)+ca)*log(real(y,dp)+ca)
      end do
    end do
    out%ey1=ex
    out%ey2=ey
    out%vy1=ex2-ex**2
    out%vy2=ey2-ey**2
    out%cov=exy-ex*ey
  end function bilgp_logmv

  real(dp) function dbinegbin(x,y,nu,p,give_log) result(v)
    integer,intent(in)::x,y
    real(dp),intent(in)::nu(3),p(3)
    logical,intent(in),optional::give_log
    logical::gl
    integer::u,umax
    real(dp)::s,c,term
    gl=.false.
    if(present(give_log)) gl=give_log
    if(x<0 .or. y<0) then
    v=merge(ninf(),0.0_dp,gl)
    return
    end if
    if(any(nu<0.0_dp) .or. any(p<=0.0_dp) .or. any(p>1.0_dp)) then
    v=qnan()
    return
    end if
    if(real_equal(nu(1),0.0_dp)) then
      if(gl) then
        v=dnbinom_prob(x,nu(2),p(2),.true.)+dnbinom_prob(y,nu(3),p(3),.true.)
      else
        v=dnbinom_prob(x,nu(2),p(2))*dnbinom_prob(y,nu(3),p(3))
      end if
      return
    end if
    umax=min(x,y)
    s=0.0_dp
    c=0.0_dp
    do u=0,umax
      term=exp(dnbinom_prob(x-u,nu(2),p(2),.true.)+dnbinom_prob(y-u,nu(3),p(3),.true.)+ &
               dnbinom_prob(u,nu(1),p(1),.true.))
      call kahan_add(s,c,term)
    end do
    if(gl) then
    v=merge(log(s),ninf(),s>0.0_dp)
    else
    v=s
    end if
  end function dbinegbin

  subroutine rbinegbin(nu,p,out)
    real(dp),intent(in)::nu(3),p(3)
    integer,intent(out)::out(2)
    integer::u
    u=rnbinom_prob(nu(1),p(1))
    out(1)=u+rnbinom_prob(nu(2),p(2))
    out(2)=u+rnbinom_prob(nu(3),p(3))
  end subroutine rbinegbin

  function rbinegbin_sample(n,nu,p) result(out)
    integer, intent(in) :: n
    real(dp), intent(in) :: nu(3),p(3)
    integer, allocatable :: out(:,:)
    integer :: i
    allocate(out(max(0,n),2))
    do i=1,n
      call rbinegbin(nu,p,out(i,:))
    end do
  end function rbinegbin_sample

  real(dp) function nb_conv_pmf(x,n0,n1,p0,p1) result(v)
    integer,intent(in)::x
    real(dp),intent(in)::n0,n1,p0,p1
    integer::u
    real(dp)::s,c,term
    if(real_equal(p0,p1)) then
    v=dnbinom_prob(x,n0+n1,p0)
    return
    end if
    s=0.0_dp
    c=0.0_dp
    do u=0,x
      term=exp(dnbinom_prob(x-u,n1,p1,.true.)+dnbinom_prob(u,n0,p0,.true.))
      call kahan_add(s,c,term)
    end do
    v=s
  end function nb_conv_pmf

  function binegbin_logmv(nu,p,const_add,tol) result(out)
    real(dp),intent(in)::nu(3),p(3)
    real(dp),intent(in),optional::const_add,tol
    type(log_moments2)::out
    real(dp)::ca,eps,pr,old,ex,ey,ex2,ey2,exy
    integer::i,j,x,y,imax,jmax
    logical::past
    ca=1.0_dp
    if(present(const_add)) ca=const_add
    eps=1.0e-14_dp
    if(present(tol)) eps=tol
    if(ca<=0.0_dp .or. eps<=0.0_dp .or. any(nu<0.0_dp) .or. any(p<=0.0_dp) .or. any(p>1.0_dp)) then
      out%ey1=qnan()
      out%ey2=qnan()
      out%vy1=qnan()
      out%vy2=qnan()
      out%cov=qnan()
      return
    end if
    ex=0.0_dp
    ex2=0.0_dp
    old=0.0_dp
    past=.false.
    imax=0
    do i=0,1000000
      pr=nb_conv_pmf(i,nu(1),nu(2),p(1),p(2))
      if(pr<old)past=.true.
      ex=ex+pr*log(real(i,dp)+ca)
      ex2=ex2+pr*log(real(i,dp)+ca)**2
      imax=i
      if(past .and. pr*log(real(i,dp)+ca)**2<eps)exit
      old=pr
    end do
    if(real_equal(nu(2),nu(3)) .and. real_equal(p(2),p(3))) then
      ey=ex
      ey2=ex2
      jmax=imax
    else
      ey=0.0_dp
      ey2=0.0_dp
      old=0.0_dp
      past=.false.
      jmax=0
      do j=0,1000000
        pr=nb_conv_pmf(j,nu(1),nu(3),p(1),p(3))
        if(pr<old)past=.true.
        ey=ey+pr*log(real(j,dp)+ca)
        ey2=ey2+pr*log(real(j,dp)+ca)**2
        jmax=j
        if(past .and. pr*log(real(j,dp)+ca)**2<eps)exit
        old=pr
      end do
    end if
    exy=0.0_dp
    do x=0,imax
      do y=0,jmax
        exy=exy+dbinegbin(x,y,nu,p)*log(real(x,dp)+ca)*log(real(y,dp)+ca)
      end do
    end do
    out%ey1=ex
    out%ey2=ey
    out%vy1=ex2-ex**2
    out%vy2=ey2-ey**2
    out%cov=exy-ex*ey
  end function binegbin_logmv
end module rmkdiscrete_bivariate

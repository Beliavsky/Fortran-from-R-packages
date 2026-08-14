! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_simulation
  use mclust_kinds, only : dp, pi_dp
  use mclust_types, only : mclust_fit
  implicit none
  private
  public :: simulate_mixture, simulate_fit, sample_multivariate_normal

contains

  subroutine simulate_fit(fit,n,x,component,status)
    type(mclust_fit),intent(in)::fit
    integer,intent(in)::n
    real(dp),allocatable,intent(out)::x(:,:)
    integer,allocatable,intent(out),optional::component(:)
    integer,intent(out),optional::status
    integer,allocatable::comp(:)
    integer::info
    allocate(x(n,fit%d),comp(n))
    call simulate_mixture(fit%pro,fit%mean,fit%sigma,x,comp,info)
    if(present(component)) call move_alloc(comp,component)
    if(present(status)) status=info
  end subroutine simulate_fit

  subroutine simulate_mixture(pro,mu,sigma,x,component,status)
    real(dp),intent(in)::pro(:),mu(:,:),sigma(:,:,:)
    real(dp),intent(out)::x(:,:)
    integer,intent(out)::component(:)
    integer,intent(out)::status
    real(dp),allocatable::cs(:),z(:)
    real(dp)::u
    integer::i,k,g,d,info
    g=size(pro); d=size(mu,1); status=0
    if(size(mu,2)/=g .or. size(sigma,1)/=d .or. size(sigma,2)/=d .or. &
       size(sigma,3)/=g .or. size(x,1)/=size(component) .or. size(x,2)/=d) then
      status=-1; return
    end if
    allocate(cs(g),z(d)); cs(1)=pro(1)
    do k=2,g; cs(k)=cs(k-1)+pro(k); end do
    if(cs(g)<=0.0_dp) then; status=-2; return; end if
    cs=cs/cs(g)
    do i=1,size(x,1)
      call random_number(u); k=1
      do while(k<g .and. u>cs(k)); k=k+1; end do
      component(i)=k
      call sample_multivariate_normal(mu(:,k),sigma(:,:,k),z,info)
      if(info/=0) then; status=10+k; return; end if
      x(i,:)=z
    end do
  end subroutine simulate_mixture

  subroutine sample_multivariate_normal(mu,sigma,x,status)
    real(dp),intent(in)::mu(:),sigma(:,:)
    real(dp),intent(out)::x(:)
    integer,intent(out)::status
    real(dp),allocatable::l(:,:),z(:)
    integer::d,i,j
    d=size(mu); allocate(l(d,d),z(d)); l=sigma
    call chol_lower(l,status); if(status/=0) return
    call normal_vector(z)
    x=mu
    do i=1,d
      do j=1,i
        x(i)=x(i)+l(i,j)*z(j)
      end do
    end do
  end subroutine sample_multivariate_normal

  subroutine normal_vector(z)
    real(dp),intent(out)::z(:)
    real(dp)::u1,u2,r
    integer::i
    i=1
    do while(i<=size(z))
      call random_number(u1); call random_number(u2)
      u1=max(u1,tiny(1.0_dp)); r=sqrt(-2.0_dp*log(u1))
      z(i)=r*cos(2.0_dp*pi_dp*u2)
      if(i+1<=size(z)) z(i+1)=r*sin(2.0_dp*pi_dp*u2)
      i=i+2
    end do
  end subroutine normal_vector

  subroutine chol_lower(a,info)
    real(dp),intent(inout)::a(:,:)
    integer,intent(out)::info
    integer::i,j,k,n
    real(dp)::s
    n=size(a,1); info=0
    do j=1,n
      s=a(j,j); do k=1,j-1; s=s-a(j,k)*a(j,k); end do
      if(s<=0.0_dp) then; info=j; return; end if
      a(j,j)=sqrt(s)
      do i=j+1,n
        s=a(i,j); do k=1,j-1; s=s-a(i,k)*a(j,k); end do
        a(i,j)=s/a(j,j)
      end do
      if(j<n) a(j,j+1:n)=0.0_dp
    end do
  end subroutine chol_lower
end module mclust_simulation

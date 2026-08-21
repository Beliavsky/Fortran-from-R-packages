! SPDX-License-Identifier: GPL-2.0-or-later
module dirichletreg_start
  use dirichletreg_kinds, only : dp
  use dirichletreg_special, only : digamma
  use dirichletreg_types, only : design_block
  use dirichletreg_common, only : common_npar
  use dirichletreg_alternative, only : alternative_predict
  use dirichletreg_linalg, only : weighted_least_squares
  use dirichletreg_optimize, only : maximize_bfgs
  implicit none
  private
  public :: common_starting_values, alternative_starting_values

contains

  subroutine common_starting_values(y,xblocks,w,start,stat)
    real(dp), intent(in) :: y(:,:),w(:)
    type(design_block), intent(in) :: xblocks(:)
    real(dp), intent(out) :: start(:)
    integer, intent(out), optional :: stat
    integer :: d,j,p,lo,hi,it,conv
    real(dp) :: f
    real(dp), allocatable :: b(:)

    if(present(stat)) stat=0
    d=size(y,2)
    if(size(xblocks)/=d .or. size(start)/=common_npar(xblocks) .or. size(w)/=size(y,1)) then
      if(present(stat)) stat=1
      start=0.0_dp
      return
    end if
    lo=1
    do j=1,d
      p=size(xblocks(j)%x,2); hi=lo+p-1
      allocate(b(2*p)); b=0.0_dp
      call maximize_bfgs(beta_obj,b,f,it,conv,iterlim=500,tol=1.0e-6_dp)
      start(lo:hi)=b(1:p)/real(d,dp)
      deallocate(b)
      lo=hi+1
    end do

  contains
    subroutine beta_obj(theta,f,g)
      real(dp), intent(in) :: theta(:)
      real(dp), intent(out) :: f,g(:)
      real(dp), allocatable :: aa(:),bb(:),ab(:),ly(:),l1y(:)
      integer :: i,v,pp
      pp=size(xblocks(j)%x,2)
      allocate(aa(size(y,1)),bb(size(y,1)),ab(size(y,1)),ly(size(y,1)),l1y(size(y,1)))
      aa=exp(matmul(xblocks(j)%x,theta(1:pp)))
      bb=exp(matmul(xblocks(j)%x,theta(pp+1:2*pp)))
      ab=aa+bb
      ly=log(y(:,j)); l1y=log(1.0_dp-y(:,j))
      f=0.0_dp; g=0.0_dp
      do i=1,size(y,1)
        f=f+w(i)*(log_gamma(ab(i))-log_gamma(aa(i))-log_gamma(bb(i))+(aa(i)-1.0_dp)*ly(i)+(bb(i)-1.0_dp)*l1y(i))
      end do
      do v=1,pp
        do i=1,size(y,1)
          g(v)=g(v)+w(i)*xblocks(j)%x(i,v)*aa(i)*(digamma(ab(i))-digamma(aa(i))+ly(i))
          g(pp+v)=g(pp+v)+w(i)*xblocks(j)%x(i,v)*bb(i)*(digamma(ab(i))-digamma(bb(i))+l1y(i))
        end do
      end do
    end subroutine beta_obj
  end subroutine common_starting_values


  subroutine alternative_starting_values(y,x,z,base,w,start,stat)
    real(dp), intent(in) :: y(:,:),x(:,:),z(:,:),w(:)
    integer, intent(in) :: base
    real(dp), intent(out) :: start(:)
    integer, intent(out), optional :: stat
    integer :: n,d,p,q,j,b,lo,hi,ierr
    real(dp) :: logphi
    real(dp), allocatable :: target(:),beta(:),mu(:,:),alpha(:,:),phi(:),dummytheta(:),gamma(:)

    if(present(stat)) stat=0
    n=size(y,1); d=size(y,2); p=size(x,2); q=size(z,2)
    if(size(x,1)/=n .or. size(z,1)/=n .or. size(w)/=n .or. base<1 .or. base>d .or. size(start)/=((d-1)*p+q)) then
      if(present(stat)) stat=1
      start=0.0_dp
      return
    end if
    allocate(target(n),beta(p),mu(n,d),alpha(n,d),phi(n),dummytheta((d-1)*p+q),gamma(q))
    start=0.0_dp; b=0
    do j=1,d
      if(j==base) cycle
      b=b+1; lo=(b-1)*p+1; hi=b*p
      target=log(y(:,j)/y(:,base))
      call weighted_least_squares(x,target,w,beta,ierr)
      if(ierr/=0) beta=0.0_dp
      start(lo:hi)=beta
    end do

    dummytheta=start
    dummytheta((d-1)*p+1:)=0.0_dp
    call alternative_predict(dummytheta,x,z,d,base,alpha,mu,phi,ierr)
    call maximize_logphi(mu,logphi)
    target=logphi
    call weighted_least_squares(z,target,w,gamma,ierr)
    if(ierr/=0) gamma=0.0_dp
    start((d-1)*p+1:)=gamma

  contains
    subroutine maximize_logphi(m,lp)
      real(dp), intent(in) :: m(:,:)
      real(dp), intent(out) :: lp
      real(dp), parameter :: gr=0.6180339887498948482_dp
      real(dp) :: a0,b0,c0,d0,fc,fd
      integer :: it2
      a0=-20.0_dp; b0=20.0_dp
      c0=b0-gr*(b0-a0); d0=a0+gr*(b0-a0)
      fc=phi_obj(c0,m); fd=phi_obj(d0,m)
      do it2=1,200
        if(abs(b0-a0)<1.0e-10_dp) exit
        if(fc>fd) then
          b0=d0; d0=c0; fd=fc; c0=b0-gr*(b0-a0); fc=phi_obj(c0,m)
        else
          a0=c0; c0=d0; fc=fd; d0=a0+gr*(b0-a0); fd=phi_obj(d0,m)
        end if
      end do
      lp=0.5_dp*(a0+b0)
    end subroutine maximize_logphi

    real(dp) function phi_obj(lp,m) result(val)
      real(dp), intent(in) :: lp,m(:,:)
      real(dp) :: ph
      integer :: ii
      ph=exp(lp); val=0.0_dp
      do ii=1,n
        val=val+w(ii)*(log_gamma(ph)-sum(log_gamma(m(ii,:)*ph))+sum((m(ii,:)*ph-1.0_dp)*log(y(ii,:))))
      end do
    end function phi_obj
  end subroutine alternative_starting_values

end module dirichletreg_start

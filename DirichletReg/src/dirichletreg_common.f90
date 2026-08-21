! SPDX-License-Identifier: GPL-2.0-or-later
module dirichletreg_common
  use dirichletreg_kinds, only : dp
  use dirichletreg_special, only : digamma, trigamma
  use dirichletreg_types, only : design_block
  implicit none
  private
  public :: common_loglik_score, common_loglik_score_hessian, common_predict, common_npar

contains

  integer function common_npar(xblocks) result(npar)
    type(design_block), intent(in) :: xblocks(:)
    integer :: j
    npar = 0
    do j = 1, size(xblocks)
      npar = npar + size(xblocks(j)%x,2)
    end do
  end function common_npar


  subroutine common_loglik_score(theta, y, xblocks, w, f, g)
    real(dp), intent(in) :: theta(:), y(:,:), w(:)
    type(design_block), intent(in) :: xblocks(:)
    real(dp), intent(out) :: f, g(:)
    real(dp), allocatable :: alpha(:,:), aplus(:), logy(:,:), da(:), dap(:)
    integer :: n, d, i, j, v, p, lo, hi
    real(dp) :: term

    n = size(y,1); d = size(y,2)
    if (size(xblocks) /= d .or. size(w) /= n .or. size(theta) /= common_npar(xblocks) .or. size(g) /= size(theta)) then
      f = -huge(1.0_dp); g = 0.0_dp; return
    end if
    allocate(alpha(n,d),aplus(n),logy(n,d),da(d),dap(n))
    logy = log(y)
    lo = 1
    do j = 1, d
      p = size(xblocks(j)%x,2); hi = lo+p-1
      if (size(xblocks(j)%x,1) /= n) then
        f = -huge(1.0_dp); g=0.0_dp; return
      end if
      alpha(:,j) = exp(matmul(xblocks(j)%x,theta(lo:hi)))
      lo = hi+1
    end do
    aplus = sum(alpha,dim=2)
    do i=1,n
      dap(i)=digamma(aplus(i))
    end do

    f=0.0_dp; g=0.0_dp
    do i=1,n
      f = f + w(i)*(log_gamma(aplus(i))-sum(log_gamma(alpha(i,:))) + sum((alpha(i,:)-1.0_dp)*logy(i,:)))
    end do
    lo=1
    do j=1,d
      p=size(xblocks(j)%x,2); hi=lo+p-1
      do i=1,n
        da(j)=digamma(alpha(i,j))
        term = w(i)*alpha(i,j)*(logy(i,j)+dap(i)-da(j))
        do v=1,p
          g(lo+v-1)=g(lo+v-1)+xblocks(j)%x(i,v)*term
        end do
      end do
      lo=hi+1
    end do
  end subroutine common_loglik_score


  subroutine common_loglik_score_hessian(theta, y, xblocks, w, f, g, h)
    real(dp), intent(in) :: theta(:), y(:,:), w(:)
    type(design_block), intent(in) :: xblocks(:)
    real(dp), intent(out) :: f, g(:), h(:,:)
    real(dp), allocatable :: alpha(:,:), aplus(:), logy(:,:), dap(:), tap(:), da(:,:), ta(:,:)
    integer, allocatable :: start(:)
    integer :: n,d,i,j,l,v,u,pj,pl,ij,il
    real(dp) :: commonterm

    call common_loglik_score(theta,y,xblocks,w,f,g)
    if (size(h,1) /= size(theta) .or. size(h,2) /= size(theta)) then
      h=0.0_dp; return
    end if
    n=size(y,1); d=size(y,2)
    allocate(alpha(n,d),aplus(n),logy(n,d),dap(n),tap(n),da(n,d),ta(n,d),start(d+1))
    logy=log(y); start(1)=1
    do j=1,d
      start(j+1)=start(j)+size(xblocks(j)%x,2)
      alpha(:,j)=exp(matmul(xblocks(j)%x,theta(start(j):start(j+1)-1)))
    end do
    aplus=sum(alpha,dim=2)
    do i=1,n
      dap(i)=digamma(aplus(i)); tap(i)=trigamma(aplus(i))
      do j=1,d
        da(i,j)=digamma(alpha(i,j)); ta(i,j)=trigamma(alpha(i,j))
      end do
    end do

    h=0.0_dp
    do j=1,d
      pj=size(xblocks(j)%x,2)
      do l=j,d
        pl=size(xblocks(l)%x,2)
        do v=1,pj
          ij=start(j)+v-1
          do u=1,pl
            il=start(l)+u-1
            h(ij,il)=0.0_dp
            if (j==l) then
              do i=1,n
                commonterm=alpha(i,j)*(logy(i,j)+dap(i)-da(i,j)+alpha(i,j)*(tap(i)-ta(i,j)))
                h(ij,il)=h(ij,il)+w(i)*xblocks(j)%x(i,v)*xblocks(j)%x(i,u)*commonterm
              end do
            else
              do i=1,n
                h(ij,il)=h(ij,il)+w(i)*xblocks(j)%x(i,v)*xblocks(l)%x(i,u)*alpha(i,j)*alpha(i,l)*tap(i)
              end do
            end if
            h(il,ij)=h(ij,il)
          end do
        end do
      end do
    end do
  end subroutine common_loglik_score_hessian


  subroutine common_predict(theta, xblocks, alpha, mu, phi, stat)
    real(dp), intent(in) :: theta(:)
    type(design_block), intent(in) :: xblocks(:)
    real(dp), intent(out) :: alpha(:,:), mu(:,:), phi(:)
    integer, intent(out), optional :: stat
    integer :: n,d,j,p,lo,hi

    if (present(stat)) stat=0
    d=size(xblocks); n=size(alpha,1)
    if (size(alpha,2)/=d .or. any(shape(mu)/=shape(alpha)) .or. size(phi)/=n .or. size(theta)/=common_npar(xblocks)) then
      if (present(stat)) stat=1
      alpha=0.0_dp; mu=0.0_dp; phi=0.0_dp; return
    end if
    lo=1
    do j=1,d
      p=size(xblocks(j)%x,2); hi=lo+p-1
      if (size(xblocks(j)%x,1)/=n) then
        if (present(stat)) stat=1
        return
      end if
      alpha(:,j)=exp(matmul(xblocks(j)%x,theta(lo:hi)))
      lo=hi+1
    end do
    phi=sum(alpha,dim=2)
    do j=1,d
      mu(:,j)=alpha(:,j)/phi
    end do
  end subroutine common_predict

end module dirichletreg_common

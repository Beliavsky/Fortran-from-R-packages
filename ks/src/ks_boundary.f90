! SPDX-License-Identifier: GPL-2.0-only
module ks_boundary
  use ks_kinds, only: dp
  use ks_normal, only: mvn_pdf
  implicit none
  private
  public :: beta_kernel2_pdf, boundary_kde_pdf, pseudo_uniform_empirical, copula_density_empirical
contains
  pure real(dp) function beta_pdf01(z,a,b) result(v)
    real(dp), intent(in) :: z,a,b
    if (z < 0.0_dp .or. z > 1.0_dp .or. a <= 0.0_dp .or. b <= 0.0_dp) then
      v=0.0_dp; return
    end if
    if (z <= 0.0_dp) then
      if (a < 1.0_dp) then; v=huge(1.0_dp)
      else if (a > 1.0_dp) then; v=0.0_dp
      else; v=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b))
      end if
      return
    end if
    if (z >= 1.0_dp) then
      if (b < 1.0_dp) then; v=huge(1.0_dp)
      else if (b > 1.0_dp) then; v=0.0_dp
      else; v=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b))
      end if
      return
    end if
    v=exp((a-1.0_dp)*log(z)+(b-1.0_dp)*log(1.0_dp-z)+log_gamma(a+b)-log_gamma(a)-log_gamma(b))
  end function

  pure real(dp) function rho_x(y,h) result(r)
    real(dp), intent(in) :: y,h
    real(dp) :: h4
    h4=h**4
    if (abs(y) <= epsilon(1.0_dp)) then
      r=1.0_dp
    else
      r=2.0_dp*h4+2.5_dp-sqrt(max(0.0_dp,4.0_dp*h4*h4+6.0_dp*h4+2.25_dp-y*y-y/(h*h)))
    end if
  end function

  pure real(dp) function beta_kernel2_pdf(center,z,h) result(v)
    real(dp), intent(in) :: center,z,h
    real(dp) :: a,b,h2,c
    if(h<=0.0_dp) then; v=0.0_dp; return; end if
    c=min(1.0_dp,max(0.0_dp,center)); h2=h*h
    if(c <= 2.0_dp*h2) then
      a=rho_x(c,h); b=(1.0_dp-c)/h2
    else if(c < 1.0_dp-2.0_dp*h2) then
      a=c/h2; b=(1.0_dp-c)/h2
    else
      a=c/h2; b=rho_x(1.0_dp-c,h)
    end if
    v=beta_pdf01(z,a,b)
  end function

  function boundary_kde_pdf(x,H,eval,xmin,xmax,w,boundary_supp) result(f)
    real(dp), intent(in) :: x(:,:),H(:,:),eval(:,:),xmin(:),xmax(:)
    real(dp), intent(in), optional :: w(:),boundary_supp
    real(dp) :: f(size(eval,1))
    real(dp), allocatable :: xs(:,:), es(:,:), Hs(:,:), wi(:), bw(:), mu(:)
    real(dp) :: bs,kern
    logical :: boundary
    integer :: n,d,m,i,j,k
    n=size(x,1); d=size(x,2); m=size(eval,1)
    if(size(H,1)/=d.or.size(H,2)/=d.or.size(eval,2)/=d.or.size(xmin)/=d.or.size(xmax)/=d) error stop 'boundary_kde_pdf: shape'
    if(any(xmax<=xmin).or.n<=0) error stop 'boundary_kde_pdf: bounds'
    allocate(xs(n,d),es(m,d),Hs(d,d),wi(n),bw(d),mu(d))
    do j=1,d
      xs(:,j)=(x(:,j)-xmin(j))/(xmax(j)-xmin(j))
      es(:,j)=(eval(:,j)-xmin(j))/(xmax(j)-xmin(j))
    end do
    Hs=H
    do j=1,d; do k=1,d; Hs(j,k)=H(j,k)/((xmax(j)-xmin(j))*(xmax(k)-xmin(k))); end do; end do
    do j=1,d
      bw(j)=sqrt(max(0.0_dp,Hs(j,j)))
    end do
    wi=1.0_dp; if(present(w)) then; if(size(w)/=n) error stop 'boundary_kde_pdf: weights'; wi=w; end if
    if(sum(wi)<=0.0_dp) error stop 'boundary_kde_pdf: weights'; wi=wi*real(n,dp)/sum(wi)
    bs=1.0_dp; if(present(boundary_supp)) bs=boundary_supp
    f=0.0_dp
    do i=1,n
      boundary=.false.
      do j=1,d
        if(abs(xs(i,j))<=bs*bw(j) .or. abs(1.0_dp-xs(i,j))<=bs*bw(j)) boundary=.true.
      end do
      if(boundary) then
        do k=1,m
          kern=1.0_dp
          do j=1,d
            kern=kern*beta_kernel2_pdf(xs(i,j),es(k,j),2.0_dp*bw(j))
          end do
          f(k)=f(k)+wi(i)*kern
        end do
      else
        mu=xs(i,:)
        do k=1,m
          f(k)=f(k)+wi(i)*mvn_pdf(es(k,:),mu,Hs)
        end do
      end if
    end do
    f=f/real(n,dp)/product(xmax-xmin)
  end function

  function pseudo_uniform_empirical(x,y) result(u)
    real(dp), intent(in) :: x(:,:),y(:,:)
    real(dp) :: u(size(y,1),size(y,2))
    integer :: i,j,k,n,d
    n=size(x,1); d=size(x,2)
    if(size(y,2)/=d.or.n<=0) error stop 'pseudo_uniform_empirical: shape'
    do j=1,d
      do i=1,size(y,1)
        k=count(x(:,j)<=y(i,j))
        u(i,j)=real(k,dp)/real(n,dp)
      end do
    end do
  end function

  function copula_density_empirical(x,H,u_eval,w) result(c)
    real(dp), intent(in) :: x(:,:),H(:,:),u_eval(:,:)
    real(dp), intent(in), optional :: w(:)
    real(dp) :: c(size(u_eval,1))
    real(dp), allocatable :: u(:,:)
    real(dp) :: lo(size(x,2)),hi(size(x,2))
    allocate(u(size(x,1),size(x,2)))
    u=pseudo_uniform_empirical(x,x)
    lo=0.0_dp; hi=1.0_dp
    c=boundary_kde_pdf(u,H,u_eval,lo,hi,w=w,boundary_supp=1.0_dp)
  end function
end module ks_boundary

! SPDX-License-Identifier: GPL-2.0-only
module mvtnorm_distributions
  use mvtnorm_kinds, only : dp, pi, log_two_pi
  use mvtnorm_linalg, only : cholesky_lower, solve_lower, jacobi_eigen, symmetrize
  use mvtnorm_random, only : seed_random, random_normals
  use mvtnorm_special, only : chi_square_quantile
  implicit none
  private
  public :: dmvnorm_one, dmvnorm, rmvnorm
  public :: dmvt_one, dmvt, rmvt_shifted, rmvt_kshirsagar

contains

  real(dp) function dmvnorm_one(x, mean, sigma, log_density, ok) result(v)
    real(dp), intent(in) :: x(:),mean(:),sigma(:,:)
    logical, intent(in), optional :: log_density
    logical, intent(out), optional :: ok
    real(dp), allocatable :: l(:,:),b(:,:),y(:,:)
    logical :: lok, llog
    character(len=256) :: message
    integer :: n,i
    real(dp) :: q,ld
    n=size(x); llog=.true.; if(present(log_density)) llog=log_density
    if(size(mean)/=n .or. size(sigma,1)/=n .or. size(sigma,2)/=n) then
      v=merge(-huge(1.0_dp),0.0_dp,llog); if(present(ok)) ok=.false.; return
    end if
    call cholesky_lower(sigma,l,lok,message)
    if(.not.lok) then
      if(maxval(abs(x-mean))<=100.0_dp*epsilon(1.0_dp)) then
        v=huge(1.0_dp)
      else
        v=merge(-huge(1.0_dp),0.0_dp,llog)
      end if
      if(present(ok)) ok=.false.; return
    end if
    allocate(b(n,1)); b(:,1)=x-mean
    call solve_lower(l,b,y,lok)
    q=sum(y(:,1)**2)
    ld=0.0_dp
    do i=1,n
      ld=ld+log(l(i,i))
    end do
    v=-0.5_dp*real(n,dp)*log_two_pi-ld-0.5_dp*q
    if(.not.llog) v=exp(v)
    if(present(ok)) ok=.true.
  end function dmvnorm_one

  function dmvnorm(x,mean,sigma,log_density) result(v)
    real(dp), intent(in) :: x(:,:),mean(:),sigma(:,:)
    logical, intent(in), optional :: log_density
    real(dp), allocatable :: v(:)
    integer :: i
    allocate(v(size(x,1)))
    do i=1,size(x,1)
      v(i)=dmvnorm_one(x(i,:),mean,sigma,log_density)
    end do
  end function dmvnorm

  function rmvnorm(n,mean,sigma,seed) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:),sigma(:,:)
    integer,intent(in),optional::seed
    real(dp),allocatable::x(:,:)
    real(dp),allocatable::l(:,:),vals(:),vecs(:,:),root(:,:),z(:)
    logical::ok
    character(len=256)::message
    integer::i,j,p
    p=size(mean); allocate(x(n,p)); x=0.0_dp
    if(present(seed)) call seed_random(seed)
    call cholesky_lower(sigma,l,ok,message,tolerance=-100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(sigma))))
    if(.not.ok) then
      call jacobi_eigen(sigma,vals,vecs,ok)
      allocate(root(p,p)); root=0.0_dp
      do j=1,p
        if(vals(j)>0.0_dp) root=root+sqrt(vals(j))*spread(vecs(:,j),2,p)*spread(vecs(:,j),1,p)
      end do
    else
      root=l
    end if
    allocate(z(p))
    do i=1,n
      call random_normals(z)
      x(i,:)=mean+matmul(root,z)
    end do
  end function rmvnorm

  real(dp) function dmvt_one(x,delta,sigma,df,log_density,ok) result(v)
    real(dp),intent(in)::x(:),delta(:),sigma(:,:),df
    logical,intent(in),optional::log_density
    logical,intent(out),optional::ok
    real(dp),allocatable::l(:,:),b(:,:),y(:,:)
    logical::lok,llog
    character(len=256)::message
    integer::p,i
    real(dp)::q,ld
    p=size(x); llog=.true.; if(present(log_density)) llog=log_density
    if(df<=0.0_dp .or. df>=huge(1.0_dp)/10.0_dp) then
      v=dmvnorm_one(x,delta,sigma,llog,lok); if(present(ok)) ok=lok; return
    end if
    call cholesky_lower(sigma,l,lok,message)
    if(.not.lok) then
      v=merge(-huge(1.0_dp),0.0_dp,llog); if(present(ok)) ok=.false.; return
    end if
    allocate(b(p,1)); b(:,1)=x-delta
    call solve_lower(l,b,y,lok); q=sum(y(:,1)**2)
    ld=0.0_dp; do i=1,p; ld=ld+log(l(i,i)); end do
    v=log_gamma(0.5_dp*(df+real(p,dp)))-log_gamma(0.5_dp*df)-ld &
      -0.5_dp*real(p,dp)*log(pi*df)-0.5_dp*(df+real(p,dp))*log(1.0_dp+q/df)
    if(.not.llog) v=exp(v)
    if(present(ok)) ok=.true.
  end function dmvt_one

  function dmvt(x,delta,sigma,df,log_density) result(v)
    real(dp),intent(in)::x(:,:),delta(:),sigma(:,:),df
    logical,intent(in),optional::log_density
    real(dp),allocatable::v(:)
    integer::i
    allocate(v(size(x,1)))
    do i=1,size(x,1)
      v(i)=dmvt_one(x(i,:),delta,sigma,df,log_density)
    end do
  end function dmvt

  function rmvt_shifted(n,sigma,df,delta,seed) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::sigma(:,:),df,delta(:)
    integer,intent(in),optional::seed
    real(dp),allocatable::x(:,:),z(:,:)
    real(dp)::u,s
    integer::i
    z=rmvnorm(n,spread(0.0_dp,1,size(delta)),sigma,seed)
    allocate(x(n,size(delta)))
    do i=1,n
      call random_number(u); u=min(1.0_dp-epsilon(1.0_dp),max(epsilon(1.0_dp),u))
      s=sqrt(chi_square_quantile(u,df)/df)
      x(i,:)=delta+z(i,:)/s
    end do
  end function rmvt_shifted

  function rmvt_kshirsagar(n,sigma,df,delta,seed) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::sigma(:,:),df,delta(:)
    integer,intent(in),optional::seed
    real(dp),allocatable::x(:,:),z(:,:)
    real(dp)::u,s
    integer::i
    z=rmvnorm(n,delta,sigma,seed)
    allocate(x(n,size(delta)))
    do i=1,n
      call random_number(u); u=min(1.0_dp-epsilon(1.0_dp),max(epsilon(1.0_dp),u))
      s=sqrt(chi_square_quantile(u,df)/df)
      x(i,:)=z(i,:)/s
    end do
  end function rmvt_kshirsagar
end module mvtnorm_distributions

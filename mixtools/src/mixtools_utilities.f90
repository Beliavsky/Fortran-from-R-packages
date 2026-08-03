! SPDX-License-Identifier: GPL-2.0-or-later
module mixtools_utilities
  use mixtools_kinds, only : dp, pi
  use mixtools_status, only : MIXTOOLS_SUCCESS, MIXTOOLS_DIMENSION_ERROR
  use mixtools_status, only : MIXTOOLS_INVALID_ARGUMENT
  use mixtools_linalg, only : inverse_spd
  implicit none
  private
  public :: wquantile, wiqr, wkde, weighted_bandwidth, matsqrt
  public :: mahalanobis_depth, ellipse_points, aug_x, ldmult, ldc
  public :: kernel_value, lambda_weights, lambda_pert, permutations
  public :: sort_real, median_real, standardize_columns, effective_sample_size
contains
  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i,j
    real(dp) :: key
    do i=2,size(x)
      key=x(i);j=i-1
      do while(j>=1)
        if(x(j)<=key)exit
        x(j+1)=x(j);j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real

  function median_real(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: v
    real(dp), allocatable :: y(:)
    integer :: n
    n=size(x);allocate(y(n));y=x;call sort_real(y)
    if(mod(n,2)==1)then
      v=y((n+1)/2)
    else
      v=0.5_dp*(y(n/2)+y(n/2+1))
    end if
  end function median_real

  function wquantile(x, weights, probability, status) result(q)
    real(dp), intent(in) :: x(:),weights(:),probability
    integer,intent(out),optional::status
    real(dp)::q
    real(dp),allocatable::xs(:),ws(:)
    integer::i,j,n
    real(dp)::tx,tw,target,cum,s
    n=size(x)
    if(size(weights)/=n.or.n==0.or.any(weights<0.0_dp).or.sum(weights)<=0.0_dp)then
      q=0.0_dp;if(present(status))status=MIXTOOLS_INVALID_ARGUMENT;return
    end if
    allocate(xs(n),ws(n));xs=x;ws=weights
    do i=2,n
      tx=xs(i);tw=ws(i);j=i-1
      do while(j>=1)
        if(xs(j)<=tx)exit
        xs(j+1)=xs(j);ws(j+1)=ws(j);j=j-1
      end do
      xs(j+1)=tx;ws(j+1)=tw
    end do
    s=sum(ws);target=min(1.0_dp,max(0.0_dp,probability))*s;cum=0.0_dp;q=xs(n)
    do i=1,n
      cum=cum+ws(i)
      if(cum>=target)then;q=xs(i);exit;end if
    end do
    if(present(status))status=MIXTOOLS_SUCCESS
  end function wquantile

  function wiqr(x,weights,status) result(v)
    real(dp),intent(in)::x(:),weights(:)
    integer,intent(out),optional::status
    real(dp)::v
    integer::s1,s2
    v=wquantile(x,weights,0.75_dp,s1)-wquantile(x,weights,0.25_dp,s2)
    if(present(status))status=max(s1,s2)
  end function wiqr

  function weighted_bandwidth(x,weights) result(h)
    real(dp),intent(in)::x(:),weights(:)
    real(dp)::h,mu,var,neff,iqr,s
    s=sum(weights)
    if(s<=0.0_dp.or.size(x)<2)then;h=1.0_dp;return;end if
    mu=sum(weights*x)/s
    var=sum(weights*(x-mu)**2)/s
    iqr=wiqr(x,weights)
    neff=effective_sample_size(weights)
    h=0.9_dp*min(sqrt(max(var,0.0_dp)),iqr/1.34_dp)*max(neff,2.0_dp)**(-0.2_dp)
    if(h<=sqrt(epsilon(1.0_dp)))h=max(sqrt(max(var,0.0_dp))*0.1_dp,1.0e-3_dp)
  end function weighted_bandwidth

  pure function effective_sample_size(weights) result(v)
    real(dp),intent(in)::weights(:)
    real(dp)::v,s
    s=sum(weights)
    if(sum(weights*weights)>0.0_dp)then
      v=s*s/sum(weights*weights)
    else
      v=0.0_dp
    end if
  end function effective_sample_size

  elemental function kernel_value(u,kernel) result(v)
    real(dp),intent(in)::u
    integer,intent(in),optional::kernel
    real(dp)::v,a
    integer::kk
    kk=1;if(present(kernel))kk=kernel;a=abs(u)
    select case(kk)
    case(1);v=exp(-0.5_dp*u*u)/sqrt(2.0_dp*pi)
    case(2);v=merge(0.75_dp*(1.0_dp-u*u),0.0_dp,a<=1.0_dp)
    case(3);v=merge(1.0_dp-a,0.0_dp,a<=1.0_dp)
    case(4);v=merge((pi/4.0_dp)*cos(pi*u/2.0_dp),0.0_dp,a<=1.0_dp)
    case(5);v=merge((3.0_dp/(4.0_dp*sqrt(5.0_dp)))*(1.0_dp-u*u/5.0_dp),0.0_dp,a<=sqrt(5.0_dp))
    case default;v=exp(-0.5_dp*u*u)/sqrt(2.0_dp*pi)
    end select
  end function kernel_value

  subroutine wkde(x,weights,grid,bandwidth,density,symmetric,kernel,status)
    real(dp),intent(in)::x(:),weights(:),grid(:),bandwidth
    real(dp),intent(out)::density(size(grid))
    logical,intent(in),optional::symmetric
    integer,intent(in),optional::kernel
    integer,intent(out),optional::status
    integer::i,j,kk
    logical::sym
    real(dp)::s,h
    sym=.false.;if(present(symmetric))sym=symmetric
    kk=1;if(present(kernel))kk=kernel
    if(size(weights)/=size(x).or.bandwidth<=0.0_dp.or.sum(weights)<=0.0_dp)then
      density=0.0_dp;if(present(status))status=MIXTOOLS_INVALID_ARGUMENT;return
    end if
    s=sum(weights);h=bandwidth
    do i=1,size(grid)
      density(i)=0.0_dp
      do j=1,size(x)
        density(i)=density(i)+weights(j)*kernel_value((grid(i)-x(j))/h,kk)
        if(sym)density(i)=density(i)+weights(j)*kernel_value((grid(i)+x(j))/h,kk)
      end do
      density(i)=density(i)/(s*h*merge(2.0_dp,1.0_dp,sym))
    end do
    if(present(status))status=MIXTOOLS_SUCCESS
  end subroutine wkde

  subroutine lambda_weights(z,x,centers,bandwidth,lambda_out,kernel,status)
    real(dp),intent(in)::z(:,:),x(:),centers(:),bandwidth
    real(dp),intent(out)::lambda_out(size(centers),size(z,2))
    integer,intent(in),optional::kernel
    integer,intent(out),optional::status
    integer::i,j,c,kk
    real(dp)::den,w
    kk=1;if(present(kernel))kk=kernel
    if(size(z,1)/=size(x).or.bandwidth<=0.0_dp)then
      lambda_out=0.0_dp;if(present(status))status=MIXTOOLS_DIMENSION_ERROR;return
    end if
    do c=1,size(z,2)
      do j=1,size(centers)
        den=0.0_dp;lambda_out(j,c)=0.0_dp
        do i=1,size(x)
          w=kernel_value((x(i)-centers(j))/bandwidth,kk)
          den=den+w;lambda_out(j,c)=lambda_out(j,c)+w*z(i,c)
        end do
        if(den>0.0_dp)lambda_out(j,c)=lambda_out(j,c)/den
      end do
    end do
    if(present(status))status=MIXTOOLS_SUCCESS
  end subroutine lambda_weights

  subroutine lambda_pert(lambda_in,perturbation,lambda_out)
    real(dp),intent(in)::lambda_in(:),perturbation(:)
    real(dp),intent(out)::lambda_out(size(lambda_in))
    lambda_out=max(0.0_dp,lambda_in+perturbation)
    if(sum(lambda_out)>0.0_dp)lambda_out=lambda_out/sum(lambda_out)
  end subroutine lambda_pert

  subroutine jacobi_eigen(a,eval,evec,status)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::eval(:),evec(:,:)
    integer,intent(out)::status
    real(dp),allocatable::b(:,:)
    real(dp)::app,aqq,apq,tau,t,c,s,bip,biq,maxoff
    integer::n,i,j,p,q,iter
    n=size(a,1)
    if(size(a,2)/=n)then;allocate(eval(0),evec(0,0));status=MIXTOOLS_DIMENSION_ERROR;return;end if
    allocate(b(n,n),eval(n),evec(n,n));b=0.5_dp*(a+transpose(a));evec=0.0_dp
    do i=1,n;evec(i,i)=1.0_dp;end do
    do iter=1,100*n*n
      maxoff=0.0_dp;p=1;q=min(2,n)
      do i=1,n-1;do j=i+1,n
        if(abs(b(i,j))>maxoff)then;maxoff=abs(b(i,j));p=i;q=j;end if
      end do;end do
      if(maxoff<1.0e-12_dp)exit
      app=b(p,p);aqq=b(q,q);apq=b(p,q);tau=(aqq-app)/(2.0_dp*apq)
      t=sign(1.0_dp,tau)/(abs(tau)+sqrt(1.0_dp+tau*tau));c=1.0_dp/sqrt(1.0_dp+t*t);s=t*c
      do i=1,n
        if(i/=p.and.i/=q)then
          bip=b(i,p);biq=b(i,q);b(i,p)=c*bip-s*biq;b(p,i)=b(i,p)
          b(i,q)=s*bip+c*biq;b(q,i)=b(i,q)
        end if
      end do
      b(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
      b(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq;b(p,q)=0.0_dp;b(q,p)=0.0_dp
      do i=1,n
        bip=evec(i,p);biq=evec(i,q);evec(i,p)=c*bip-s*biq;evec(i,q)=s*bip+c*biq
      end do
    end do
    do i=1,n;eval(i)=b(i,i);end do
    status=MIXTOOLS_SUCCESS
  end subroutine jacobi_eigen

  subroutine matsqrt(a,sqrta,status)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::sqrta(:,:)
    integer,intent(out)::status
    real(dp),allocatable::eval(:),evec(:,:),d(:,:)
    integer::i,n
    call jacobi_eigen(a,eval,evec,status);if(status/=MIXTOOLS_SUCCESS)then;allocate(sqrta(0,0));return;end if
    n=size(eval);allocate(d(n,n));d=0.0_dp
    do i=1,n;d(i,i)=sqrt(max(0.0_dp,eval(i)));end do
    sqrta=matmul(evec,matmul(d,transpose(evec)))
  end subroutine matsqrt

  subroutine mahalanobis_depth(points,x,covariance,depth,status)
    real(dp),intent(in)::points(:,:),x(:,:),covariance(:,:)
    real(dp),intent(out)::depth(size(points,1))
    integer,intent(out)::status
    real(dp),allocatable::inv(:,:),mu(:),d(:)
    integer::i,p
    p=size(x,2)
    if(size(points,2)/=p)then;depth=0.0_dp;status=MIXTOOLS_DIMENSION_ERROR;return;end if
    call inverse_spd(covariance,inv,status);if(status/=MIXTOOLS_SUCCESS)then;depth=0.0_dp;return;end if
    allocate(mu(p),d(p));mu=sum(x,dim=1)/real(size(x,1),dp)
    do i=1,size(points,1);d=points(i,:)-mu;depth(i)=1.0_dp/(1.0_dp+dot_product(d,matmul(inv,d)));end do
  end subroutine mahalanobis_depth

  subroutine ellipse_points(mu,sigma,alpha,npoints,points,status)
    real(dp),intent(in)::mu(2),sigma(2,2),alpha
    integer,intent(in)::npoints
    real(dp),allocatable,intent(out)::points(:,:)
    integer,intent(out)::status
    real(dp),allocatable::root(:,:)
    real(dp)::radius,theta
    integer::i
    call matsqrt(sigma,root,status);if(status/=MIXTOOLS_SUCCESS)then;allocate(points(0,0));return;end if
    radius=sqrt(-2.0_dp*log(max(alpha,tiny(1.0_dp))))
    allocate(points(npoints,2))
    do i=1,npoints
      theta=2.0_dp*pi*real(i-1,dp)/real(max(1,npoints-1),dp)
      points(i,:)=mu+radius*matmul(root,[cos(theta),sin(theta)])
    end do
  end subroutine ellipse_points

  subroutine aug_x(x,cp_locs,cp,delta,xaug,status)
    real(dp),intent(in)::x(:,:),cp(:),delta(:)
    integer,intent(in)::cp_locs(:)
    real(dp),allocatable,intent(out)::xaug(:,:)
    integer,intent(out)::status
    integer::n,p,j
    n=size(x,1);p=size(x,2)
    if(size(cp)/=size(cp_locs).or.size(delta)/=size(cp_locs).or.any(cp_locs<1).or.any(cp_locs>p))then
      allocate(xaug(0,0));status=MIXTOOLS_DIMENSION_ERROR;return
    end if
    allocate(xaug(n,p+2*size(cp)));xaug(:,1:p)=x
    do j=1,size(cp)
      xaug(:,p+2*j-1)=max(0.0_dp,x(:,cp_locs(j))-cp(j))
      xaug(:,p+2*j)=-delta(j)*merge(1.0_dp,0.0_dp,x(:,cp_locs(j))>cp(j))
    end do
    status=MIXTOOLS_SUCCESS
  end subroutine aug_x

  function ldmult(y,theta) result(v)
    real(dp),intent(in)::y(:),theta(:)
    real(dp)::v
    v=log_gamma(sum(y)+1.0_dp)-sum(log_gamma(y+1.0_dp))+sum(y*log(max(theta,tiny(1.0_dp))))
  end function ldmult

  subroutine ldc(data,classes,scores,status)
    real(dp),intent(in)::data(:,:)
    integer,intent(in)::classes(:)
    real(dp),intent(out)::scores(size(data,1))
    integer,intent(out)::status
    real(dp),allocatable::m1(:),m2(:),cov(:,:),inv(:,:)
    integer::n1,n2,i,p
    p=size(data,2)
    if(size(classes)/=size(data,1).or.maxval(classes)>2.or.minval(classes)<1)then
      scores=0.0_dp;status=MIXTOOLS_DIMENSION_ERROR;return
    end if
    n1=count(classes==1);n2=count(classes==2);allocate(m1(p),m2(p),cov(p,p));m1=0.0_dp;m2=0.0_dp
    do i=1,size(classes);if(classes(i)==1)m1=m1+data(i,:);if(classes(i)==2)m2=m2+data(i,:);end do
    m1=m1/real(n1,dp);m2=m2/real(n2,dp);cov=0.0_dp
    do i=1,size(classes)
      if(classes(i)==1)then;cov=cov+outer(data(i,:)-m1);else;cov=cov+outer(data(i,:)-m2);end if
    end do
    cov=cov/real(max(1,size(classes)-2),dp);call inverse_spd(cov,inv,status)
    if(status/=MIXTOOLS_SUCCESS)then;scores=0.0_dp;return;end if
    scores=matmul(data,matmul(inv,m1-m2))
  contains
    pure function outer(v) result(a)
      real(dp),intent(in)::v(:);real(dp)::a(size(v),size(v));integer::ii,jj
      do ii=1,size(v);do jj=1,size(v);a(ii,jj)=v(ii)*v(jj);end do;end do
    end function outer
  end subroutine ldc

  recursive subroutine permutations(n,r,out,status)
    integer,intent(in)::n,r
    integer,allocatable,intent(out)::out(:,:)
    integer,intent(out)::status
    integer::count,row
    integer,allocatable::current(:)
    logical,allocatable::used(:)
    if(r<0.or.r>n)then;allocate(out(0,0));status=MIXTOOLS_INVALID_ARGUMENT;return;end if
    count=1
    if(r>0)then
      count=1
      do row=0,r-1;count=count*(n-row);end do
    end if
    allocate(out(count,r),current(r),used(n));used=.false.;row=0
    call fill(1)
    status=MIXTOOLS_SUCCESS
  contains
    recursive subroutine fill(pos)
      integer,intent(in)::pos
      integer::j
      if(pos>r)then;row=row+1;out(row,:)=current;return;end if
      do j=1,n
        if(.not.used(j))then;used(j)=.true.;current(pos)=j;call fill(pos+1);used(j)=.false.;end if
      end do
    end subroutine fill
  end subroutine permutations

  subroutine standardize_columns(x,z,center,scale)
    real(dp),intent(in)::x(:,:)
    real(dp),intent(out)::z(size(x,1),size(x,2)),center(size(x,2)),scale(size(x,2))
    integer::j
    do j=1,size(x,2)
      center(j)=sum(x(:,j))/real(size(x,1),dp)
      scale(j)=sqrt(sum((x(:,j)-center(j))**2)/real(max(1,size(x,1)-1),dp))
      if(scale(j)<=tiny(1.0_dp))scale(j)=1.0_dp
      z(:,j)=(x(:,j)-center(j))/scale(j)
    end do
  end subroutine standardize_columns
end module mixtools_utilities

module mlr_preprocess
  use mlr_kinds, only : dp
  use mlr_types, only : scaler_model
  use mlr_utils, only : mean_dp, median_dp
  use mlr_rng, only : rng_state, rng_integer, rng_normal
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  implicit none
  private
  public :: fit_standardizer, apply_standardizer, standardize
  public :: impute_mean, impute_median, impute_constant, cap_large_values, one_hot_encode
  public :: downsample_classes, oversample_classes
contains
  subroutine fit_standardizer(x, model)
    real(dp), intent(in) :: x(:,:)
    type(scaler_model), intent(out) :: model
    integer :: j
    allocate(model%center(size(x,2)), model%scale(size(x,2)))
    do j=1,size(x,2)
      model%center(j)=mean_dp(x(:,j))
      if(size(x,1)>1)then
        model%scale(j)=sqrt(sum((x(:,j)-model%center(j))**2)/real(size(x,1)-1,dp))
      else
        model%scale(j)=1.0_dp
      end if
      if(model%scale(j)<=sqrt(tiny(1.0_dp)))model%scale(j)=1.0_dp
    end do
  end subroutine

  subroutine apply_standardizer(x, model, z)
    real(dp), intent(in) :: x(:,:)
    type(scaler_model), intent(in) :: model
    real(dp), allocatable, intent(out) :: z(:,:)
    integer :: j
    allocate(z(size(x,1),size(x,2)))
    do j=1,size(x,2); z(:,j)=(x(:,j)-model%center(j))/model%scale(j); end do
  end subroutine

  subroutine standardize(x,z,model)
    real(dp),intent(in)::x(:,:); real(dp),allocatable,intent(out)::z(:,:)
    type(scaler_model),intent(out)::model
    call fit_standardizer(x,model); call apply_standardizer(x,model,z)
  end subroutine

  subroutine impute_mean(x)
    real(dp),intent(inout)::x(:,:)
    integer::i,j,n; real(dp)::s,m
    do j=1,size(x,2); s=0.0_dp;n=0
      do i=1,size(x,1); if(.not.ieee_is_nan(x(i,j)))then;s=s+x(i,j);n=n+1;end if;end do
      if(n==0)cycle; m=s/real(n,dp)
      do i=1,size(x,1);if(ieee_is_nan(x(i,j)))x(i,j)=m;end do
    end do
  end subroutine

  subroutine impute_median(x)
    real(dp),intent(inout)::x(:,:)
    real(dp),allocatable::tmp(:); integer::i,j,n
    do j=1,size(x,2)
      n=count(.not.ieee_is_nan(x(:,j))); if(n==0)cycle
      allocate(tmp(n)); n=0
      do i=1,size(x,1);if(.not.ieee_is_nan(x(i,j)))then;n=n+1;tmp(n)=x(i,j);end if;end do
      do i=1,size(x,1);if(ieee_is_nan(x(i,j)))x(i,j)=median_dp(tmp);end do
      deallocate(tmp)
    end do
  end subroutine

  subroutine impute_constant(x,value)
    real(dp),intent(inout)::x(:,:); real(dp),intent(in)::value
    integer::i,j
    do j=1,size(x,2);do i=1,size(x,1);if(ieee_is_nan(x(i,j)))x(i,j)=value;end do;end do
  end subroutine

  subroutine cap_large_values(x,lower,upper)
    real(dp),intent(inout)::x(:,:); real(dp),intent(in)::lower,upper
    x=max(lower,min(upper,x))
  end subroutine

  subroutine one_hot_encode(codes, nlevels, out)
    integer,intent(in)::codes(:,:),nlevels(:)
    real(dp),allocatable,intent(out)::out(:,:)
    integer::i,j,k,pos,p
    if(size(codes,2)/=size(nlevels))error stop "one_hot_encode: level mismatch"
    p=sum(nlevels);allocate(out(size(codes,1),p));out=0.0_dp;pos=0
    do j=1,size(codes,2)
      do i=1,size(codes,1)
        k=codes(i,j);if(k>=1.and.k<=nlevels(j))out(i,pos+k)=1.0_dp
      end do
      pos=pos+nlevels(j)
    end do
  end subroutine

  subroutine downsample_classes(x,y,target_per_class,rng,xout,yout)
    real(dp),intent(in)::x(:,:); integer,intent(in)::y(:),target_per_class
    type(rng_state),intent(inout)::rng
    real(dp),allocatable,intent(out)::xout(:,:);integer,allocatable,intent(out)::yout(:)
    integer::c,nclass,nout,i,j,nc,sel; integer,allocatable::pool(:),take(:)
    nclass=maxval(y);nout=0
    do c=1,nclass;nout=nout+min(target_per_class,count(y==c));end do
    allocate(xout(nout,size(x,2)),yout(nout));j=0
    do c=1,nclass
      nc=count(y==c);allocate(pool(nc));i=0
      do sel=1,size(y);if(y(sel)==c)then;i=i+1;pool(i)=sel;end if;end do
      allocate(take(nc));take=pool
      do i=1,min(target_per_class,nc)
        sel=rng_integer(rng,i,nc); call swap_int(take(i),take(sel))
        j=j+1;xout(j,:)=x(take(i),:);yout(j)=c
      end do
      deallocate(pool,take)
    end do
  contains
    subroutine swap_int(a,b);integer,intent(inout)::a,b;integer::t;t=a;a=b;b=t;end subroutine
  end subroutine

  subroutine oversample_classes(x,y,target_per_class,rng,xout,yout)
    real(dp),intent(in)::x(:,:); integer,intent(in)::y(:),target_per_class
    type(rng_state),intent(inout)::rng
    real(dp),allocatable,intent(out)::xout(:,:);integer,allocatable,intent(out)::yout(:)
    integer::nclass,c,nout,i,j,nc,sel;integer,allocatable::pool(:)
    nclass=maxval(y);nout=0
    do c=1,nclass;nout=nout+max(target_per_class,count(y==c));end do
    allocate(xout(nout,size(x,2)),yout(nout));j=0
    do c=1,nclass
      nc=count(y==c);allocate(pool(nc));i=0
      do sel=1,size(y);if(y(sel)==c)then;i=i+1;pool(i)=sel;end if;end do
      do i=1,max(target_per_class,nc)
        if(i<=nc)then;sel=pool(i);else;sel=pool(rng_integer(rng,1,nc));end if
        j=j+1;xout(j,:)=x(sel,:);yout(j)=c
      end do
      deallocate(pool)
    end do
  end subroutine
end module mlr_preprocess

! SPDX-License-Identifier: GPL-2.0-only
module kernlab_kernels
  use kernlab_kinds, only: dp, KL_SUCCESS, KL_INVALID_ARGUMENT
  use kernlab_types
  implicit none
  private
  public :: rbfdot, laplacedot, besseldot, polydot, tanhdot, vanilladot
  public :: anovadot, splinedot, stringdot, fourierdot
  public :: kernel_value, kernel_matrix, kernel_mult, kernel_pol, kernel_fast
  public :: string_kernel_value, string_kernel_matrix

contains

  pure function rbfdot(sigma) result(k)
    real(dp), intent(in), optional :: sigma
    type(kernel_spec) :: k
    k%kind=KERNEL_RBF; if(present(sigma)) k%sigma=sigma
  end function rbfdot
  pure function laplacedot(sigma) result(k)
    real(dp), intent(in), optional :: sigma
    type(kernel_spec) :: k
    k%kind=KERNEL_LAPLACE; if(present(sigma)) k%sigma=sigma
  end function laplacedot
  pure function besseldot(sigma,order,degree) result(k)
    real(dp), intent(in), optional :: sigma
    integer, intent(in), optional :: order,degree
    type(kernel_spec) :: k
    k%kind=KERNEL_BESSEL; if(present(sigma)) k%sigma=sigma
    if(present(order)) k%order=order; if(present(degree)) k%degree=degree
  end function besseldot
  pure function polydot(degree,scale,offset) result(k)
    integer,intent(in),optional :: degree
    real(dp),intent(in),optional :: scale,offset
    type(kernel_spec)::k
    k%kind=KERNEL_POLY; if(present(degree)) k%degree=degree
    if(present(scale)) k%scale=scale; if(present(offset)) k%offset=offset
  end function polydot
  pure function tanhdot(scale,offset) result(k)
    real(dp),intent(in),optional::scale,offset
    type(kernel_spec)::k
    k%kind=KERNEL_TANH; if(present(scale)) k%scale=scale; if(present(offset)) k%offset=offset
  end function tanhdot
  pure function vanilladot() result(k)
    type(kernel_spec)::k; k%kind=KERNEL_LINEAR
  end function vanilladot
  pure function anovadot(sigma,degree) result(k)
    real(dp),intent(in),optional::sigma
    integer,intent(in),optional::degree
    type(kernel_spec)::k
    k%kind=KERNEL_ANOVA; if(present(sigma)) k%sigma=sigma; if(present(degree)) k%degree=degree
  end function anovadot
  pure function splinedot() result(k)
    type(kernel_spec)::k; k%kind=KERNEL_SPLINE
  end function splinedot
  pure function fourierdot(sigma) result(k)
    real(dp),intent(in),optional::sigma
    type(kernel_spec)::k; k%kind=KERNEL_FOURIER; if(present(sigma)) k%sigma=sigma
  end function fourierdot
  pure function stringdot(length,lambda,normalized) result(k)
    integer,intent(in),optional::length
    real(dp),intent(in),optional::lambda
    logical,intent(in),optional::normalized
    type(kernel_spec)::k
    k%kind=KERNEL_STRING; if(present(length)) k%string_length=length
    if(present(lambda)) k%lambda=lambda; if(present(normalized)) k%normalized=normalized
  end function stringdot

  pure real(dp) function kernel_value(kernel,x,y)
    type(kernel_spec),intent(in)::kernel
    real(dp),intent(in)::x(:),y(:)
    real(dp)::d2,r,minv,v,lim,bkt
    integer::j
    if(size(x)/=size(y)) then; kernel_value=0.0_dp; return; end if
    select case(kernel%kind)
    case(KERNEL_RBF)
      d2=sum((x-y)**2); kernel_value=exp(-kernel%sigma*d2)
    case(KERNEL_LAPLACE)
      d2=sum((x-y)**2); kernel_value=exp(-kernel%sigma*sqrt(max(0.0_dp,d2)))
    case(KERNEL_BESSEL)
      bkt=kernel%sigma*sqrt(max(0.0_dp,sum((x-y)**2)))
      lim=1.0_dp/(gamma(real(kernel%order+1,dp))*2.0_dp**kernel%order)
      if(bkt<1.0e-4_dp) then; v=lim
      else; v=bessel_jn(max(0,kernel%order),bkt)*bkt**(-kernel%order); end if
      kernel_value=(v/lim)**kernel%degree
    case(KERNEL_POLY)
      kernel_value=(kernel%scale*dot_product(x,y)+kernel%offset)**kernel%degree
    case(KERNEL_TANH)
      kernel_value=tanh(kernel%scale*dot_product(x,y)+kernel%offset)
    case(KERNEL_LINEAR)
      kernel_value=dot_product(x,y)
    case(KERNEL_ANOVA)
      kernel_value=sum(exp(-kernel%sigma*(x-y)**2))**kernel%degree
    case(KERNEL_SPLINE)
      r=1.0_dp
      do j=1,size(x)
        minv=min(x(j),y(j))
        r=r*(1.0_dp+x(j)*y(j)*(1.0_dp+minv)-0.5_dp*(x(j)+y(j))*minv**2+minv**3/3.0_dp)
      end do
      kernel_value=r
    case(KERNEL_FOURIER)
      r=1.0_dp
      do j=1,size(x)
        r=r*((1.0_dp-kernel%sigma**2)/2.0_dp*(1.0_dp-2.0_dp*kernel%sigma*cos(x(j)-y(j))+kernel%sigma**2))
      end do
      kernel_value=r
    case default
      kernel_value=0.0_dp
    end select
  end function kernel_value

  subroutine kernel_matrix(kernel,x,k,status,y)
    type(kernel_spec),intent(in)::kernel
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::k(:,:)
    integer,intent(out)::status
    real(dp),intent(in),optional::y(:,:)
    integer::i,j,n,m
    status=KL_INVALID_ARGUMENT; allocate(k(0,0))
    if(present(y)) then
      if(size(x,2)/=size(y,2)) return
      n=size(x,1);m=size(y,1); deallocate(k);allocate(k(n,m))
      do i=1,n;do j=1,m;k(i,j)=kernel_value(kernel,x(i,:),y(j,:));end do;end do
    else
      n=size(x,1); deallocate(k);allocate(k(n,n))
      do i=1,n
        do j=1,i
          k(i,j)=kernel_value(kernel,x(i,:),x(j,:)); k(j,i)=k(i,j)
        end do
      end do
    end if
    status=KL_SUCCESS
  end subroutine kernel_matrix

  subroutine kernel_mult(kernel,x,z,result,status,y,blocksize)
    type(kernel_spec),intent(in)::kernel
    real(dp),intent(in)::x(:,:),z(:,:)
    real(dp),allocatable,intent(out)::result(:,:)
    integer,intent(out)::status
    real(dp),intent(in),optional::y(:,:)
    integer,intent(in),optional::blocksize
    real(dp),allocatable::k(:,:)
    integer::m, bs
    bs=128; if(present(blocksize)) bs=max(1,blocksize)
    if(present(y)) then; call kernel_matrix(kernel,x,k,status,y)
    else; call kernel_matrix(kernel,x,k,status); end if
    if(status/=KL_SUCCESS) then; allocate(result(0,0));return;end if
    m=size(k,2); if(size(z,1)/=m) then;status=KL_INVALID_ARGUMENT;allocate(result(0,0));return;end if
    allocate(result(size(k,1),size(z,2))); result=matmul(k,z)
    if(bs<0) result=result
  end subroutine kernel_mult

  subroutine kernel_pol(kernel,x,z,value,status,y,k0)
    type(kernel_spec),intent(in)::kernel
    real(dp),intent(in)::x(:,:),z(:)
    real(dp),intent(out)::value
    integer,intent(out)::status
    real(dp),intent(in),optional::y(:,:),k0
    real(dp),allocatable::v(:,:),zz(:,:)
    allocate(zz(size(z),1));zz(:,1)=z
    if(present(y)) then; call kernel_mult(kernel,x,zz,v,status,y)
    else; call kernel_mult(kernel,x,zz,v,status);end if
    if(status/=KL_SUCCESS) then;value=0.0_dp;return;end if
    value=dot_product(z,v(:,1)); if(present(k0)) value=value+k0
  end subroutine kernel_pol

  subroutine kernel_fast(kernel,x,y,values,status,dota)
    type(kernel_spec),intent(in)::kernel
    real(dp),intent(in)::x(:,:),y(:)
    real(dp),allocatable,intent(out)::values(:)
    integer,intent(out)::status
    real(dp),intent(in),optional::dota(:)
    integer::i
    if(present(dota)) then
      if(size(dota)<0) continue
    end if
    status=KL_INVALID_ARGUMENT;allocate(values(0));if(size(x,2)/=size(y))return
    deallocate(values);allocate(values(size(x,1)))
    do i=1,size(x,1);values(i)=kernel_value(kernel,x(i,:),y);end do
    status=KL_SUCCESS
  end subroutine kernel_fast

  pure real(dp) function raw_spectrum(a,b,l)
    character(len=*),intent(in)::a,b
    integer,intent(in)::l
    integer::i,j,na,nb
    raw_spectrum=0.0_dp;na=len_trim(a);nb=len_trim(b)
    if(l<1.or.na<l.or.nb<l)return
    do i=1,na-l+1;do j=1,nb-l+1
      if(a(i:i+l-1)==b(j:j+l-1))raw_spectrum=raw_spectrum+1.0_dp
    end do;end do
  end function raw_spectrum

  pure real(dp) function string_kernel_value(kernel,a,b)
    type(kernel_spec),intent(in)::kernel
    character(len=*),intent(in)::a,b
    real(dp)::kab,kaa,kbb
    kab=raw_spectrum(a,b,kernel%string_length)
    if(kernel%normalized)then
      kaa=raw_spectrum(a,a,kernel%string_length);kbb=raw_spectrum(b,b,kernel%string_length)
      if(kaa>0.0_dp.and.kbb>0.0_dp)kab=kab/sqrt(kaa*kbb)
    end if
    string_kernel_value=kab
  end function string_kernel_value

  subroutine string_kernel_matrix(kernel,x,k,status,y)
    type(kernel_spec),intent(in)::kernel
    character(len=*),intent(in)::x(:)
    real(dp),allocatable,intent(out)::k(:,:)
    integer,intent(out)::status
    character(len=*),intent(in),optional::y(:)
    integer::i,j,n,m
    status=KL_INVALID_ARGUMENT;allocate(k(0,0));if(kernel%kind/=KERNEL_STRING)return
    n=size(x)
    if(present(y))then
      m=size(y);deallocate(k);allocate(k(n,m));do i=1,n;do j=1,m;k(i,j)=string_kernel_value(kernel,x(i),y(j));end do;end do
    else
      deallocate(k);allocate(k(n,n));do i=1,n;do j=1,i;k(i,j)=string_kernel_value(kernel,x(i),x(j));k(j,i)=k(i,j);end do;end do
    end if
    status=KL_SUCCESS
  end subroutine string_kernel_matrix

end module kernlab_kernels

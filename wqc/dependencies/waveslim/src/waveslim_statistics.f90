! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
module waveslim_statistics
  use waveslim_kinds, only : dp, pi
  use waveslim_types, only : wavelet_transform, real_vector, test_result
  use waveslim_math, only : normal_quantile, chi_square_cdf, autocovariance, &
    autocorrelation, cross_correlation, mean_value, variance_value, fft_complex
  use waveslim_linalg, only : symmetric_eigen
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: wave_variance, wave_covariance, wave_correlation
  public :: spin_covariance, spin_correlation, my_acf, my_ccf
  public :: wave_variance_2d, sine_taper, dpss_taper, periodogram
  public :: rotcumvar, wavelet_filter, squared_gain

  type, public :: interval_vector
    real(dp), allocatable :: estimate(:), lower(:), upper(:)
  end type interval_vector
contains
  function wave_variance(wt,ci_type,p) result(ans)
    type(wavelet_transform),intent(in)::wt
    character(len=*),intent(in),optional::ci_type
    real(dp),intent(in),optional::p
    type(interval_vector)::ans
    character(len=20)::kind
    real(dp)::alpha,z,df
    integer::j,n,levels
    real(dp),allocatable::x(:)
    kind='eta3'
    if(present(ci_type))kind=trim(ci_type)
    alpha=0.025_dp
    if(present(p))alpha=p
    levels=wt%levels()+1
    allocate(ans%estimate(levels),ans%lower(levels),ans%upper(levels))
    z=normal_quantile(1.0_dp-alpha)
    do j=1,levels
      if(j<=wt%levels())then
      x=pack(wt%detail(j)%values,ieee_is_finite(wt%detail(j)%values))
      else
      x=pack(wt%smooth,ieee_is_finite(wt%smooth))
      end if
      n=size(x)
      if(n==0)then
      ans%estimate(j)=0
      ans%lower(j)=0
      ans%upper(j)=0
      cycle
      end if
      ans%estimate(j)=sum(x*x)/real(n,dp)
      select case(kind)
      case('gaussian')
        ans%lower(j)=ans%estimate(j)-z*sqrt(max(2.0_dp*ans%estimate(j)**2/real(n,dp),0.0_dp))
        ans%upper(j)=ans%estimate(j)+z*sqrt(max(2.0_dp*ans%estimate(j)**2/real(n,dp),0.0_dp))
      case('nongaussian')
        ans%lower(j)=ans%estimate(j)-z*sqrt(max(variance_value(x*x)/real(n,dp),0.0_dp))
        ans%upper(j)=ans%estimate(j)+z*sqrt(max(variance_value(x*x)/real(n,dp),0.0_dp))
      case default
        df=max(real(n,dp)/real(2**min(j,wt%levels()),dp),1.0_dp)
        ans%lower(j)=df*ans%estimate(j)/chi_square_quantile(1.0_dp-alpha,df)
        ans%upper(j)=df*ans%estimate(j)/chi_square_quantile(alpha,df)
      end select
    end do
  end function wave_variance

  function wave_covariance(x,y,p) result(ans)
    type(wavelet_transform),intent(in)::x,y
    real(dp),intent(in),optional::p
    type(interval_vector)::ans
    integer::j,n,levels
    real(dp)::alpha,z,se
    real(dp),allocatable::a(:),b(:),a0(:),b0(:)
    alpha=0.025_dp
    if(present(p))alpha=p
    z=normal_quantile(1.0_dp-alpha)
    levels=min(x%levels(),y%levels())+1
    allocate(ans%estimate(levels),ans%lower(levels),ans%upper(levels))
    do j=1,levels
      if(j<=levels-1)then
        a0=x%detail(j)%values
        b0=y%detail(j)%values
      else
        a0=x%smooth
        b0=y%smooth
      end if
      n=min(size(a0),size(b0))
      a=pack(a0(1:n),ieee_is_finite(a0(1:n)).and.ieee_is_finite(b0(1:n)))
      b=pack(b0(1:n),ieee_is_finite(a0(1:n)).and.ieee_is_finite(b0(1:n)))
      n=min(size(a),size(b))
      if(n==0)then
      ans%estimate(j)=0
      ans%lower(j)=0
      ans%upper(j)=0
      cycle
      end if
      ans%estimate(j)=sum(a(1:n)*b(1:n))/real(n,dp)
      se=sqrt(max(variance_value(a(1:n)*b(1:n))/real(n,dp),0.0_dp))
      ans%lower(j)=ans%estimate(j)-z*se
      ans%upper(j)=ans%estimate(j)+z*se
    end do
  end function wave_covariance

  function wave_correlation(x,y,n_original,p) result(ans)
    type(wavelet_transform),intent(in)::x,y
    integer,intent(in),optional::n_original
    real(dp),intent(in),optional::p
    type(interval_vector)::ans
    integer::j,n,levels,n0
    real(dp)::prob,z,r,s
    real(dp),allocatable::a(:),b(:)
    prob=0.975_dp
    if(present(p))prob=p
    z=normal_quantile(prob)
    n0=x%original_length
    if(present(n_original))n0=n_original
    levels=min(x%levels(),y%levels())+1
    allocate(ans%estimate(levels),ans%lower(levels),ans%upper(levels))
    do j=1,levels
      if(j<=levels-1)then
      a=x%detail(j)%values
      b=y%detail(j)%values
      else
      a=x%smooth
      b=y%smooth
      end if
      n=min(size(a),size(b))
      if(n==0)cycle
      r=sum(a(1:n)*b(1:n))/sqrt(max(sum(a(1:n)**2)*sum(b(1:n)**2),tiny(1.0_dp)))
      r=max(-0.999999_dp,min(0.999999_dp,r))
      ans%estimate(j)=r
      n=max(4,n0/2**min(j,levels-1))
      s=z/sqrt(real(n-3,dp))
      ans%lower(j)=tanh(atanh(r)-s)
      ans%upper(j)=tanh(atanh(r)+s)
    end do
  end function wave_correlation

  function spin_covariance(x,y,lag_max) result(c)
    real(dp),intent(in)::x(:),y(:)
    integer,intent(in),optional::lag_max
    real(dp),allocatable::c(:)
    real(dp),allocatable::a(:),b(:)
    integer::lmax,l,n
    a=pack(x,ieee_is_finite(x))
    b=pack(y,ieee_is_finite(y))
    n=min(size(a),size(b))
    lmax=n-1
    if(present(lag_max))lmax=min(lmax,lag_max)
    allocate(c(2*lmax+1))
    do l=-lmax,lmax
    c(l+lmax+1)=lag_product(a(1:n),b(1:n),l)
    end do
  end function spin_covariance

  function spin_correlation(x,y,lag_max) result(c)
    real(dp),intent(in)::x(:),y(:)
    integer,intent(in),optional::lag_max
    real(dp),allocatable::c(:)
    real(dp)::den
    c=spin_covariance(x,y,lag_max)
    den=sqrt(max(mean_square(x)*mean_square(y),tiny(1.0_dp)))
    c=c/den
  end function spin_correlation

  pure function lag_product(x,y,lag) result(v)
    real(dp),intent(in)::x(:),y(:)
    integer,intent(in)::lag
    real(dp)::v
    integer::i,n
    n=min(size(x),size(y))
    v=0.0_dp
    if(lag>=0)then
    do i=1,n-lag
    v=v+x(i)*y(i+lag)
    end do
    else
    do i=1,n+lag
    v=v+x(i-lag)*y(i)
    end do
    end if
    v=v/real(n,dp)
  end function lag_product

  function mean_square(x) result(v)
    real(dp),intent(in)::x(:)
    real(dp)::v
    real(dp),allocatable::a(:)
    a=pack(x,ieee_is_finite(x))
    if(size(a)>0)then
    v=sum(a*a)/real(size(a),dp)
    else
    v=0.0_dp
    end if
  end function mean_square

  function my_acf(x,max_lag) result(v)
    real(dp),intent(in)::x(:)
    integer,intent(in),optional::max_lag
    real(dp),allocatable::v(:)
    integer::l,lmax
    lmax=size(x)-1
    if(present(max_lag))lmax=min(lmax,max_lag)
    allocate(v(0:lmax))
    do l=0,lmax
    v(l)=autocovariance(x,l)
    end do
  end function my_acf

  function my_ccf(x,y,max_lag) result(v)
    real(dp),intent(in)::x(:),y(:)
    integer,intent(in),optional::max_lag
    real(dp),allocatable::v(:)
    integer::l,lmax
    lmax=min(size(x),size(y))-1
    if(present(max_lag))lmax=min(lmax,max_lag)
    allocate(v(-lmax:lmax))
    do l=-lmax,lmax
    v(l)=cross_correlation(x,y,l)
    end do
  end function my_ccf

  function periodogram(z) result(p)
    real(dp),intent(in)::z(:)
    real(dp),allocatable::p(:)
    complex(dp),allocatable::x(:)
    integer::n
    n=size(z)
    allocate(x(n))
    x=cmplx(z,0.0_dp,dp)
    call fft_complex(x)
    allocate(p(n/2+1))
    p=abs(x(1:n/2+1))**2/real(n,dp)
  end function periodogram

  function sine_taper(n,k) result(tapers)
    integer,intent(in)::n,k
    real(dp),allocatable::tapers(:,:)
    integer::i,j
    allocate(tapers(n,k))
    do j=1,k
    do i=1,n
    tapers(i,j)=sqrt(2.0_dp/real(n+1,dp))*sin(pi*real(j*i,dp)/real(n+1,dp))
    end do
    end do
  end function sine_taper

  function dpss_taper(n,k,nw,eigenvalues) result(tapers)
    integer,intent(in)::n,k
    real(dp),intent(in),optional::nw
    real(dp),allocatable,intent(out),optional::eigenvalues(:)
    real(dp),allocatable::tapers(:,:),a(:,:),eval(:),evec(:,:)
    real(dp)::band
    integer::i,j,info,kk
    band=4.0_dp
    if(present(nw))band=nw
    band=band/real(n,dp)
    allocate(a(n,n))
    do i=1,n
    do j=1,n
    if(i==j)then
    a(i,j)=2.0_dp*band
    else
    a(i,j)=sin(2.0_dp*pi*band*real(i-j,dp))/(pi*real(i-j,dp))
    end if
    end do
    end do
    call symmetric_eigen(a,eval,evec,info)
    kk=min(k,n)
    allocate(tapers(n,kk))
    tapers=evec(:,1:kk)
    if(present(eigenvalues))eigenvalues=eval(1:kk)
  end function dpss_taper

  function wave_variance_2d(wt,p) result(ans)
    use waveslim_types,only:wavelet_transform_2d
    type(wavelet_transform_2d),intent(in)::wt
    real(dp),intent(in),optional::p
    real(dp),allocatable::ans(:,:,:)
    integer::j,n
    real(dp)::alpha,z
    alpha=0.025_dp
    if(present(p))alpha=p
    z=normal_quantile(1.0_dp-alpha)
    allocate(ans(size(wt%level),3,3))
    do j=1,size(wt%level)
    n=size(wt%level(j)%lh)
    ans(j,1,1)=sum(wt%level(j)%lh**2)/real(n,dp)
    ans(j,2,1)=sum(wt%level(j)%hl**2)/real(n,dp)
    ans(j,3,1)=sum(wt%level(j)%hh**2)/real(n,dp)
    ans(j,:,2)=max(0.0_dp,ans(j,:,1)-z*sqrt(2.0_dp/real(n,dp))*ans(j,:,1))
    ans(j,:,3)=ans(j,:,1)+z*sqrt(2.0_dp/real(n,dp))*ans(j,:,1)
    end do
  end function wave_variance_2d

  function rotcumvar(x) result(ans)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable::ans(:,:)
    real(dp),allocatable::cov(:,:),eval(:),evec(:,:)
    integer::n,p,i,j,info
    n=size(x,1)
    p=size(x,2)
    allocate(cov(p,p))
    do i=1,p
    do j=1,p
    cov(i,j)=sum((x(:,i)-mean_value(x(:,i)))*(x(:,j)-mean_value(x(:,j))))/real(max(1,n-1),dp)
    end do
    end do
    call symmetric_eigen(cov,eval,evec,info)
    allocate(ans(p,2))
    ans(:,1)=eval/sum(eval)
    do i=1,p
    ans(i,2)=sum(ans(1:i,1))
    end do
  end function rotcumvar

  function wavelet_filter(wf, filter_sequence) result(coefficients)
    use waveslim_filters, only : wave_filter
    use waveslim_status, only : status_type
    use waveslim_types, only : wavelet_filter_type
    character(len=*), intent(in) :: wf
    character(len=*), intent(in), optional :: filter_sequence
    real(dp), allocatable :: coefficients(:)
    type(wavelet_filter_type) :: filter
    type(status_type) :: status
    character(len=:), allocatable :: sequence
    real(dp), allocatable :: work(:), selected(:)
    integer :: stage, count

    sequence = 'L'
    if (present(filter_sequence)) sequence = trim(filter_sequence)
    filter = wave_filter(wf,status)
    if (.not. status%ok()) then
      allocate(coefficients(0))
      return
    end if
    allocate(work(1))
    work = 1.0_dp
    count = len_trim(sequence)
    do stage = 1, count
      select case (sequence(count-stage+1:count-stage+1))
      case ('H','h')
        selected = filter%hpf
      case ('L','l')
        selected = filter%lpf
      case default
        allocate(coefficients(0))
        return
      end select
      work = cascade_filter(selected,work,stage-1)
    end do
    coefficients = work
  end function wavelet_filter

  function cascade_filter(filter, x, stage) result(y)
    real(dp), intent(in) :: filter(:), x(:)
    integer, intent(in) :: stage
    real(dp), allocatable :: y(:), padded(:), work(:)
    integer :: filter_length, input_length, offset, padded_length
    integer :: output_length, stride, i

    filter_length = size(filter)
    input_length = size(x)
    stride = 2**stage
    offset = (filter_length-1)*stride
    padded_length = offset-filter_length+2
    output_length = 2*offset-filter_length+2
    if (input_length > padded_length) then
      allocate(y(0))
      return
    end if
    allocate(padded(padded_length))
    padded = 0.0_dp
    padded(1:input_length) = x
    allocate(work(2*offset+padded_length))
    work = 0.0_dp
    work(offset+1:offset+padded_length) = padded
    allocate(y(output_length))
    y = 0.0_dp
    do i = 1, filter_length
      y = y + filter(filter_length-i+1)* &
        work(1+(i-1)*stride:output_length+(i-1)*stride)
    end do
  end function cascade_filter

  function squared_gain(wf,filter_sequence,n) result(gain)
    character(len=*), intent(in) :: wf
    character(len=*), intent(in), optional :: filter_sequence
    integer, intent(in), optional :: n
    real(dp), allocatable :: gain(:)
    real(dp), allocatable :: coefficients(:)
    complex(dp), allocatable :: transformed(:)
    integer :: length_fft
    character(len=:), allocatable :: sequence

    length_fft = 512
    if (present(n)) length_fft = n
    sequence = 'L'
    if (present(filter_sequence)) sequence = trim(filter_sequence)
    coefficients = wavelet_filter(wf,sequence)
    if (size(coefficients) == 0 .or. size(coefficients) > length_fft .or. &
        iand(length_fft,length_fft-1) /= 0) then
      allocate(gain(0))
      return
    end if
    allocate(transformed(length_fft))
    transformed = cmplx(0.0_dp,0.0_dp,dp)
    transformed(1:size(coefficients)) = cmplx(coefficients,0.0_dp,dp)
    call fft_complex(transformed)
    gain = abs(transformed(1:length_fft/2+1))**2
  end function squared_gain

  function chi_square_quantile(p,df) result(x)
    real(dp),intent(in)::p,df
    real(dp)::x,lo,hi,mid
    integer::i
    lo=0.0_dp
    hi=max(10.0_dp,df+10.0_dp*sqrt(2.0_dp*df))
    do while(chi_square_cdf(hi,df)<p)
    hi=2*hi
    end do
    do i=1,100
    mid=(lo+hi)/2
    if(chi_square_cdf(mid,df)<p)then
    lo=mid
    else
    hi=mid
    end if
    end do
    x=(lo+hi)/2
  end function chi_square_quantile
end module waveslim_statistics

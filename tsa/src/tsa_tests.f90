! SPDX-License-Identifier: GPL-2.0-or-later
module tsa_tests
  use tsa_kinds, only : dp
  use tsa_types, only : tsa_test_result, outlier_result
  use tsa_utils, only : variance_n, mean_value, sd_n, normal_quantile, arma_to_ma, build_lag_matrix
  use tsa_statistics, only : autocorrelation
  use tsa_arma, only : ar_ols_fit
  use tseries_linalg, only : least_squares
  use tseries_special, only : chi_square_cdf, f_cdf
  implicit none
  private
  public :: lb_test, mcleod_li_test, keenan_test, tsay_test, detect_io, detect_ao

contains

  function lb_test(residuals,lag,nparm,ljung_box) result(res)
    real(dp),intent(in)::residuals(:)
    integer,intent(in)::lag,nparm
    logical,intent(in),optional::ljung_box
    type(tsa_test_result)::res
    real(dp),allocatable::ac(:)
    logical::lb
    integer::k,n
    lb=.true.
    if(present(ljung_box))lb=ljung_box
    n=size(residuals)
    res%lag=lag
    res%df=lag-nparm
    if(lag<=nparm .or. lag>=n .or. res%df<=0)then
    res%status=1
    res%p_value=1.0_dp
    return
    end if
    call autocorrelation(residuals,lag,ac)
    if(lb)then
      res%statistic=0.0_dp
      do k=1,lag
      res%statistic=res%statistic+ac(k)**2/real(n-k,dp)
      end do
      res%statistic=real(n*(n+2),dp)*res%statistic
    else
      res%statistic=real(n,dp)*sum(ac**2)
    end if
    res%p_value=1.0_dp-chi_square_cdf(res%statistic,real(res%df,dp))
  end function lb_test

  subroutine mcleod_li_test(residuals,lag_max,p_values)
    real(dp),intent(in)::residuals(:)
    integer,intent(in)::lag_max
    real(dp),allocatable,intent(out)::p_values(:)
    type(tsa_test_result)::r
    integer::i
    allocate(p_values(lag_max))
    do i=1,lag_max
    r=lb_test(residuals**2,i,0,.true.)
    p_values(i)=r%p_value
    end do
  end subroutine mcleod_li_test

  function keenan_test(x, order) result(res)
    real(dp), intent(in) :: x(:)
    integer, intent(in), optional :: order
    type(tsa_test_result) :: res
    integer :: m, n, status
    real(dp), allocatable :: y(:), xx(:,:), b(:), r1(:), fit1(:)
    real(dp), allocatable :: b2(:), r2(:), z(:,:), b3(:), r3(:)

    n = size(x)
    if (present(order)) then
      m = order
    else
      m = choose_order(x)
    end if
    res%order = m
    if (m < 1 .or. n <= 2*m + 3) then
      res%status = 1
      return
    end if

    call build_lag_matrix(x, m, m+1, y, xx, .false.)
    allocate(b(size(xx,2)), r1(size(y)))
    call least_squares(xx, y, b, residuals=r1, status=status)
    if (status /= 0) then
      res%status = 2
      return
    end if
    fit1 = (y-r1)**2

    allocate(b2(size(xx,2)), r2(size(y)))
    call least_squares(xx, fit1, b2, residuals=r2, status=status)
    if (status /= 0) then
      res%status = 3
      return
    end if

    allocate(z(size(r2),1), b3(1), r3(size(r2)))
    z(:,1) = r2
    call least_squares(z, r1, b3, residuals=r3, status=status)
    if (status /= 0) then
      res%status = 4
      return
    end if
    res%statistic = (sum(r1**2)-sum(r3**2)) / &
      max(sum(r3**2)/real(size(r3)-1,dp), tiny(1.0_dp))
    res%statistic = res%statistic * real(n-2*m-2,dp) / real(n-m-1,dp)
    res%df = 1
    res%p_value = 1.0_dp - f_cdf(res%statistic, 1.0_dp, real(n-2*m-2,dp))
  contains
    integer function choose_order(v) result(p)
      use tsa_types, only : ar_fit_result
      real(dp), intent(in) :: v(:)
      type(ar_fit_result) :: af
      af = ar_ols_fit(v, min(10,size(v)/5), select_aic=.true.)
      p = max(1, af%order)
    end function choose_order
  end function keenan_test

  function tsay_test(x,order) result(res)
    real(dp),intent(in)::x(:)
    integer,intent(in),optional::order
    type(tsa_test_result)::res
    integer::m,n,nobs,p1,p2,status,i,j,k,col,start
    real(dp),allocatable::y(:),x1(:,:),x2(:,:),b1(:),b2(:),r1(:),r2(:)
    n=size(x)
    if(present(order))then
    m=order
    else
    m=choose_order(x)
    end if
    res%order=m
    if(m<1 .or. n<=m+3)then
    res%status=1
    return
    end if
    start=m+1
    nobs=n-start+1
    p1=m+1
    p2=m+m*(m+1)/2+1
    allocate(y(nobs),x1(nobs,p1),x2(nobs,p2))
    y=x(start:)
    x1(:,1)=1.0_dp
    x2(:,1)=1.0_dp
    do i=1,nobs
      do j=1,m
      x1(i,1+j)=x(start+i-1-j)
      x2(i,1+j)=x(start+i-1-j)
      end do
      col=1+m
      do j=1,m
      do k=1,j
      col=col+1
      x2(i,col)=x(start+i-1-j)*x(start+i-1-k)
      end do
      end do
    end do
    allocate(b1(p1), b2(p2), r1(nobs), r2(nobs))
    call least_squares(x1, y, b1, residuals=r1, status=status)
    if (status /= 0) then
      res%status = 2
      return
    end if
    call least_squares(x2,y,b2,residuals=r2,status=status)
    if(status/=0)then
    res%status=3
    return
    end if
    res%df=p2-p1
    res%statistic=((sum(r1**2)-sum(r2**2))/real(res%df,dp))/(sum(r2**2)/real(nobs-p2,dp))
    res%p_value=1.0_dp-f_cdf(res%statistic,real(res%df,dp),real(nobs-p2,dp))
  contains
    integer function choose_order(v) result(p)
      use tsa_types,only:ar_fit_result
      real(dp),intent(in)::v(:)
      type(ar_fit_result)::af
      af=ar_ols_fit(v,min(10,size(v)/5),select_aic=.true.)
      p=max(1,af%order)
    end function choose_order
  end function tsay_test

  function detect_io(residuals,sigma2,alpha,robust) result(out)
    real(dp),intent(in)::residuals(:),sigma2
    real(dp),intent(in),optional::alpha
    logical,intent(in),optional::robust
    type(outlier_result)::out
    real(dp)::a,sigma,cut,pi
    logical::rb
    real(dp),allocatable::stat(:)
    integer::i,nfound
    a=0.05_dp
    if(present(alpha))a=alpha
    rb=.true.
    if(present(robust))rb=robust
    pi=acos(-1.0_dp)
    if(rb)then
    sigma=sqrt(pi/2.0_dp)*sum(abs(residuals))/real(size(residuals),dp)
    else
    sigma=sqrt(max(sigma2,tiny(1.0_dp)))
    end if
    allocate(stat(size(residuals)))
    stat = residuals/max(sigma,tiny(1.0_dp))
    cut = normal_quantile(1.0_dp-a/(2.0_dp*real(size(stat),dp)))
    out%cutoff = cut
    nfound=count(abs(stat)>cut)
    allocate(out%index(nfound),out%statistic(nfound))
    nfound=0
    do i=1,size(stat)
    if(abs(stat(i))>cut)then
    nfound=nfound+1
    out%index(nfound)=i
    out%statistic(nfound)=stat(i)
    end if
    end do
  end function detect_io

  function detect_ao(residuals,ar,ma,sigma2,alpha,robust) result(out)
    real(dp),intent(in)::residuals(:),ar(:),ma(:),sigma2
    real(dp),intent(in),optional::alpha
    logical,intent(in),optional::robust
    type(outlier_result)::out
    real(dp)::a,sigma,cut,pi
    logical::rb
    real(dp),allocatable::piwt(:),piseq(:),rho2(:),omega(:),stat(:)
    integer::n,t,j,nfound
    n=size(residuals)
    allocate(piwt(max(0,n-1)))
    if(n>1)call arma_to_ma(-ma,-ar,n-1,piwt)
    allocate(piseq(n))
    piseq(1)=1.0_dp
    if(n>1)piseq(2:)=piwt
    allocate(rho2(n))
    rho2=0.0_dp
    rho2(1)=1.0_dp/(piseq(1)**2)
    do t=2,n
    rho2(t)=1.0_dp/sum(piseq(1:t)**2)
    end do
    allocate(omega(n),stat(n))
    omega=0.0_dp
    ! Equivalent to filtering the reversed residuals by pi-weights and reversing back.
    do t=1,n
      do j=1,n-t+1
      omega(t)=omega(t)+piseq(j)*residuals(t+j-1)
      end do
      omega(t)=omega(t)*rho2(n-t+1)
    end do
    pi=acos(-1.0_dp)
    a=0.05_dp
    if(present(alpha))a=alpha
    rb=.true.
    if(present(robust))rb=robust
    if(rb)then
    sigma=sqrt(pi/2.0_dp)*sum(abs(residuals))/real(n,dp)
    else
    sigma=sqrt(max(sigma2,tiny(1.0_dp)))
    end if
    do t=1,n
    stat(t)=omega(t)/max(sigma*sqrt(rho2(n-t+1)),tiny(1.0_dp))
    end do
    cut = normal_quantile(1.0_dp-a/(2.0_dp*real(n,dp)))
    out%cutoff = cut
    nfound = count(abs(stat)>cut)
    allocate(out%index(nfound), out%statistic(nfound))
    nfound = 0
    do t=1,n
    if(abs(stat(t))>cut)then
    nfound=nfound+1
    out%index(nfound)=t
    out%statistic(nfound)=stat(t)
    end if
    end do
  end function detect_ao
end module tsa_tests

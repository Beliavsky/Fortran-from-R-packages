! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2016 Marius Hofert, Kurt Hornik and Alexander J. McNeil
module qrmtools_tests
  use qrmtools_kinds, only : dp
  use qrmtools_types, only : test_result
  use qrmtools_stats, only : covariance_matrix, invert_matrix, chi_square_cdf, beta_cdf, normal_cdf, sort_increasing
  implicit none
  private
  public :: mahalanobis_squared, maha2_test, mardia_test
contains
  function mahalanobis_squared(x) result(distances)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: distances(:)
    real(dp), allocatable :: covariance(:,:),inverse(:,:),means(:),diff(:)
    logical :: ok
    integer :: i
    means=sum(x,dim=1)/real(size(x,1),dp); covariance=covariance_matrix(x)
    call invert_matrix(covariance,inverse,ok); allocate(distances(size(x,1)))
    if(.not.ok) then; distances=huge(1.0_dp); return; end if
    allocate(diff(size(x,2)))
    do i=1,size(x,1); diff=x(i,:)-means; distances(i)=dot_product(diff,matmul(inverse,diff)); end do
  end function mahalanobis_squared

  function maha2_test(x,use_beta) result(output)
    real(dp), intent(in) :: x(:,:)
    logical, intent(in), optional :: use_beta
    type(test_result) :: output
    real(dp), allocatable :: d2(:),u(:)
    real(dp) :: dstat,p,expected
    integer :: i,n,d
    logical :: beta
    n=size(x,1); d=size(x,2); beta=.false.; if(present(use_beta))beta=use_beta
    if(n<=d+1) then; output%message='Too few observations for the Mahalanobis test.'; return; end if
    d2=mahalanobis_squared(x); call sort_increasing(d2); allocate(u(n)); dstat=0.0_dp
    do i=1,n
      if(beta) then
        p=beta_cdf(d2(i)*real(n,dp)/real((n-1)*(n-1),dp),0.5_dp*real(d,dp),0.5_dp*real(n-d-1,dp))
      else
        p=chi_square_cdf(d2(i),real(d,dp))
      end if
      expected=real(i,dp)/real(n,dp); dstat=max(dstat,abs(p-expected)); u(i)=p
    end do
    output%statistic=dstat; output%p_value=ks_pvalue(dstat,n); output%distances=d2; output%ok=.true.
  end function maha2_test

  function mardia_test(x,skewness) result(output)
    real(dp), intent(in) :: x(:,:)
    logical, intent(in), optional :: skewness
    type(test_result) :: output
    real(dp), allocatable :: covariance(:,:),inverse(:,:),centered(:,:),angles(:,:),d2(:)
    real(dp) :: b,t
    logical :: use_skew,ok
    integer :: n,d,i,df
    n=size(x,1); d=size(x,2); use_skew=.false.; if(present(skewness))use_skew=skewness
    covariance=covariance_matrix(x); call invert_matrix(covariance,inverse,ok)
    if(.not.ok) then; output%message='Covariance matrix is singular.'; return; end if
    allocate(centered(n,d)); centered=x-spread(sum(x,dim=1)/real(n,dp),1,n)
    if(use_skew) then
      allocate(angles(n,n)); angles=matmul(matmul(centered,inverse),transpose(centered))
      b=sum(angles**3)/real(n*n,dp); t=real(n,dp)*b/6.0_dp
      df=d*(d+1)*(d+2)/6; output%statistic=t; output%p_value=1.0_dp-chi_square_cdf(t,real(df,dp))
    else
      allocate(d2(n)); do i=1,n; d2(i)=dot_product(centered(i,:),matmul(inverse,centered(i,:))); end do
      b=sum(d2**2)/real(n,dp); t=(b-real(d*(d+2),dp))/sqrt(8.0_dp*real(d*(d+2),dp)/real(n,dp))
      output%statistic=t; output%p_value=2.0_dp*(1.0_dp-normal_cdf(abs(t)))
    end if
    output%ok=.true.
  end function mardia_test

  pure real(dp) function ks_pvalue(d,n) result(value)
    real(dp), intent(in) :: d
    integer, intent(in) :: n
    real(dp) :: lambda,sumv,term
    integer :: k
    lambda=(sqrt(real(n,dp))+0.12_dp+0.11_dp/sqrt(real(n,dp)))*d; sumv=0.0_dp
    do k=1,100
      term=2.0_dp*(-1.0_dp)**(k-1)*exp(-2.0_dp*real(k*k,dp)*lambda*lambda); sumv=sumv+term
      if(abs(term)<1.0e-14_dp)exit
    end do
    value=min(max(sumv,0.0_dp),1.0_dp)
  end function ks_pvalue
end module qrmtools_tests

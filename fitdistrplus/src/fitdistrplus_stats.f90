! SPDX-License-Identifier: GPL-2.0-or-later
module fitdistrplus_stats
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fitdistrplus_kinds, only : dp
  use fitdistrplus_types
  use fitdistrplus_math, only : type7_quantile, chi_square_survival
  implicit none
  private

  public :: descdist, gofstat, detectbound, parameter_quantiles

contains

  subroutine descdist(data,result,unbiased)
    real(dp),intent(in)::data(:)
    type(descriptive_result),intent(out)::result
    logical,intent(in),optional::unbiased
    real(dp)::m,m2,m3,m4,nr
    logical::ub
    integer::n
    ub=.true.;if(present(unbiased))ub=unbiased
    n=size(data);result%n=n;result%unbiased=ub
    if(n<4 .or. any(.not.ieee_is_finite(data)))then
      result%status=fit_invalid_argument;return
    end if
    nr=real(n,dp);m=sum(data)/nr
    m2=sum((data-m)**2)/nr;m3=sum((data-m)**3)/nr;m4=sum((data-m)**4)/nr
    result%minimum=minval(data);result%maximum=maxval(data)
    result%median=type7_quantile(data,0.5_dp);result%mean=m
    if(ub)then
      result%standard_deviation=sqrt(sum((data-m)**2)/real(n-1,dp))
      result%skewness=sqrt(nr*(nr-1.0_dp))/(nr-2.0_dp)*m3/max(m2**1.5_dp,tiny(1.0_dp))
      result%kurtosis=(nr-1.0_dp)/((nr-2.0_dp)*(nr-3.0_dp))* &
        ((nr+1.0_dp)*m4/max(m2*m2,tiny(1.0_dp))-3.0_dp*(nr-1.0_dp))+3.0_dp
    else
      result%standard_deviation=sqrt(m2)
      result%skewness=m3/max(m2**1.5_dp,tiny(1.0_dp))
      result%kurtosis=m4/max(m2*m2,tiny(1.0_dp))
    end if
    result%status=fit_success
  end subroutine descdist

  subroutine gofstat(data,dist,fit,result,breaks)
    real(dp),intent(in)::data(:)
    type(distribution_model),intent(in)::dist
    type(fit_result),intent(in)::fit
    type(gof_result),intent(out)::result
    real(dp),intent(in),optional::breaks(:)
    real(dp),allocatable::s(:),p(:),boundaries(:)
    real(dp)::expected,observed,cdf_lo,cdf_hi
    integer::n,i,j,k,cells

    n=size(data)
    if(n<2 .or. .not.associated(dist%cdf) .or. .not.allocated(fit%estimate))then
      result%status=fit_invalid_argument;return
    end if
    s=data;call sort_real_local(s);allocate(p(n))
    do i=1,n;p(i)=dist%cdf(s(i),fit%estimate);end do
    if(any(p<=0.0_dp) .or. any(p>=1.0_dp))then
      result%status=fit_numerical_error;return
    end if
    result%ks=maxval(max(abs(p-[(real(i,dp)/real(n,dp),i=1,n)]), &
      abs(p-[(real(i-1,dp)/real(n,dp),i=1,n)])))
    result%cvm=1.0_dp/(12.0_dp*real(n,dp))+ &
      sum((p-[(real(2*i-1,dp)/(2.0_dp*real(n,dp)),i=1,n)])**2)
    result%ad=-real(n,dp)-sum([(real(2*i-1,dp)*(log(p(i))+ &
      log(1.0_dp-p(n+1-i))),i=1,n)])/real(n,dp)
    result%aic=fit%aic;result%bic=fit%bic

    if(present(breaks))then
      if(size(breaks)<1)then;result%status=fit_invalid_argument;return;end if
      boundaries=breaks;call sort_real_local(boundaries);cells=size(boundaries)+1
      result%chi_square=0.0_dp
      do k=1,cells
        if(k==1)then
          cdf_lo=0.0_dp
        else
          cdf_lo=dist%cdf(boundaries(k-1),fit%estimate)
        end if
        if(k==cells)then
          cdf_hi=1.0_dp
        else
          cdf_hi=dist%cdf(boundaries(k),fit%estimate)
        end if
        expected=real(n,dp)*max(cdf_hi-cdf_lo,0.0_dp)
        observed=0.0_dp
        do j=1,n
          if(k==1)then
            if(data(j)<=boundaries(1))observed=observed+1.0_dp
          else if(k==cells)then
            if(data(j)>boundaries(cells-1))observed=observed+1.0_dp
          else
            if(data(j)>boundaries(k-1) .and. data(j)<=boundaries(k))observed=observed+1.0_dp
          end if
        end do
        if(expected>tiny(1.0_dp))result%chi_square=result%chi_square+(observed-expected)**2/expected
      end do
      result%chi_square_df=cells-1-size(fit%estimate)
      result%chi_square_pvalue=chi_square_survival(result%chi_square,result%chi_square_df)
    end if
    result%status=fit_success
  end subroutine gofstat

  subroutine detectbound(dist,lower,upper,status)
    type(distribution_model),intent(in)::dist
    real(dp),allocatable,intent(out)::lower(:),upper(:)
    integer,intent(out)::status
    if(dist%npar<=0 .or. .not.allocated(dist%default_lower) .or. &
       .not.allocated(dist%default_upper))then
      allocate(lower(0),upper(0));status=fit_not_supported;return
    end if
    lower=dist%default_lower;upper=dist%default_upper;status=fit_success
  end subroutine detectbound

  subroutine parameter_quantiles(draws,probs,quantiles,status)
    real(dp),intent(in)::draws(:,:),probs(:)
    real(dp),allocatable,intent(out)::quantiles(:,:)
    integer,intent(out)::status
    integer::i,j
    if(size(draws,1)==0 .or. any(probs<0.0_dp) .or. any(probs>1.0_dp))then
      allocate(quantiles(0,0));status=fit_invalid_argument;return
    end if
    allocate(quantiles(size(probs),size(draws,2)))
    do j=1,size(draws,2)
      do i=1,size(probs)
        quantiles(i,j)=type7_quantile(draws(:,j),probs(i))
      end do
    end do
    status=fit_success
  end subroutine parameter_quantiles

  subroutine sort_real_local(x)
    real(dp),intent(inout)::x(:)
    real(dp)::key
    integer::i,j
    do i=2,size(x);key=x(i);j=i-1
      do while(j>=1)
        if(x(j)<=key)exit
        x(j+1)=x(j);j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real_local

end module fitdistrplus_stats

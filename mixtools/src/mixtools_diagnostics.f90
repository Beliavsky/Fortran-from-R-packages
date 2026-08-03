! SPDX-License-Identifier: GPL-2.0-or-later
module mixtools_diagnostics
  use mixtools_kinds, only : dp
  use mixtools_status
  use mixtools_types
  use mixtools_parametric, only : normalmix_em
  use mixtools_regression, only : regmix_em, regmix_em_mixed
  use mixtools_utilities, only : wkde
  implicit none
  private
  public :: component_cdf, fdr_from_posterior, density_from_semiparametric
  public :: integrated_squared_error, test_equality_normal, test_equality_regression
  public :: test_equality_mixed
contains
  subroutine component_cdf(data,weights,grid,cdf,status)
    real(dp),intent(in)::data(:),weights(:,:),grid(:)
    real(dp),intent(out)::cdf(size(grid),size(weights,2))
    integer,intent(out)::status
    integer::i,j,h
    real(dp)::den
    if(size(weights,1)/=size(data))then;cdf=0.0_dp;status=MIXTOOLS_DIMENSION_ERROR;return;end if
    do j=1,size(weights,2)
      den=sum(weights(:,j));do i=1,size(grid);cdf(i,j)=0.0_dp
        do h=1,size(data);if(data(h)<=grid(i))cdf(i,j)=cdf(i,j)+weights(h,j);end do
        if(den>0.0_dp)cdf(i,j)=cdf(i,j)/den
      end do
    end do
    status=MIXTOOLS_SUCCESS
  end subroutine component_cdf

  subroutine fdr_from_posterior(posterior,null_component,order,fdr,status)
    real(dp),intent(in)::posterior(:,:)
    integer,intent(in)::null_component
    integer,intent(out)::order(size(posterior,1))
    real(dp),intent(out)::fdr(size(posterior,1))
    integer,intent(out)::status
    integer::i,j,n,tmp
    real(dp)::cum
    n=size(posterior,1)
    if(null_component<1.or.null_component>size(posterior,2))then
      order=0;fdr=1.0_dp;status=MIXTOOLS_INVALID_ARGUMENT;return
    end if
    do i=1,n;order(i)=i;end do
    do i=2,n
      tmp=order(i);j=i-1
      do while(j>=1)
        if(posterior(order(j),null_component)<=posterior(tmp,null_component))exit
        order(j+1)=order(j);j=j-1
      end do
      order(j+1)=tmp
    end do
    cum=0.0_dp
    do i=1,n;cum=cum+posterior(order(i),null_component);fdr(i)=cum/real(i,dp);end do
    status=MIXTOOLS_SUCCESS
  end subroutine fdr_from_posterior

  subroutine density_from_semiparametric(result,component,grid,density,status)
    type(semiparametric_result),intent(in)::result
    integer,intent(in)::component
    real(dp),intent(in)::grid(:)
    real(dp),intent(out)::density(size(grid))
    integer,intent(out)::status
    integer::i,j,m
    real(dp)::t
    if(component<1.or.component>size(result%density,1).or.size(result%grid)<2)then
      density=0.0_dp;status=MIXTOOLS_INVALID_ARGUMENT;return
    end if
    m=size(result%grid)
    do i=1,size(grid)
      if(grid(i)<=result%grid(1))then;density(i)=result%density(component,1)
      else if(grid(i)>=result%grid(m))then;density(i)=result%density(component,m)
      else
        do j=1,m-1;if(grid(i)<=result%grid(j+1))exit;end do
        t=(grid(i)-result%grid(j))/(result%grid(j+1)-result%grid(j))
        density(i)=(1.0_dp-t)*result%density(component,j)+t*result%density(component,j+1)
      end if
    end do
    status=MIXTOOLS_SUCCESS
  end subroutine density_from_semiparametric

  function integrated_squared_error(grid,estimate,truth) result(value)
    real(dp),intent(in)::grid(:),estimate(:),truth(:)
    real(dp)::value
    integer::i
    value=0.0_dp
    if(size(grid)/=size(estimate).or.size(grid)/=size(truth))then;value=huge(1.0_dp);return;end if
    do i=1,size(grid)-1
      value=value+0.5_dp*(grid(i+1)-grid(i))*((estimate(i)-truth(i))**2+(estimate(i+1)-truth(i+1))**2)
    end do
  end function integrated_squared_error

  subroutine test_equality_normal(x,k,which_test,statistic,p_value,control)
    real(dp),intent(in)::x(:)
    integer,intent(in)::k,which_test
    real(dp),intent(out)::statistic,p_value
    type(em_control),intent(in),optional::control
    type(mixture_result)::nullfit,altfit
    integer::df
    select case(which_test)
    case(1)
      call normalmix_em(x,k,nullfit,control,common_mean=.true.)
      df=k-1
    case(2)
      call normalmix_em(x,k,nullfit,control,common_sigma=.true.)
      df=k-1
    case default
      call normalmix_em(x,k,nullfit,control,common_mean=.true.,common_sigma=.true.)
      df=2*(k-1)
    end select
    call normalmix_em(x,k,altfit,control)
    statistic=max(0.0_dp,2.0_dp*(altfit%loglik-nullfit%loglik));call chisq_tail(statistic,df,p_value)
  end subroutine test_equality_normal

  subroutine test_equality_regression(y,x,k,which_test,statistic,p_value,control)
    real(dp),intent(in)::y(:),x(:,:)
    integer,intent(in)::k,which_test
    real(dp),intent(out)::statistic,p_value
    type(em_control),intent(in),optional::control
    type(regression_mixture_result)::nullfit,altfit
    integer::df
    select case(which_test)
    case(1)
      call regmix_em(y,x,k,nullfit,control,.true.,common_beta=.true.)
      df=(k-1)*(size(x,2)+1)
    case(2)
      call regmix_em(y,x,k,nullfit,control,.true.,common_sigma=.true.)
      df=k-1
    case default
      call regmix_em(y,x,k,nullfit,control,.true.,common_beta=.true.,common_sigma=.true.)
      df=(k-1)*(size(x,2)+2)
    end select
    call regmix_em(y,x,k,altfit,control,.true.)
    statistic=max(0.0_dp,2.0_dp*(altfit%loglik-nullfit%loglik));call chisq_tail(statistic,df,p_value)
  end subroutine test_equality_regression

  subroutine test_equality_mixed(y,x,groups,k,statistic,p_value,control)
    real(dp),intent(in)::y(:),x(:,:)
    integer,intent(in)::groups(:),k
    real(dp),intent(out)::statistic,p_value
    type(em_control),intent(in),optional::control
    type(regression_mixture_result)::mixed,ordinary
    real(dp),allocatable::re(:,:)
    integer::df
    call regmix_em_mixed(y,x,groups,k,mixed,re,control,.true.)
    call regmix_em(y,x,k,ordinary,control,.true.)
    statistic=max(0.0_dp,2.0_dp*(mixed%loglik-ordinary%loglik));df=maxval(groups)*k
    call chisq_tail(statistic,df,p_value)
  end subroutine test_equality_mixed

  subroutine chisq_tail(x,df,p)
    real(dp),intent(in)::x
    integer,intent(in)::df
    real(dp),intent(out)::p
    real(dp)::z
    if(df<=0)then;p=1.0_dp;return;end if
    z=((max(x,0.0_dp)/real(df,dp))**(1.0_dp/3.0_dp)-(1.0_dp-2.0_dp/(9.0_dp*real(df,dp)))) &
      /sqrt(2.0_dp/(9.0_dp*real(df,dp)))
    p=0.5_dp*erfc(z/sqrt(2.0_dp))
  end subroutine chisq_tail
end module mixtools_diagnostics

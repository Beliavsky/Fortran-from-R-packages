! SPDX-License-Identifier: GPL-2.0-or-later
module fitdistrplus_bootstrap
  use fitdistrplus_kinds, only : dp
  use fitdistrplus_types
  use fitdistrplus_math, only : seed_rng, type7_quantile
  use fitdistrplus_fit, only : fitdist, mledist_censored
  implicit none
  private

  public :: bootdist, bootdistcens, turnbull_npmle, cdf_bootstrap_band

contains

  subroutine bootdist(data,dist,fit,nsim,result,parametric,seed,control,orders,probs,gof,phidiv,power)
    real(dp),intent(in)::data(:)
    type(distribution_model),intent(in)::dist
    type(fit_result),intent(in)::fit
    integer,intent(in)::nsim
    type(bootstrap_result),intent(out)::result
    logical,intent(in),optional::parametric
    integer,intent(in),optional::seed
    type(fit_control),intent(in),optional::control
    integer,intent(in),optional::orders(:),gof,phidiv
    real(dp),intent(in),optional::probs(:),power
    logical::parboot
    real(dp),allocatable::sample(:)
    type(fit_result)::refit
    integer::b,i,index,method
    real(dp)::u

    parboot=.true.;if(present(parametric))parboot=parametric
    if(nsim<1 .or. size(data)<2 .or. .not.allocated(fit%estimate))then
      result%status=fit_invalid_argument;return
    end if
    if(parboot .and. .not.associated(dist%random_value))then
      result%status=fit_not_supported;return
    end if
    method=method_code(fit%method)
    if(method==0)then;result%status=fit_not_supported;return;end if
    call seed_rng(seed)
    allocate(result%estimates(nsim,size(fit%estimate)),result%convergence(nsim),sample(size(data)))
    result%estimates=0.0_dp;result%convergence=fit_invalid_argument
    do b=1,nsim
      if(parboot)then
        do i=1,size(data);sample(i)=dist%random_value(fit%estimate);end do
      else
        do i=1,size(data)
          call random_number(u);index=1+min(int(u*real(size(data),dp)),size(data)-1)
          sample(i)=data(index)
        end do
      end if
      call fitdist(sample,dist,method,fit%estimate,refit,control,orders=orders,probs=probs, &
        gof=gof,phidiv=phidiv,power=power)
      result%convergence(b)=refit%convergence
      if(refit%convergence==fit_success)result%estimates(b,:)=refit%estimate
    end do
    result%successful=count(result%convergence==fit_success);result%nsim=nsim
    allocate(result%mean(size(fit%estimate)),result%standard_deviation(size(fit%estimate)))
    if(result%successful>0)then
      do i=1,size(fit%estimate)
        result%mean(i)=sum(result%estimates(:,i),mask=result%convergence==fit_success)/real(result%successful,dp)
        if(result%successful>1)then
          result%standard_deviation(i)=sqrt(sum((result%estimates(:,i)-result%mean(i))**2, &
            mask=result%convergence==fit_success)/real(result%successful-1,dp))
        else
          result%standard_deviation(i)=0.0_dp
        end if
      end do
      result%status=fit_success
    else
      result%mean=0.0_dp;result%standard_deviation=0.0_dp;result%status=fit_no_convergence
    end if
  end subroutine bootdist

  subroutine bootdistcens(sample,dist,fit,nsim,result,seed,control)
    type(censored_sample),intent(in)::sample
    type(distribution_model),intent(in)::dist
    type(fit_result),intent(in)::fit
    integer,intent(in)::nsim
    type(bootstrap_result),intent(out)::result
    integer,intent(in),optional::seed
    type(fit_control),intent(in),optional::control
    type(censored_sample)::draw
    type(fit_result)::refit
    integer::b,i,index,n
    real(dp)::u

    n=size(sample%left)
    if(n<2 .or. size(sample%right)/=n .or. nsim<1 .or. .not.allocated(fit%estimate))then
      result%status=fit_invalid_argument;return
    end if
    call seed_rng(seed);allocate(draw%left(n),draw%right(n))
    allocate(result%estimates(nsim,size(fit%estimate)),result%convergence(nsim))
    result%estimates=0.0_dp;result%convergence=fit_invalid_argument
    do b=1,nsim
      do i=1,n
        call random_number(u);index=1+min(int(u*real(n,dp)),n-1)
        draw%left(i)=sample%left(index);draw%right(i)=sample%right(index)
      end do
      call mledist_censored(draw,dist,fit%estimate,refit,control)
      result%convergence(b)=refit%convergence
      if(refit%convergence==fit_success)result%estimates(b,:)=refit%estimate
    end do
    result%successful=count(result%convergence==fit_success);result%nsim=nsim
    call summarize_bootstrap(result)
  end subroutine bootdistcens

  subroutine turnbull_npmle(sample,result,max_iterations,tolerance)
    type(censored_sample),intent(in)::sample
    type(npmle_result),intent(out)::result
    integer,intent(in),optional::max_iterations
    real(dp),intent(in),optional::tolerance
    real(dp),allocatable::points(:),unique_points(:),prob(:),newprob(:),denom(:)
    logical,allocatable::incidence(:,:)
    real(dp)::tol,diff
    integer::n,m,i,j,k,it,maxit

    n=size(sample%left);maxit=10000;if(present(max_iterations))maxit=max_iterations
    tol=1.0e-10_dp;if(present(tolerance))tol=tolerance
    if(n<1 .or. size(sample%right)/=n .or. any(sample%left>sample%right))then
      result%status=fit_invalid_argument;return
    end if
    allocate(points(2*n));k=0
    do i=1,n
      if(abs(sample%left(i))<0.25_dp*huge(1.0_dp))then;k=k+1;points(k)=sample%left(i);end if
      if(abs(sample%right(i))<0.25_dp*huge(1.0_dp))then;k=k+1;points(k)=sample%right(i);end if
    end do
    if(k==0)then;result%status=fit_invalid_argument;return;end if
    call sort_local(points(:k));allocate(unique_points(k));m=0
    do i=1,k
      if(m==0)then
        m=1;unique_points(m)=points(i)
      else if(abs(points(i)-unique_points(m))>tol*max(1.0_dp,abs(points(i))))then
        m=m+1;unique_points(m)=points(i)
      end if
    end do
    unique_points=unique_points(:m)
    allocate(incidence(n,m),prob(m),newprob(m),denom(n));incidence=.false.
    do i=1,n
      do j=1,m
        incidence(i,j)=point_in_interval(unique_points(j),sample%left(i),sample%right(i),tol)
      end do
      if(.not.any(incidence(i,:)))then;result%status=fit_numerical_error;return;end if
    end do
    prob=1.0_dp/real(m,dp)
    do it=1,maxit
      do i=1,n;denom(i)=sum(prob,mask=incidence(i,:));end do
      if(any(denom<=tiny(1.0_dp)))then;result%status=fit_numerical_error;return;end if
      newprob=0.0_dp
      do j=1,m
        do i=1,n
          if(incidence(i,j))newprob(j)=newprob(j)+prob(j)/denom(i)
        end do
        newprob(j)=newprob(j)/real(n,dp)
      end do
      newprob=newprob/sum(newprob);diff=maxval(abs(newprob-prob));prob=newprob
      if(diff<=tol)exit
    end do
    allocate(result%interval_left(m),result%interval_right(m),result%probability(m))
    result%interval_left=unique_points;result%interval_right=unique_points;result%probability=prob
    result%iterations=min(it,maxit)
    result%status=merge(fit_success,fit_no_convergence,it<=maxit)
  end subroutine turnbull_npmle

  subroutine cdf_bootstrap_band(grid,dist,bootstrap,probabilities,band,status)
    real(dp),intent(in)::grid(:),probabilities(:)
    type(distribution_model),intent(in)::dist
    type(bootstrap_result),intent(in)::bootstrap
    real(dp),allocatable,intent(out)::band(:,:)
    integer,intent(out)::status
    real(dp),allocatable::values(:)
    integer::i,j,k,ns
    if(.not.associated(dist%cdf) .or. size(grid)==0 .or. &
       any(probabilities<0.0_dp) .or. any(probabilities>1.0_dp) .or. &
       .not.allocated(bootstrap%estimates))then
      allocate(band(0,0));status=fit_invalid_argument;return
    end if
    ns=count(bootstrap%convergence==fit_success)
    if(ns<1)then;allocate(band(0,0));status=fit_no_convergence;return;end if
    allocate(band(size(grid),size(probabilities)),values(ns))
    do i=1,size(grid)
      k=0
      do j=1,bootstrap%nsim
        if(bootstrap%convergence(j)==fit_success)then
          k=k+1;values(k)=dist%cdf(grid(i),bootstrap%estimates(j,:))
        end if
      end do
      do j=1,size(probabilities)
        band(i,j)=type7_quantile(values,probabilities(j))
      end do
    end do
    status=fit_success
  end subroutine cdf_bootstrap_band

  subroutine summarize_bootstrap(result)
    type(bootstrap_result),intent(inout)::result
    integer::j,p
    p=size(result%estimates,2)
    allocate(result%mean(p),result%standard_deviation(p))
    if(result%successful<1)then
      result%mean=0.0_dp;result%standard_deviation=0.0_dp;result%status=fit_no_convergence;return
    end if
    do j=1,p
      result%mean(j)=sum(result%estimates(:,j),mask=result%convergence==fit_success)/real(result%successful,dp)
      if(result%successful>1)then
        result%standard_deviation(j)=sqrt(sum((result%estimates(:,j)-result%mean(j))**2, &
          mask=result%convergence==fit_success)/real(result%successful-1,dp))
      else
        result%standard_deviation(j)=0.0_dp
      end if
    end do
    result%status=fit_success
  end subroutine summarize_bootstrap

  integer function method_code(name) result(code)
    character(len=*),intent(in)::name
    select case(trim(name))
    case("mle");code=method_mle
    case("mme");code=method_mme
    case("qme");code=method_qme
    case("mge");code=method_mge
    case("mse");code=method_mse
    case default;code=0
    end select
  end function method_code

  logical function point_in_interval(x,left,right,tol) result(value)
    real(dp),intent(in)::x,left,right,tol
    logical::leftinf,rightinf,exact
    leftinf=abs(left)>=0.25_dp*huge(1.0_dp);rightinf=abs(right)>=0.25_dp*huge(1.0_dp)
    exact=.not.leftinf .and. .not.rightinf .and. abs(left-right)<=tol*max(1.0_dp,abs(left),abs(right))
    if(exact)then
      value=abs(x-left)<=tol*max(1.0_dp,abs(left))
    else if(leftinf)then
      value=x<=right+tol*max(1.0_dp,abs(right))
    else if(rightinf)then
      value=x>=left-tol*max(1.0_dp,abs(left))
    else
      value=x>=left-tol*max(1.0_dp,abs(left)) .and. x<=right+tol*max(1.0_dp,abs(right))
    end if
  end function point_in_interval

  subroutine sort_local(x)
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
  end subroutine sort_local

end module fitdistrplus_bootstrap

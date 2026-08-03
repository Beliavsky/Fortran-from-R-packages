! SPDX-License-Identifier: GPL-3.0-only
module imputefin_ar1_gaussian
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
  use imputefin_kinds, only : dp
  use imputefin_types, only : ar1_options, imputation_options, ar1_fit_result, imputation_result, &
       impute_ok, impute_invalid_input, impute_insufficient_data, impute_singular, impute_not_converged
  use imputefin_rng, only : rng_state, rng_seed
  use imputefin_linalg, only : mvn_sample
  use imputefin_math, only : normal_cdf
  use imputefin_missing, only : missing_blocks, is_inner_na, any_inner_na, collect_indices, trim_observed_range
  implicit none
  private
  public :: fit_ar1_gaussian, impute_ar1_gaussian, impute_rolling_ar1_gaussian
  public :: conditional_gaussian_moments

  interface fit_ar1_gaussian
    module procedure fit_ar1_gaussian_vector
    module procedure fit_ar1_gaussian_matrix
  end interface
  interface impute_ar1_gaussian
    module procedure impute_ar1_gaussian_vector
    module procedure impute_ar1_gaussian_matrix
  end interface
  interface impute_rolling_ar1_gaussian
    module procedure impute_rolling_ar1_gaussian_vector
    module procedure impute_rolling_ar1_gaussian_matrix
  end interface
contains
  pure logical function close_rel(a,b,tol)
    real(dp), intent(in) :: a,b,tol
    close_rel = abs(a-b) <= tol*max(1.0_dp,0.5_dp*(abs(a)+abs(b)))
  end function close_rel

  subroutine set_error(res,status,message)
    type(ar1_fit_result), intent(inout) :: res
    integer, intent(in) :: status
    character(*), intent(in) :: message
    res%status=status; res%message=message; res%converged=.false.
  end subroutine set_error

  subroutine complete_fit(y,opt,res)
    real(dp), intent(in) :: y(:)
    type(ar1_options), intent(in) :: opt
    type(ar1_fit_result), intent(inout) :: res
    real(dp) :: sy2,sy1,sy22,sy11,sy21,den
    integer :: n
    n=size(y)
    sy2=sum(y(2:n)); sy1=sum(y(1:n-1))
    sy22=sum(y(2:n)**2); sy11=sum(y(1:n-1)**2); sy21=sum(y(2:n)*y(1:n-1))
    if(.not.opt%random_walk .and. .not.opt%zero_mean) then
      den=sy11-sy1*sy1/real(n-1,dp)
      if(abs(den)<=epsilon(1.0_dp)) then
        call set_error(res,impute_singular,'constant or singular time series'); return
      end if
      res%phi1=(sy21-sy2*sy1/real(n-1,dp))/den
      res%phi0=(sy2-res%phi1*sy1)/real(n-1,dp)
    else if(opt%random_walk .and. .not.opt%zero_mean) then
      res%phi1=1.0_dp; res%phi0=(sy2-sy1)/real(n-1,dp)
    else if(.not.opt%random_walk .and. opt%zero_mean) then
      if(abs(sy11)<=epsilon(1.0_dp)) then
        call set_error(res,impute_singular,'constant or singular time series'); return
      end if
      res%phi1=sy21/sy11; res%phi0=0.0_dp
    else
      res%phi1=1.0_dp; res%phi0=0.0_dp
    end if
    res%sigma2=(sy22+real(n-1,dp)*res%phi0**2+res%phi1**2*sy11-2.0_dp*res%phi0*sy2-&
         2.0_dp*res%phi1*sy21+2.0_dp*res%phi0*res%phi1*sy1)/real(n-1,dp)
    if(res%sigma2<=tiny(1.0_dp)) then
      call set_error(res,impute_singular,'constant time series'); return
    end if
    res%iterations=1; res%converged=.true.; res%status=impute_ok; res%message='ok'
  end subroutine complete_fit

  subroutine consecutive_pairs(y,x1,x2,status)
    real(dp), intent(in) :: y(:)
    real(dp), allocatable, intent(out) :: x1(:),x2(:)
    integer, intent(out) :: status
    integer :: i,k,m
    m=0
    do i=2,size(y)
      if(.not.ieee_is_nan(y(i-1)) .and. .not.ieee_is_nan(y(i))) m=m+1
    end do
    allocate(x1(m),x2(m)); k=0
    do i=2,size(y)
      if(.not.ieee_is_nan(y(i-1)) .and. .not.ieee_is_nan(y(i))) then
        k=k+1; x1(k)=y(i-1); x2(k)=y(i)
      end if
    end do
    status=merge(0,1,m>=2)
  end subroutine consecutive_pairs

  subroutine heuristic_fit(y,opt,res)
    real(dp), intent(in) :: y(:)
    type(ar1_options), intent(in) :: opt
    type(ar1_fit_result), intent(inout) :: res
    real(dp), allocatable :: x1(:),x2(:)
    integer :: st
    call consecutive_pairs(y,x1,x2,st)
    if(st/=0) then; call set_error(res,impute_insufficient_data,'too few consecutive observed pairs'); return; end if
    call complete_fit_pairs(x1,x2,opt,res)
  end subroutine heuristic_fit

  subroutine complete_fit_pairs(x1,x2,opt,res)
    real(dp), intent(in) :: x1(:),x2(:)
    type(ar1_options), intent(in) :: opt
    type(ar1_fit_result), intent(inout) :: res
    real(dp) :: s1,s2,s11,s22,s12,den
    integer :: m
    m=size(x1); s1=sum(x1); s2=sum(x2); s11=sum(x1*x1); s22=sum(x2*x2); s12=sum(x1*x2)
    if(.not.opt%random_walk .and. .not.opt%zero_mean) then
      den=s11-s1*s1/real(m,dp)
      if(abs(den)<=epsilon(1.0_dp)) then; call set_error(res,impute_singular,'singular observed pairs'); return; end if
      res%phi1=(s12-s2*s1/real(m,dp))/den; res%phi0=(s2-res%phi1*s1)/real(m,dp)
    else if(opt%random_walk .and. .not.opt%zero_mean) then
      res%phi1=1.0_dp; res%phi0=(s2-s1)/real(m,dp)
    else if(.not.opt%random_walk .and. opt%zero_mean) then
      if(abs(s11)<=epsilon(1.0_dp)) then; call set_error(res,impute_singular,'singular observed pairs'); return; end if
      res%phi1=s12/s11; res%phi0=0.0_dp
    else
      res%phi1=1.0_dp; res%phi0=0.0_dp
    end if
    res%sigma2=sum((x2-res%phi0-res%phi1*x1)**2)/real(m,dp)
    if(res%sigma2<=tiny(1.0_dp)) then; call set_error(res,impute_singular,'constant time series'); return; end if
    res%converged=.true.; res%status=impute_ok; res%iterations=1; res%message='ok'
  end subroutine complete_fit_pairs

  subroutine conditional_gaussian_moments(y,phi0,phi1,sigma2,mean_y,cov_y,status)
    real(dp), intent(in) :: y(:),phi0,phi1,sigma2
    real(dp), allocatable, intent(out) :: mean_y(:),cov_y(:,:)
    integer, intent(out) :: status
    type(missing_blocks) :: blocks
    integer :: b,m,i,j,idx,n
    real(dp), allocatable :: covobs(:,:), lastcol(:), diagv(:), muobs(:), blockcov(:,:)
    real(dp) :: sumphi, denom
    n=size(y); allocate(mean_y(n),cov_y(n,n)); mean_y=y; cov_y=0.0_dp; status=0
    call blocks%build(y)
    do b=1,blocks%n_block
      m=blocks%length(b)
      allocate(covobs(m+1,m+1),lastcol(m+1),diagv(m+1),muobs(m+1),blockcov(m,m))
      covobs=0.0_dp; diagv(1)=1.0_dp
      do i=2,m+1; diagv(i)=diagv(i-1)*phi1*phi1+1.0_dp; end do
      do i=1,m+1
        covobs(i,i)=diagv(i)
        do j=i+1,m+1
          covobs(i,j)=diagv(i)*phi1**(j-i); covobs(j,i)=covobs(i,j)
        end do
      end do
      do i=1,m+1
        sumphi=0.0_dp
        do j=0,i-1; sumphi=sumphi+phi1**j; end do
        muobs(i)=sumphi*phi0+phi1**i*y(blocks%first(b)-1)
      end do
      lastcol=covobs(:,m+1); denom=lastcol(m+1)
      if(denom<=tiny(1.0_dp)) then; status=1; return; end if
      do i=1,m
        idx=blocks%first(b)+i-1
        mean_y(idx)=muobs(i)+lastcol(i)/denom*(y(blocks%last(b)+1)-muobs(m+1))
      end do
      blockcov=sigma2*(covobs(1:m,1:m)-spread(lastcol(1:m),2,m)*spread(lastcol(1:m),1,m)/denom)
      cov_y(blocks%first(b):blocks%last(b),blocks%first(b):blocks%last(b))=blockcov
      deallocate(covobs,lastcol,diagv,muobs,blockcov)
    end do
  end subroutine conditional_gaussian_moments

  subroutine find_outliers_gaussian(y,fit,threshold,idx)
    real(dp), intent(in) :: y(:),threshold
    type(ar1_fit_result), intent(in) :: fit
    integer, allocatable, intent(out) :: idx(:)
    real(dp), allocatable :: work(:)
    logical, allocatable :: mark(:)
    integer, allocatable :: obs(:)
    integer :: i,d,j
    real(dp) :: mu,delta,sump
    allocate(work(size(y)),mark(size(y))); work=y; mark=.false.
    call collect_indices(.not.[(ieee_is_nan(work(i)),i=1,size(work))],obs)
    do i=2,size(obs)
      d=obs(i)-obs(i-1); sump=0.0_dp
      do j=0,d-1
        sump=sump+fit%phi1**j
      end do
      mu=sump*fit%phi0+fit%phi1**(obs(i)-obs(i-1))*work(obs(i-1))
      delta=abs(work(obs(i))-mu)
      if(normal_cdf(-delta/sqrt(fit%sigma2))<threshold) then
        mark(obs(i))=.true.; work(obs(i))=mu
      end if
    end do
    call collect_indices(mark,idx)
  end subroutine find_outliers_gaussian

  recursive subroutine fit_ar1_gaussian_vector(y,res,options,return_iterates,return_conditional)
    real(dp), intent(in) :: y(:)
    type(ar1_fit_result), intent(out) :: res
    type(ar1_options), intent(in), optional :: options
    logical, intent(in), optional :: return_iterates,return_conditional
    type(ar1_options) :: opt
    real(dp), allocatable :: work(:),trim(:),mean_y(:),cov_y(:,:),p0(:),p1(:),s2(:)
    logical, allocatable :: missmask(:)
    integer :: first,last,st,k,k_used,n,nobs
    real(dp) :: sy2,sy1,sy22,sy11,sy21,np0,np1,ns2,den
    logical :: keep_iter,keep_cond
    opt=ar1_options(); if(present(options)) opt=options
    keep_iter=.false.; if(present(return_iterates)) keep_iter=return_iterates
    keep_cond=.false.; if(present(return_conditional)) keep_cond=return_conditional
    res=ar1_fit_result()
    if(size(y)<5 .or. count(.not.[(ieee_is_nan(y(k)),k=1,size(y))])<5) then
      call set_error(res,impute_insufficient_data,'at least five observed values are required'); return
    end if
    if(opt%tol<=0.0_dp .or. opt%maxiter<1) then
      call set_error(res,impute_invalid_input,'invalid tolerance or iteration limit'); return
    end if
    allocate(work(size(y))); work=y
    if(opt%remove_outliers) then
      call fit_ar1_gaussian_vector(work,res,ar1_options(random_walk=opt%random_walk,zero_mean=opt%zero_mean,&
           remove_outliers=.false.,tol=opt%tol,maxiter=opt%maxiter))
      if(res%status/=impute_ok) return
      call find_outliers_gaussian(work,res,opt%outlier_prob_th,res%index_outliers)
      if(allocated(res%index_outliers)) work(res%index_outliers)=ieee_value(work(1),ieee_quiet_nan)
    end if
    call trim_observed_range(work,first,last,st)
    if(st/=0) then; call set_error(res,impute_insufficient_data,'no observed values'); return; end if
    trim=work(first:last); missmask=is_inner_na(trim); call collect_indices(missmask,res%index_miss)
    if(allocated(res%index_miss)) res%index_miss=res%index_miss+first-1
    if(.not.any(missmask)) then
      call complete_fit(trim,opt,res)
      return
    end if
    call heuristic_fit(trim,opt,res); if(res%status/=impute_ok) return
    n=size(trim); allocate(p0(opt%maxiter+1),p1(opt%maxiter+1),s2(opt%maxiter+1))
    p0(1)=res%phi0; p1(1)=res%phi1; s2(1)=res%sigma2; res%converged=.false.
    do k=1,opt%maxiter
      call conditional_gaussian_moments(trim,p0(k),p1(k),s2(k),mean_y,cov_y,st)
      if(st/=0) then; call set_error(res,impute_singular,'conditional covariance failure'); return; end if
      sy2=sum(mean_y(2:n)); sy1=sum(mean_y(1:n-1))
      sy22=sum(mean_y(2:n)**2+[(cov_y(k,k),k=2,n)])
      sy11=sum(mean_y(1:n-1)**2+[(cov_y(k,k),k=1,n-1)])
      sy21=sum(mean_y(2:n)*mean_y(1:n-1)+[(cov_y(k,k+1),k=1,n-1)])
      nobs=n-1
      if(.not.opt%random_walk .and. .not.opt%zero_mean) then
        den=sy11-sy1*sy1/real(nobs,dp)
        if(abs(den)<=epsilon(1.0_dp)) then; call set_error(res,impute_singular,'singular EM update'); return; end if
        np1=(sy21-sy2*sy1/real(nobs,dp))/den; np0=(sy2-np1*sy1)/real(nobs,dp)
      else if(opt%random_walk .and. .not.opt%zero_mean) then
        np1=1.0_dp; np0=(sy2-sy1)/real(nobs,dp)
      else if(.not.opt%random_walk .and. opt%zero_mean) then
        np1=sy21/sy11; np0=0.0_dp
      else
        np1=1.0_dp; np0=0.0_dp
      end if
      ns2=(sy22+real(nobs,dp)*np0*np0+np1*np1*sy11-2.0_dp*np0*sy2-2.0_dp*np1*sy21+&
           2.0_dp*np0*np1*sy1)/real(nobs,dp)
      ns2=max(ns2,100.0_dp*tiny(1.0_dp))
      p0(k+1)=np0; p1(k+1)=np1; s2(k+1)=ns2
      if(close_rel(p0(k+1),p0(k),opt%tol).and.close_rel(p1(k+1),p1(k),opt%tol).and.&
           close_rel(s2(k+1),s2(k),opt%tol)) then; res%converged=.true.; exit; end if
      if(allocated(mean_y)) deallocate(mean_y,cov_y)
    end do
    k_used=min(k,opt%maxiter)
    res%iterations=k_used
    res%phi0=p0(k_used+1)
    res%phi1=p1(k_used+1)
    res%sigma2=s2(k_used+1)
    res%status=merge(impute_ok,impute_not_converged,res%converged)
    if (res%converged) then
      res%message='ok'
    else
      res%message='maximum iterations reached'
    end if
    if(keep_iter) then
      res%phi0_iterates=p0(1:k_used+1); res%phi1_iterates=p1(1:k_used+1); res%sigma2_iterates=s2(1:k_used+1)
    end if
    if(keep_cond) then
      call conditional_gaussian_moments(trim,res%phi0,res%phi1,res%sigma2,res%cond_mean,res%cond_cov,st)
    end if
  end subroutine fit_ar1_gaussian_vector

  subroutine fit_ar1_gaussian_matrix(y,res,options,return_iterates,return_conditional)
    real(dp), intent(in) :: y(:,:)
    type(ar1_fit_result), allocatable, intent(out) :: res(:)
    type(ar1_options), intent(in), optional :: options
    logical, intent(in), optional :: return_iterates,return_conditional
    integer :: j
    allocate(res(size(y,2)))
    do j=1,size(y,2)
      call fit_ar1_gaussian_vector(y(:,j),res(j),options,return_iterates,return_conditional)
    end do
  end subroutine fit_ar1_gaussian_matrix

  subroutine sample_gaussian_missing(y,fit,state,out,status)
    real(dp), intent(in) :: y(:)
    type(ar1_fit_result), intent(in) :: fit
    type(rng_state), intent(inout) :: state
    real(dp), intent(out) :: out(:)
    integer, intent(out) :: status
    real(dp), allocatable :: mean_y(:),cov_y(:,:),mu(:),cov(:,:),draw(:)
    integer, allocatable :: miss(:)
    logical, allocatable :: mask(:)
    integer :: st
    out=y; mask=is_inner_na(y); call collect_indices(mask,miss); status=0
    if(size(miss)==0) return
    call conditional_gaussian_moments(y,fit%phi0,fit%phi1,fit%sigma2,mean_y,cov_y,st)
    if(st/=0) then; status=1; return; end if
    mu=mean_y(miss); cov=cov_y(miss,miss); allocate(draw(size(miss)))
    call mvn_sample(mu,cov,state,draw,st)
    if(st/=0) then; status=1; return; end if
    out(miss)=draw
  end subroutine sample_gaussian_missing

  subroutine impute_ar1_gaussian_vector(y,result,options,impute_options_in)
    real(dp), intent(in) :: y(:)
    type(imputation_result), intent(out) :: result
    type(ar1_options), intent(in), optional :: options
    type(imputation_options), intent(in), optional :: impute_options_in
    type(ar1_options) :: opt
    type(imputation_options) :: iopt
    type(rng_state) :: state
    real(dp), allocatable :: work(:)
    integer :: s,st
    opt=ar1_options(); if(present(options)) opt=options
    iopt=imputation_options(); if(present(impute_options_in)) iopt=impute_options_in
    allocate(result%fits(1),result%values(size(y),1,iopt%n_samples)); work=y
    call fit_ar1_gaussian_vector(work,result%fits(1),opt,return_conditional=.true.)
    if(result%fits(1)%status/=impute_ok .and. result%fits(1)%status/=impute_not_converged) then
      result%status=result%fits(1)%status; result%message=result%fits(1)%message; return
    end if
    if(opt%remove_outliers .and. allocated(result%fits(1)%index_outliers)) then
      work(result%fits(1)%index_outliers)=ieee_value(work(1),ieee_quiet_nan)
    end if
    call rng_seed(state,iopt%seed)
    do s=1,iopt%n_samples
      call sample_gaussian_missing(work,result%fits(1),state,result%values(:,1,s),st)
      if(st/=0) then; result%status=impute_singular; result%message='sampling failure'; return; end if
    end do
    result%status=impute_ok; result%message='ok'
  end subroutine impute_ar1_gaussian_vector

  subroutine impute_ar1_gaussian_matrix(y,result,options,impute_options_in)
    real(dp), intent(in) :: y(:,:)
    type(imputation_result), intent(out) :: result
    type(ar1_options), intent(in), optional :: options
    type(imputation_options), intent(in), optional :: impute_options_in
    type(imputation_options) :: iopt
    type(imputation_result), allocatable :: one
    integer :: j,s
    iopt=imputation_options(); if(present(impute_options_in)) iopt=impute_options_in
    allocate(one)
    allocate(result%fits(size(y,2)),result%values(size(y,1),size(y,2),iopt%n_samples))
    do j=1,size(y,2)
      call impute_ar1_gaussian_vector(y(:,j),one,options,iopt)
      result%fits(j)=one%fits(1)
      do s=1,iopt%n_samples; result%values(:,j,s)=one%values(:,1,s); end do
      if(one%status/=impute_ok) then; result%status=one%status; result%message=one%message; return; end if
      iopt%seed=iopt%seed+104729_8
    end do
    result%status=impute_ok; result%message='ok'
  end subroutine impute_ar1_gaussian_matrix


  subroutine fill_linear_inner(y,out)
    real(dp), intent(in) :: y(:)
    real(dp), intent(out) :: out(:)
    integer :: i,left,right
    real(dp) :: frac
    out=y
    do i=1,size(y)
      if(.not.ieee_is_nan(y(i))) cycle
      left=i-1
      do while(left>=1)
        if(.not.ieee_is_nan(y(left))) exit
        left=left-1
      end do
      right=i+1
      do while(right<=size(y))
        if(.not.ieee_is_nan(y(right))) exit
        right=right+1
      end do
      if(left>=1 .and. right<=size(y)) then
        frac=real(i-left,dp)/real(right-left,dp)
        out(i)=(1.0_dp-frac)*y(left)+frac*y(right)
      end if
    end do
  end subroutine fill_linear_inner

  subroutine impute_rolling_ar1_gaussian_vector(y,y_imputed,options,rolling_window,seed,status)
    real(dp), intent(in) :: y(:)
    real(dp), allocatable, intent(out) :: y_imputed(:)
    type(ar1_options), intent(in), optional :: options
    integer, intent(in), optional :: rolling_window
    integer(kind=8), intent(in), optional :: seed
    integer, intent(out), optional :: status
    integer :: first,last,st,w,start,finish
    type(imputation_result), allocatable :: tmp
    type(imputation_options) :: io
    allocate(tmp)
    y_imputed=y; w=252; if(present(rolling_window)) w=rolling_window
    io%n_samples=1; if(present(seed)) io%seed=seed
    call trim_observed_range(y,first,last,st)
    if(st/=0) then; if(present(status)) status=impute_insufficient_data; return; end if
    start=first
    do while(start<=last)
      finish=min(last,start+w-1)
      do while(finish<last .and. ieee_is_nan(y(finish))); finish=finish+1; end do
      if(.not.ieee_is_nan(y(start)) .and. .not.ieee_is_nan(y(finish)) .and. finish-start+1>=5) then
        call impute_ar1_gaussian_vector(y(start:finish),tmp,options,io)
        if(tmp%status==impute_ok) then
          y_imputed(start:finish)=tmp%values(:,1,1)
        else
          call fill_linear_inner(y(start:finish),y_imputed(start:finish))
        end if
      end if
      start=finish+1; io%seed=io%seed+65537_8
    end do
    if(present(status)) status=impute_ok
  end subroutine impute_rolling_ar1_gaussian_vector

  subroutine impute_rolling_ar1_gaussian_matrix(y,y_imputed,options,rolling_window,seed,status)
    real(dp), intent(in) :: y(:,:)
    real(dp), allocatable, intent(out) :: y_imputed(:,:)
    type(ar1_options), intent(in), optional :: options
    integer, intent(in), optional :: rolling_window
    integer(kind=8), intent(in), optional :: seed
    integer, intent(out), optional :: status
    integer :: j,st
    integer(kind=8) :: sd
    real(dp), allocatable :: col(:)
    allocate(y_imputed(size(y,1),size(y,2))); sd=5489_8; if(present(seed)) sd=seed
    do j=1,size(y,2)
      call impute_rolling_ar1_gaussian_vector(y(:,j),col,options,rolling_window,sd,st)
      y_imputed(:,j)=col; if(st/=impute_ok) then; if(present(status)) status=st; return; end if
      sd=sd+104729_8
    end do
    if(present(status)) status=impute_ok
  end subroutine impute_rolling_ar1_gaussian_matrix
end module imputefin_ar1_gaussian

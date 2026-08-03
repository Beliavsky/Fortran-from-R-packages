! SPDX-License-Identifier: GPL-3.0-only
module imputefin_ar1_t
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
  use imputefin_kinds, only : dp
  use imputefin_types, only : ar1_options, imputation_options, ar1_fit_result, imputation_result, &
       impute_ok, impute_invalid_input, impute_insufficient_data, impute_singular, impute_not_converged
  use imputefin_rng, only : rng_state, rng_seed, rng_gamma
  use imputefin_linalg, only : mvn_sample
  use imputefin_math, only : student_t_cdf, log_gamma_dp, digamma_dp
  use imputefin_missing, only : missing_blocks, is_inner_na, collect_indices, trim_observed_range
  use imputefin_ar1_gaussian, only : fit_ar1_gaussian, conditional_gaussian_moments
  implicit none
  private
  public :: fit_ar1_t, impute_ar1_t

  interface fit_ar1_t
    module procedure fit_ar1_t_vector
    module procedure fit_ar1_t_matrix
  end interface
  interface impute_ar1_t
    module procedure impute_ar1_t_vector
    module procedure impute_ar1_t_matrix
  end interface
contains
  pure logical function close_rel(a,b,tol)
    real(dp), intent(in) :: a,b,tol
    close_rel=abs(a-b)<=tol*max(1.0_dp,0.5_dp*(abs(a)+abs(b)))
  end function close_rel

  subroutine set_error(res,status,message)
    type(ar1_fit_result), intent(inout) :: res
    integer, intent(in) :: status
    character(*), intent(in) :: message
    res%status=status; res%message=message; res%converged=.false.
  end subroutine set_error

  subroutine consecutive_pairs(y,x1,x2,status)
    real(dp), intent(in) :: y(:)
    real(dp), allocatable, intent(out) :: x1(:),x2(:)
    integer, intent(out) :: status
    integer :: i,k,m
    m=0
    do i=2,size(y)
      if(.not.ieee_is_nan(y(i-1)).and..not.ieee_is_nan(y(i))) m=m+1
    end do
    allocate(x1(m),x2(m)); k=0
    do i=2,size(y)
      if(.not.ieee_is_nan(y(i-1)).and..not.ieee_is_nan(y(i))) then
        k=k+1; x1(k)=y(i-1); x2(k)=y(i)
      end if
    end do
    status=merge(0,1,m>=2)
  end subroutine consecutive_pairs

  function nu_objective(nu,q,dim) result(v)
    real(dp), intent(in) :: nu,q(:)
    integer, intent(in) :: dim
    real(dp) :: v
    if(nu<=0.05_dp) then; v=huge(1.0_dp); return; end if
    v=sum(0.5_dp*(nu+real(dim,dp))*log(nu+q)-log_gamma_dp(0.5_dp*(nu+real(dim,dp)))+&
         log_gamma_dp(0.5_dp*nu)-0.5_dp*nu*log(nu))
  end function nu_objective

  subroutine optimize_nu(q,dim,nu)
    real(dp), intent(in) :: q(:)
    integer, intent(in) :: dim
    real(dp), intent(out) :: nu
    real(dp) :: a,b,x1,x2,f1,f2,gr
    integer :: it
    a=0.2_dp; b=200.0_dp; gr=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
    x1=b-gr*(b-a); x2=a+gr*(b-a); f1=nu_objective(x1,q,dim); f2=nu_objective(x2,q,dim)
    do it=1,150
      if(abs(b-a)<=1.0e-7_dp*(1.0_dp+abs(x1)+abs(x2))) exit
      if(f1>f2) then
        a=x1; x1=x2; f1=f2; x2=a+gr*(b-a); f2=nu_objective(x2,q,dim)
      else
        b=x2; x2=x1; f2=f1; x1=b-gr*(b-a); f1=nu_objective(x1,q,dim)
      end if
    end do
    nu=merge(x1,x2,f1<=f2)
  end subroutine optimize_nu

  subroutine fit_pairs_t(x1,x2,opt,res,keep_iter)
    real(dp), intent(in) :: x1(:),x2(:)
    type(ar1_options), intent(in) :: opt
    type(ar1_fit_result), intent(inout) :: res
    logical, intent(in) :: keep_iter
    real(dp), allocatable :: w(:),q(:),p0(:),p1(:),s2(:),nuv(:)
    real(dp) :: sw,sw1,sw2,sw11,sw12,np0,np1,ns2,den,nnu
    integer :: k,k_used,m
    type(ar1_fit_result) :: gres
    m=size(x1)
    if(m<2) then; call set_error(res,impute_insufficient_data,'too few consecutive pairs'); return; end if
    call fit_pairs_gaussian(x1,x2,opt,gres)
    if(gres%status/=impute_ok) then; res=gres; return; end if
    allocate(w(m),q(m),p0(opt%maxiter+1),p1(opt%maxiter+1),s2(opt%maxiter+1),nuv(opt%maxiter+1))
    p0(1)=gres%phi0; p1(1)=gres%phi1; s2(1)=gres%sigma2; nuv(1)=3.0_dp
    res%converged=.false.
    do k=1,opt%maxiter
      q=(x2-p0(k)-p1(k)*x1)**2/max(s2(k),tiny(1.0_dp))
      w=(nuv(k)+1.0_dp)/(nuv(k)+q)
      sw=sum(w); sw1=sum(w*x1); sw2=sum(w*x2); sw11=sum(w*x1*x1); sw12=sum(w*x1*x2)
      if(.not.opt%random_walk.and..not.opt%zero_mean) then
        den=sw*sw11-sw1*sw1
        if(abs(den)<=epsilon(1.0_dp)) then; call set_error(res,impute_singular,'singular weighted update'); return; end if
        np1=(sw*sw12-sw2*sw1)/den; np0=(sw2-np1*sw1)/sw
      else if(opt%random_walk.and..not.opt%zero_mean) then
        np1=1.0_dp; np0=(sw2-sw1)/sw
      else if(.not.opt%random_walk.and.opt%zero_mean) then
        np1=sw12/sw11; np0=0.0_dp
      else
        np1=1.0_dp; np0=0.0_dp
      end if
      ns2=sum(w*(x2-np0-np1*x1)**2)/real(m,dp); ns2=max(ns2,100.0_dp*tiny(1.0_dp))
      q=(x2-np0-np1*x1)**2/ns2; call optimize_nu(q,1,nnu)
      p0(k+1)=np0; p1(k+1)=np1; s2(k+1)=ns2; nuv(k+1)=nnu
      if(close_rel(p0(k+1),p0(k),opt%tol).and.close_rel(p1(k+1),p1(k),opt%tol).and.&
           close_rel(s2(k+1),s2(k),opt%tol).and.close_rel(nuv(k+1),nuv(k),sqrt(opt%tol))) then
        res%converged=.true.; exit
      end if
    end do
    k_used=min(k,opt%maxiter)
    res%phi0=p0(k_used+1)
    res%phi1=p1(k_used+1)
    res%sigma2=s2(k_used+1)
    res%nu=nuv(k_used+1)
    res%iterations=k_used
    res%status=merge(impute_ok,impute_not_converged,res%converged)
    if (res%converged) then
      res%message='ok'
    else
      res%message='maximum iterations reached'
    end if
    if(keep_iter) then
      res%phi0_iterates=p0(1:k_used+1); res%phi1_iterates=p1(1:k_used+1)
      res%sigma2_iterates=s2(1:k_used+1); res%nu_iterates=nuv(1:k_used+1)
    end if
  end subroutine fit_pairs_t

  subroutine fit_pairs_gaussian(x1,x2,opt,res)
    real(dp), intent(in) :: x1(:),x2(:)
    type(ar1_options), intent(in) :: opt
    type(ar1_fit_result), intent(out) :: res
    real(dp) :: s1,s2,s11,s22,s12,den
    integer :: m
    res=ar1_fit_result(); m=size(x1)
    s1=sum(x1);s2=sum(x2);s11=sum(x1*x1);s22=sum(x2*x2);s12=sum(x1*x2)
    if(.not.opt%random_walk.and..not.opt%zero_mean) then
      den=s11-s1*s1/real(m,dp)
      if(abs(den)<=epsilon(1.0_dp)) then
        call set_error(res,impute_singular,'singular pairs')
        return
      end if
      res%phi1=(s12-s2*s1/real(m,dp))/den;res%phi0=(s2-res%phi1*s1)/real(m,dp)
    else if(opt%random_walk.and..not.opt%zero_mean) then
      res%phi1=1.0_dp;res%phi0=(s2-s1)/real(m,dp)
    else if(.not.opt%random_walk.and.opt%zero_mean) then
      res%phi1=s12/s11;res%phi0=0.0_dp
    else
      res%phi1=1.0_dp;res%phi0=0.0_dp
    end if
    res%sigma2=sum((x2-res%phi0-res%phi1*x1)**2)/real(m,dp)
    if(res%sigma2<=tiny(1.0_dp)) then;call set_error(res,impute_singular,'constant series');return;end if
    res%status=impute_ok;res%converged=.true.;res%message='ok'
  end subroutine fit_pairs_gaussian

  subroutine find_outliers_t(y,fit,threshold,idx)
    real(dp), intent(in) :: y(:),threshold
    type(ar1_fit_result), intent(in) :: fit
    integer, allocatable, intent(out) :: idx(:)
    real(dp), allocatable :: work(:)
    logical, allocatable :: mark(:),obs_mask(:)
    integer, allocatable :: obs(:)
    integer :: i,j,gap
    real(dp) :: mu,delta,sump
    allocate(work(size(y)),mark(size(y)),obs_mask(size(y)));work=y;mark=.false.
    do i=1,size(y);obs_mask(i)=.not.ieee_is_nan(work(i));end do
    call collect_indices(obs_mask,obs)
    do i=2,size(obs)
      gap=obs(i)-obs(i-1);sump=0.0_dp
      do j=0,gap-1;sump=sump+fit%phi1**j;end do
      mu=sump*fit%phi0+fit%phi1**gap*work(obs(i-1));delta=abs(work(obs(i))-mu)
      if(student_t_cdf(-delta/sqrt(fit%sigma2),max(2.0_dp,anint(fit%nu)))<threshold) then
        mark(obs(i))=.true.;work(obs(i))=mu
      end if
    end do
    call collect_indices(mark,idx)
  end subroutine find_outliers_t

  subroutine gibbs_ar1_step(y_template,y_current,fit,state,y_new,tau,status)
    real(dp), intent(in) :: y_template(:),y_current(:)
    type(ar1_fit_result), intent(in) :: fit
    type(rng_state), intent(inout) :: state
    real(dp), intent(out) :: y_new(:),tau(:)
    integer, intent(out) :: status
    type(missing_blocks) :: blocks
    real(dp), allocatable :: cov(:,:),muobs(:),lastcol(:),draw(:),blockcov(:,:),condmean(:)
    real(dp) :: residual,q,var,sumphi,denom
    integer :: i,j,b,m,t,st
    y_new=y_current; tau=1.0_dp; status=0
    do i=2,size(y_new)
      residual=y_new(i)-fit%phi0-fit%phi1*y_new(i-1)
      q=residual*residual/max(fit%sigma2,tiny(1.0_dp))
      tau(i)=rng_gamma(state,0.5_dp*(fit%nu+1.0_dp),0.5_dp*(fit%nu+q))
      tau(i)=max(tau(i),1.0e-12_dp)
    end do
    call blocks%build(y_template)
    do b=1,blocks%n_block
      m=blocks%length(b)
      allocate(cov(m+1,m+1),muobs(m+1),lastcol(m+1),draw(m),blockcov(m,m),condmean(m))
      cov=0.0_dp
      do i=1,m+1
        t=blocks%first(b)+i-1;var=fit%sigma2/tau(t)
        if(i==1) then;cov(i,i)=var;else;cov(i,i)=fit%phi1**2*cov(i-1,i-1)+var;end if
        do j=1,i-1;cov(j,i)=cov(j,j)*fit%phi1**(i-j);cov(i,j)=cov(j,i);end do
        sumphi=0.0_dp;do j=0,i-1;sumphi=sumphi+fit%phi1**j;end do
        muobs(i)=sumphi*fit%phi0+fit%phi1**i*y_new(blocks%first(b)-1)
      end do
      lastcol=cov(:,m+1);denom=cov(m+1,m+1)
      condmean=muobs(1:m)+lastcol(1:m)/denom*(y_template(blocks%last(b)+1)-muobs(m+1))
      blockcov=cov(1:m,1:m)-spread(lastcol(1:m),2,m)*spread(lastcol(1:m),1,m)/denom
      call mvn_sample(condmean,blockcov,state,draw,st)
      if(st/=0) then;status=1;return;end if
      y_new(blocks%first(b):blocks%last(b))=draw
      deallocate(cov,muobs,lastcol,draw,blockcov,condmean)
    end do
  end subroutine gibbs_ar1_step

  recursive subroutine fit_ar1_t_vector(y,res,options,return_iterates,return_conditional)
    real(dp), intent(in) :: y(:)
    type(ar1_fit_result), intent(out) :: res
    type(ar1_options), intent(in), optional :: options
    logical, intent(in), optional :: return_iterates,return_conditional
    type(ar1_options) :: opt,noout
    logical :: keep_iter,keep_cond
    real(dp), allocatable :: work(:),trim(:),x1(:),x2(:),chains(:,:),next(:),tau(:)
    real(dp), allocatable :: p0(:),p1(:),s2(:),nuv(:),savg(:),snew(:),q(:)
    type(ar1_fit_result) :: gres
    type(rng_state) :: state
    integer :: first,last,st,j,k,k_used,thin,n
    real(dp) :: gamma,sw,sw1,sw2,sw11,sw12,np0,np1,ns2,den,nnu
    opt=ar1_options();if(present(options))opt=options
    keep_iter=.false.;if(present(return_iterates))keep_iter=return_iterates
    keep_cond=.false.;if(present(return_conditional))keep_cond=return_conditional
    res=ar1_fit_result()
    if(size(y)<5.or.count(.not.[(ieee_is_nan(y(k)),k=1,size(y))])<5) then
      call set_error(res,impute_insufficient_data,'at least five observed values are required');return
    end if
    work=y
    if(opt%remove_outliers) then
      noout=opt;noout%remove_outliers=.false.
      call fit_ar1_t_vector(work,res,noout)
      if(res%status/=impute_ok.and.res%status/=impute_not_converged)return
      call find_outliers_t(work,res,opt%outlier_prob_th,res%index_outliers)
      if(allocated(res%index_outliers))work(res%index_outliers)=ieee_value(work(1),ieee_quiet_nan)
    end if
    call trim_observed_range(work,first,last,st)
    if(st/=0)then
      call set_error(res,impute_insufficient_data,'no observed values')
      return
    end if
    trim=work(first:last);call collect_indices(is_inner_na(trim),res%index_miss)
    if(allocated(res%index_miss))res%index_miss=res%index_miss+first-1
    call consecutive_pairs(trim,x1,x2,st)
    if(st/=0)then
      call set_error(res,impute_insufficient_data,'too few observed pairs')
      return
    end if
    if(.not.any(is_inner_na(trim)).or.opt%fast_and_heuristic) then
      call fit_pairs_t(x1,x2,opt,res,keep_iter)
      if(keep_cond.and.any(is_inner_na(trim))) then
        call fit_ar1_gaussian(trim,gres,ar1_options(random_walk=opt%random_walk,zero_mean=opt%zero_mean),return_conditional=.true.)
        if(allocated(gres%cond_mean))res%cond_mean=gres%cond_mean
      end if
      return
    end if
    call fit_ar1_gaussian(trim,gres,ar1_options(random_walk=opt%random_walk,zero_mean=opt%zero_mean),return_conditional=.true.)
    if(gres%status/=impute_ok.and.gres%status/=impute_not_converged)then;res=gres;return;end if
    n=size(trim);allocate(chains(n,opt%n_chain),next(n),tau(n),p0(opt%maxiter+1),p1(opt%maxiter+1),&
         s2(opt%maxiter+1),nuv(opt%maxiter+1),savg(7),snew(7),q(n-1))
    do j=1,opt%n_chain;chains(:,j)=gres%cond_mean;end do
    p0(1)=gres%phi0;p1(1)=gres%phi1;s2(1)=gres%sigma2;nuv(1)=3.0_dp;savg=0.0_dp
    call rng_seed(state,7919_8);res%converged=.false.
    do k=1,opt%maxiter
      snew=0.0_dp
      do j=1,opt%n_chain
        do thin=1,opt%n_thin
          call gibbs_ar1_step(trim,chains(:,j),ar1_fit_result(phi0=p0(k),phi1=p1(k),sigma2=s2(k),nu=nuv(k)),&
               state,next,tau,st)
          if(st/=0)then;call set_error(res,impute_singular,'Student-t Gibbs step failed');return;end if
          chains(:,j)=next
        end do
        snew(1)=snew(1)+sum(log(tau(2:n))-tau(2:n))
        snew(2)=snew(2)+sum(tau(2:n)*next(2:n)**2)
        snew(3)=snew(3)+sum(tau(2:n))
        snew(4)=snew(4)+sum(tau(2:n)*next(1:n-1)**2)
        snew(5)=snew(5)+sum(tau(2:n)*next(2:n))
        snew(6)=snew(6)+sum(tau(2:n)*next(2:n)*next(1:n-1))
        snew(7)=snew(7)+sum(tau(2:n)*next(1:n-1))
      end do
      snew=snew/real(opt%n_chain,dp)
      if(k<=opt%saem_burn)then
        gamma=1.0_dp
      else
        gamma=1.0_dp/real(k-opt%saem_burn,dp)
      end if
      savg=savg+gamma*(snew-savg);sw=savg(3);sw2=savg(5);sw1=savg(7);sw11=savg(4);sw12=savg(6)
      if(.not.opt%random_walk.and..not.opt%zero_mean)then
        den=sw*sw11-sw1*sw1
        if(abs(den)<=epsilon(1.0_dp))then
          call set_error(res,impute_singular,'singular SAEM update')
          return
        end if
        np1=(sw*sw12-sw2*sw1)/den;np0=(sw2-np1*sw1)/sw
      else if(opt%random_walk.and..not.opt%zero_mean)then
        np1=1.0_dp;np0=(sw2-sw1)/sw
      else if(.not.opt%random_walk.and.opt%zero_mean)then
        np1=sw12/sw11;np0=0.0_dp
      else
        np1=1.0_dp;np0=0.0_dp
      end if
      ns2=(savg(2)+np0*np0*savg(3)+np1*np1*savg(4)-2*np0*savg(5)-2*np1*savg(6)+2*np0*np1*savg(7))/real(n-1,dp)
      ns2=max(ns2,100.0_dp*tiny(1.0_dp))
      q=(chains(2:n,1)-np0-np1*chains(1:n-1,1))**2/ns2;call optimize_nu(q,1,nnu)
      p0(k+1)=np0;p1(k+1)=np1;s2(k+1)=ns2;nuv(k+1)=nnu
      if(k>opt%saem_burn.and.close_rel(np0,p0(k),opt%tol).and.close_rel(np1,p1(k),opt%tol).and.&
           close_rel(ns2,s2(k),opt%tol).and.close_rel(nnu,nuv(k),sqrt(opt%tol)))then;res%converged=.true.;exit;end if
    end do
    k_used=min(k,opt%maxiter)
    res%phi0=p0(k_used+1)
    res%phi1=p1(k_used+1)
    res%sigma2=s2(k_used+1)
    res%nu=nuv(k_used+1)
    res%iterations=k_used
    res%status=merge(impute_ok,impute_not_converged,res%converged)
    if (res%converged) then
      res%message='ok'
    else
      res%message='maximum iterations reached'
    end if
    if(keep_iter)then
      res%phi0_iterates=p0(1:k_used+1)
      res%phi1_iterates=p1(1:k_used+1)
      res%sigma2_iterates=s2(1:k_used+1)
      res%nu_iterates=nuv(1:k_used+1)
    end if
    if(keep_cond)res%cond_mean=gres%cond_mean
  end subroutine fit_ar1_t_vector

  subroutine fit_ar1_t_matrix(y,res,options,return_iterates,return_conditional)
    real(dp),intent(in)::y(:,:)
    type(ar1_fit_result),allocatable,intent(out)::res(:)
    type(ar1_options),intent(in),optional::options
    logical,intent(in),optional::return_iterates,return_conditional
    integer::j
    allocate(res(size(y,2)))
    do j=1,size(y,2);call fit_ar1_t_vector(y(:,j),res(j),options,return_iterates,return_conditional);end do
  end subroutine fit_ar1_t_matrix

  subroutine impute_ar1_t_vector(y,result,options,impute_options_in)
    real(dp),intent(in)::y(:)
    type(imputation_result),intent(out)::result
    type(ar1_options),intent(in),optional::options
    type(imputation_options),intent(in),optional::impute_options_in
    type(ar1_options)::opt
    type(imputation_options)::io
    type(rng_state)::state
    real(dp),allocatable::work(:),current(:),next(:),tau(:),mean_y(:),cov_y(:,:)
    integer::first,last,st,it,s
    opt=ar1_options();if(present(options))opt=options
    io=imputation_options();if(present(impute_options_in))io=impute_options_in
    allocate(result%fits(1),result%values(size(y),1,io%n_samples));work=y
    call fit_ar1_t_vector(work,result%fits(1),opt,return_conditional=.true.)
    if(result%fits(1)%status/=impute_ok .and. &
       result%fits(1)%status/=impute_not_converged)then
      result%status=result%fits(1)%status
      result%message=result%fits(1)%message
      return
    end if
    if(opt%remove_outliers .and. allocated(result%fits(1)%index_outliers)) then
      work(result%fits(1)%index_outliers)=ieee_value(work(1),ieee_quiet_nan)
    end if
    call trim_observed_range(work,first,last,st);current=work(first:last)
    if(any(is_inner_na(current)))then
      call conditional_gaussian_moments(current,result%fits(1)%phi0,result%fits(1)%phi1,result%fits(1)%sigma2,mean_y,cov_y,st)
      where(ieee_is_nan(current))current=mean_y
    end if
    allocate(next(size(current)),tau(size(current)));call rng_seed(state,io%seed)
    do it=1,io%n_burn
      call gibbs_ar1_step(work(first:last),current,result%fits(1),state,next,tau,st)
      if(st/=0)then
        result%status=impute_singular
        result%message='Gibbs imputation failed'
        return
      end if
      current=next
    end do
    do s=1,io%n_samples
      do it=1,io%n_thin;call gibbs_ar1_step(work(first:last),current,result%fits(1),state,next,tau,st);current=next;end do
      result%values(:,1,s)=work;result%values(first:last,1,s)=current
    end do
    result%status=impute_ok;result%message='ok'
  end subroutine impute_ar1_t_vector

  subroutine impute_ar1_t_matrix(y,result,options,impute_options_in)
    real(dp),intent(in)::y(:,:)
    type(imputation_result),intent(out)::result
    type(ar1_options),intent(in),optional::options
    type(imputation_options),intent(in),optional::impute_options_in
    type(imputation_options)::io
    type(imputation_result), allocatable :: one
    integer::j,s
    io=imputation_options();if(present(impute_options_in))io=impute_options_in
    allocate(one)
    allocate(result%fits(size(y,2)),result%values(size(y,1),size(y,2),io%n_samples))
    do j=1,size(y,2)
      call impute_ar1_t_vector(y(:,j),one,options,io);result%fits(j)=one%fits(1)
      do s=1,io%n_samples;result%values(:,j,s)=one%values(:,1,s);end do
      if(one%status/=impute_ok)then;result%status=one%status;result%message=one%message;return;end if
      io%seed=io%seed+104729_8
    end do
    result%status=impute_ok;result%message='ok'
  end subroutine impute_ar1_t_matrix
end module imputefin_ar1_t

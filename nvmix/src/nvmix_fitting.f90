! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2018 Marius Hofert, Erik Hintz and Christiane Lemieux
module nvmix_fitting
  use nvmix_kinds, only : dp
  use nvmix_types
  use nvmix_core
  use nvmix_distributions, only : fit_norm,dstudent_copula,dgrouped_student_copula
  use nvmix_linalg, only : sample_mean_covariance,quadratic_form_spd,covariance_to_correlation
  use nvmix_special, only : student_quantile,chi_square_quantile
  use nvmix_gamma_mix, only : qgammamix
  implicit none
  private
  public :: fit_student,fit_nvmix,fit_student_copula,fit_grouped_student_copula,qqplot_maha
contains
  function fit_student(x,estimate_df,df_fixed,max_iterations) result(fit)
    real(dp), intent(in) :: x(:,:)
    logical, intent(in), optional :: estimate_df
    real(dp), intent(in), optional :: df_fixed
    integer, intent(in), optional :: max_iterations
    type(fit_result) :: fit
    logical :: est,ok
    real(dp) :: a,b,c,d,fc,fd,bestdf
    integer :: i,maxit
    est=.true.; if(present(estimate_df))est=estimate_df
    maxit=80; if(present(max_iterations))maxit=max_iterations
    if(size(x,1)<3 .or. size(x,2)<1)then; fit%ok=.false.; fit%message='insufficient data'; return; end if
    if(.not.est)then
      bestdf=10.0_dp; if(present(df_fixed))bestdf=df_fixed
      call student_profile(x,bestdf,maxit,fit%loc,fit%scale,fit%log_likelihood,ok)
      fit%mixing_parameter=[bestdf]; fit%converged=ok; fit%iterations=maxit
    else
      a=1.01_dp; b=200.0_dp; c=b-(b-a)/1.618033988749895_dp; d=a+(b-a)/1.618033988749895_dp
      fc=negative_profile(x,c,maxit); fd=negative_profile(x,d,maxit)
      do i=1,60
        if(fc<fd)then; b=d; d=c; fd=fc; c=b-(b-a)/1.618033988749895_dp; fc=negative_profile(x,c,maxit)
        else; a=c; c=d; fc=fd; d=a+(b-a)/1.618033988749895_dp; fd=negative_profile(x,d,maxit); end if
      end do
      bestdf=0.5_dp*(a+b)
      call student_profile(x,bestdf,maxit,fit%loc,fit%scale,fit%log_likelihood,ok)
      fit%mixing_parameter=[bestdf]; fit%converged=ok; fit%iterations=60
    end if
    fit%ok=fit%converged
    if(.not.fit%ok)fit%message='Student fitting failed'
    fit%aic=-2.0_dp*fit%log_likelihood+2.0_dp*real(size(x,2)+size(x,2)*(size(x,2)+1)/2+merge(1,0,est),dp)
    fit%bic=-2.0_dp*fit%log_likelihood+log(real(size(x,1),dp))*&
      real(size(x,2)+size(x,2)*(size(x,2)+1)/2+merge(1,0,est),dp)
  end function

  real(dp) function negative_profile(x,df,maxit) result(v)
    real(dp), intent(in) :: x(:,:),df
    integer, intent(in) :: maxit
    real(dp), allocatable :: loc(:),scale(:,:)
    real(dp) :: ll
    logical :: ok
    call student_profile(x,df,maxit,loc,scale,ll,ok)
    if(ok)then; v=-ll; else; v=huge(1.0_dp); end if
  end function

  subroutine student_profile(x,df,maxit,loc,scale,ll,ok)
    real(dp), intent(in) :: x(:,:),df
    integer, intent(in) :: maxit
    real(dp), allocatable, intent(out) :: loc(:),scale(:,:)
    real(dp), intent(out) :: ll
    logical, intent(out) :: ok
    real(dp), allocatable :: newloc(:),newscale(:,:),delta(:),weights(:)
    real(dp) :: q,sw,change
    type(nvmix_model) :: model
    integer :: n,p,i,j,k,iter
    n=size(x,1); p=size(x,2); call sample_mean_covariance(x,loc,scale,ok); if(.not.ok)return
    allocate(newloc(p),newscale(p,p),delta(p),weights(n))
    do iter=1,maxit
      do i=1,n
        q=mahalanobis_squared(x(i,:),loc,scale,ok); if(.not.ok)return
        weights(i)=(df+real(p,dp))/(df+q)
      end do
      sw=sum(weights); newloc=0.0_dp
      do i=1,n; newloc=newloc+weights(i)*x(i,:); end do; newloc=newloc/sw
      newscale=0.0_dp
      do i=1,n
        delta=x(i,:)-newloc
        do j=1,p; do k=1,j; newscale(j,k)=newscale(j,k)+weights(i)*delta(j)*delta(k); end do; end do
      end do
      newscale=newscale/real(n,dp)
      do j=1,p; do k=1,j-1; newscale(k,j)=newscale(j,k); end do; end do
      change=max(maxval(abs(newloc-loc)),maxval(abs(newscale-scale)))
      loc=newloc; scale=newscale
      if(change<1.0e-9_dp)exit
    end do
    model=make_nvmix_model(loc,scale,mix_inverse_gamma,df); ll=0.0_dp
    do i=1,n; ll=ll+nvmix_logpdf(x(i,:),model); end do
    ok=ll>-huge(1.0_dp)/2.0_dp
  end subroutine

  function fit_nvmix(x,family,estimate_parameter,parameter,control) result(fit)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: family
    logical, intent(in), optional :: estimate_parameter
    real(dp), intent(in), optional :: parameter
    type(integration_control), intent(in), optional :: control
    type(fit_result) :: fit
    type(fit_result) :: nf
    type(nvmix_model) :: model
    type(integration_control) :: ctrl
    logical :: est
    real(dp) :: par,best,ll,bestll
    integer :: i,j
    ctrl=integration_control(samples=2048,batches=8,seed=5489); if(present(control))ctrl=control
    if(family==mix_constant)then; fit=fit_norm(x); fit%mixing_parameter=[1.0_dp]; return; end if
    if(family==mix_inverse_gamma)then
      est=.true.; if(present(estimate_parameter))est=estimate_parameter
      if(present(parameter))then; fit=fit_student(x,est,parameter); else; fit=fit_student(x,est); end if
      return
    end if
    nf=fit_norm(x); if(.not.nf%ok)then; fit=nf; return; end if
    est=.true.; if(present(estimate_parameter))est=estimate_parameter
    par=3.0_dp; if(present(parameter))par=parameter
    best=par; bestll=-huge(1.0_dp)
    if(est)then
      do j=1,50
        if(family==mix_pareto)then; par=0.55_dp+0.25_dp*real(j-1,dp)
        else; par=0.2_dp+0.3_dp*real(j-1,dp); end if
        model=make_nvmix_model(nf%loc,nf%scale,family,par); ll=0.0_dp
        do i=1,size(x,1); ll=ll+nvmix_logpdf(x(i,:),model,ctrl); end do
        if(ll>bestll)then; bestll=ll; best=par; end if
      end do
    else
      model=make_nvmix_model(nf%loc,nf%scale,family,par); bestll=0.0_dp
      do i=1,size(x,1); bestll=bestll+nvmix_logpdf(x(i,:),model,ctrl); end do; best=par
    end if
    fit%loc=nf%loc; fit%scale=nf%scale; fit%mixing_parameter=[best]; fit%log_likelihood=bestll
    fit%converged=.true.; fit%iterations=merge(50,0,est)
    fit%aic=-2.0_dp*bestll+2.0_dp*real(size(x,2)+size(x,2)*(size(x,2)+1)/2+merge(1,0,est),dp)
    fit%bic=-2.0_dp*bestll+log(real(size(x,1),dp))*&
      real(size(x,2)+size(x,2)*(size(x,2)+1)/2+merge(1,0,est),dp)
  end function

  function fit_student_copula(u) result(fit)
    real(dp), intent(in) :: u(:,:)
    type(fit_result) :: fit
    real(dp), allocatable :: x(:,:),mean(:),cov(:,:),corr(:,:),sd(:)
    real(dp) :: df,ll,bestll,bestdf
    logical :: ok
    integer :: n,d,i,j,k
    n=size(u,1); d=size(u,2); allocate(x(n,d)); bestll=-huge(1.0_dp); bestdf=10.0_dp
    do k=1,50
      df=1.2_dp+0.8_dp*real(k-1,dp)
      do i=1,n; do j=1,d; x(i,j)=student_quantile(min(1.0_dp-epsilon(1.0_dp),max(epsilon(1.0_dp),u(i,j))),df); end do; end do
      call sample_mean_covariance(x,mean,cov,ok); if(.not.ok)cycle
      call covariance_to_correlation(cov,corr,sd,ok); if(.not.ok)cycle
      ll=0.0_dp; do i=1,n; ll=ll+dstudent_copula(u(i,:),df,corr,.true.); end do
      if(ll>bestll)then; bestll=ll; bestdf=df; fit%scale=corr; end if
    end do
    fit%loc=[(0.0_dp,j=1,d)]; fit%mixing_parameter=[bestdf]; fit%log_likelihood=bestll
    fit%converged=bestll>-huge(1.0_dp)/2.0_dp; fit%ok=fit%converged
  end function

  function fit_grouped_student_copula(u,groupings) result(fit)
    real(dp), intent(in) :: u(:,:)
    integer, intent(in) :: groupings(:)
    type(fit_result) :: fit
    real(dp), allocatable :: x(:,:),mean(:),cov(:,:),corr(:,:),sd(:),dfs(:)
    real(dp) :: m2,m4,kurt
    logical :: ok
    integer :: n,d,g,i,j
    n=size(u,1); d=size(u,2); g=maxval(groupings); allocate(x(n,d),dfs(g)); dfs=10.0_dp
    do j=1,d
      do i=1,n; x(i,j)=student_quantile(min(1.0_dp-epsilon(1.0_dp),max(epsilon(1.0_dp),u(i,j))),10.0_dp); end do
      m2=sum((x(:,j)-sum(x(:,j))/real(n,dp))**2)/real(n,dp)
      m4=sum((x(:,j)-sum(x(:,j))/real(n,dp))**4)/real(n,dp)
      kurt=m4/max(m2*m2,tiny(1.0_dp))
      if(kurt>3.05_dp)dfs(groupings(j))=max(4.1_dp,min(100.0_dp,4.0_dp+6.0_dp/(kurt-3.0_dp)))
    end do
    do j=1,d; do i=1,n; x(i,j)=student_quantile(u(i,j),dfs(groupings(j))); end do; end do
    call sample_mean_covariance(x,mean,cov,ok); call covariance_to_correlation(cov,corr,sd,ok)
    fit%loc=[(0.0_dp,j=1,d)]; fit%scale=corr; fit%mixing_parameter=dfs; fit%log_likelihood=0.0_dp
    do i=1,n; fit%log_likelihood=fit%log_likelihood+dgrouped_student_copula(u(i,:),groupings,dfs,corr,log_density=.true.); end do
    fit%ok=ok; fit%converged=ok
  end function

  function qqplot_maha(x,model,control) result(qr)
    real(dp), intent(in) :: x(:,:)
    type(nvmix_model), intent(in) :: model
    type(integration_control), intent(in), optional :: control
    type(qq_result) :: qr
    real(dp), allocatable :: values(:)
    real(dp) :: p
    logical :: ok
    integer :: n,i
    n=size(x,1); allocate(values(n),qr%observed(n),qr%theoretical(n))
    do i=1,n; values(i)=mahalanobis_squared(x(i,:),model%loc,model%scale,ok); if(.not.ok)then; qr%ok=.false.; return; end if; end do
    call sort_values(values); qr%observed=values
    do i=1,n
      p=(real(i,dp)-0.5_dp)/real(n,dp)
      if(model%groups()==1)then
        if(present(control))then
          qr%theoretical(i)=qgammamix(p,model%dimension(),model%mix_family(1),model%mix_parameter(1),control)
        else
          qr%theoretical(i)=qgammamix(p,model%dimension(),model%mix_family(1),model%mix_parameter(1))
        end if
      else
        qr%theoretical(i)=chi_square_quantile(p,real(model%dimension(),dp))
      end if
    end do
  end function
  subroutine sort_values(x)
    real(dp), intent(inout) :: x(:)
    real(dp) :: key
    integer :: i,j
    do i=2,size(x); key=x(i); j=i-1; do while(j>=1); if(x(j)<=key)exit; x(j+1)=x(j); j=j-1; end do; x(j+1)=key; end do
  end subroutine
end module nvmix_fitting

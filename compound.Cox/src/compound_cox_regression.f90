! SPDX-License-Identifier: GPL-2.0-only
module compound_cox_regression
  use compound_cox_kinds, only : dp
  use compound_cox_types, only : depend_cox_result, depend_cv_result, compound_result
  use compound_cox_math, only : bfgs_minimize, invert_matrix, chisq_cdf, scalar_second_derivative
  use numderiv, only : hessian, hessian_options, deriv_options, nd_success
  use survival_cox, only : coxph_fit
  use survival_types, only : coxph_result, concordance_result
  use survival_stats, only : concordance_right
  implicit none
  private
  public :: depend_cox_reg, depend_cox_reg_cv, cindex_cv, compound_reg
  public :: cox_loglik_breslow
contains
  real(dp) function cox_loglik_breslow(time,status,x,beta) result(ll)
    real(dp),intent(in)::time(:),x(:,:),beta(:)
    integer,intent(in)::status(:)
    real(dp),allocatable::eta(:)
    real(dp)::den
    integer::n,i,j
    n=size(time)
    allocate(eta(n))
    eta=matmul(x,beta)
    ll=0
    do i=1,n
      if(status(i)==0)cycle
      den=0
      do j=1,n
      if(time(j)>=time(i))den=den+exp(min(eta(j),700.0_dp))
      end do
      ll=ll+eta(i)-log(max(den,tiny(1.0_dp)))
    end do
  end function cox_loglik_breslow

  real(dp) function univ_sum_loglik(time,status,x,beta) result(ll)
    real(dp),intent(in)::time(:),x(:,:),beta(:)
    integer,intent(in)::status(:)
    integer::j
    real(dp)::xx(size(time),1),bb(1)
    ll=0
    do j=1,size(x,2)
    xx(:,1)=x(:,j)
    bb(1)=beta(j)
    ll=ll+cox_loglik_breslow(time,status,xx,bb)
    end do
  end function univ_sum_loglik

  subroutine subset_rows(time,status,x,keep,t2,d2,x2)
    real(dp),intent(in)::time(:),x(:,:)
    integer,intent(in)::status(:)
    logical,intent(in)::keep(:)
    real(dp),allocatable,intent(out)::t2(:),x2(:,:)
    integer,allocatable,intent(out)::d2(:)
    integer::m,i,k
    m=count(keep)
    allocate(t2(m),d2(m),x2(m,size(x,2)))
    k=0
    do i=1,size(time)
    if(keep(i))then
    k=k+1
    t2(k)=time(i)
    d2(k)=status(i)
    x2(k,:)=x(i,:)
    end if
    end do
  end subroutine subset_rows

  subroutine fit_shrinkage(time,status,x,a,beta,fval,converged)
    real(dp),intent(in)::time(:),x(:,:),a
    integer,intent(in)::status(:)
    real(dp),allocatable,intent(out)::beta(:)
    real(dp),intent(out)::fval
    logical,intent(out)::converged
    real(dp)::aa
    allocate(beta(size(x,2)))
    beta=0
    aa=min(a,0.99_dp)
    call bfgs_minimize(obj,beta,fval,converged,maxiter=250,tol=2e-6_dp)
  contains
    real(dp) function obj(b) result(v)
      real(dp),intent(in)::b(:)
      v=-(aa*cox_loglik_breslow(time,status,x,b)+(1-aa)*univ_sum_loglik(time,status,x,b))
    end function obj
  end subroutine fit_shrinkage

  real(dp) function cv_a_value(time,status,x,kfold,a) result(cv)
    real(dp),intent(in)::time(:),x(:,:),a
    integer,intent(in)::status(:),kfold
    logical,allocatable::keep(:)
    real(dp),allocatable::tt(:),xx(:,:),b(:)
    integer,allocatable::dd(:)
    integer::n,k,lo,hi
    real(dp)::fv
    logical::conv
    n=size(time)
    allocate(keep(n))
    cv=0
    do k=1,kfold
      lo=(k-1)*n/kfold+1
      hi=k*n/kfold
      keep=.true.
      keep(lo:hi)=.false.
      call subset_rows(time,status,x,keep,tt,dd,xx)
      call fit_shrinkage(tt,dd,xx,a,b,fv,conv)
      cv=cv+cox_loglik_breslow(time,status,x,b)-cox_loglik_breslow(tt,dd,xx,b)
    end do
  end function cv_a_value

  subroutine compound_reg(time,status,x,res,kfold,delta_a,a0,with_variance)
    real(dp),intent(in)::time(:),x(:,:)
    integer,intent(in)::status(:)
    type(compound_result),intent(out)::res
    integer,intent(in),optional::kfold
    real(dp),intent(in),optional::delta_a,a0
    logical,intent(in),optional::with_variance
    integer::k,p,n,istat
    real(dp)::da,aold,anew,cvold,cvnew,fv,ez
    logical::dovar,conv,ok
    real(dp),allocatable::h(:,:),hinv(:,:),vhat(:,:),hdot(:),ahat(:,:),sig(:,:),bb(:)
    type(deriv_options)::hopts
    k=5
    if(present(kfold))k=kfold
    da=0.025_dp
    if(present(delta_a))da=delta_a
    aold=0
    if(present(a0))aold=a0
    dovar=.false.
    if(present(with_variance))dovar=with_variance
    n=size(time)
    p=size(x,2)
    cvold=cv_a_value(time,status,x,k,aold)
    do
      anew=aold+da
      if(anew>=1)exit
      cvnew=cv_a_value(time,status,x,k,anew)
      if(cvnew<=cvold)exit
      cvold=cvnew
      aold=anew
    end do
    res%a=min(aold,1.0_dp)
    call fit_shrinkage(time,status,x,res%a,bb,fv,conv)
    res%converged=conv
    allocate(res%beta(p))
    res%beta=bb
    if(.not.dovar)return
    hopts=hessian_options()
    call hessian(obj_final,bb,h,hopts,status=istat)
    if(istat/=nd_success)then
    allocate(h(p,p))
    h=0
    end if
    allocate(vhat(p,p))
    vhat=h/real(n,dp)
    call invert_matrix(vhat,hinv,ok)
    if(.not.ok)then
    hinv=0
    end if
    allocate(hdot(p))
    hdot=0
    call calculate_hdot(time,status,x,bb,hdot)
    ez=-scalar_second_derivative(cvscalar,res%a)/real(n,dp)
    if(abs(ez)<1e-10_dp)ez=sign(1e-10_dp,ez+1e-30_dp)
    allocate(ahat(p,p))
    ahat=identity(p)+matmul(reshape(matmul(hinv,hdot),[p,1]),reshape(hdot,[1,p]))/ez
    allocate(sig(p,p))
    sig=matmul(ahat,matmul(hinv,transpose(ahat)))
    allocate(res%se(p),res%lower95(p),res%upper95(p),res%sigma(p,p),res%v(p,p),res%h_dot(p))
    res%se=sqrt(max(diagonal(sig)/real(n,dp),0.0_dp))
    res%lower95=bb-1.96_dp*res%se
    res%upper95=bb+1.96_dp*res%se
    res%sigma=sig
    res%v=vhat
    res%h_dot=hdot
    res%hessian_cv=scalar_second_derivative(cvscalar,res%a)
  contains
    real(dp) function obj_final(b) result(v)
      real(dp),intent(in)::b(:)
      v=-(min(res%a,0.99_dp)*cox_loglik_breslow(time,status,x,b)+(1-min(res%a,0.99_dp))*univ_sum_loglik(time,status,x,b))
    end function obj_final
    real(dp) function cvscalar(a) result(v)
      real(dp),intent(in)::a
      v=cv_a_value(time,status,x,k,a)
    end function cvscalar
  end subroutine compound_reg

  subroutine calculate_hdot(time,status,x,b,hdot)
    real(dp),intent(in)::time(:),x(:,:),b(:)
    integer,intent(in)::status(:)
    real(dp),intent(out)::hdot(:)
    integer::n,p,i,j,k
    real(dp)::s0
    real(dp),allocatable::s1(:),s0v(:),s1v(:),eta(:)
    n=size(time)
    p=size(x,2)
    allocate(s1(p),s0v(p),s1v(p),eta(n))
    hdot=0
    do i=1,n
      if(status(i)==0)cycle
      s0=0
      s1=0
      s0v=0
      s1v=0
      do j=1,n
        if(time(j)<time(i))cycle
        eta(j)=dot_product(x(j,:),b)
        s0=s0+exp(min(eta(j),700.0_dp))
        s1=s1+x(j,:)*exp(min(eta(j),700.0_dp))
        do k=1,p
        s0v(k)=s0v(k)+exp(min(x(j,k)*b(k),700.0_dp))
        s1v(k)=s1v(k)+x(j,k)*exp(min(x(j,k)*b(k),700.0_dp))
        end do
      end do
      hdot=hdot-s1/s0+s1v/s0v
    end do
    hdot=hdot/real(n,dp)
  end subroutine calculate_hdot

  pure function identity(n) result(a)
    integer,intent(in)::n
    real(dp)::a(n,n)
    integer::i
    a=0
    do i=1,n
    a(i,i)=1
    end do
  end function identity
  pure function diagonal(a) result(d)
    real(dp),intent(in)::a(:,:)
    real(dp)::d(min(size(a,1),size(a,2)))
    integer::i
    do i=1,size(d)
    d(i)=a(i,i)
    end do
  end function diagonal

  subroutine sort_triplet(time,status,x,ts,ds,xs)
    real(dp),intent(in)::time(:),x(:)
    integer,intent(in)::status(:)
    real(dp),allocatable,intent(out)::ts(:),xs(:)
    integer,allocatable,intent(out)::ds(:)
    integer::n,i,j
    real(dp)::tv,xv
    integer::dv
    n=size(time)
    allocate(ts(n),xs(n),ds(n))
    ts=time
    xs=x
    ds=status
    do i=2,n
    tv=ts(i)
    xv=xs(i)
    dv=ds(i)
    j=i-1
    do while(j>=1)
    if(ts(j)<=tv)exit
    ts(j+1)=ts(j)
    xs(j+1)=xs(j)
    ds(j+1)=ds(j)
    j=j-1
    end do
    ts(j+1)=tv
    xs(j+1)=xv
    ds(j+1)=dv
    end do
  end subroutine sort_triplet

  real(dp) function dep_objective(par,ds,xs,alpha) result(v)
    real(dp),intent(in)::par(:),xs(:),alpha
    integer,intent(in)::ds(:)
    integer::n,n1,i,k1,k2
    real(dp)::b1,b2,r1,r2,dr1,dr2,le1,le2,lden,ll
    n=size(ds)
    n1=sum(ds)
    b1=par(n+1)
    b2=par(n+2)
    r1=0
    r2=0
    k1=0
    k2=n1
    ll=0
    do i=1,n
      dr1=0
      dr2=0
      if(ds(i)==1)then
      k1=k1+1
      dr1=exp(min(par(k1),700.0_dp))
      r1=r1+dr1
      else
      k2=k2+1
      dr2=exp(min(par(k2),700.0_dp))
      r2=r2+dr2
      end if
      le1=alpha*r1*exp(min(b1*xs(i),700.0_dp))
      le2=alpha*r2*exp(min(b2*xs(i),700.0_dp))
      lden=log_exp_sum_minus_one(le1,le2)
      if(ds(i)==1)ll=ll+b1*xs(i)+(le1-lden)+log(dr1)
      if(ds(i)==0)ll=ll+b2*xs(i)+(le2-lden)+log(dr2)
      ll=ll-lden/alpha
    end do
    v=-ll
  end function dep_objective

  pure real(dp) function log_exp_sum_minus_one(a,b) result(v)
    real(dp),intent(in)::a,b
    real(dp)::m,s
    m=max(a,b)
    s=exp(a-m)+exp(b-m)-exp(-m)
    v=m+log(max(s,tiny(1.0_dp)))
  end function log_exp_sum_minus_one

  subroutine depend_cox_reg(time,status,x,alpha,res,with_variance,censor_reg,baseline)
    real(dp),intent(in)::time(:),x(:),alpha
    integer,intent(in)::status(:)
    type(depend_cox_result),intent(out)::res
    logical,intent(in),optional::with_variance,censor_reg,baseline
    real(dp),allocatable::ts(:),xs(:),par(:),hh(:,:),hinv(:,:)
    integer,allocatable::ds(:)
    integer::n,n1,i,k1,k2,istat
    real(dp)::aa,fv,r1,r2
    logical::dovar,docen,dobase,conv,ok
    type(deriv_options)::hopts
    aa=max(alpha,0.01_dp)
    dovar=.true.
    if(present(with_variance))dovar=with_variance
    docen=.false.
    if(present(censor_reg))docen=censor_reg
    dobase=.false.
    if(present(baseline))dobase=baseline
    call sort_triplet(time,status,x,ts,ds,xs)
    n=size(time)
    n1=sum(ds)
    allocate(par(n+2))
    par(1:n)=log(1.0_dp/real(n,dp))
    par(n+1:n+2)=0
    call bfgs_minimize(obj,par,fv,conv,maxiter=500,tol=2e-6_dp)
    res%converged=conv
    res%beta=par(n+1)
    res%beta_censor=par(n+2)
    if(dovar)then
      hopts=hessian_options()
      call hessian(obj,par,hh,hopts,status=istat)
      if(istat==nd_success)then
      call invert_matrix(hh,hinv,ok)
      if(ok)then
      res%se=sqrt(max(hinv(n+1,n+1),0.0_dp))
      res%z=res%beta/res%se
      res%p=1-chisq_cdf(res%z**2,1.0_dp)
      res%se_censor=sqrt(max(hinv(n+2,n+2),0.0_dp))
      res%z_censor=res%beta_censor/res%se_censor
      res%p_censor=1-chisq_cdf(res%z_censor**2,1.0_dp)
      end if
      end if
    end if
    if(dobase)then
      allocate(res%baseline(n),res%censor_baseline(n))
      res%baseline=0
      res%censor_baseline=0
      r1=0
      r2=0
      k1=0
      k2=n1
      do i=1,n
      if(ds(i)==1)then
      k1=k1+1
      r1=r1+exp(par(k1))
      else
      k2=k2+1
      r2=r2+exp(par(k2))
      end if
      res%baseline(i)=r1
      res%censor_baseline(i)=r2
      end do
    end if
  contains
    real(dp) function obj(q) result(v)
      real(dp),intent(in)::q(:)
      v=dep_objective(q,ds,xs,aa)
    end function obj
  end subroutine depend_cox_reg

  real(dp) function cindex_cv(time,status,x,alpha,kfold) result(ci)
    real(dp),intent(in)::time(:),x(:,:),alpha
    integer,intent(in)::status(:),kfold
    integer::n,p,k,j,lo,hi,i
    logical,allocatable::keep(:)
    real(dp),allocatable::tt(:),xx(:,:),risk(:)
    integer,allocatable::dd(:)
    type(depend_cox_result)::dr
    type(concordance_result)::cc
    n=size(time)
    p=size(x,2)
    allocate(keep(n),risk(n))
    risk=0
    do k=1,kfold
    lo=(k-1)*n/kfold+1
    hi=k*n/kfold
    keep=.true.
    keep(lo:hi)=.false.
    call subset_rows(time,status,x,keep,tt,dd,xx)
      do j=1,p
      call depend_cox_reg(tt,dd,xx(:,j),alpha,dr,with_variance=.false.)
      do i=lo,hi
      risk(i)=risk(i)+x(i,j)*dr%beta
      end do
      end do
    end do
    call concordance_right(time,status,risk,cc)
    ci=cc%cindex
  end function cindex_cv

  subroutine depend_cox_reg_cv(time,status,x,res,kfold,ngrid)
    real(dp),intent(in)::time(:),x(:,:)
    integer,intent(in)::status(:)
    type(depend_cv_result),intent(out)::res
    integer,intent(in),optional::kfold,ngrid
    integer::k,g,i,j,p
    real(dp)::tau,alpha,ci,best
    type(depend_cox_result)::dr
    k=5
    if(present(kfold))k=kfold
    g=20
    if(present(ngrid))g=ngrid
    best=-1
    res%alpha=0
    do i=1,g
    tau=0.0001_dp+(0.9_dp-0.0001_dp)*real(i-1,dp)/real(max(1,g-1),dp)
    alpha=2*tau/(1-tau)
    ci=cindex_cv(time,status,x,alpha,k)
    if(ci>best)then
    best=ci
    res%alpha=alpha
    end if
    end do
    p=size(x,2)
    allocate(res%beta(p),res%se(p),res%z(p),res%p(p))
    do j=1,p
    call depend_cox_reg(time,status,x(:,j),res%alpha,dr)
    res%beta(j)=dr%beta
    res%se(j)=dr%se
    end do
    res%z=res%beta/res%se
    do j=1,p
    res%p(j)=1-chisq_cdf(res%z(j)**2,1.0_dp)
    end do
    res%c_index=best
  end subroutine depend_cox_reg_cv
end module compound_cox_regression

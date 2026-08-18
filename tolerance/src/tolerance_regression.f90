! SPDX-License-Identifier: GPL-2.0-or-later
module tolerance_regression
  use tolerance_kinds, only : dp
  use tolerance_types, only : regression_band
  use tolerance_math, only : solve_linear, invert_spd, noncentral_t_quantile, &
       noncentral_chisq_quantile, chisq_quantile, normal_quantile
  use tolerance_nonparametric, only : nptol_int
  use tolerance_optimize, only : nelder_mead
  use tolerance_types, only : tolerance_interval
  use tolerance_normal, only : k_factor
  implicit none
  private
  public :: linear_regression_fit, regtol_int, npregtol_int, nlregtol_int
  public :: anova_group_tol

  abstract interface
    subroutine nonlinear_model(x,beta,yhat)
      import dp
      real(dp),intent(in)::x(:,:),beta(:)
      real(dp),intent(out)::yhat(:)
    end subroutine nonlinear_model
  end interface
contains

  subroutine linear_regression_fit(x,y,beta,mse,df,cov_unscaled,info)
    real(dp),intent(in)::x(:,:),y(:)
    real(dp),intent(out)::beta(:),mse,cov_unscaled(:,:)
    integer,intent(out)::df,info
    real(dp),allocatable::xtx(:,:),xty(:),res(:)
    integer::n,p
    n=size(x,1);p=size(x,2);allocate(xtx(p,p),xty(p),res(n))
    xtx=matmul(transpose(x),x);xty=matmul(transpose(x),y)
    call solve_linear(xtx,xty,beta,info);if(info/=0)then;mse=huge(1.0_dp);df=0;return;end if
    call invert_spd(xtx,cov_unscaled,info);res=y-matmul(x,beta);df=n-p
    if(df>0)then;mse=dot_product(res,res)/real(df,dp);else;mse=huge(1.0_dp);end if
  end subroutine linear_regression_fit

  function regtol_int(x,y,xpred,alpha,p,side) result(band)
    real(dp),intent(in)::x(:,:),y(:),xpred(:,:)
    real(dp),intent(in),optional::alpha,p
    integer,intent(in),optional::side
    type(regression_band)::band
    real(dp),allocatable::beta(:),cov(:,:),fit(:)
    real(dp)::a,pp,mse,h,nstar,delta,k
    integer::n,pars,df,info,i,sd
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(p))pp=p;sd=1;if(present(side))sd=side
    n=size(x,1);pars=size(x,2);allocate(beta(pars),cov(pars,pars),fit(size(xpred,1)))
    call linear_regression_fit(x,y,beta,mse,df,cov,info)
    fit=matmul(xpred,beta);allocate(band%fit(size(fit)),band%lower(size(fit)),band%upper(size(fit)));band%fit=fit
    do i=1,size(fit)
      h=dot_product(xpred(i,:),matmul(cov,xpred(i,:)))
      if(h<=0.0_dp)then;k=huge(1.0_dp)
      else
        nstar=1.0_dp/h
        if(sd==1)then
          delta=sqrt(nstar)*normal_quantile(pp)
          k=noncentral_t_quantile(1.0_dp-a,real(n-pars,dp),delta)/sqrt(nstar)
        else
          k=sqrt(real(df,dp)*noncentral_chisq_quantile(pp,1.0_dp,1.0_dp/nstar)/ &
               chisq_quantile(a,real(df,dp)))
        end if
      end if
      band%lower(i)=fit(i)-sqrt(mse)*k;band%upper(i)=fit(i)+sqrt(mse)*k
    end do
  end function regtol_int

  function npregtol_int(y,yhat,alpha,p,side,method,lower_bound,upper_bound) result(band)
    real(dp),intent(in)::y(:),yhat(:)
    real(dp),intent(in),optional::alpha,p,lower_bound,upper_bound
    integer,intent(in),optional::side
    character(len=*),intent(in),optional::method
    type(regression_band)::band
    type(tolerance_interval)::ti
    real(dp),allocatable::res(:)
    real(dp)::a,pp
    integer::sd
    character(len=8)::meth
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(p))pp=p;sd=1;if(present(side))sd=side
    meth='WILKS';if(present(method))meth=adjustl(method)
    allocate(res(size(y)));res=y-yhat;ti=nptol_int(res,a,pp,sd,meth)
    allocate(band%fit(size(y)),band%lower(size(y)),band%upper(size(y)))
    band%fit=yhat;band%lower=yhat+ti%lower;band%upper=yhat+ti%upper
    if(present(lower_bound))then;band%lower=max(band%lower,lower_bound);band%upper=max(band%upper,lower_bound);end if
    if(present(upper_bound))then;band%lower=min(band%lower,upper_bound);band%upper=min(band%upper,upper_bound);end if
  end function npregtol_int

  subroutine nlregtol_int(model,x,y,beta,xpred,band,alpha,p,side,max_iter,info)
    procedure(nonlinear_model)::model
    real(dp),intent(in)::x(:,:),y(:),xpred(:,:)
    real(dp),intent(inout)::beta(:)
    type(regression_band),intent(out)::band
    real(dp),intent(in),optional::alpha,p
    integer,intent(in),optional::side,max_iter
    integer,intent(out),optional::info
    real(dp),allocatable::yh(:),yp(:),jac(:,:),jtj(:,:),jinv(:,:),bp(:),bm(:)
    real(dp),allocatable::yp_work_plus(:),yp_work_minus(:),pgrad(:)
    real(dp)::a,pp,mse,h,nstar,delta,k,step
    integer::n,np,nb,df,i,j,mi,inf
    n=size(x,1);np=size(xpred,1);nb=size(beta);a=0.05_dp;if(present(alpha))a=alpha
    pp=0.99_dp;if(present(p))pp=p;mi=1000;if(present(max_iter))mi=max_iter
    allocate(yh(n),yp(np),jac(n,nb),jtj(nb,nb),jinv(nb,nb),bp(nb),bm(nb))
    allocate(yp_work_plus(n),yp_work_minus(n))
    call nelder_mead(nll,beta,step=0.1_dp,tol=1.0e-10_dp,max_iter=mi)
    call model(x,beta,yh);df=n-nb;mse=sum((y-yh)**2)/real(max(df,1),dp)
    do j=1,nb
      step=1.0e-5_dp*max(1.0_dp,abs(beta(j)));bp=beta;bm=beta;bp(j)=bp(j)+step;bm(j)=bm(j)-step
      call model(x,bp,yp_work_plus);call model(x,bm,yp_work_minus)
      jac(:,j)=(yp_work_plus-yp_work_minus)/(2.0_dp*step)
    end do
    jtj=matmul(transpose(jac),jac);call invert_spd(jtj,jinv,inf);if(present(info))info=inf
    call model(xpred,beta,yp);allocate(band%fit(np),band%lower(np),band%upper(np));band%fit=yp
    do i=1,np
      call prediction_gradient(model,xpred(i:i,:),beta,pgrad)
      h=dot_product(pgrad,matmul(jinv,pgrad));nstar=1.0_dp/max(h,tiny(1.0_dp))
      if(present(side))then;j=side;else;j=1;end if
      if(j==1)then
        delta=sqrt(nstar)*normal_quantile(pp);k=noncentral_t_quantile(1.0_dp-a,real(df,dp),delta)/sqrt(nstar)
      else
        k=sqrt(real(df,dp)*noncentral_chisq_quantile(pp,1.0_dp,1.0_dp/nstar)/chisq_quantile(a,real(df,dp)))
      end if
      band%lower(i)=yp(i)-sqrt(mse)*k;band%upper(i)=yp(i)+sqrt(mse)*k
    end do
  contains
    real(dp) function nll(b) result(v)
      real(dp),intent(in)::b(:)
      call model(x,b,yh);v=sum((y-yh)**2)
    end function nll
  end subroutine nlregtol_int

  subroutine prediction_gradient(model,xrow,beta,g)
    procedure(nonlinear_model)::model
    real(dp),intent(in)::xrow(:,:),beta(:)
    real(dp),allocatable,intent(out)::g(:)
    real(dp),allocatable::bp(:),bm(:),fp(:),fm(:)
    real(dp)::h
    integer::j
    allocate(g(size(beta)),bp(size(beta)),bm(size(beta)),fp(1),fm(1))
    do j=1,size(beta)
      h=1.0e-5_dp*max(1.0_dp,abs(beta(j)));bp=beta;bm=beta;bp(j)=bp(j)+h;bm(j)=bm(j)-h
      call model(xrow,bp,fp);call model(xrow,bm,fm);g(j)=(fp(1)-fm(1))/(2.0_dp*h)
    end do
  end subroutine prediction_gradient

  subroutine anova_group_tol(values,groups,ngroups,mse,error_df,lower,upper,kvals,alpha,p,side,method,m)
    real(dp),intent(in)::values(:),mse
    integer,intent(in)::groups(:),ngroups,error_df
    real(dp),intent(out)::lower(ngroups),upper(ngroups),kvals(ngroups)
    real(dp),intent(in),optional::alpha,p
    integer,intent(in),optional::side,m
    character(len=*),intent(in),optional::method
    real(dp)::a,pp,mu
    integer::g,n,sd,mm
    character(len=8)::meth
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(p))pp=p;sd=1;if(present(side))sd=side
    mm=50;if(present(m))mm=m;meth='HE';if(present(method))meth=adjustl(method)
    do g=1,ngroups
      n=count(groups==g);if(n>0)then;mu=sum(values,mask=groups==g)/real(n,dp);else;mu=0.0_dp;end if
      kvals(g)=k_factor(n,a,pp,sd,trim(meth),f=error_df,m=mm);lower(g)=mu-kvals(g)*sqrt(mse);upper(g)=mu+kvals(g)*sqrt(mse)
    end do
  end subroutine anova_group_tol
end module tolerance_regression

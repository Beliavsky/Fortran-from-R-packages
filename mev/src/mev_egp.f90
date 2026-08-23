module mev_egp
  use mev_kinds, only: dp
  use mev_distributions, only: degp, qegp
  use mev_math, only: pattern_minimize, finite_diff_hessian, inverse_matrix, mean_real, variance_real
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  type, public :: egp_fit_result
    real(dp) :: kappa=1.0_dp, scale=1.0_dp, shape=0.0_dp
    real(dp) :: loglik=-huge(1.0_dp)
    real(dp) :: vcov(3,3)=0.0_dp
    integer :: convergence=1, nobs=0
  end type egp_fit_result
  public :: egp_ll, egp_fit, egp_retlev
contains
  real(dp) function egp_ll(xdat,thresh,kappa,scale,shape,model) result(ll)
    real(dp),intent(in)::xdat(:),thresh,kappa,scale,shape
    character(len=*),intent(in)::model
    integer::i
    real(dp)::z,d
    ll=0.0_dp
    if(kappa<=0.0_dp.or.scale<=0.0_dp.or.shape<=-1.0_dp)then
      ll=-huge(1.0_dp);return
    end if
    do i=1,size(xdat)
      if(xdat(i)>thresh)then
        z=xdat(i)-thresh
        if(shape<0.0_dp.and.1.0_dp+shape*z/scale<=0.0_dp)then
          ll=-huge(1.0_dp);return
        end if
        d=degp(z,scale,shape,kappa,model,.true.)
        if(.not.ieee_is_finite(d))then;ll=-huge(1.0_dp);return;end if
        ll=ll+d
      end if
    end do
  end function egp_ll

  subroutine egp_fit(xdat,thresh,model,fit,start,info)
    real(dp),intent(in)::xdat(:),thresh
    character(len=*),intent(in)::model
    type(egp_fit_result),intent(out)::fit
    real(dp),intent(in),optional::start(3)
    integer,intent(out),optional::info
    real(dp)::x(3),fval,mu,var,h(3,3),hinv(3,3),jac(3,3)
    real(dp)::par0(3)
    integer::ier,nexc
    nexc=count(xdat>thresh);fit%nobs=nexc
    if(nexc<5)then;fit%convergence=2;if(present(info))info=2;return;end if
    if(present(start))then
      par0=start
    else
      mu=sum(pack(xdat-thresh,xdat>thresh))/real(nexc,dp)
      var=sum(pack((xdat-thresh-mu)**2,xdat>thresh))/real(max(1,nexc-1),dp)
      par0(1)=1.01_dp
      if(var>0.0_dp)then
        par0(3)=max(-0.8_dp,min(0.8_dp,0.5_dp*(1.0_dp-mu*mu/var)))
      else
        par0(3)=0.05_dp
      end if
      par0(2)=max(tiny(1.0_dp),mu*(1.0_dp-par0(3)))
    end if
    x=[log(max(par0(1),1.0e-6_dp)),log(max(par0(2),1.0e-8_dp)),log(max(par0(3)+1.0_dp,1.0e-6_dp))]
    call pattern_minimize(obj,x,fval,ier,maxiter=1800,tol=1.0e-8_dp,initial_step=0.2_dp)
    fit%kappa=exp(x(1));fit%scale=exp(x(2));fit%shape=exp(x(3))-1.0_dp
    fit%loglik=-fval;fit%convergence=ier
    call finite_diff_hessian(obj,x,h,1.0e-4_dp)
    call inverse_matrix(h,hinv,ier)
    if(ier==0)then
      jac=0.0_dp;jac(1,1)=fit%kappa;jac(2,2)=fit%scale;jac(3,3)=fit%shape+1.0_dp
      fit%vcov=matmul(jac,matmul(hinv,transpose(jac)))
    end if
    if(present(info))info=fit%convergence
  contains
    function obj(z) result(v)
      real(dp),intent(in)::z(:)
      real(dp)::v,ka,sc,sh,xmax
      ka=exp(z(1));sc=exp(z(2));sh=exp(z(3))-1.0_dp
      xmax=maxval(pack(xdat-thresh,xdat>thresh))
      if(sh<0.0_dp.and.sc+sh*xmax<=0.0_dp)then
        v=1.0e100_dp+1.0e8_dp*abs(sc+sh*xmax);return
      end if
      v=-egp_ll(xdat,thresh,ka,sc,sh,model)
      if(.not.ieee_is_finite(v))v=1.0e100_dp
    end function obj
  end subroutine egp_fit

  pure real(dp) function egp_retlev(thresh,rate,p,kappa,scale,shape,model) result(q)
    real(dp),intent(in)::thresh,rate,p,kappa,scale,shape
    character(len=*),intent(in)::model
    if(p<=0.0_dp)then
      q=huge(1.0_dp)
    else if(p>=rate)then
      q=thresh
    else
      q=thresh+qegp(1.0_dp-p/rate,scale,shape,kappa,model)
    end if
  end function egp_retlev
end module mev_egp

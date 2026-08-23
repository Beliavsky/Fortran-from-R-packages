module mev_stein
  use mev_kinds, only: dp
  use mev_math, only: pattern_minimize, finite_diff_hessian, inverse_matrix, variance_real, sort_descending
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: wgpd_fit_result, stein_weights, stein_gp_lik, fit_wgpd

  type :: wgpd_fit_result
    real(dp) :: estimate(2) = 0.0_dp
    real(dp) :: std_error(2) = 0.0_dp
    real(dp) :: vcov(2,2) = 0.0_dp
    real(dp) :: threshold = 0.0_dp
    real(dp) :: nllh = huge(1.0_dp)
    integer :: nat = 0
    real(dp) :: pat = 0.0_dp
    integer :: convergence = 1
    real(dp), allocatable :: weights(:)
    real(dp), allocatable :: exceedances(:)
  end type wgpd_fit_result

contains

  subroutine stein_weights(n,weights,gamma)
    integer, intent(in) :: n
    real(dp), allocatable, intent(out) :: weights(:)
    real(dp), intent(in), optional :: gamma
    real(dp) :: g
    integer :: j
    g=1.0_dp; if(present(gamma)) g=gamma
    if(n<1 .or. g<=0.0_dp) then
      allocate(weights(0)); return
    end if
    allocate(weights(n))
    do j=1,n
      weights(j)=(g+1.0_dp)/g*(1.0_dp-(real(j-1,dp)/real(n,dp))**g)
    end do
  end subroutine stein_weights

  real(dp) function stein_gp_lik(pars,xdat,weights,sorted) result(ll)
    real(dp), intent(in) :: pars(2),xdat(:),weights(:)
    logical, intent(in), optional :: sorted
    real(dp), allocatable :: x(:),xs(:),lp(:)
    real(dp) :: scale,shape
    integer :: n,j
    logical :: srt
    ll=-huge(1.0_dp)
    n=size(xdat); if(n<1 .or. size(weights)/=n) return
    scale=pars(1); shape=pars(2)
    if(scale<=0.0_dp .or. shape< -1.0_dp .or. any(xdat<=0.0_dp)) return
    srt=.false.; if(present(sorted)) srt=sorted
    allocate(x(n))
    if(srt) then
      x=xdat
    else
      allocate(xs(n)); call sort_descending(xdat,xs); x=xs
    end if
    if(shape<0.0_dp .and. 1.0_dp+shape*x(1)/scale<=0.0_dp) return
    if(abs(shape)>1.0e-5_dp) then
      allocate(lp(n)); lp=log(1.0_dp+shape*x/scale)
      ll=-sum(weights)*log(scale)
      do j=1,n
        ll=ll-weights(j)*(real(j,dp)/shape+1.0_dp)*lp(j)
      end do
      do j=1,n-1
        ll=ll+weights(j)*real(j,dp)*lp(j+1)/shape
      end do
    else
      ll=sum(weights*(-log(scale)-x/scale))
    end if
  end function stein_gp_lik

  subroutine fit_wgpd(xdat,result,threshold,gamma,start)
    real(dp), intent(in) :: xdat(:)
    type(wgpd_fit_result), intent(out) :: result
    real(dp), intent(in), optional :: threshold,gamma,start(2)
    real(dp), allocatable :: exc(:),xs(:),xnorm(:),w(:)
    real(dp) :: th,sc,x(2),fval,pars(2),h(2,2),hinv(2,2),det
    integer :: n,j,idx,info
    result=wgpd_fit_result()
    th=0.0_dp; if(present(threshold)) th=threshold
    n=count(xdat>th)
    if(n<3) return
    allocate(exc(n)); idx=0
    do j=1,size(xdat)
      if(xdat(j)>th) then; idx=idx+1; exc(idx)=xdat(j)-th; end if
    end do
    allocate(xs(n)); call sort_descending(exc,xs); exc=xs
    sc=sqrt(variance_real(exc)); if(sc<=0.0_dp) sc=max(sum(exc)/real(n,dp),1.0_dp)
    allocate(xnorm(n)); xnorm=exc/sc
    call stein_weights(n,w,gamma)
    if(present(start)) then
      pars=[max(start(1)/sc,1.0e-6_dp),min(1.999_dp,max(-0.999_dp,start(2)))]
    else
      pars=[1.0_dp,0.1_dp]
    end if
    x(1)=log(pars(1))
    x(2)=log((pars(2)+1.0_dp)/(2.0_dp-pars(2)))
    call pattern_minimize(obj_trans,x,fval,info,2000,1.0e-8_dp,0.2_dp)
    call decode(x,pars)
    result%estimate=[pars(1)*sc,pars(2)]
    result%threshold=th; result%nat=n; result%pat=real(n,dp)/real(size(xdat),dp)
    result%nllh=-stein_gp_lik(result%estimate,exc,w,.true.)
    result%convergence=info
    allocate(result%weights(n),result%exceedances(n)); result%weights=w; result%exceedances=exc
    call finite_diff_hessian(obj_actual,pars,h,1.0e-4_dp)
    call inverse_matrix(h,hinv,info)
    if(info==0 .and. all(ieee_is_finite(hinv))) then
      result%vcov=hinv
      result%vcov(1,:)=result%vcov(1,:)*sc
      result%vcov(:,1)=result%vcov(:,1)*sc
      det=hinv(1,1); if(det>0.0_dp) result%std_error(1)=sqrt(det)*sc
      det=hinv(2,2); if(det>0.0_dp) result%std_error(2)=sqrt(det)
    else
      result%vcov=ieee_value(0.0_dp,ieee_quiet_nan)
      result%std_error=ieee_value(0.0_dp,ieee_quiet_nan)
    end if
  contains
    real(dp) function obj_trans(z) result(v)
      real(dp), intent(in) :: z(:)
      real(dp) :: p(2)
      call decode(z,p)
      if(p(2)<0.0_dp .and. p(1)+p(2)*xnorm(1)<=1.0e-8_dp) then
        v=1.0e12_dp+1.0e8_dp*abs(p(1)+p(2)*xnorm(1)); return
      end if
      v=-stein_gp_lik(p,xnorm,w,.true.)
      if(.not.ieee_is_finite(v)) v=1.0e12_dp
    end function obj_trans
    real(dp) function obj_actual(p) result(v)
      real(dp), intent(in) :: p(:)
      if(p(1)<=0.0_dp .or. p(2)<=-1.0_dp .or. p(2)>=2.0_dp) then
        v=1.0e12_dp; return
      end if
      if(p(2)<0.0_dp .and. p(1)+p(2)*xnorm(1)<=0.0_dp) then
        v=1.0e12_dp; return
      end if
      v=-stein_gp_lik(p(1:2),xnorm,w,.true.)
    end function obj_actual
    subroutine decode(z,p)
      real(dp), intent(in) :: z(:)
      real(dp), intent(out) :: p(2)
      real(dp) :: s
      p(1)=exp(min(50.0_dp,max(-50.0_dp,z(1))))
      s=1.0_dp/(1.0_dp+exp(-min(50.0_dp,max(-50.0_dp,z(2)))))
      p(2)=-1.0_dp+3.0_dp*s
    end subroutine decode
  end subroutine fit_wgpd

end module mev_stein

module mev_bias
  use mev_kinds, only: dp
  use mev_univariate, only: gpd_score, gpd_infomat, gev_score, gev_infomat
  use nleqslv_fortran, only: nleq_options, nleq_result, solve_nleqslv, NLEQ_BROYDEN, NLEQ_DBLDOG, NLEQ_NONE
  use mev_distributions, only: qgev
  use mev_math, only: inverse_matrix
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: bias_correction_result, gpd_bias, gpd_fscore, gpd_bcor
  public :: gev_bias, gev_fscore, gev_bcor

  type :: bias_correction_result
    real(dp), allocatable :: estimate(:)
    real(dp), allocatable :: residual(:)
    integer :: convergence = 1
    integer :: termcd = 0
    integer :: iterations = 0
  end type bias_correction_result

contains

  pure subroutine gpd_bias(par,n,bias)
    real(dp), intent(in) :: par(2)
    integer, intent(in) :: n
    real(dp), intent(out) :: bias(2)
    real(dp) :: den
    if(n<1 .or. par(1)<=0.0_dp .or. par(2)<=-1.0_dp/3.0_dp) then
      bias=ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    den=real(n,dp)*(1.0_dp+3.0_dp*par(2))
    bias(1)=par(1)*(3.0_dp+5.0_dp*par(2)+4.0_dp*par(2)**2)/den
    bias(2)=-(1.0_dp+par(2))*(3.0_dp+par(2))/den
  end subroutine gpd_bias

  subroutine gpd_fscore(par,dat,score,expected)
    real(dp), intent(in) :: par(2),dat(:)
    real(dp), intent(out) :: score(2)
    logical, intent(in), optional :: expected
    real(dp) :: raw(2), info(2,2), b(2)
    logical :: ex
    ex=.false.; if(present(expected)) ex=expected
    call gpd_score(par,dat,raw)
    call gpd_bias(par,size(dat),b)
    call gpd_infomat(par,dat,info,ex)
    if(any(.not.ieee_is_finite(raw)) .or. any(.not.ieee_is_finite(info)) .or. any(.not.ieee_is_finite(b))) then
      score=huge(1.0_dp)/1000.0_dp
    else
      score=raw-matmul(info,b)
    end if
  end subroutine gpd_fscore

  subroutine gpd_bcor(par,dat,result,correction,expected)
    real(dp), intent(in) :: par(2),dat(:)
    type(bias_correction_result), intent(out) :: result
    character(len=*), intent(in), optional :: correction
    logical, intent(in), optional :: expected
    type(nleq_options) :: opt
    type(nleq_result) :: sol
    real(dp) :: start(2), maxdat
    character(len=16) :: corr
    logical :: ex

    result=bias_correction_result()
    allocate(result%estimate(2),result%residual(2))
    result%estimate=ieee_value(0.0_dp,ieee_quiet_nan)
    result%residual=ieee_value(0.0_dp,ieee_quiet_nan)
    if(size(dat)<2 .or. any(dat<0.0_dp) .or. par(1)<=0.0_dp) return
    corr='subtract'; if(present(correction)) corr=trim(adjustl(correction))
    ex=.false.; if(present(expected)) ex=expected
    maxdat=maxval(dat)
    start=par; start(2)=max(start(2),-0.25_dp)
    if(.not.valid_par(start,maxdat)) then
      start(1)=max(maxval(dat)*0.75_dp,epsilon(1.0_dp))
      start(2)=max(0.05_dp,start(2))
    end if
    opt=nleq_options(); opt%method=NLEQ_BROYDEN; opt%global=NLEQ_DBLDOG
    opt%maxit=1000; opt%xtol=1.0e-10_dp; opt%ftol=1.0e-9_dp
    if(corr=='firth') then
      call solve_nleqslv(start,firth_fn,sol,opt)
    else
      call solve_nleqslv(start,subtract_fn,sol,opt)
    end if
    if(sol%termcd/=1 .and. maxval(abs(sol%fvec))>1.0e-6_dp) then
      opt%global=NLEQ_NONE
      if(corr=='firth') then
        call solve_nleqslv(start,firth_fn,sol,opt)
      else
        call solve_nleqslv(start,subtract_fn,sol,opt)
      end if
    end if
    result%termcd=sol%termcd; result%iterations=sol%iter
    if(valid_par(sol%x,maxdat) .and. maxval(abs(sol%fvec))<=1.0e-5_dp) then
      result%estimate=sol%x; result%residual=sol%fvec; result%convergence=0
    else
      result%estimate=sol%x; result%residual=sol%fvec
    end if
  contains
    subroutine subtract_fn(x,f)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f(:)
      real(dp) :: b(2)
      if(.not.valid_par(x,maxdat)) then
        call invalid_residual(x,f); return
      end if
      call gpd_bias(x(1:2),size(dat),b)
      f=x(1:2)-par+b
      if(any(.not.ieee_is_finite(f))) call invalid_residual(x,f)
    end subroutine subtract_fn

    subroutine firth_fn(x,f)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f(:)
      if(.not.valid_par(x,maxdat)) then
        call invalid_residual(x,f); return
      end if
      call gpd_fscore(x(1:2),dat,f,ex)
      if(any(.not.ieee_is_finite(f))) call invalid_residual(x,f)
    end subroutine firth_fn

    subroutine invalid_residual(x,f)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f(:)
      real(dp) :: pscale,pshape,psupp
      pscale=max(0.0_dp,1.0e-8_dp-x(1))
      pshape=max(0.0_dp,-1.0_dp/3.0_dp+1.0e-8_dp-x(2))
      psupp=max(0.0_dp,-(x(1)+x(2)*maxdat)+1.0e-8_dp)
      f(1)=1.0e3_dp*(1.0_dp+pscale+psupp)
      f(2)=1.0e3_dp*(1.0_dp+pshape+psupp)
    end subroutine invalid_residual
  end subroutine gpd_bcor


  subroutine gev_bias(par,n,bias,nquad)
    !! Cox-Snell first-order bias for the GEV MLE.
    !!
    !! The upstream package evaluates a very large closed symbolic expression.
    !! Here the same defining cumulants are evaluated deterministically.  The
    !! expected-information derivative is differentiated directly and the
    !! third log-likelihood cumulant is integrated over the exact GEV quantile
    !! transform.  A power transform clusters quadrature nodes near the upper
    !! probability endpoint, which is important for negative shape values.
    real(dp), intent(in) :: par(3)
    integer, intent(in) :: n
    real(dp), intent(out) :: bias(3)
    integer, intent(in), optional :: nquad
    real(dp) :: pm(3),pp(3),bm(3),bp(3),a
    integer :: nq
    bias=ieee_value(0.0_dp,ieee_quiet_nan)
    if(n<1 .or. par(2)<=0.0_dp .or. par(3)<=-1.0_dp/3.0_dp .or. par(3)>=1.0_dp) return
    nq=96; if(present(nquad)) nq=max(32,nquad)
    if(abs(par(3))<2.5e-3_dp) then
      pm=par;pp=par;pm(3)=-2.5e-3_dp;pp(3)=2.5e-3_dp
      call gev_bias_core(pm,n,bm,nq)
      call gev_bias_core(pp,n,bp,nq)
      a=(par(3)+2.5e-3_dp)/5.0e-3_dp
      bias=(1.0_dp-a)*bm+a*bp
    else
      call gev_bias_core(par,n,bias,nq)
    end if
  end subroutine gev_bias

  subroutine gev_bias_core(par,n,bias,nq)
    real(dp),intent(in)::par(3)
    integer,intent(in)::n,nq
    real(dp),intent(out)::bias(3)
    integer::i,t,ier,it
    real(dp),allocatable::nodes(:),weights(:)
    real(dp)::v,u,wgt,x,hp(3,3),hm(3,3),third(3,3),k3(3,3,3)
    real(dp)::info(3,3),iinv(3,3),ip(3,3),im(3,3),dk(3,3,3),rhs(3),pt(3),h
    logical::okp,okm
    real(dp)::dat1(1)
    bias=ieee_value(0.0_dp,ieee_quiet_nan);k3=0.0_dp
    allocate(nodes(nq),weights(nq));call gauss_legendre_unit(nq,nodes,weights)
    do i=1,nq
      v=nodes(i);u=1.0_dp-v**3
      wgt=weights(i)*3.0_dp*v*v
      x=qgev(u,loc=par(1),scale=par(2),shape=par(3))
      if(.not.ieee_is_finite(x)) cycle
      do t=1,3
        h=gev_fd_step(par,t);okp=.false.;okm=.false.
        do it=1,12
          pt=par;pt(t)=par(t)+h;call gev_score_jacobian(pt,x,hp,okp)
          pt=par;pt(t)=par(t)-h;call gev_score_jacobian(pt,x,hm,okm)
          if(okp.and.okm) exit
          h=0.5_dp*h
        end do
        if(.not.(okp.and.okm)) return
        third=(hp-hm)/(2.0_dp*h)
        k3(:,:,t)=k3(:,:,t)+wgt*third
      end do
    end do
    dat1(1)=par(1);call gev_infomat(par,dat1,info,.true.,1)
    do t=1,3
      h=gev_fd_step(par,t)
      pt=par;pt(t)=par(t)+h;call gev_infomat(pt,dat1,ip,.true.,1)
      pt=par;pt(t)=par(t)-h;call gev_infomat(pt,dat1,im,.true.,1)
      ! kappa_rs = E(l_rs) = -I_rs
      dk(:,:,t)=(-ip+im)/(2.0_dp*h)
    end do
    call inverse_matrix(info,iinv,ier);if(ier/=0.or.any(.not.ieee_is_finite(iinv))) return
    rhs=0.0_dp
    do t=1,3
      rhs=rhs+matmul(dk(:,:,t)-0.5_dp*k3(:,:,t),iinv(:,t))
    end do
    bias=matmul(iinv,rhs)/real(n,dp)
  end subroutine gev_bias_core

  subroutine gev_fscore(par,dat,score,expected,nquad)
    real(dp), intent(in) :: par(3),dat(:)
    real(dp), intent(out) :: score(3)
    logical, intent(in), optional :: expected
    integer, intent(in), optional :: nquad
    real(dp) :: raw(3),info(3,3),b(3)
    logical :: ex
    integer :: nq
    ex=.false.; if(present(expected)) ex=expected
    nq=64; if(present(nquad)) nq=nquad
    call gev_score(par,dat,raw)
    call gev_bias(par,size(dat),b,nq)
    call gev_infomat(par,dat,info,ex)
    if(any(.not.ieee_is_finite(raw)) .or. any(.not.ieee_is_finite(info)) .or. any(.not.ieee_is_finite(b))) then
      score=huge(1.0_dp)/1000.0_dp
    else
      score=raw-matmul(info,b)
    end if
  end subroutine gev_fscore

  subroutine gev_bcor(par,dat,result,correction,expected,nquad)
    real(dp), intent(in) :: par(3),dat(:)
    type(bias_correction_result), intent(out) :: result
    character(len=*), intent(in), optional :: correction
    logical, intent(in), optional :: expected
    integer, intent(in), optional :: nquad
    type(nleq_options) :: opt
    type(nleq_result) :: sol
    real(dp) :: start(3),xmin,xmax
    character(len=16) :: corr
    logical :: ex
    integer :: nq

    result=bias_correction_result()
    allocate(result%estimate(3),result%residual(3))
    result%estimate=ieee_value(0.0_dp,ieee_quiet_nan)
    result%residual=ieee_value(0.0_dp,ieee_quiet_nan)
    if(size(dat)<3 .or. par(2)<=0.0_dp) return
    corr='subtract'; if(present(correction)) corr=trim(adjustl(correction))
    ex=.false.; if(present(expected)) ex=expected
    nq=64; if(present(nquad)) nq=max(24,nquad)
    xmin=minval(dat);xmax=maxval(dat)
    start=par; start(3)=max(-0.25_dp,min(0.9_dp,start(3)))
    if(.not.valid_gev_par(start,xmin,xmax)) then
      start(1)=0.5_dp*(xmin+xmax)
      start(2)=max(0.25_dp*(xmax-xmin),epsilon(1.0_dp))
      start(3)=0.05_dp
    end if
    opt=nleq_options();opt%method=NLEQ_BROYDEN;opt%global=NLEQ_DBLDOG
    opt%maxit=1000;opt%xtol=1.0e-9_dp;opt%ftol=2.0e-8_dp
    if(corr=='firth') then
      call solve_nleqslv(start,firth_fn,sol,opt)
    else
      call solve_nleqslv(start,subtract_fn,sol,opt)
    end if
    if(sol%termcd/=1 .and. maxval(abs(sol%fvec))>2.0e-5_dp) then
      opt%global=NLEQ_NONE
      if(corr=='firth') then
        call solve_nleqslv(start,firth_fn,sol,opt)
      else
        call solve_nleqslv(start,subtract_fn,sol,opt)
      end if
    end if
    result%termcd=sol%termcd;result%iterations=sol%iter
    result%estimate=sol%x;result%residual=sol%fvec
    if(valid_gev_par(sol%x,xmin,xmax) .and. maxval(abs(sol%fvec))<=5.0e-4_dp) result%convergence=0
  contains
    subroutine subtract_fn(x,f)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::f(:)
      real(dp)::b(3)
      if(.not.valid_gev_par(x,xmin,xmax)) then
        call invalid_gev_residual(x,f);return
      end if
      call gev_bias(x(1:3),size(dat),b,nq)
      f=x(1:3)-par+b
      if(any(.not.ieee_is_finite(f))) call invalid_gev_residual(x,f)
    end subroutine subtract_fn
    subroutine firth_fn(x,f)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::f(:)
      if(.not.valid_gev_par(x,xmin,xmax)) then
        call invalid_gev_residual(x,f);return
      end if
      call gev_fscore(x(1:3),dat,f,ex,nq)
      if(any(.not.ieee_is_finite(f))) call invalid_gev_residual(x,f)
    end subroutine firth_fn
    subroutine invalid_gev_residual(x,f)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::f(:)
      real(dp)::pscale,pshape,psupp
      pscale=max(0.0_dp,1.0e-8_dp-x(2))
      pshape=max(0.0_dp,-1.0_dp/3.0_dp+1.0e-7_dp-x(3))+max(0.0_dp,x(3)-0.999_dp)
      psupp=max(0.0_dp,-min(1.0_dp+x(3)*(xmin-x(1))/max(x(2),1.0e-12_dp), &
                                1.0_dp+x(3)*(xmax-x(1))/max(x(2),1.0e-12_dp))+1.0e-8_dp)
      f=1.0e3_dp*(1.0_dp+pscale+pshape+psupp)
    end subroutine invalid_gev_residual
  end subroutine gev_bcor

  subroutine gev_score_point(par,x,score,ok)
    real(dp),intent(in)::par(3),x
    real(dp),intent(out)::score(3)
    logical,intent(out)::ok
    real(dp)::d(1)
    d(1)=x;call gev_score(par,d,score)
    ok=all(ieee_is_finite(score))
  end subroutine gev_score_point

  subroutine gev_score_jacobian(par,x,jac,ok)
    real(dp),intent(in)::par(3),x
    real(dp),intent(out)::jac(3,3)
    logical,intent(out)::ok
    real(dp)::pp(3),pm(3),sp(3),sm(3),h
    logical::op,om
    integer::j,it
    ok=.true.;jac=0.0_dp
    do j=1,3
      h=gev_fd_step(par,j)
      op=.false.;om=.false.
      do it=1,10
        pp=par;pm=par;pp(j)=pp(j)+h;pm(j)=pm(j)-h
        call gev_score_point(pp,x,sp,op);call gev_score_point(pm,x,sm,om)
        if(op.and.om) exit
        h=0.5_dp*h
      end do
      if(.not.(op.and.om)) then;ok=.false.;return;end if
      jac(:,j)=(sp-sm)/(2.0_dp*h)
    end do
  end subroutine gev_score_jacobian

  pure real(dp) function gev_fd_step(par,j) result(h)
    real(dp),intent(in)::par(3)
    integer,intent(in)::j
    select case(j)
    case(1);h=5.0e-5_dp*max(par(2),max(abs(par(1)),1.0_dp))
    case(2);h=5.0e-5_dp*max(par(2),1.0_dp)
    case default;h=1.0e-4_dp*max(abs(par(3)),1.0_dp)
    end select
  end function gev_fd_step

  subroutine gauss_legendre_unit(n,x,w)
    integer,intent(in)::n
    real(dp),intent(out)::x(n),w(n)
    integer::i,j,m
    real(dp)::z,z1,p1,p2,p3,pp,eps
    m=(n+1)/2;eps=5.0e-15_dp
    do i=1,m
      z=cos(acos(-1.0_dp)*(real(i,dp)-0.25_dp)/(real(n,dp)+0.5_dp))
      do
        p1=1.0_dp;p2=0.0_dp
        do j=1,n
          p3=p2;p2=p1;p1=((2.0_dp*real(j,dp)-1.0_dp)*z*p2-(real(j,dp)-1.0_dp)*p3)/real(j,dp)
        end do
        pp=real(n,dp)*(z*p1-p2)/(z*z-1.0_dp)
        z1=z;z=z1-p1/pp
        if(abs(z-z1)<=eps) exit
      end do
      x(i)=0.5_dp*(1.0_dp-z);x(n+1-i)=0.5_dp*(1.0_dp+z)
      w(i)=1.0_dp/((1.0_dp-z*z)*pp*pp);w(n+1-i)=w(i)
    end do
  end subroutine gauss_legendre_unit

  pure logical function valid_gev_par(x,xmin,xmax)
    real(dp),intent(in)::x(:),xmin,xmax
    real(dp)::a,b
    if(size(x)<3) then;valid_gev_par=.false.;return;end if
    if(x(2)<=0.0_dp .or. x(3)<=-1.0_dp/3.0_dp .or. x(3)>=1.0_dp) then
      valid_gev_par=.false.;return
    end if
    a=1.0_dp+x(3)*(xmin-x(1))/x(2);b=1.0_dp+x(3)*(xmax-x(1))/x(2)
    valid_gev_par=(a>0.0_dp.and.b>0.0_dp)
  end function valid_gev_par

  pure logical function valid_par(x,maxdat)
    real(dp), intent(in) :: x(:),maxdat
    valid_par=size(x)>=2 .and. x(1)>0.0_dp .and. x(2)>-1.0_dp/3.0_dp .and. x(1)+x(2)*maxdat>0.0_dp
  end function valid_par

end module mev_bias

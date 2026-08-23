module skewhyperbolic_fit
  use skewhyp_math, only: dp
  use skewhyperbolic_distribution, only: dskewhyp
  implicit none
  private
  type, public :: skewhyp_fit_result
    real(dp) :: param(4)=0.0_dp
    real(dp) :: loglik=-huge(1.0_dp)
    integer :: iterations=0
    integer :: convergence=1
  end type
  public :: skewhyp_fit, skewhyp_fit_start
contains
  subroutine skewhyp_fit_start(x,param)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::param(4)
    real(dp)::m,s2,s,sk
    integer::n
    n=size(x);m=sum(x)/real(n,dp)
    if(n>1)then;s2=sum((x-m)**2)/real(n-1,dp);else;s2=1.0_dp;endif
    s=sqrt(max(s2,1.0e-8_dp))
    sk=sum((x-m)**3)/real(max(n,1),dp)/max(s,1.0e-8_dp)**3
    param(4)=10.0_dp
    param(3)=sign(1.0_dp,sk)*min(abs(sk),2.0_dp)*0.15_dp/max(s,1.0e-8_dp)
    param(2)=s*sqrt(max(param(4)-2.0_dp,1.0_dp))
    param(1)=m-param(3)*param(2)**2/(param(4)-2.0_dp)
  end subroutine skewhyp_fit_start

  subroutine skewhyp_fit(x,result,start,maxiter,tol)
    real(dp),intent(in)::x(:)
    type(skewhyp_fit_result),intent(out)::result
    real(dp),intent(in),optional::start(4),tol
    integer,intent(in),optional::maxiter
    real(dp)::z(4),best,trial(4),f,step(4),eps,p0(4)
    integer::iter,j,sgn,mx
    if(present(start))then;p0=start;else;call skewhyp_fit_start(x,p0);endif
    z=[p0(1),log(max(p0(2),1.0e-8_dp)),p0(3),log(max(p0(4),1.0e-6_dp))]
    step=[max(0.2_dp,0.2_dp*abs(z(1))+0.1_dp),0.25_dp,max(0.1_dp,0.25_dp*abs(z(3))+0.05_dp),0.25_dp]
    eps=1.0e-6_dp;if(present(tol))eps=tol
    mx=800;if(present(maxiter))mx=maxiter
    best=nll(z)
    do iter=1,mx
      do j=1,4
        do sgn=-1,1,2
          trial=z;trial(j)=trial(j)+real(sgn,dp)*step(j)
          f=nll(trial)
          if(f<best)then;z=trial;best=f;endif
        enddo
      enddo
      if(maxval(step)<eps)exit
      step=step*0.92_dp
    enddo
    result%param=[z(1),exp(z(2)),z(3),exp(z(4))]
    result%loglik=-best;result%iterations=iter
    if(maxval(step)<10.0_dp*eps)then;result%convergence=0;else;result%convergence=1;endif
  contains
    function nll(v) result(val)
      real(dp),intent(in)::v(4)
      real(dp)::val,p(4),ld
      integer::i
      p=[v(1),exp(v(2)),v(3),exp(v(4))]
      if(p(4)<0.05_dp .or. p(2)<1.0e-8_dp)then;val=huge(1.0_dp);return;endif
      val=0.0_dp
      do i=1,size(x)
        ld=dskewhyp(x(i),p(1),p(2),p(3),p(4),log_density=.true.)
        if(.not.(ld>-huge(1.0_dp)))then;val=huge(1.0_dp);return;endif
        val=val-ld
      enddo
    end function nll
  end subroutine skewhyp_fit
end module skewhyperbolic_fit

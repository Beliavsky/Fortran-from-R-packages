module trawl_fit
  use trawl_kinds, only : dp
  use trawl_types, only : trawl_fit_result,poisson_fit_result,nb_fit_result, &
    trawl_ok,trawl_invalid_argument,trawl_optimization_failed
  use trawl_statistics, only : sample_mean,sample_variance,empirical_acf
  use trawl_functions, only : acf_exp,acf_dexp,acf_supig,acf_lm
  use trawl_optimize, only : differential_evolution
  implicit none
  private
  public :: fit_exptrawl,fit_supigtrawl,fit_lmtrawl,fit_dexptrawl
  public :: fit_marginal_poisson,fit_marginal_nb
contains
  function fit_marginal_poisson(x,lm) result(res)
    real(dp),intent(in)::x(:),lm
    type(poisson_fit_result)::res
    if(size(x)==0 .or. lm==0.0_dp) then
      res%status=trawl_invalid_argument; return
    end if
    res%v=sample_mean(x)/lm
  end function

  function fit_marginal_nb(x,lm) result(res)
    real(dp),intent(in)::x(:),lm
    type(nb_fit_result)::res
    real(dp)::mu,vx
    if(size(x)<2 .or. lm==0.0_dp) then
      res%status=trawl_invalid_argument; return
    end if
    mu=sample_mean(x); vx=sample_variance(x)
    if(vx==0.0_dp) then; res%status=trawl_invalid_argument; return; end if
    res%theta=1.0_dp-mu/vx
    if(res%theta==0.0_dp) then; res%status=trawl_invalid_argument; return; end if
    res%m=mu/lm*(1.0_dp-res%theta)/res%theta
    res%a=res%theta/(1.0_dp-res%theta)
  end function

  function fit_exptrawl(x,delta_t) result(res)
    real(dp),intent(in)::x(:)
    real(dp),intent(in),optional::delta_t
    type(trawl_fit_result)::res
    real(dp),allocatable::ac(:)
    real(dp)::dt
    dt=1.0_dp; if(present(delta_t)) dt=delta_t
    if(size(x)<2 .or. dt<=0.0_dp) then; res%status=trawl_invalid_argument; return; end if
    call empirical_acf(x,1,ac)
    if(size(ac)<1 .or. ac(1)<=0.0_dp) then; res%status=trawl_invalid_argument; return; end if
    res%lambda1=-log(ac(1))/dt
    if(res%lambda1==0.0_dp) then
      res%lm=huge(1.0_dp)
    else
      res%lm=1.0_dp/res%lambda1
    end if
  end function

  function fit_supigtrawl(x,delta_t,gmm_lag,itermax) result(res)
    real(dp),intent(in)::x(:)
    real(dp),intent(in),optional::delta_t
    integer,intent(in),optional::gmm_lag,itermax
    type(trawl_fit_result)::res
    real(dp),allocatable::ac(:)
    real(dp)::dt,best(2),val,lo(2),hi(2)
    integer::lag,it,s
    dt=1.0_dp; if(present(delta_t)) dt=delta_t
    lag=5; if(present(gmm_lag)) lag=gmm_lag
    it=1000; if(present(itermax)) it=itermax
    if(size(x)<=lag .or. lag<1 .or. dt<=0.0_dp) then; res%status=trawl_invalid_argument; return; end if
    call empirical_acf(x,lag,ac); lo=[0.0_dp,0.0_dp]; hi=[100.0_dp,100.0_dp]
    call differential_evolution(obj,lo,hi,best,val,itermax=it,status=s)
    res%delta=best(1); res%gamma=best(2); res%objective=val
    if(res%delta==0.0_dp) then; res%lm=huge(1.0_dp); else; res%lm=res%gamma/res%delta; end if
    if(s/=0) res%status=trawl_optimization_failed
  contains
    function obj(y) result(f)
      real(dp),intent(in)::y(:); real(dp)::f; integer::i
      if(y(1)<0.0_dp .or. y(2)<=sqrt(tiny(1.0_dp))) then; f=huge(1.0_dp); return; end if
      f=0.0_dp
      do i=1,lag; f=f+(ac(i)-acf_supig(dt*real(i,dp),y(1),y(2)))**2; end do
    end function
  end function

  function fit_lmtrawl(x,delta_t,gmm_lag,itermax) result(res)
    real(dp),intent(in)::x(:)
    real(dp),intent(in),optional::delta_t
    integer,intent(in),optional::gmm_lag,itermax
    type(trawl_fit_result)::res
    real(dp),allocatable::ac(:)
    real(dp)::dt,best(2),val,lo(2),hi(2)
    integer::lag,it,s
    dt=1.0_dp; if(present(delta_t)) dt=delta_t
    lag=5; if(present(gmm_lag)) lag=gmm_lag
    it=1000; if(present(itermax)) it=itermax
    if(size(x)<=lag .or. lag<1 .or. dt<=0.0_dp) then; res%status=trawl_invalid_argument; return; end if
    call empirical_acf(x,lag,ac); lo=[0.0_dp,0.0_dp]; hi=[100.0_dp,100.0_dp]
    call differential_evolution(obj,lo,hi,best,val,itermax=it,status=s)
    res%alpha=best(1); res%h=best(2); res%objective=val
    if(abs(res%h-1.0_dp)<=epsilon(1.0_dp)) then
      res%lm=sign(huge(1.0_dp),res%h-1.0_dp)
    else
      res%lm=res%alpha/(res%h-1.0_dp)
    end if
    if(s/=0) res%status=trawl_optimization_failed
  contains
    function obj(y) result(f)
      real(dp),intent(in)::y(:); real(dp)::f; integer::i
      if(y(1)<=sqrt(tiny(1.0_dp))) then; f=huge(1.0_dp); return; end if
      f=0.0_dp
      do i=1,lag; f=f+(ac(i)-acf_lm(dt*real(i,dp),y(1),y(2)))**2; end do
    end function
  end function

  function fit_dexptrawl(x,delta_t,gmm_lag,itermax,preserve_upstream_delta_bug) result(res)
    real(dp),intent(in)::x(:)
    real(dp),intent(in),optional::delta_t
    integer,intent(in),optional::gmm_lag,itermax
    logical,intent(in),optional::preserve_upstream_delta_bug
    type(trawl_fit_result)::res
    real(dp),allocatable::ac(:)
    real(dp)::dt,fit_dt,best(3),val,lo(3),hi(3)
    integer::lag,it,s
    logical::bug
    dt=1.0_dp; if(present(delta_t)) dt=delta_t
    bug=.true.; if(present(preserve_upstream_delta_bug)) bug=preserve_upstream_delta_bug
    fit_dt=merge(1.0_dp,dt,bug)
    lag=5; if(present(gmm_lag)) lag=gmm_lag
    it=1000; if(present(itermax)) it=itermax
    if(size(x)<=lag .or. lag<1 .or. dt<=0.0_dp) then; res%status=trawl_invalid_argument; return; end if
    call empirical_acf(x,lag,ac); lo=[0.0_dp,0.0_dp,0.0_dp]; hi=[0.5_dp,100.0_dp,100.0_dp]
    call differential_evolution(obj,lo,hi,best,val,itermax=it,status=s)
    res%w=best(1); res%lambda1=best(2); res%lambda2=best(3); res%objective=val
    if(res%lambda1<=0.0_dp .or. res%lambda2<=0.0_dp) then
      res%lm=huge(1.0_dp)
    else
      res%lm=res%w/res%lambda1+(1.0_dp-res%w)/res%lambda2
    end if
    if(s/=0) res%status=trawl_optimization_failed
  contains
    function obj(y) result(f)
      real(dp),intent(in)::y(:); real(dp)::f; integer::i
      if(y(2)<=sqrt(tiny(1.0_dp)) .or. y(3)<=sqrt(tiny(1.0_dp))) then; f=huge(1.0_dp); return; end if
      f=0.0_dp
      do i=1,lag; f=f+(ac(i)-acf_dexp(fit_dt*real(i,dp),y(1),y(2),y(3)))**2; end do
    end function
  end function
end module trawl_fit

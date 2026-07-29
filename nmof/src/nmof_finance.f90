! SPDX-License-Identifier: GPL-3.0-only
module nmof_finance
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use nmof_kinds, only: dp, pi
   use nmof_math, only: normal_cdf, normal_pdf, brent_root, log_factorial
   use nmof_utilities, only: xw_gauss
   use nmof_types, only: option_result, bond_return_result, quadrature_rule, &
                         nmof_ok, nmof_invalid_input, nmof_numerical_failure
   implicit none
   private
   public :: vanilla_bond, yield_to_maturity, yield_to_maturity_curve, bond_duration, bond_convexity, approximate_bond_return
   public :: vanilla_option_european, vanilla_option_american, vanilla_option_implied_vol
   public :: put_call_parity, barrier_option_european, european_call_tree, european_call_binomial_expectation
   public :: cf_heston, cf_bsm, cf_bates, cf_merton, cf_variance_gamma
   public :: call_cf, call_heston_cf, call_merton
   public :: bund_future, bund_future_implied_rate, xt_contract_value, xt_tick_value
   public :: cf_callback

   abstract interface
      function cf_callback(omega, spot, tau, r, q, context) result(phi)
         import dp
         complex(dp), intent(in) :: omega
         real(dp), intent(in) :: spot, tau, r, q
         class(*), intent(in), optional :: context
         complex(dp) :: phi
      end function cf_callback
   end interface
contains
   pure real(dp) function vanilla_bond(cashflows,times,discount_factors,yields) result(price)
      real(dp),intent(in)::cashflows(:)
      real(dp),intent(in),optional::times(:),discount_factors(:),yields(:)
      real(dp),allocatable::tm(:),df(:)
      integer::i,n
      n=size(cashflows); allocate(tm(n),df(n)); tm=[(real(i,dp),i=1,n)]
      if(present(times)) tm=times
      if(present(discount_factors)) then
         df=discount_factors
      else if(present(yields)) then
         df=1.0_dp/(1.0_dp+yields)**tm
      else
         df=1.0_dp
      end if
      price=dot_product(cashflows,df)
   end function vanilla_bond

   function yield_to_maturity(cashflows,times,y0,tol,maxiter,offset,status) result(y)
      real(dp),intent(in)::cashflows(:),times(:)
      real(dp),intent(in),optional::y0,tol,offset
      integer,intent(in),optional::maxiter
      integer,intent(out),optional::status
      real(dp)::y,eps,off,y1,g,dg,dr
      integer::i,mit
      y=0.05_dp; if(present(y0)) y=y0
      eps=1.0e-5_dp; if(present(tol)) eps=tol
      off=0.0_dp; if(present(offset)) off=offset
      mit=1000; if(present(maxiter)) mit=maxiter
      if(size(cashflows)/=size(times).or.any(times<0.0_dp).or.size(times)<1) then
         y=ieee_value(y,ieee_quiet_nan); if(present(status)) status=nmof_invalid_input; return
      end if
      do i=1,mit
         y1=1.0_dp+off+y
         if(y1<=0.0_dp) then; y=0.5_dp*(y-off); cycle; end if
         g=sum(cashflows/y1**times)
         dg=sum(times*cashflows*y1**(times-1.0_dp))
         if(abs(dg)<=tiny(1.0_dp)) exit
         dr=g/dg; y=y+dr
         if(abs(dr)<eps) then; if(present(status)) status=nmof_ok; return; end if
      end do
      if(present(status)) status=nmof_numerical_failure
   end function yield_to_maturity

   function yield_to_maturity_curve(cashflows,times,offsets,y0,tol,maxiter,status) result(y)
      real(dp),intent(in)::cashflows(:),times(:),offsets(:)
      real(dp),intent(in),optional::y0,tol
      integer,intent(in),optional::maxiter
      integer,intent(out),optional::status
      real(dp)::y,eps,y1(size(cashflows)),g,dg,dr
      integer::i,mit
      y=0.05_dp; if(present(y0)) y=y0
      eps=1.0e-5_dp; if(present(tol)) eps=tol
      mit=1000; if(present(maxiter)) mit=maxiter
      if(size(cashflows)/=size(times).or.(size(offsets)/=1.and.size(offsets)/=size(times)).or. &
         any(times<0.0_dp).or.size(times)<1) then
         y=ieee_value(y,ieee_quiet_nan); if(present(status)) status=nmof_invalid_input; return
      end if
      do i=1,mit
         if(size(offsets)==1) then
            y1=1.0_dp+offsets(1)+y
         else
            y1=1.0_dp+offsets+y
         end if
         if(any(y1<=0.0_dp)) then; y=0.5_dp*y; cycle; end if
         g=sum(cashflows/y1**times)
         dg=sum(times*cashflows*y1**(times-1.0_dp))
         if(abs(dg)<=tiny(1.0_dp)) exit
         dr=g/dg; y=y+dr
         if(abs(dr)<eps) then; if(present(status)) status=nmof_ok; return; end if
      end do
      if(present(status)) status=nmof_numerical_failure
   end function yield_to_maturity_curve

   pure real(dp) function bond_duration(cashflows,times,yield_value,modified,raw) result(value)
      real(dp),intent(in)::cashflows(:),times(:),yield_value
      logical,intent(in),optional::modified,raw
      real(dp)::y1,denom
      real(dp)::dcf(size(cashflows))
      logical::modi,rawv
      modi=.true.; if(present(modified)) modi=modified; rawv=.false.; if(present(raw)) rawv=raw
      y1=1.0_dp+yield_value; dcf=cashflows/y1**times
      if(rawv) then; value=-sum(times*dcf)/y1
      else
         denom=sum(dcf); value=sum(times*dcf)/denom; if(modi) value=value/y1
      end if
   end function bond_duration

   pure real(dp) function bond_convexity(cashflows,times,yield_value,raw) result(value)
      real(dp),intent(in)::cashflows(:),times(:),yield_value
      logical,intent(in),optional::raw
      real(dp)::y1,dcf(size(cashflows)); logical::r
      r=.false.; if(present(raw)) r=raw; y1=1.0_dp+yield_value; dcf=cashflows/y1**times
      value=sum((times+1.0_dp)*times*dcf)/(y1*y1); if(.not.r) value=value/sum(dcf)
   end function bond_convexity

   function approximate_bond_return(yields,tm,coupons_per_year,scale,pad,use_pad,status) result(ans)
      real(dp),intent(in)::yields(:),tm
      integer,intent(in),optional::coupons_per_year
      real(dp),intent(in),optional::scale,pad
      logical,intent(in),optional::use_pad
      integer,intent(out),optional::status
      type(bond_return_result)::ans
      integer::nper,n,i,start
      real(dp)::sc,dy
      logical::up
      nper=2; if(present(coupons_per_year)) nper=coupons_per_year
      sc=1.0_dp/250.0_dp; if(present(scale)) sc=scale; up=.false.; if(present(use_pad)) up=use_pad
      n=size(yields)
      if(n<2.or.any(abs(yields)<=tiny(1.0_dp)).or.nper<1) then; ans%status=nmof_invalid_input; if(present(status)) status=ans%status; return; end if
      allocate(ans%duration(n),ans%convexity(n),ans%returns(merge(n,n-1,up)))
      ans%duration=1.0_dp/yields*(1.0_dp-1.0_dp/(1.0_dp+yields/real(nper,dp))**(real(nper,dp)*tm))
      ans%convexity=2.0_dp/(yields*yields)*(1.0_dp-1.0_dp/(1.0_dp+yields/real(nper,dp))**(real(nper,dp)*tm))- &
         (2.0_dp*tm)/(yields*(1.0_dp+yields/real(nper,dp))**(real(nper,dp)*tm+1.0_dp))
      start=1
      if(up) then; ans%returns(1)=0.0_dp; if(present(pad)) ans%returns(1)=pad; start=2; end if
      do i=1,n-1
         dy=yields(i+1)-yields(i)
         ans%returns(start+i-1)=(yields(i)+1.0_dp)**sc-dy*ans%duration(i)+0.5_dp*ans%convexity(i)*dy*dy-1.0_dp
      end do
      ans%status=nmof_ok; if(present(status)) status=ans%status
   end function approximate_bond_return

   function vanilla_option_european(spot,strike,tau,r,q,variance,option_type,dividend_times,dividends,greeks,status) result(ans)
      real(dp),intent(in)::spot,strike,tau,r,q,variance
      character(len=*),intent(in),optional::option_type
      real(dp),intent(in),optional::dividend_times(:),dividends(:)
      logical,intent(in),optional::greeks
      integer,intent(out),optional::status
      type(option_result)::ans
      real(dp)::sadj,d1,d2,n1,n2,exq,exr,sexq,xexr,sgn,pvdiv
      logical::g
      character(len=8)::typ
      g=.true.; if(present(greeks)) g=greeks; typ='call'; if(present(option_type)) typ=lowercase(option_type)
      if(tau<=0.0_dp.or.variance<=0.0_dp.or.spot<=0.0_dp.or.strike<=0.0_dp) then
         ans%status=nmof_invalid_input; if(present(status)) status=ans%status; return
      end if
      pvdiv=present_value_dividends(r,tau,q,dividend_times,dividends,ans%status)
      if(ans%status/=nmof_ok) then; if(present(status)) status=ans%status; return; end if
      sadj=spot-pvdiv; if(sadj<=0.0_dp) then; ans%status=nmof_invalid_input; if(present(status)) status=ans%status; return; end if
      sgn=merge(1.0_dp,-1.0_dp,index(typ,'call')==1); exq=exp(-q*tau); exr=exp(-r*tau)
      sexq=sadj*exq; xexr=strike*exr
      d1=(log(sadj/strike)+(r-q+0.5_dp*variance)*tau)/sqrt(variance*tau); d2=d1-sqrt(variance*tau)
      n1=normal_cdf(sgn*d1); n2=normal_cdf(sgn*d2)
      ans%value=sgn*(sexq*n1-xexr*n2)
      if(g) then
         ans%delta=sgn*exq*n1; ans%theta=-sexq*normal_pdf(d1)*sqrt(variance)/(2.0_dp*sqrt(tau))- &
            sgn*(-q*sexq*n1+r*xexr*n2)
         ans%rho=sgn*tau*xexr*n2; ans%rho_div=-sgn*tau*sexq*n1
         ans%gamma=exq*normal_pdf(sgn*d1)/(sadj*sqrt(variance*tau))
         ans%vega=sexq*normal_pdf(sgn*d1)*sqrt(tau)
         ans%dvega_dspot=ans%vega/sadj*(1.0_dp-d1/sqrt(variance*tau))
         ans%dvega_dvol=ans%vega*d1*d2/sqrt(variance)
      end if
      ans%status=nmof_ok; if(present(status)) status=ans%status
   end function vanilla_option_european

   function vanilla_option_american(spot,strike,tau,r,q,variance,option_type,steps,dividend_times,dividends,greeks,status) result(ans)
      real(dp),intent(in)::spot,strike,tau,r,q,variance
      character(len=*),intent(in),optional::option_type
      integer,intent(in),optional::steps
      real(dp),intent(in),optional::dividend_times(:),dividends(:)
      logical,intent(in),optional::greeks
      integer,intent(out),optional::status
      type(option_result)::ans
      real(dp),allocatable::w(:),wnew(:),si(:)
      real(dp)::sadj,dt,u,d,p,v1,v2,sgn,pvdiv,tval,delta0,gamma0,theta0
      integer::m,i,j
      logical::g
      character(len=8)::typ
      m=101; if(present(steps)) m=steps; g=.true.; if(present(greeks)) g=greeks
      typ='call'; if(present(option_type)) typ=lowercase(option_type)
      if(tau<=0.0_dp.or.variance<=0.0_dp.or.m<1) then; ans%status=nmof_invalid_input; if(present(status)) status=ans%status; return; end if
      pvdiv=present_value_dividends(r,tau,q,dividend_times,dividends,ans%status)
      if(ans%status/=nmof_ok) then; if(present(status)) status=ans%status; return; end if
      sadj=spot-pvdiv; dt=tau/real(m,dp); u=exp(sqrt(variance*dt)); d=1.0_dp/u
      p=(exp((r-q)*dt)-d)/(u-d); v1=p*exp(-r*dt); v2=(1.0_dp-p)*exp(-r*dt)
      sgn=merge(1.0_dp,-1.0_dp,index(typ,'call')==1)
      allocate(w(m+1)); do j=0,m; w(j+1)=max(sgn*(sadj*d**real(m-j,dp)*u**real(j,dp)-strike),0.0_dp); end do
      delta0=ieee_value(0.0_dp,ieee_quiet_nan)
      gamma0=ieee_value(0.0_dp,ieee_quiet_nan)
      theta0=ieee_value(0.0_dp,ieee_quiet_nan)
      do i=m,1,-1
         allocate(wnew(i),si(i)); tval=real(i-1,dp)*dt
         do j=0,i-1
            si(j+1)=sadj*d**real(i-1-j,dp)*u**real(j,dp)
            wnew(j+1)=max(sgn*(si(j+1)+future_dividend_value(r,tval,tau,dividend_times,dividends)-strike), &
                          v1*w(j+2)+v2*w(j+1),0.0_dp)
         end do
         if(g) then
            if(i==2) delta0=(wnew(2)-wnew(1))/(si(2)-si(1))
            if(i==3) then
               gamma0=((wnew(3)-wnew(2))/(si(3)-si(2))-(wnew(2)-wnew(1))/(si(2)-si(1)))/(0.5_dp*(si(3)-si(1)))
               theta0=wnew(2)
            end if
            if(i==1) theta0=(theta0-wnew(1))/(2.0_dp*dt)
         end if
         call move_alloc(wnew,w); deallocate(si)
      end do
      ans%value=w(1); ans%delta=delta0; ans%gamma=gamma0; ans%theta=theta0
      ans%vega=ieee_value(0.0_dp,ieee_quiet_nan); ans%rho=ans%vega; ans%rho_div=ans%vega
      ans%status=nmof_ok; if(present(status)) status=ans%status
   end function vanilla_option_american

   function vanilla_option_implied_vol(exercise,price,spot,strike,tau,r,q,option_type,steps,dividend_times,dividends,lower,upper,tol,status) result(vol)
      character(len=*),intent(in)::exercise
      real(dp),intent(in)::price,spot,strike,tau,r,q
      character(len=*),intent(in),optional::option_type
      integer,intent(in),optional::steps
      real(dp),intent(in),optional::dividend_times(:),dividends(:),lower,upper,tol
      integer,intent(out),optional::status
      real(dp)::vol,lo,hi,eps
      integer::info
      lo=1.0e-5_dp; hi=2.0_dp; eps=epsilon(1.0_dp)**0.25_dp
      if(present(lower)) lo=lower; if(present(upper)) hi=upper; if(present(tol)) eps=tol
      call brent_root(price_difference,lo,hi,vol,info,eps,1000)
      if(present(status)) status=merge(nmof_ok,nmof_numerical_failure,info==0)
   contains
      function price_difference(vv) result(f)
         real(dp),intent(in)::vv; real(dp)::f; type(option_result)::o
         if(index(lowercase(exercise),'european')==1) then
            o=vanilla_option_european(spot,strike,tau,r,q,vv*vv,option_type,dividend_times,dividends,.false.)
         else
            o=vanilla_option_american(spot,strike,tau,r,q,vv*vv,option_type,steps,dividend_times,dividends,.false.)
         end if
         f=price-o%value
      end function price_difference
   end function vanilla_option_implied_vol

   function put_call_parity(what,call_value,put_value,spot,strike,tau,r,q,dividend_times,dividends,status) result(value)
      character(len=*),intent(in)::what
      real(dp),intent(in)::call_value,put_value,spot,strike,tau,r,q
      real(dp),intent(in),optional::dividend_times(:),dividends(:)
      integer,intent(out),optional::status
      real(dp)::value,pvdiv
      integer::st
      pvdiv=present_value_dividends(r,tau,q,dividend_times,dividends,st)
      if(st/=nmof_ok) then; value=ieee_value(0.0_dp,ieee_quiet_nan); if(present(status)) status=st; return; end if
      if(index(lowercase(what),'call')==1) then
         value=put_value+spot*exp(-q*tau)-pvdiv-strike*exp(-r*tau)
      else
         value=call_value+strike*exp(-r*tau)-spot*exp(-q*tau)+pvdiv
      end if
      if(present(status)) status=nmof_ok
   end function put_call_parity

   function barrier_option_european(spot,strike,barrier,tau,r,q,variance,option_type,barrier_type,rebate,status) result(value)
      real(dp),intent(in)::spot,strike,barrier,tau,r,q,variance
      character(len=*),intent(in),optional::option_type,barrier_type
      real(dp),intent(in),optional::rebate
      integer,intent(out),optional::status
      real(dp)::value,kreb,mu,lambda,z,x1,x2,y1,y2,a,b,c,d,e,f,nu,phi,svt
      character(len=16)::typ,bt
      typ='call'; bt='downin'; if(present(option_type)) typ=lowercase(option_type); if(present(barrier_type)) bt=lowercase(barrier_type)
      kreb=0.0_dp; if(present(rebate)) kreb=rebate
      if(spot<=0.or.strike<=0.or.barrier<=0.or.tau<=0.or.variance<=0) then; value=0; if(present(status)) status=nmof_invalid_input; return; end if
      if(index(typ,'call')==1) then; phi=1.0_dp; else; phi=-1.0_dp; end if
      if(index(bt,'down')==1) then; nu=1.0_dp; else; nu=-1.0_dp; end if
      mu=(r-q-variance/2.0_dp)/variance; lambda=sqrt(mu*mu+2.0_dp*r/variance); svt=sqrt(variance*tau)
      z=log(barrier/spot)/svt+lambda*svt; x1=log(spot/strike)/svt+(1.0_dp+mu)*svt
      x2=log(spot/barrier)/svt+(1.0_dp+mu)*svt; y1=log(barrier*barrier/(spot*strike))/svt+(1.0_dp+mu)*svt
      y2=log(barrier/spot)/svt+(1.0_dp+mu)*svt
      a=phi*spot*exp(-q*tau)*normal_cdf(phi*x1)-phi*strike*exp(-r*tau)*normal_cdf(phi*x1-phi*svt)
      b=phi*spot*exp(-q*tau)*normal_cdf(phi*x2)-phi*strike*exp(-r*tau)*normal_cdf(phi*x2-phi*svt)
      c=phi*spot*exp(-q*tau)*(barrier/spot)**(2*(mu+1))*normal_cdf(nu*y1)- &
        phi*strike*exp(-r*tau)*(barrier/spot)**(2*mu)*normal_cdf(nu*y1-nu*svt)
      d=phi*spot*exp(-q*tau)*(barrier/spot)**(2*(mu+1))*normal_cdf(nu*y2)- &
        phi*strike*exp(-r*tau)*(barrier/spot)**(2*mu)*normal_cdf(nu*y2-nu*svt)
      e=kreb*exp(-r*tau)*(normal_cdf(nu*x2-nu*svt)-(barrier/spot)**(2*mu)*normal_cdf(nu*y2-nu*svt))
      f=kreb*((barrier/spot)**(mu+lambda)*normal_cdf(nu*z)+(barrier/spot)**(mu-lambda)*normal_cdf(nu*z-2*nu*lambda*svt))
      if(index(typ,'call')==1) then
         select case(trim(bt))
         case('downin'); if(strike>barrier) then; value=c+e; else; value=a-b+d+e; end if
         case('upin'); if(strike>barrier) then; value=a+e; else; value=b-c+d+e; end if
         case('downout'); if(strike>barrier) then; value=a-c+f; else; value=b-d+f; end if
         case('upout'); if(strike>barrier) then; value=f; else; value=a-b+c-d+f; end if
         case default; value=0.0_dp
         end select
      else
         select case(trim(bt))
         case('downin'); if(strike>barrier) then; value=b-c+d+e; else; value=a+e; end if
         case('upin'); if(strike>barrier) then; value=a-b+d+e; else; value=c+e; end if
         case('downout'); if(strike>barrier) then; value=a-b+c-d+f; else; value=f; end if
         case('upout'); if(strike>barrier) then; value=b-d+f; else; value=a-c+f; end if
         case default; value=0.0_dp
         end select
      end if
      if(present(status)) status=nmof_ok
   end function barrier_option_european

   pure real(dp) function european_call_tree(s0,strike,r,tau,sigma,steps) result(value)
      real(dp),intent(in)::s0,strike,r,tau,sigma
      integer,intent(in),optional::steps
      real(dp),allocatable::c(:),next(:)
      real(dp)::dt,disc,u,d,p
      integer::m,i,j
      m=101; if(present(steps)) m=steps; dt=tau/real(m,dp); disc=exp(-r*dt); u=exp(sigma*sqrt(dt)); d=1/u
      p=(exp(r*dt)-d)/(u-d); allocate(c(m+1)); do j=0,m; c(j+1)=max(s0*d**real(m-j,dp)*u**real(j,dp)-strike,0.0_dp); end do
      do i=m,1,-1; allocate(next(i)); do j=1,i; next(j)=disc*(p*c(j+1)+(1-p)*c(j)); end do; call move_alloc(next,c); end do
      value=c(1)
   end function european_call_tree

   pure real(dp) function european_call_binomial_expectation(s0,strike,r,tau,sigma,steps) result(value)
      real(dp),intent(in)::s0,strike,r,tau,sigma
      integer,intent(in),optional::steps
      real(dp)::dt,u,d,p,logprob,pay
      integer::m,j
      m=101; if(present(steps)) m=steps; dt=tau/real(m,dp); u=exp(sigma*sqrt(dt)); d=1/u; p=(exp(r*dt)-d)/(u-d)
      value=0.0_dp
      do j=0,m
         pay=max(s0*d**real(m-j,dp)*u**real(j,dp)-strike,0.0_dp)
         logprob=log_factorial(m)-log_factorial(j)-log_factorial(m-j)+real(j,dp)*log(p)+real(m-j,dp)*log(1-p)
         value=value+exp(logprob)*pay
      end do
      value=exp(-r*tau)*value
   end function european_call_binomial_expectation

   pure complex(dp) function cf_heston(omega,spot,tau,r,q,v0,vlong,rho,kappa,sigma) result(phi)
      complex(dp),intent(in)::omega
      real(dp),intent(in)::spot,tau,r,q,v0,vlong,rho,kappa,sigma
      complex(dp)::ii,d,g,c1,c2,c3
      real(dp)::sig
      ii=cmplx(0.0_dp,1.0_dp,dp); sig=max(sigma,1.0e-8_dp)
      d=sqrt((rho*sig*ii*omega-kappa)**2+sig*sig*(ii*omega+omega*omega))
      g=(kappa-rho*sig*ii*omega-d)/(kappa-rho*sig*ii*omega+d)
      c1=ii*omega*(log(spot)+(r-q)*tau)
      c2=vlong*kappa/(sig*sig)*((kappa-rho*sig*ii*omega-d)*tau-2.0_dp*log((1.0_dp-g*exp(-d*tau))/(1.0_dp-g)))
      c3=v0/(sig*sig)*(kappa-rho*sig*ii*omega-d)*(1.0_dp-exp(-d*tau))/(1.0_dp-g*exp(-d*tau))
      phi=exp(c1+c2+c3)
   end function cf_heston

   pure complex(dp) function cf_bsm(omega,spot,tau,r,q,variance) result(phi)
      complex(dp),intent(in)::omega
      real(dp),intent(in)::spot,tau,r,q,variance
      complex(dp)::ii; ii=cmplx(0.0_dp,1.0_dp,dp)
      phi=exp(ii*omega*log(spot)+ii*tau*(r-q)*omega-0.5_dp*tau*variance*(ii*omega+omega*omega))
   end function cf_bsm

   pure complex(dp) function cf_bates(omega,spot,tau,r,q,v0,vlong,rho,kappa,sigma,lambda,mu_j,v_j) result(phi)
      complex(dp),intent(in)::omega
      real(dp),intent(in)::spot,tau,r,q,v0,vlong,rho,kappa,sigma,lambda,mu_j,v_j
      complex(dp)::base,ii,omii,jump
      ii=cmplx(0.0_dp,1.0_dp,dp); omii=omega*ii
      base=cf_heston(omega,spot,tau,r,q,v0,vlong,rho,kappa,max(sigma,1.0e-4_dp))
      jump=-lambda*mu_j*omii*tau+lambda*tau*((1.0_dp+mu_j)**omii*exp(v_j*(omii/2.0_dp)*(omii-1.0_dp))-1.0_dp)
      phi=base*exp(jump)
   end function cf_bates

   pure complex(dp) function cf_merton(omega,spot,tau,r,q,variance,lambda,mu_j,v_j) result(phi)
      complex(dp),intent(in)::omega
      real(dp),intent(in)::spot,tau,r,q,variance,lambda,mu_j,v_j
      complex(dp)::ii,omii,c1,c2
      ii=cmplx(0.0_dp,1.0_dp,dp); omii=omega*ii
      c1=omii*log(spot)+omii*tau*(r-q-0.5_dp*variance-lambda*mu_j)-0.5_dp*omega*omega*variance*tau
      c2=lambda*tau*(exp(omii*log(1.0_dp+mu_j)-0.5_dp*omii*v_j-0.5_dp*v_j*omega*omega)-1.0_dp)
      phi=exp(c1+c2)
   end function cf_merton

   pure complex(dp) function cf_variance_gamma(omega,spot,tau,r,q,nu,theta,sigma) result(phi)
      complex(dp),intent(in)::omega
      real(dp),intent(in)::spot,tau,r,q,nu,theta,sigma
      complex(dp)::ii,omii,temp
      real(dp)::w,sigma2
      ii=cmplx(0.0_dp,1.0_dp,dp); omii=omega*ii; sigma2=sigma*sigma
      w=log(1.0_dp-theta*nu-0.5_dp*nu*sigma2)/nu
      temp=exp(omii*log(spot)+omii*(r-q+w)*tau)
      phi=temp/(1.0_dp-omii*theta*nu+0.5_dp*sigma2*nu*omega*omega)**(tau/nu)
   end function cf_variance_gamma

   function call_cf(cf,spot,strike,tau,r,q,context,implied_vol,n_quad,status) result(value)
      procedure(cf_callback)::cf
      real(dp),intent(in)::spot,strike,tau,r,q
      class(*),intent(in),optional::context
      real(dp),intent(out),optional::implied_vol
      integer,intent(in),optional::n_quad
      integer,intent(out),optional::status
      real(dp)::value,p1,p2,omega,t,jac,integ1,integ2,vol
      complex(dp)::ii,phi1,phi2
      type(quadrature_rule)::rule
      integer::i,nq,st
      nq=512; if(present(n_quad)) nq=n_quad; ii=cmplx(0.0_dp,1.0_dp,dp); rule=xw_gauss(nq,'legendre')
      if(rule%status/=nmof_ok) then; value=0.0_dp; if(present(status)) status=rule%status; return; end if
      call transform_rule(rule%nodes,rule%weights,-1.0_dp,1.0_dp,0.0_dp,1.0_dp)
      integ1=0.0_dp; integ2=0.0_dp
      do i=1,nq
         t=rule%nodes(i); omega=t/(1.0_dp-t); jac=1.0_dp/(1.0_dp-t)**2
         if(present(context)) then
            phi1=cf(cmplx(omega,-1.0_dp,dp),spot,tau,r,q,context); phi2=cf(cmplx(omega,0.0_dp,dp),spot,tau,r,q,context)
         else
            phi1=cf(cmplx(omega,-1.0_dp,dp),spot,tau,r,q); phi2=cf(cmplx(omega,0.0_dp,dp),spot,tau,r,q)
         end if
         integ1=integ1+rule%weights(i)*real(exp(-ii*log(strike)*omega)*phi1/(ii*omega*spot*exp((r-q)*tau)),dp)*jac
         integ2=integ2+rule%weights(i)*real(exp(-ii*log(strike)*omega)*phi2/(ii*omega),dp)*jac
      end do
      p1=0.5_dp+integ1/pi; p2=0.5_dp+integ2/pi; value=exp(-q*tau)*spot*p1-exp(-r*tau)*strike*p2
      if(present(implied_vol)) then
         vol=vanilla_option_implied_vol('european',value,spot,strike,tau,r,q,'call',lower=1e-5_dp,upper=2.0_dp,status=st)
         implied_vol=vol
      end if
      if(present(status)) status=nmof_ok
   contains
      subroutine transform_rule(nodes,weights,oldmin,oldmax,newmin,newmax)
         real(dp),intent(inout)::nodes(:),weights(:); real(dp),intent(in)::oldmin,oldmax,newmin,newmax
         real(dp)::nr,orr; nr=newmax-newmin; orr=oldmax-oldmin
         nodes=(nr*nodes+newmin*oldmax-newmax*oldmin)/orr; weights=weights*nr/orr
      end subroutine transform_rule
   end function call_cf

   function call_heston_cf(spot,strike,tau,r,q,v0,vlong,rho,kappa,sigma,implied_vol,n_quad,status) result(value)
      real(dp),intent(in)::spot,strike,tau,r,q,v0,vlong,rho,kappa,sigma
      real(dp),intent(out),optional::implied_vol
      integer,intent(in),optional::n_quad
      integer,intent(out),optional::status
      real(dp)::value
      type heston_context
         real(dp)::v0,vlong,rho,kappa,sigma
      end type
      type(heston_context)::ctx
      ctx=heston_context(v0,vlong,rho,kappa,max(sigma,0.01_dp))
      value=call_cf(wrapper,spot,strike,tau,r,q,ctx,implied_vol,n_quad,status)
   contains
      function wrapper(omega,s,t,rr,qq,context) result(phi)
         complex(dp),intent(in)::omega; real(dp),intent(in)::s,t,rr,qq; class(*),intent(in),optional::context
         complex(dp)::phi
         select type(context); type is(heston_context); phi=cf_heston(omega,s,t,rr,qq,context%v0,context%vlong,context%rho,context%kappa,context%sigma)
         class default; phi=cmplx(0.0_dp,0.0_dp,dp); end select
      end function wrapper
   end function call_heston_cf

   function call_merton(spot,strike,tau,r,q,variance,lambda,mu_j,v_j,n_jumps,implied_vol,status) result(value)
      real(dp),intent(in)::spot,strike,tau,r,q,variance,lambda,mu_j,v_j
      integer,intent(in)::n_jumps
      real(dp),intent(out),optional::implied_vol
      integer,intent(out),optional::status
      real(dp)::value,lambda2,vn,rn,weight,vol
      type(option_result)::opt
      integer::n,st
      lambda2=lambda*(1.0_dp+mu_j); value=0.0_dp
      do n=0,n_jumps
         vn=variance+real(n,dp)*v_j/tau; rn=r-lambda*mu_j+real(n,dp)*log(1.0_dp+mu_j)/tau
         weight=exp(-lambda2*tau+real(n,dp)*log(max(lambda2*tau,tiny(1.0_dp)))-log_factorial(n))
         if(n==0.and.abs(lambda2*tau)<=tiny(1.0_dp)) weight=1.0_dp
         opt=vanilla_option_european(spot,strike,tau,rn,q,vn,'call',greeks=.false.); value=value+weight*opt%value
      end do
      if(present(implied_vol)) then; vol=vanilla_option_implied_vol('european',value,spot,strike,tau,r,q,'call',status=st); implied_vol=vol; end if
      if(present(status)) status=nmof_ok
   end function call_merton

   pure real(dp) function bund_future(clean,coupon,trade_day,expiry_day,last_coupon_day,r,conversion_factor) result(value)
      real(dp),intent(in)::clean,coupon,r,conversion_factor
      integer,intent(in)::trade_day,expiry_day,last_coupon_day
      real(dp)::days
      days=real(expiry_day-trade_day,dp)
      value=(clean+(clean+coupon*real(trade_day-last_coupon_day,dp)/365.0_dp)*r*days/360.0_dp-coupon*days/365.0_dp)/conversion_factor
   end function bund_future

   pure real(dp) function bund_future_implied_rate(future,clean,coupon,trade_day,expiry_day,last_coupon_day,conversion_factor) result(value)
      real(dp),intent(in)::future,clean,coupon,conversion_factor
      integer,intent(in)::trade_day,expiry_day,last_coupon_day
      real(dp)::days
      days=real(expiry_day-trade_day,dp)
      value=(future*conversion_factor-clean+coupon*days/365.0_dp)/(clean+coupon*real(trade_day-last_coupon_day,dp)/365.0_dp)/(days/360.0_dp)
   end function bund_future_implied_rate

   pure real(dp) function xt_contract_value(quoted_price,coupon,do_round) result(value)
      real(dp),intent(in)::quoted_price,coupon
      logical,intent(in),optional::do_round
      real(dp)::yield_value,i,v,c
      logical::roundv
      roundv=.true.; if(present(do_round)) roundv=do_round; yield_value=100.0_dp-quoted_price; i=yield_value/200.0_dp
      v=1.0_dp/(1.0_dp+i); if(roundv) v=round_to(v,8); c=coupon/2.0_dp
      if(abs(i)<tiny(1.0_dp)) then; value=1000.0_dp*(20.0_dp*c+100.0_dp)
      else
         value=1000.0_dp*((c*(1.0_dp-v**20)/i)+100.0_dp*v**20)
         if(roundv) value=1000.0_dp*(round_to(c*(1.0_dp-v**20)/i,8)+100.0_dp*round_to(v**20,8))
      end if
   end function xt_contract_value

   pure real(dp) function xt_tick_value(quoted_price,coupon,do_round) result(value)
      real(dp),intent(in)::quoted_price,coupon; logical,intent(in),optional::do_round
      value=(xt_contract_value(quoted_price+0.01_dp,coupon,do_round)-xt_contract_value(quoted_price-0.01_dp,coupon,do_round))*50.0_dp
   end function xt_tick_value

   function present_value_dividends(r,tau,q,times,divs,status) result(pv)
      real(dp),intent(in)::r,tau,q
      real(dp),intent(in),optional::times(:),divs(:)
      integer,intent(out)::status
      real(dp)::pv
      integer::i
      pv=0.0_dp; status=nmof_ok
      if(present(divs)) then
         if(.not.present(times).or.size(times)/=size(divs)) then; status=nmof_invalid_input; return; end if
         if(abs(q)>tiny(1.0_dp).and.any(abs(divs)>tiny(1.0_dp))) then; status=nmof_invalid_input; return; end if
         do i=1,size(divs); if(times(i)<=tau.and.times(i)>0.0_dp) pv=pv+exp(-r*times(i))*divs(i); end do
      end if
   end function present_value_dividends

   function future_dividend_value(r,current_time,tau,times,divs) result(pv)
      real(dp),intent(in)::r,current_time,tau
      real(dp),intent(in),optional::times(:),divs(:)
      real(dp)::pv
      integer::i
      pv=0.0_dp
      if(present(times).and.present(divs)) then
         do i=1,min(size(times),size(divs)); if(times(i)>current_time.and.times(i)<=tau) pv=pv+divs(i)*exp(-r*(times(i)-current_time)); end do
      end if
   end function future_dividend_value

   pure elemental real(dp) function round_to(x,digits) result(y)
      real(dp),intent(in)::x; integer,intent(in)::digits; real(dp)::s
      s=10.0_dp**digits; y=anint(x*s)/s
   end function round_to

   pure function lowercase(s) result(t)
      character(len=*),intent(in)::s; character(len=len(s))::t; integer::i,c
      do i=1,len(s); c=iachar(s(i:i)); if(c>=65.and.c<=90) then; t(i:i)=achar(c+32); else; t(i:i)=s(i:i); end if; end do
   end function lowercase
end module nmof_finance

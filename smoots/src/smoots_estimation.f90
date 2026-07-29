! SPDX-License-Identifier: GPL-3.0-only
module smoots_estimation
   use smoots_kinds, only : dp
   use smoots_status, only : sm_ok, sm_invalid_input, sm_iteration_limit
   use smoots_types
   use smoots_stats, only : factorial_real, normal_quantile, mean_value
   use smoots_linalg, only : least_squares_normal
   use smoots_smoothing
   use smoots_arma, only : estimate_cf0_ar, estimate_cf0_ma, estimate_cf0_arma
   implicit none
   private
   public :: tsmooth, msmooth, dsmooth, fixed_gsmooth, fixed_knsmooth
   public :: conf_bounds
contains
   subroutine fixed_gsmooth(y,v,p,mu,b,bb,result)
      real(dp),intent(in)::y(:),b
      integer,intent(in)::v,p,mu,bb
      type(smooth_result),intent(out)::result
      integer::status
      result%n=size(y);result%v=v;result%p=p;result%mu=mu;result%bb=bb
      result%b0=b;result%method=sm_method_lpr
      allocate(result%original(size(y)));result%original=y
      call local_polynomial_smooth(y,v,p,mu,b,bb,result%estimate,result%weights,status)
      result%status=status
      if(v==0)then;allocate(result%residuals(size(y)));result%residuals=y-result%estimate;end if
   end subroutine fixed_gsmooth

   subroutine fixed_knsmooth(y,mu,b,bb,result)
      real(dp),intent(in)::y(:),b
      integer,intent(in)::mu,bb
      type(smooth_result),intent(out)::result
      integer::status
      result%n=size(y);result%v=0;result%p=1;result%mu=mu;result%bb=bb
      result%b0=b;result%method=sm_method_kernel
      allocate(result%original(size(y)));result%original=y
      call kernel_smooth(y,mu,b,bb,result%estimate,status)
      result%status=status
      allocate(result%residuals(size(y)));result%residuals=y-result%estimate
   end subroutine fixed_knsmooth

   subroutine tsmooth(y,p,mu,cf_method,inflation,b_start,enlarged_variance,bb,boundary_cut,method,result)
      real(dp),intent(in)::y(:),b_start,boundary_cut
      integer,intent(in)::p,mu,cf_method,inflation,bb,method
      logical,intent(in)::enlarged_variance
      type(smooth_result),intent(out)::result
      real(dp)::rp,muk,c1,c2,c3,bd,bv,bold,bold1,bopt,i2,cf0,expo1,expo2,minb
      real(dp),allocatable::yed(:),ye_cf(:),dummy_weights(:,:)
      type(arma_model)::arma
      integer::n,k,pd,n1,i,status,l0
      logical::running

      n=size(y);result%n=n;result%p=p;result%mu=mu;result%v=0;result%bb=bb
      result%method=method;result%cf_method=cf_method;result%inflation=inflation
      result%b_start=b_start;result%boundary_cut=boundary_cut
      result%enlarged_variance_bandwidth=enlarged_variance
      if(n<10.or.(p/=1.and.p/=3).or.mu<0.or.mu>3.or.b_start<=0.0_dp.or. &
         boundary_cut<0.0_dp.or.boundary_cut>=0.5_dp)then
         result%status=sm_invalid_input;return
      end if
      allocate(result%original(n));result%original=y
      k=p+1;pd=p+2;n1=int(real(n,dp)*boundary_cut)
      call trend_kernel_constants(p,mu,boundary_cut,rp,muk,status)
      if(status/=sm_ok.or.abs(muk)<=tiny(1.0_dp))then;result%status=sm_invalid_input;return;end if
      c1=factorial_real(k)**2/(2.0_dp*real(k,dp))
      c2=(1.0_dp-2.0_dp*boundary_cut)*rp/(muk*muk)
      call trend_exponents(p,expo1,expo2);minb=real(n,dp)**expo2
      allocate(result%bandwidth_steps(40));result%bandwidth_steps=0.0_dp
      running=.true.;bopt=min(0.49_dp,max(minb,b_start));bold1=bopt
      do i=1,40
         if(.not.running)exit
         if(i==1)then;bold=b_start;else;bold1=bold;bold=bopt;end if
         bd=min(0.49_dp,inflation_bandwidth(p,inflation,bold))
         call local_polynomial_smooth(y,k,pd,mu,bd,bb,yed,dummy_weights,status)
         if(status/=sm_ok)then;result%status=status;return;end if
         i2=sum(yed(n1+1:n-n1)**2)/real(n-2*n1,dp)
         if(enlarged_variance)then;bv=variance_bandwidth(p,mu,bold);else;bv=bold;end if
         bv=min(0.49_dp,bv)
         call local_polynomial_smooth(y,0,p,mu,bv,bb,ye_cf,dummy_weights,status)
         if(status/=sm_ok)then;result%status=status;return;end if
         call calculate_cf0(y-ye_cf,cf_method,cf0,l0,arma,status)
         if(status/=sm_ok.and.status/=sm_iteration_limit)then;result%status=status;return;end if
         c3=max(cf0,epsilon(1.0_dp))/max(i2,epsilon(1.0_dp))
         bopt=(c1*c2*c3)**expo1*real(n,dp)**(-expo1)
         bopt=min(0.49_dp,max(minb,bopt));result%bandwidth_steps(i)=bopt
         if(i>2.and.abs(bold-bopt)/max(bopt,epsilon(1.0_dp))<1.0_dp/real(n,dp))running=.false.
         if(i>3.and.abs(bold1-bopt)/max(bopt,epsilon(1.0_dp))<1.0_dp/real(n,dp))then
            bopt=0.5_dp*(bold+bopt);result%bandwidth_steps(i)=bopt;running=.false.
         end if
      end do
      result%niterations=count(result%bandwidth_steps>0.0_dp)
      if(result%niterations<40)result%bandwidth_steps=result%bandwidth_steps(1:result%niterations)
      result%b0=bopt;result%cf0=cf0;result%curvature_integral=i2;result%lag_window=l0
      result%ar_order=arma%p;result%ma_order=arma%q
      if(method==sm_method_lpr)then
         call local_polynomial_smooth(y,0,p,mu,bopt,bb,result%estimate,result%weights,status)
      else
         call kernel_smooth(y,mu,bopt,bb,result%estimate,status)
      end if
      result%status=status
      if(status==sm_ok)then;allocate(result%residuals(n));result%residuals=y-result%estimate;end if
   end subroutine tsmooth

   subroutine msmooth(y,result,p,mu,b_start,algorithm,method)
      real(dp),intent(in)::y(:)
      type(smooth_result),intent(out)::result
      integer,intent(in),optional::p,mu,method
      real(dp),intent(in),optional::b_start
      character(len=*),intent(in),optional::algorithm
      integer::pp,mm,meth,cfm,infl
      real(dp)::bs
      logical::enlarged
      character(len=8)::alg
      pp=1;if(present(p))pp=p;mm=1;if(present(mu))mm=mu
      bs=0.15_dp;if(present(b_start))bs=b_start
      meth=sm_method_lpr;if(present(method))meth=method
      if(meth==sm_method_kernel)pp=1
      if(present(algorithm))then;alg=adjustl(algorithm);else;if(pp==1)then;alg='A';else;alg='B';end if;end if
      call algorithm_settings(trim(alg),cfm,infl,enlarged)
      call tsmooth(y,pp,mm,cfm,infl,bs,enlarged,1,0.05_dp,meth,result)
   end subroutine msmooth

   subroutine algorithm_settings(alg,cf_method,inflation,enlarged)
      character(len=*),intent(in)::alg
      integer,intent(out)::cf_method,inflation
      logical,intent(out)::enlarged
      enlarged=.false.;inflation=sm_infl_opt;cf_method=sm_cf_np
      select case(trim(alg))
      case('A');cf_method=sm_cf_np;inflation=sm_infl_opt;enlarged=.true.
      case('B');cf_method=sm_cf_np;inflation=sm_infl_naive;enlarged=.true.
      case('N');cf_method=sm_cf_np;inflation=sm_infl_naive
      case('NA');cf_method=sm_cf_ar;inflation=sm_infl_naive
      case('O');cf_method=sm_cf_np;inflation=sm_infl_opt
      case('OA');cf_method=sm_cf_ar;inflation=sm_infl_opt
      case('OM');cf_method=sm_cf_ma;inflation=sm_infl_opt
      case('NM');cf_method=sm_cf_ma;inflation=sm_infl_naive
      case('OAM');cf_method=sm_cf_arma;inflation=sm_infl_opt
      case('NAM');cf_method=sm_cf_arma;inflation=sm_infl_naive
      case default;cf_method=sm_cf_np;inflation=sm_infl_opt;enlarged=.true.
      end select
   end subroutine algorithm_settings

   subroutine calculate_cf0(x,cf_method,cf0,l0,arma,status)
      real(dp),intent(in)::x(:)
      integer,intent(in)::cf_method
      real(dp),intent(out)::cf0
      integer,intent(out)::l0,status
      type(arma_model),intent(out)::arma
      integer::lg
      l0=0
      select case(cf_method)
      case(sm_cf_np);call lag_window_variance(x,cf0,l0,lg,status)
      case(sm_cf_ar);call estimate_cf0_ar(x,cf0,arma,status)
      case(sm_cf_ma);call estimate_cf0_ma(x,cf0,arma,status)
      case(sm_cf_arma);call estimate_cf0_arma(x,cf0,arma,status)
      case default;status=sm_invalid_input;cf0=0.0_dp
      end select
   end subroutine calculate_cf0

   subroutine dsmooth(y,result,d,mu,pilot_p,pilot_b_start,b_start)
      real(dp),intent(in)::y(:)
      type(smooth_result),intent(out)::result
      integer,intent(in),optional::d,mu,pilot_p
      real(dp),intent(in),optional::pilot_b_start,b_start
      integer::dd,mm,pp,n,p,k,pd,n1,i,status,infl
      real(dp)::bps,bs,cf0,rp,muk,c1,c2,c3,bd,bopt,bold,bold1,i2,expo1,expo2,minb
      real(dp),allocatable::yed(:),dw(:,:)
      type(smooth_result)::pilot
      logical::running
      dd=1;if(present(d))dd=d;mm=1;if(present(mu))mm=mu;pp=1;if(present(pilot_p))pp=pilot_p
      bps=0.15_dp;if(present(pilot_b_start))bps=pilot_b_start
      bs=0.15_dp;if(present(b_start))bs=b_start
      if(pp==1)then;call msmooth(y,pilot,p=1,mu=mm,b_start=bps,algorithm='A');else;call msmooth(y,pilot,p=3,mu=mm,b_start=bps,algorithm='B');end if
      if(pilot%status/=sm_ok)then;result%status=pilot%status;return;end if
      cf0=pilot%cf0;n=size(y);p=dd+1;k=p+1;pd=p+2;n1=int(0.05_dp*real(n,dp))
      infl=merge(sm_infl_naive,sm_infl_variance,dd==1)
      call derivative_kernel_constants(dd,mm,rp,muk,status)
      c1=factorial_real(k)**2*real(2*dd+1,dp)/(2.0_dp*real(k-dd,dp))
      c2=0.9_dp*rp/(muk*muk);call derivative_exponents(dd,expo1,expo2);minb=real(n,dp)**expo2
      allocate(result%bandwidth_steps(40));result%bandwidth_steps=0.0_dp
      running=.true.;bopt=min(0.49_dp,max(minb,bs));bold1=bopt
      do i=1,40
         if(.not.running)exit
         if(i==1)then;bold=bs;else;bold1=bold;bold=bopt;end if
         bd=min(0.49_dp,derivative_inflation_bandwidth(dd,infl,bold))
         call local_polynomial_smooth(y,k,pd,mm,bd,1,yed,dw,status)
         i2=sum(yed(n1+1:n-n1)**2)/real(n-2*n1,dp)
         c3=max(cf0,epsilon(1.0_dp))/max(i2,epsilon(1.0_dp))
         bopt=(c1*c2*c3)**expo1*real(n,dp)**(-expo1)
         bopt=min(0.49_dp,max(minb,bopt));result%bandwidth_steps(i)=bopt
         if(i>2.and.abs(bold-bopt)/bopt<1.0_dp/real(n,dp))running=.false.
         if(i>3.and.abs(bold1-bopt)/bopt<1.0_dp/real(n,dp))then;bopt=0.5_dp*(bold+bopt);running=.false.;end if
      end do
      result%n=n;result%v=dd;result%p=p;result%mu=mm;result%bb=1;result%b0=bopt
      result%b_start=bs;result%pilot_b_start=bps;result%pilot_p=pp;result%cf0=cf0
      result%curvature_integral=i2;result%inflation=infl;result%cf_method=sm_cf_np
      result%enlarged_variance_bandwidth=.true.;result%niterations=count(result%bandwidth_steps>0.0_dp)
      if(result%niterations<40)result%bandwidth_steps=result%bandwidth_steps(1:result%niterations)
      allocate(result%original(n));result%original=y
      call local_polynomial_smooth(y,dd,p,mm,bopt,1,result%estimate,result%weights,status)
      result%status=status
   end subroutine dsmooth

   subroutine conf_bounds(object,confidence,result,parametric_order)
      type(smooth_result),intent(in)::object
      real(dp),intent(in),optional::confidence
      type(confidence_result),intent(out)::result
      integer,intent(in),optional::parametric_order
      real(dp)::alpha,bub,bv,cf,z,bpilotub
      real(dp),allocatable::ye_cf(:),weights_cf(:,:),ub_weights(:,:),rx_sub(:),rx(:),xmat(:,:),beta(:)
      type(arma_model)::arma
      type(smooth_result)::pilot_model
      integer::n,v,pp,pilot_pp,po,status,l0,i,l,boundary,nw
      alpha=0.95_dp;if(present(confidence))alpha=confidence
      po=1;if(present(parametric_order))po=parametric_order
      n=object%n;v=object%v;pp=object%p;pilot_pp=object%pilot_p
      if(v>0)po=0
      result%confidence_level=alpha;result%derivative_order=v
      if(n<3.or.alpha<=0.0_dp.or.alpha>=1.0_dp)then;result%status=sm_invalid_input;return;end if
      bub=object%b0**((2.0_dp*real(pp+1,dp)+1.0_dp)/(2.0_dp*real(pp+1,dp)))
      call local_polynomial_smooth(object%original,v,pp,object%mu,bub,object%bb,result%estimate,ub_weights,status)
      if(status/=sm_ok)then;result%status=status;return;end if
      if(v==0)then
         if(object%enlarged_variance_bandwidth)then;bv=min(0.49_dp,variance_bandwidth(pp,object%mu,bub));else;bv=bub;end if
         call local_polynomial_smooth(object%original,0,pp,object%mu,bv,object%bb,ye_cf,weights_cf,status)
      else
         if(pilot_pp==3)then
            call msmooth(object%original,pilot_model,p=3,mu=object%mu,b_start=object%pilot_b_start,algorithm='B')
         else
            pilot_pp=1
            call msmooth(object%original,pilot_model,p=1,mu=object%mu,b_start=object%pilot_b_start,algorithm='A')
         end if
         bpilotub=pilot_model%b0**((2.0_dp*real(pilot_pp+1,dp)+1.0_dp)/(2.0_dp*real(pilot_pp+1,dp)))
         bv=min(0.49_dp,variance_bandwidth(pilot_pp,object%mu,bpilotub))
         call local_polynomial_smooth(object%original,0,pilot_pp,object%mu,bv,1,ye_cf,weights_cf,status)
      end if
      call calculate_cf0(object%original-ye_cf,object%cf_method,cf,l0,arma,status)
      nw=size(ub_weights,1);allocate(rx_sub(nw));rx_sub=sum(ub_weights**2,dim=2)
      boundary=(nw-1)/2;allocate(rx(n))
      if(boundary>0)then
         rx(1:boundary)=rx_sub(1:boundary);rx(boundary+1:n-boundary)=rx_sub(boundary+1);rx(n-boundary+1:n)=rx_sub(boundary+2:nw)
      else;rx=rx_sub(1);end if
      z=normal_quantile(alpha+(1.0_dp-alpha)/2.0_dp)
      allocate(result%lower(n),result%upper(n));result%lower=result%estimate-z*sqrt(max(0.0_dp,cf*rx));result%upper=result%estimate+z*sqrt(max(0.0_dp,cf*rx))
      allocate(result%parametric(n));result%parametric=0.0_dp
      if(v==0.and.po>=0)then
         allocate(xmat(n,po+1),beta(po+1));xmat(:,1)=1.0_dp
         do i=1,po;xmat(:,i+1)=[(real(l,dp)**i,l=1,n)];end do
         call least_squares_normal(xmat,object%original,beta,status,1.0e-12_dp)
         result%parametric=matmul(xmat,beta)
      else if(v==1.and.po==0)then
         allocate(xmat(n,2),beta(2));xmat(:,1)=1.0_dp;xmat(:,2)=[(real(i,dp)/real(n,dp),i=1,n)]
         call least_squares_normal(xmat,object%original,beta,status,1.0e-12_dp);result%parametric=beta(2)
      end if
      result%unbiased_bandwidth=bub;result%status=sm_ok
   end subroutine conf_bounds
end module smoots_estimation

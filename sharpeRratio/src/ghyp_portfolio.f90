! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from ghyp 1.6.5 by Marc Weibel, David Luethi, and Henriette-Elise Breymann.
module ghyp_portfolio
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ghyp_kinds, only : dp
   use ghyp_model, only : ghyp_model_type, ghyp_moments, moments_result, transform_ghyp
   use ghyp_distribution, only : qghyp
   use ghyp_risk, only : esghyp
   use ghyp_linalg, only : solve_spd, solve_linear
   use ghyp_optimize, only : nelder_mead
   implicit none
   private

   type, public :: portfolio_result
      real(dp), allocatable :: weights(:)
      real(dp) :: expected_return = 0.0_dp
      real(dp) :: standard_deviation = 0.0_dp
      real(dp) :: risk = 0.0_dp
      type(ghyp_model_type) :: distribution
      integer :: iterations = 0
      logical :: converged = .false.
      logical :: ok = .false.
      character(len=160) :: message = ''
   end type portfolio_result

   public :: portfolio_optimize

   type(ghyp_model_type), save :: active_model
   real(dp), save :: active_level = 0.95_dp
   real(dp), save :: active_target = 0.0_dp
   real(dp), save :: active_risk_free = 0.0_dp
   integer, save :: active_risk = 1
   integer, save :: active_problem = 1
   logical, save :: active_loss = .true.

contains

   function portfolio_optimize(model, risk_measure, problem, level, target_return, &
      risk_free, loss, max_iter) result(result)
      type(ghyp_model_type), intent(in) :: model
      character(len=*), intent(in), optional :: risk_measure, problem
      real(dp), intent(in), optional :: level, target_return, risk_free
      logical, intent(in), optional :: loss
      integer, intent(in), optional :: max_iter
      type(portfolio_result) :: result
      type(moments_result) :: mom
      real(dp), allocatable :: invone(:), invmu(:), ones(:), rhs(:), solution(:)
      real(dp), allocatable :: start(:), optimum(:), a(:,:), kkt(:,:)
      real(dp) :: denom, obj, tr, rf
      integer :: d, rm, pb, maxit
      logical :: ok, use_loss, conv
      character(len=:), allocatable :: rkey, pkey

      if (.not. model%ok .or. model%dimension() < 2) then
         result%message = 'a valid multivariate model is required'
         return
      end if
      d = model%dimension()
      mom = ghyp_moments(model)
      if (.not. mom%ok) then
         result%message = 'model moments are unavailable'
         return
      end if
      rkey = 'sd'; if (present(risk_measure)) rkey = lower_string(trim(risk_measure))
      pkey = 'minimum.risk'; if (present(problem)) pkey = lower_string(trim(problem))
      select case(rkey)
      case('sd','standard.deviation','standard-deviation'); rm=1
      case('value.at.risk','var','value-at-risk'); rm=2
      case('expected.shortfall','es'); rm=3
      case default; result%message='unknown risk measure'; return
      end select
      select case(pkey)
      case('minimum.risk','minimum','min'); pb=1
      case('tangency'); pb=2
      case('target.return','target'); pb=3
      case default; result%message='unknown optimization problem'; return
      end select
      tr=0.0_dp;rf=0.0_dp;use_loss=.true.
      if(present(target_return))tr=target_return
      if(present(risk_free))rf=risk_free
      if(present(loss))use_loss=loss
      if(pb==3.and..not.present(target_return))then
         result%message='target return is required';return
      end if
      if(pb==2.and..not.present(risk_free))then
         result%message='risk-free rate is required';return
      end if
      allocate(result%weights(d),ones(d))
      ones=1.0_dp

      if(rm==1.or.model%is_symmetric())then
         select case(pb)
         case(1)
            call solve_spd(mom%covariance,ones,invone,ok)
            if(.not.ok)then;result%message='singular covariance matrix';return;end if
            denom=sum(invone);result%weights=invone/denom
         case(2)
            call solve_spd(mom%covariance,mom%mean-rf,invmu,ok)
            if(.not.ok.or.abs(sum(invmu))<=tiny(1.0_dp))then
               result%message='tangency portfolio is undefined';return
            end if
            result%weights=invmu/sum(invmu)
         case(3)
            allocate(kkt(d+2,d+2),rhs(d+2))
            kkt=0.0_dp
            kkt(1:d,1:d)=mom%covariance
            kkt(1:d,d+1)=-ones;kkt(d+1,1:d)=ones
            kkt(1:d,d+2)=-mom%mean;kkt(d+2,1:d)=mom%mean
            rhs=0.0_dp;rhs(d+1)=1.0_dp;rhs(d+2)=tr
            call solve_linear(kkt,rhs,solution,ok)
            if(.not.ok)then;result%message='target-return system is singular';return;end if
            result%weights=solution(1:d)
         end select
         result%converged=.true.
      else
         active_model=model;active_level=0.95_dp;if(present(level))active_level=level
         active_target=tr;active_risk_free=rf;active_risk=rm;active_problem=pb
         active_loss=use_loss
         allocate(start(d-1));start=1.0_dp/real(d,dp)
         maxit=1000;if(present(max_iter))maxit=max_iter
         call nelder_mead(portfolio_objective,start,optimum,obj,conv,result%iterations, &
            maxit,1.0e-8_dp)
         result%weights(2:d)=optimum
         result%weights(1)=1.0_dp-sum(optimum)
         result%converged=conv
      end if
      allocate(a(1,d));a(1,:)=result%weights
      result%distribution=transform_ghyp(model,a)
      if(.not.result%distribution%ok)then
         result%message='portfolio transformation failed';return
      end if
      result%expected_return=dot_product(result%weights,mom%mean)
      result%standard_deviation=sqrt(max(0.0_dp, &
         dot_product(result%weights,matmul(mom%covariance,result%weights))))
      select case(rm)
      case(1);result%risk=result%standard_deviation
      case(2);result%risk=qghyp(active_level_or(level),result%distribution)
      case(3);result%risk=esghyp(active_level_or(level),result%distribution,use_loss)
      end select
      result%ok=ieee_is_finite(result%risk)
      if(.not.result%ok)result%message='risk evaluation failed'
   end function portfolio_optimize

   function portfolio_objective(pars) result(value)
      real(dp),intent(in)::pars(:)
      real(dp)::value,risk,ret,penalty
      real(dp),allocatable::weights(:),a(:,:)
      type(ghyp_model_type)::ptf
      integer::d
      value=huge(1.0_dp)/100.0_dp
      d=active_model%dimension();allocate(weights(d),a(1,d))
      weights(2:d)=pars;weights(1)=1.0_dp-sum(pars);a(1,:)=weights
      ptf=transform_ghyp(active_model,a)
      if(.not.ptf%ok)then;value=huge(1.0_dp)/100.0_dp;return;end if
      if(active_risk==2)then
         risk=qghyp(active_level,ptf)
      else
         risk=esghyp(active_level,ptf,active_loss)
      end if
      if(.not.ieee_is_finite(risk))then;value=huge(1.0_dp)/100.0_dp;return;end if
      ret=ghyp_portfolio_mean(active_model,weights)
      penalty=0.0_dp
      select case(active_problem)
      case(1);value=risk
      case(2);value=-(ret-active_risk_free)/max(abs(risk),1.0e-12_dp)
      case(3)
         penalty=1.0e5_dp*(ret-active_target)**2
         value=risk+penalty
      end select
   end function portfolio_objective

   function ghyp_portfolio_mean(model,weights) result(value)
      type(ghyp_model_type),intent(in)::model
      real(dp),intent(in)::weights(:)
      real(dp)::value
      type(moments_result)::mom
      mom=ghyp_moments(model)
      if(mom%ok)then;value=dot_product(weights,mom%mean);else;value=huge(1.0_dp);end if
   end function ghyp_portfolio_mean

   pure function active_level_or(level) result(value)
      real(dp),intent(in),optional::level
      real(dp)::value
      value=0.95_dp;if(present(level))value=level
   end function active_level_or

   pure function lower_string(s) result(out)
      character(len=*),intent(in)::s
      character(len=len(s))::out
      integer::i,c
      do i=1,len(s)
         c=iachar(s(i:i));if(c>=iachar('A').and.c<=iachar('Z'))then
            out(i:i)=achar(c+32)
         else;out(i:i)=s(i:i);end if
      end do
   end function lower_string

end module ghyp_portfolio

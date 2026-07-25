! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module msgarch_simulation
   use msgarch_kinds, only : dp
   use msgarch_types, only : msgarch_spec, regime_state, simulation_result, filter_result
   use msgarch_models, only : stationary_distribution, initialize_state, update_state, spec_valid
   use msgarch_distributions, only : random_innovation
   use msgarch_filter, only : hamilton_filter, terminal_states
   use msgarch_rng, only : sample_discrete
   implicit none
   private
   public :: simulate_msgarch, simulate_ahead, forecast_volatility, forecast_mean
   public :: unconditional_volatility
contains
   function simulate_msgarch(spec,n,nsim) result(result)
      type(msgarch_spec),intent(in)::spec
      integer,intent(in)::n,nsim
      type(simulation_result)::result
      type(regime_state),allocatable::state(:)
      real(dp),allocatable::initial_probability(:)
      integer::i,t,k,current
      if(n<1.or.nsim<1)return
      if(.not.spec_valid(spec))return
      allocate(result%draw(nsim,n),result%state(nsim,n),result%conditional_sd(nsim,n,spec%k))
      initial_probability=stationary_distribution(spec%transition)
      allocate(state(spec%k))
      do i=1,nsim
         do k=1,spec%k;call initialize_state(spec%regime(k),state(k));end do
         current=sample_discrete(initial_probability);result%state(i,1)=current
         do k=1,spec%k;result%conditional_sd(i,1,k)=sqrt(state(k)%h);end do
         result%draw(i,1)=sqrt(state(current)%h)*random_innovation(spec%regime(current)%distribution, &
            spec%regime(current)%shape,spec%regime(current)%skew)
         do t=2,n
            current=sample_discrete(spec%transition(current,:));result%state(i,t)=current
            do k=1,spec%k
               call update_state(spec%regime(k),state(k),result%draw(i,t-1))
               result%conditional_sd(i,t,k)=sqrt(state(k)%h)
            end do
            result%draw(i,t)=sqrt(state(current)%h)*random_innovation(spec%regime(current)%distribution, &
               spec%regime(current)%shape,spec%regime(current)%skew)
         end do
      end do
   end function simulate_msgarch

   function simulate_ahead(spec,y,n_ahead,nsim) result(result)
      type(msgarch_spec),intent(in)::spec
      real(dp),intent(in)::y(:)
      integer,intent(in)::n_ahead,nsim
      type(simulation_result)::result
      type(filter_result)::filtered
      type(regime_state),allocatable::terminal(:),state(:)
      integer::i,t,k,current
      if(n_ahead<1.or.nsim<1.or.size(y)<2)return
      if(.not.spec_valid(spec))return
      filtered=hamilton_filter(spec,y);call terminal_states(spec,y,terminal)
      allocate(result%draw(nsim,n_ahead),result%state(nsim,n_ahead),result%conditional_sd(nsim,n_ahead,spec%k))
      allocate(state(spec%k))
      do i=1,nsim
         state=terminal
         current=sample_discrete(filtered%next_probability);result%state(i,1)=current
         do k=1,spec%k;result%conditional_sd(i,1,k)=sqrt(state(k)%h);end do
         result%draw(i,1)=sqrt(state(current)%h)*random_innovation(spec%regime(current)%distribution, &
            spec%regime(current)%shape,spec%regime(current)%skew)
         do t=2,n_ahead
            current=sample_discrete(spec%transition(current,:));result%state(i,t)=current
            do k=1,spec%k
               call update_state(spec%regime(k),state(k),result%draw(i,t-1))
               result%conditional_sd(i,t,k)=sqrt(state(k)%h)
            end do
            result%draw(i,t)=sqrt(state(current)%h)*random_innovation(spec%regime(current)%distribution, &
               spec%regime(current)%shape,spec%regime(current)%skew)
         end do
      end do
   end function simulate_ahead

   function forecast_volatility(spec,y,n_ahead,nsim,cumulative) result(volatility)
      type(msgarch_spec),intent(in)::spec
      real(dp),intent(in)::y(:)
      integer,intent(in)::n_ahead,nsim
      logical,intent(in),optional::cumulative
      real(dp),allocatable::volatility(:)
      type(simulation_result)::simulation
      real(dp),allocatable::path(:,:)
      logical::do_cumulative
      integer::t
      do_cumulative=.false.;if(present(cumulative))do_cumulative=cumulative
      simulation=simulate_ahead(spec,y,n_ahead,nsim)
      allocate(path(nsim,n_ahead));path=simulation%draw
      if(do_cumulative)then
         do t=2,n_ahead;path(:,t)=path(:,t-1)+path(:,t);end do
      end if
      allocate(volatility(n_ahead))
      do t=1,n_ahead
         volatility(t)=sqrt(sum(path(:,t)**2)/real(nsim,dp)-(sum(path(:,t))/real(nsim,dp))**2)
      end do
   end function forecast_volatility

   function unconditional_volatility(spec,nsim,n_ahead,burn) result(value)
      type(msgarch_spec),intent(in)::spec
      integer,intent(in),optional::nsim,n_ahead,burn
      real(dp)::value
      type(simulation_result)::simulation
      integer::m,horizon,discard,t
      real(dp)::mean_t,var_t,total
      m=250;if(present(nsim))m=nsim
      horizon=5000;if(present(n_ahead))horizon=n_ahead
      discard=1000;if(present(burn))discard=burn
      if(m<2.or.horizon<2.or.discard<0.or.discard>=horizon)error stop 'unconditional_volatility: invalid controls'
      simulation=simulate_msgarch(spec,horizon,m);total=0.0_dp
      do t=discard+1,horizon
         mean_t=sum(simulation%draw(:,t))/real(m,dp)
         var_t=sum((simulation%draw(:,t)-mean_t)**2)/real(m-1,dp)
         total=total+sqrt(max(var_t,0.0_dp))
      end do
      value=total/real(horizon-discard,dp)
   end function unconditional_volatility

   function forecast_mean(spec,y,n_ahead,nsim,cumulative) result(mean_forecast)
      type(msgarch_spec),intent(in)::spec
      real(dp),intent(in)::y(:)
      integer,intent(in)::n_ahead,nsim
      logical,intent(in),optional::cumulative
      real(dp),allocatable::mean_forecast(:)
      type(simulation_result)::simulation
      real(dp),allocatable::path(:,:)
      logical::do_cumulative
      integer::t
      do_cumulative=.false.;if(present(cumulative))do_cumulative=cumulative
      simulation=simulate_ahead(spec,y,n_ahead,nsim);allocate(path(nsim,n_ahead));path=simulation%draw
      if(do_cumulative)then;do t=2,n_ahead;path(:,t)=path(:,t-1)+path(:,t);end do;end if
      allocate(mean_forecast(n_ahead));mean_forecast=sum(path,dim=1)/real(nsim,dp)
   end function forecast_mean
end module msgarch_simulation

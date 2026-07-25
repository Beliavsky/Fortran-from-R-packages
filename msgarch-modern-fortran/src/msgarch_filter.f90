! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module msgarch_filter
   use msgarch_kinds, only : dp
   use msgarch_types, only : msgarch_spec, regime_state, filter_result
   use msgarch_models, only : spec_valid, stationary_distribution, initialize_state, update_state, conditional_log_density
   use msgarch_distributions, only : innovation_pdf, innovation_cdf
   implicit none
   private
   public :: hamilton_filter, predictive_pdf, predictive_cdf, pit_values
   public :: conditional_volatility, terminal_states, in_sample_pdf, in_sample_cdf
contains
   function hamilton_filter(spec,y) result(result)
      type(msgarch_spec),intent(in)::spec
      real(dp),intent(in)::y(:)
      type(filter_result)::result
      type(regime_state),allocatable::state(:)
      real(dp),allocatable::emission(:),work(:),scale(:),delta(:),ratio(:),score(:,:),best(:)
      integer,allocatable::back(:,:)
      real(dp)::maximum,total
      integer::t,k,j,n
      n=size(y)
      if(n<2)return
      if(.not.spec_valid(spec))return
      allocate(result%log_density(n,spec%k),result%variance(n,spec%k),result%predicted(n,spec%k), &
         result%filtered(n,spec%k),result%smoothed(n,spec%k),result%viterbi(n),result%next_probability(spec%k))
      allocate(state(spec%k),emission(spec%k),work(spec%k),scale(n),delta(spec%k),ratio(spec%k), &
         score(spec%k,n),best(spec%k),back(spec%k,n))
      result%log_density=0.0_dp;result%variance=0.0_dp;result%predicted=0.0_dp;result%filtered=0.0_dp
      result%smoothed=0.0_dp;score=-huge(1.0_dp);back=1;scale=1.0_dp
      delta=stationary_distribution(spec%transition)
      result%predicted(1,:)=delta;result%filtered(1,:)=delta
      do k=1,spec%k
         call initialize_state(spec%regime(k),state(k));result%variance(1,k)=state(k)%h
         score(k,1)=log(max(delta(k),tiny(1.0_dp)))
      end do
      result%loglik=0.0_dp
      do t=2,n
         result%predicted(t,:)=matmul(result%filtered(t-1,:),spec%transition)
         do k=1,spec%k
            call update_state(spec%regime(k),state(k),y(t-1))
            result%variance(t,k)=state(k)%h
            result%log_density(t,k)=conditional_log_density(spec%regime(k),state(k),y(t))
         end do
         maximum=maxval(result%log_density(t,:))
         emission=exp(result%log_density(t,:)-maximum)
         work=result%predicted(t,:)*emission
         total=sum(work)
         if(total<=tiny(1.0_dp))then
            result%loglik=-huge(1.0_dp);return
         end if
         result%filtered(t,:)=work/total
         scale(t)=total
         result%loglik=result%loglik+maximum+log(total)
         do k=1,spec%k
            best(k)=-huge(1.0_dp);back(k,t)=1
            do j=1,spec%k
               if(score(j,t-1)+log(max(spec%transition(j,k),tiny(1.0_dp)))>best(k))then
                  best(k)=score(j,t-1)+log(max(spec%transition(j,k),tiny(1.0_dp)))
                  back(k,t)=j
               end if
            end do
            score(k,t)=best(k)+result%log_density(t,k)
         end do
      end do
      result%next_probability=matmul(result%filtered(n,:),spec%transition)
      result%smoothed(n,:)=result%filtered(n,:)
      do t=n-1,1,-1
         ratio=result%smoothed(t+1,:)/max(result%predicted(t+1,:),tiny(1.0_dp))
         result%smoothed(t,:)=result%filtered(t,:)*matmul(spec%transition,ratio)
         total=sum(result%smoothed(t,:));if(total>0.0_dp)result%smoothed(t,:)=result%smoothed(t,:)/total
      end do
      result%viterbi(n)=maxloc(score(:,n),dim=1)
      do t=n,2,-1;result%viterbi(t-1)=back(result%viterbi(t),t);end do
   end function hamilton_filter

   subroutine terminal_states(spec,y,state)
      type(msgarch_spec),intent(in)::spec
      real(dp),intent(in)::y(:)
      type(regime_state),allocatable,intent(out)::state(:)
      integer::t,k
      allocate(state(spec%k))
      do k=1,spec%k;call initialize_state(spec%regime(k),state(k));end do
      do t=1,size(y)
         do k=1,spec%k;call update_state(spec%regime(k),state(k),y(t));end do
      end do
   end subroutine terminal_states

   function predictive_pdf(spec,y,x,log_value) result(value)
      type(msgarch_spec),intent(in)::spec
      real(dp),intent(in)::y(:),x
      logical,intent(in),optional::log_value
      real(dp)::value,z
      type(filter_result)::filtered
      type(regime_state),allocatable::state(:)
      integer::k
      filtered=hamilton_filter(spec,y);call terminal_states(spec,y,state);value=0.0_dp
      do k=1,spec%k
         z=x/sqrt(state(k)%h)
         value=value+filtered%next_probability(k)*innovation_pdf(z,spec%regime(k)%distribution, &
            spec%regime(k)%shape,spec%regime(k)%skew)/sqrt(state(k)%h)
      end do
      if(present(log_value))then;if(log_value)value=log(max(value,tiny(1.0_dp)));end if
   end function predictive_pdf

   function predictive_cdf(spec,y,x,log_value) result(value)
      type(msgarch_spec),intent(in)::spec
      real(dp),intent(in)::y(:),x
      logical,intent(in),optional::log_value
      real(dp)::value,z
      type(filter_result)::filtered
      type(regime_state),allocatable::state(:)
      integer::k
      filtered=hamilton_filter(spec,y);call terminal_states(spec,y,state);value=0.0_dp
      do k=1,spec%k
         z=x/sqrt(state(k)%h)
         value=value+filtered%next_probability(k)*innovation_cdf(z,spec%regime(k)%distribution, &
            spec%regime(k)%shape,spec%regime(k)%skew)
      end do
      value=max(0.0_dp,min(1.0_dp,value))
      if(present(log_value))then;if(log_value)value=log(max(value,tiny(1.0_dp)));end if
   end function predictive_cdf

   function pit_values(spec,y) result(pit)
      type(msgarch_spec),intent(in)::spec
      real(dp),intent(in)::y(:)
      real(dp),allocatable::pit(:)
      type(filter_result)::filtered
      integer::t,k
      real(dp)::z
      filtered=hamilton_filter(spec,y);allocate(pit(size(y)));pit(1)=0.5_dp
      do t=2,size(y)
         pit(t)=0.0_dp
         do k=1,spec%k
            z=y(t)/sqrt(filtered%variance(t,k))
            pit(t)=pit(t)+filtered%predicted(t,k)*innovation_cdf(z,spec%regime(k)%distribution, &
               spec%regime(k)%shape,spec%regime(k)%skew)
         end do
      end do
   end function pit_values

   function conditional_volatility(spec,y,probability_type) result(volatility)
      type(msgarch_spec),intent(in)::spec
      real(dp),intent(in)::y(:)
      character(len=*),intent(in),optional::probability_type
      real(dp),allocatable::volatility(:)
      type(filter_result)::filtered
      integer::t
      logical::use_predicted
      filtered=hamilton_filter(spec,y);allocate(volatility(size(y)));use_predicted=.false.
      if(present(probability_type))use_predicted=trim(probability_type)=='predicted'
      do t=1,size(y)
         if(use_predicted)then
            volatility(t)=sqrt(sum(filtered%predicted(t,:)*filtered%variance(t,:)))
         else
            volatility(t)=sqrt(sum(filtered%filtered(t,:)*filtered%variance(t,:)))
         end if
      end do
   end function conditional_volatility

   function in_sample_pdf(spec,y,x) result(value)
      type(msgarch_spec),intent(in)::spec
      real(dp),intent(in)::y(:),x(:)
      real(dp),allocatable::value(:,:)
      type(filter_result)::filtered
      integer::t,k,j
      real(dp)::z
      filtered=hamilton_filter(spec,y);allocate(value(size(y),size(x)));value=0.0_dp
      do t=1,size(y)
         do j=1,size(x)
            do k=1,spec%k
               z=x(j)/sqrt(filtered%variance(t,k))
               value(t,j)=value(t,j)+filtered%predicted(t,k)*innovation_pdf(z,spec%regime(k)%distribution, &
                  spec%regime(k)%shape,spec%regime(k)%skew)/sqrt(filtered%variance(t,k))
            end do
         end do
      end do
   end function in_sample_pdf

   function in_sample_cdf(spec,y,x) result(value)
      type(msgarch_spec),intent(in)::spec
      real(dp),intent(in)::y(:),x(:)
      real(dp),allocatable::value(:,:)
      type(filter_result)::filtered
      integer::t,k,j
      real(dp)::z
      filtered=hamilton_filter(spec,y);allocate(value(size(y),size(x)));value=0.0_dp
      do t=1,size(y)
         do j=1,size(x)
            do k=1,spec%k
               z=x(j)/sqrt(filtered%variance(t,k))
               value(t,j)=value(t,j)+filtered%predicted(t,k)*innovation_cdf(z,spec%regime(k)%distribution, &
                  spec%regime(k)%shape,spec%regime(k)%skew)
            end do
         end do
      end do
   end function in_sample_cdf
end module msgarch_filter

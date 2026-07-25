! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module msgarch_forecast
   use msgarch_kinds, only : dp
   use msgarch_types, only : msgarch_spec, filter_result
   use msgarch_filter, only : hamilton_filter, predictive_pdf, predictive_cdf
   use msgarch_simulation, only : simulate_ahead
   use msgarch_types, only : simulation_result
   use msgarch_special, only : normal_pdf
   implicit none
   private
   public :: transition_matrix_power, forecast_state_probabilities
   public :: predictive_distribution_forecast
contains
   function transition_matrix_power(transition,n) result(power)
      real(dp),intent(in)::transition(:,:)
      integer,intent(in)::n
      real(dp),allocatable::power(:,:),base(:,:),temp(:,:)
      integer::exponent,i,k
      k=size(transition,1);allocate(power(k,k),base(k,k),temp(k,k));power=0.0_dp
      do i=1,k;power(i,i)=1.0_dp;end do
      base=transition;exponent=max(n,0)
      do while(exponent>0)
         if(mod(exponent,2)==1)then;temp=matmul(power,base);power=temp;end if
         exponent=exponent/2
         if(exponent>0)then;temp=matmul(base,base);base=temp;end if
      end do
   end function transition_matrix_power

   subroutine predictive_distribution_forecast(spec,y,x,n_ahead,nsim,density,cdf,simulations,cumulative)
      type(msgarch_spec),intent(in)::spec
      real(dp),intent(in)::y(:),x(:)
      integer,intent(in)::n_ahead,nsim
      real(dp),allocatable,intent(out)::density(:,:),cdf(:,:)
      real(dp),allocatable,intent(out),optional::simulations(:,:)
      logical,intent(in),optional::cumulative
      type(simulation_result)::simulation
      real(dp),allocatable::path(:,:)
      real(dp)::bandwidth,z
      logical::do_cumulative
      integer::h,j,i
      if(n_ahead<1.or.nsim<2)error stop 'predictive_distribution_forecast: invalid controls'
      do_cumulative=.false.;if(present(cumulative))do_cumulative=cumulative
      allocate(density(n_ahead,size(x)),cdf(n_ahead,size(x)));density=0.0_dp;cdf=0.0_dp
      do j=1,size(x)
         density(1,j)=predictive_pdf(spec,y,x(j))
         cdf(1,j)=predictive_cdf(spec,y,x(j))
      end do
      simulation=simulate_ahead(spec,y,n_ahead,nsim);path=simulation%draw
      if(do_cumulative)then
         do h=2,n_ahead;path(:,h)=path(:,h-1)+path(:,h);end do
      end if
      do h=2,n_ahead
         bandwidth=kernel_bandwidth(path(:,h))
         do j=1,size(x)
            do i=1,nsim
               z=(x(j)-path(i,h))/bandwidth
               density(h,j)=density(h,j)+normal_pdf(z)/bandwidth
            end do
            density(h,j)=density(h,j)/real(nsim,dp)
            cdf(h,j)=real(count(path(:,h)<=x(j)),dp)/real(nsim,dp)
         end do
      end do
      if(present(simulations))then
         allocate(simulations(size(path,1),size(path,2)));simulations=path
      end if
   end subroutine predictive_distribution_forecast

   function kernel_bandwidth(values) result(bandwidth)
      real(dp),intent(in)::values(:)
      real(dp)::bandwidth,mean_value,sd
      integer::n
      n=size(values);mean_value=sum(values)/real(n,dp)
      sd=sqrt(sum((values-mean_value)**2)/real(max(n-1,1),dp))
      bandwidth=max(1.0e-8_dp,1.06_dp*max(sd,1.0e-8_dp)*real(n,dp)**(-0.2_dp))
   end function kernel_bandwidth

   function forecast_state_probabilities(spec,y,n_ahead) result(probability)
      type(msgarch_spec),intent(in)::spec
      real(dp),intent(in)::y(:)
      integer,intent(in)::n_ahead
      real(dp),allocatable::probability(:,:)
      type(filter_result)::filtered
      integer::h
      filtered=hamilton_filter(spec,y);allocate(probability(n_ahead,spec%k))
      if(n_ahead<1)return
      probability(1,:)=filtered%next_probability
      do h=2,n_ahead;probability(h,:)=matmul(probability(h-1,:),spec%transition);end do
   end function forecast_state_probabilities
end module msgarch_forecast

! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module msgarch_risk
   use msgarch_kinds, only : dp
   use msgarch_types, only : msgarch_spec, risk_result, simulation_result, filter_result
   use msgarch_models, only : unpack_parameters
   use msgarch_simulation, only : simulate_ahead
   use msgarch_filter, only : hamilton_filter
   use msgarch_distributions, only : innovation_cdf
   implicit none
   private
   public :: risk_forecast, risk_in_sample, posterior_risk_forecast, empirical_quantile
contains
   function empirical_quantile(values,p) result(q)
      real(dp),intent(in)::values(:),p
      real(dp)::q,h,w
      real(dp),allocatable::sorted(:)
      integer::n,i,j
      sorted=values;n=size(sorted)
      do i=2,n
         h=sorted(i);j=i-1
         do while(j>=1)
            if(sorted(j)<=h)exit
            sorted(j+1)=sorted(j);j=j-1
         end do
         sorted(j+1)=h
      end do
      if(p<=0.0_dp)then;q=sorted(1);return;end if
      if(p>=1.0_dp)then;q=sorted(n);return;end if
      h=1.0_dp+(real(n-1,dp))*p;i=floor(h);w=h-real(i,dp)
      if(i>=n)then;q=sorted(n);else;q=(1.0_dp-w)*sorted(i)+w*sorted(i+1);end if
   end function empirical_quantile

   function risk_forecast(spec,y,alpha,n_ahead,nsim,cumulative,return_simulations) result(result)
      type(msgarch_spec),intent(in)::spec
      real(dp),intent(in)::y(:),alpha(:)
      integer,intent(in)::n_ahead,nsim
      logical,intent(in),optional::cumulative,return_simulations
      type(risk_result)::result
      type(simulation_result)::simulation
      real(dp),allocatable::path(:,:)
      real(dp)::q
      logical::do_cumulative,keep
      integer::h,j,nbelow
      do_cumulative=.false.;keep=.false.
      if(present(cumulative))do_cumulative=cumulative
      if(present(return_simulations))keep=return_simulations
      simulation=simulate_ahead(spec,y,n_ahead,nsim);path=simulation%draw
      if(do_cumulative)then;do h=2,n_ahead;path(:,h)=path(:,h-1)+path(:,h);end do;end if
      allocate(result%var(n_ahead,size(alpha)),result%es(n_ahead,size(alpha)))
      do h=1,n_ahead
         do j=1,size(alpha)
            q=empirical_quantile(path(:,h),alpha(j));result%var(h,j)=q
            nbelow=count(path(:,h)<=q)
            if(nbelow>0)then;result%es(h,j)=sum(path(:,h),mask=path(:,h)<=q)/real(nbelow,dp);else;result%es(h,j)=q;end if
         end do
      end do
      if (keep) then
         allocate(result%simulations(size(path,1),size(path,2)))
         result%simulations = path
      else
         allocate(result%simulations(0,0))
      end if
   end function risk_forecast


   function posterior_risk_forecast(template,draws,y,alpha,n_ahead,nsim_per_draw,cumulative) result(result)
      type(msgarch_spec),intent(in)::template
      real(dp),intent(in)::draws(:,:),y(:),alpha(:)
      integer,intent(in)::n_ahead,nsim_per_draw
      logical,intent(in),optional::cumulative
      type(risk_result)::result
      type(msgarch_spec)::spec
      type(simulation_result)::simulation
      real(dp),allocatable::path(:,:)
      real(dp)::q
      integer::i,h,j,offset,nvalid,nbelow,total_n
      logical::valid,do_cumulative
      do_cumulative=.false.;if(present(cumulative))do_cumulative=cumulative
      nvalid=0
      do i=1,size(draws,1)
         call unpack_parameters(template,draws(i,:),spec,valid)
         if(valid)nvalid=nvalid+1
      end do
      total_n=nvalid*nsim_per_draw
      allocate(path(total_n,n_ahead));offset=0
      do i=1,size(draws,1)
         call unpack_parameters(template,draws(i,:),spec,valid)
         if(.not.valid)cycle
         simulation=simulate_ahead(spec,y,n_ahead,nsim_per_draw)
         path(offset+1:offset+nsim_per_draw,:)=simulation%draw
         offset=offset+nsim_per_draw
      end do
      if(do_cumulative)then;do h=2,n_ahead;path(:,h)=path(:,h-1)+path(:,h);end do;end if
      allocate(result%var(n_ahead,size(alpha)),result%es(n_ahead,size(alpha)),result%simulations(total_n,n_ahead))
      result%simulations=path
      do h=1,n_ahead
         do j=1,size(alpha)
            q=empirical_quantile(path(:,h),alpha(j));result%var(h,j)=q
            nbelow=count(path(:,h)<=q)
            if(nbelow>0)then;result%es(h,j)=sum(path(:,h),mask=path(:,h)<=q)/real(nbelow,dp);else;result%es(h,j)=q;end if
         end do
      end do
   end function posterior_risk_forecast

   function risk_in_sample(spec,y,alpha) result(result)
      type(msgarch_spec),intent(in)::spec
      real(dp),intent(in)::y(:),alpha(:)
      type(risk_result)::result
      type(filter_result)::filtered
      real(dp)::lo,hi,mid,cdf,q,integral,x,dx
      integer::t,j,k,it,m
      filtered=hamilton_filter(spec,y)
      allocate(result%var(size(y),size(alpha)),result%es(size(y),size(alpha)))
      result%var(1,:)=0.0_dp;result%es(1,:)=0.0_dp
      do t=2,size(y)
         do j=1,size(alpha)
            lo=-20.0_dp*sqrt(maxval(filtered%variance(t,:)));hi=0.0_dp
            do it=1,90
               mid=0.5_dp*(lo+hi);cdf=0.0_dp
               do k=1,spec%k
                  cdf=cdf+filtered%predicted(t,k)*innovation_cdf(mid/sqrt(filtered%variance(t,k)), &
                     spec%regime(k)%distribution,spec%regime(k)%shape,spec%regime(k)%skew)
               end do
               if(cdf<alpha(j))then;lo=mid;else;hi=mid;end if
            end do
            q=0.5_dp*(lo+hi);result%var(t,j)=q
            m=1000;dx=(q-(-20.0_dp*sqrt(maxval(filtered%variance(t,:)))))/real(m,dp);integral=0.0_dp
            do it=0,m
               x=-20.0_dp*sqrt(maxval(filtered%variance(t,:)))+real(it,dp)*dx;cdf=0.0_dp
               do k=1,spec%k
                  cdf=cdf+filtered%predicted(t,k)*innovation_cdf(x/sqrt(filtered%variance(t,k)), &
                     spec%regime(k)%distribution,spec%regime(k)%shape,spec%regime(k)%skew)
               end do
               if(it==0.or.it==m)then;integral=integral+cdf;else if(mod(it,2)==0)then;integral=integral+2.0_dp*cdf;else;integral=integral+4.0_dp*cdf;end if
            end do
            integral=integral*dx/3.0_dp
            result%es(t,j)=q-integral/max(alpha(j),tiny(1.0_dp))
         end do
      end do
   end function risk_in_sample
end module msgarch_risk

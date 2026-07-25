! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module msgarch_posterior
   use msgarch_kinds, only : dp
   use msgarch_types, only : msgarch_spec, posterior_state_result, filter_result
   use msgarch_models, only : unpack_parameters
   use msgarch_filter, only : hamilton_filter, predictive_pdf, predictive_cdf, conditional_volatility, pit_values
   use msgarch_simulation, only : unconditional_volatility
   implicit none
   private
   public :: posterior_mean_spec, posterior_state_probabilities
   public :: posterior_predictive_pdf, posterior_predictive_cdf
   public :: posterior_volatility, posterior_pit, posterior_unconditional_volatility
contains
   function posterior_mean_spec(template,draws) result(spec)
      type(msgarch_spec),intent(in)::template
      real(dp),intent(in)::draws(:,:)
      type(msgarch_spec)::spec
      real(dp),allocatable::mean_theta(:)
      logical::valid
      mean_theta=sum(draws,dim=1)/real(size(draws,1),dp)
      call unpack_parameters(template,mean_theta,spec,valid)
      if(.not.valid)spec=template
   end function posterior_mean_spec

   function posterior_state_probabilities(template,y,draws) result(result)
      type(msgarch_spec),intent(in)::template
      real(dp),intent(in)::y(:),draws(:,:)
      type(posterior_state_result)::result
      type(msgarch_spec)::spec
      type(filter_result)::filtered
      logical::valid
      integer::i,nvalid
      allocate(result%predicted(size(y),template%k),result%filtered(size(y),template%k),result%smoothed(size(y),template%k))
      result%predicted=0.0_dp;result%filtered=0.0_dp;result%smoothed=0.0_dp;nvalid=0
      do i=1,size(draws,1)
         call unpack_parameters(template,draws(i,:),spec,valid)
         if(.not.valid)cycle
         filtered=hamilton_filter(spec,y)
         result%predicted=result%predicted+filtered%predicted
         result%filtered=result%filtered+filtered%filtered
         result%smoothed=result%smoothed+filtered%smoothed
         nvalid=nvalid+1
      end do
      if(nvalid>0)then
         result%predicted=result%predicted/real(nvalid,dp)
         result%filtered=result%filtered/real(nvalid,dp)
         result%smoothed=result%smoothed/real(nvalid,dp)
      end if
   end function posterior_state_probabilities

   function posterior_predictive_pdf(template,y,draws,x) result(value)
      type(msgarch_spec),intent(in)::template
      real(dp),intent(in)::y(:),draws(:,:),x
      real(dp)::value
      type(msgarch_spec)::spec
      logical::valid
      integer::i,nvalid
      value=0.0_dp;nvalid=0
      do i=1,size(draws,1)
         call unpack_parameters(template,draws(i,:),spec,valid)
         if(valid)then;value=value+predictive_pdf(spec,y,x);nvalid=nvalid+1;end if
      end do
      if(nvalid>0)value=value/real(nvalid,dp)
   end function posterior_predictive_pdf

   function posterior_predictive_cdf(template,y,draws,x) result(value)
      type(msgarch_spec),intent(in)::template
      real(dp),intent(in)::y(:),draws(:,:),x
      real(dp)::value
      type(msgarch_spec)::spec
      logical::valid
      integer::i,nvalid
      value=0.0_dp;nvalid=0
      do i=1,size(draws,1)
         call unpack_parameters(template,draws(i,:),spec,valid)
         if(valid)then;value=value+predictive_cdf(spec,y,x);nvalid=nvalid+1;end if
      end do
      if(nvalid>0)value=value/real(nvalid,dp)
   end function posterior_predictive_cdf

   function posterior_volatility(template,y,draws) result(volatility)
      type(msgarch_spec),intent(in)::template
      real(dp),intent(in)::y(:),draws(:,:)
      real(dp),allocatable::volatility(:),current(:)
      type(msgarch_spec)::spec
      logical::valid
      integer::i,nvalid
      allocate(volatility(size(y)));volatility=0.0_dp;nvalid=0
      do i=1,size(draws,1)
         call unpack_parameters(template,draws(i,:),spec,valid)
         if(.not.valid)cycle
         current=conditional_volatility(spec,y);volatility=volatility+current;nvalid=nvalid+1
      end do
      if(nvalid>0)volatility=volatility/real(nvalid,dp)
   end function posterior_volatility

   function posterior_unconditional_volatility(template,draws,nsim,n_ahead,burn) result(value)
      type(msgarch_spec),intent(in)::template
      real(dp),intent(in)::draws(:,:)
      integer,intent(in),optional::nsim,n_ahead,burn
      real(dp)::value
      type(msgarch_spec)::spec
      logical::valid
      integer::i,nvalid,m,horizon,discard
      m=250;if(present(nsim))m=nsim
      horizon=5000;if(present(n_ahead))horizon=n_ahead
      discard=1000;if(present(burn))discard=burn
      value=0.0_dp;nvalid=0
      do i=1,size(draws,1)
         call unpack_parameters(template,draws(i,:),spec,valid)
         if(.not.valid)cycle
         value=value+unconditional_volatility(spec,m,horizon,discard);nvalid=nvalid+1
      end do
      if(nvalid>0)value=value/real(nvalid,dp)
   end function posterior_unconditional_volatility

   function posterior_pit(template,y,draws) result(pit)
      type(msgarch_spec),intent(in)::template
      real(dp),intent(in)::y(:),draws(:,:)
      real(dp),allocatable::pit(:),current(:)
      type(msgarch_spec)::spec
      logical::valid
      integer::i,nvalid
      allocate(pit(size(y)));pit=0.0_dp;nvalid=0
      do i=1,size(draws,1)
         call unpack_parameters(template,draws(i,:),spec,valid)
         if(.not.valid)cycle
         current=pit_values(spec,y);pit=pit+current;nvalid=nvalid+1
      end do
      if(nvalid>0)pit=pit/real(nvalid,dp)
   end function posterior_pit
end module msgarch_posterior

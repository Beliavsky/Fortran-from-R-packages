! SPDX-License-Identifier: GPL-3.0-only
! Adapters to the translated coda and quantreg dependencies supplied with this translation.
module mcmcpack_dependencies
   use mcmcpack_kinds, only : dp
   use mcmcpack_samplers, only : mcmc_result
   use coda, only : mcmc_chain,make_mcmc
   use quantreg, only : rq_result,rq_fit_fnb
   use mcmc, only : mcmc_scale, metrop_result, scale_constant, scale_diagonal, scale_full, metrop, set_mcmc_seed
   implicit none
   private
   public :: as_coda_chain, quantreg_start
   public :: mcmc_scale, metrop_result, scale_constant, scale_diagonal, scale_full, metrop, set_mcmc_seed
contains
   function as_coda_chain(result,start,thin) result(chain)
      type(mcmc_result),intent(in)::result
      integer,intent(in),optional::start,thin
      type(mcmc_chain)::chain
      integer::s,t
      s=1;t=1;if(present(start))s=start;if(present(thin))t=thin
      if(.not.allocated(result%draws))then
         chain=make_mcmc(reshape([0.0_dp],[1,1]),start=s,thin=t)
      else
         chain=make_mcmc(result%draws,start=s,thin=t)
      end if
   end function as_coda_chain

   subroutine quantreg_start(x,y,tau,beta,status)
      real(dp),intent(in)::x(:,:),y(:),tau
      real(dp),intent(out)::beta(size(x,2))
      integer,intent(out)::status
      type(rq_result)::fit
      call rq_fit_fnb(x,y,tau,fit)
      status=fit%info
      if(status==0)then;beta=fit%coefficients;else;beta=0.0_dp;end if
   end subroutine quantreg_start
end module mcmcpack_dependencies

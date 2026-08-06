! SPDX-License-Identifier: GPL-2.0-only
module tsgarch_profile_module
  use ghyp_kinds, only : dp
  use tsgarch_types
  use tsgarch_fit_module, only : unpack_parameters, garch_loglikelihood
  implicit none
  private
  public :: profile_likelihood
contains
  function profile_likelihood(y,fit,parameter_index,grid,vreg) result(out)
    real(dp),intent(in)::y(:),grid(:)
    type(garch_fit),intent(in)::fit
    integer,intent(in)::parameter_index
    real(dp),intent(in),optional::vreg(:,:)
    type(profile_result)::out
    real(dp),allocatable::theta(:)
    type(garch_parameters)::par
    integer::i,status
    if(parameter_index<1.or.parameter_index>fit%npars.or.size(grid)<1)then
    out%message='invalid parameter index or grid'
    return
    end if
    out%parameter_name=fit%parameter_names(parameter_index)
    allocate(out%grid(size(grid)))
    out%grid=grid
    allocate(out%log_likelihood(size(grid)))
    do i=1,size(grid)
      theta=fit%packed_parameters
      theta(parameter_index)=grid(i)
      call unpack_parameters(fit%spec,theta,fit%parameters,par,status)
      if(status/=tsg_success)then
      out%log_likelihood(i)=-huge(1.0_dp)
      cycle
      end if
      if(present(vreg))then
      out%log_likelihood(i)=garch_loglikelihood(y,fit%spec,par,vreg)
      else
      out%log_likelihood(i)=garch_loglikelihood(y,fit%spec,par)
      end if
    end do
    out%status=tsg_success
    out%message='ok'
  end function profile_likelihood
end module tsgarch_profile_module

! High-level robust and quasi-likelihood iteration based on R locfit.robust/locfit.quasi.
! GPL-2-or-later; see LICENSE_NOTICE and upstream/locfit-R.
module locfit_robust
  use locfit_kinds, only : dp
  use locfit_constants, only : lf_ok
  use locfit_core, only : locfit_options, locfit_result, locfit_fit
  use locfit_diagnostics, only : median_value
  implicit none
  private
  public :: locfit_robust_fit, locfit_quasi_fit
contains
  subroutine locfit_robust_fit(x,y,x_eval,result,options,iterations,base,censored,scale,style)
    real(dp),intent(in)::x(:,:),y(:),x_eval(:,:)
    type(locfit_result),intent(out)::result
    type(locfit_options),intent(in),optional::options
    integer,intent(in),optional::iterations
    real(dp),intent(in),optional::base(:),scale(:)
    logical,intent(in),optional::censored(:)
    integer,intent(in),optional::style(:)
    type(locfit_options)::opts
    type(locfit_result)::at_data
    real(dp),allocatable::rw(:),res(:),ab(:)
    real(dp)::s
    integer::it,niter
    opts=locfit_options();if(present(options))opts=options
    niter=3;if(present(iterations))niter=max(0,iterations)
    allocate(rw(size(y)),res(size(y)),ab(size(y)));rw=1.0_dp
    do it=0,niter
      call locfit_fit(x,y,x,at_data,opts,prior_weights=rw,base=base,censored=censored,scale=scale,style=style)
      if(it<niter)then
        res=y-at_data%fit
        ab=abs(res);s=median_value(ab)
        if(s<=sqrt(tiny(1.0_dp)))exit
        rw=max(1.0_dp-(res/(6.0_dp*s))**2,0.0_dp)**2
      end if
    end do
    call locfit_fit(x,y,x_eval,result,opts,prior_weights=rw,base=base,censored=censored,scale=scale,style=style)
  end subroutine locfit_robust_fit

  subroutine locfit_quasi_fit(x,y,x_eval,result,variance,options,iterations,weights,base,censored,scale,style)
    real(dp),intent(in)::x(:,:),y(:),x_eval(:,:)
    type(locfit_result),intent(out)::result
    interface
      pure function variance(mu) result(v)
        import dp
        real(dp),intent(in)::mu
        real(dp)::v
      end function variance
    end interface
    type(locfit_options),intent(in),optional::options
    integer,intent(in),optional::iterations
    real(dp),intent(in),optional::weights(:),base(:),scale(:)
    logical,intent(in),optional::censored(:)
    integer,intent(in),optional::style(:)
    type(locfit_options)::opts
    type(locfit_result)::at_data
    real(dp),allocatable::w0(:),qw(:)
    real(dp)::vv
    integer::it,i,niter
    opts=locfit_options();if(present(options))opts=options
    niter=3;if(present(iterations))niter=max(0,iterations)
    allocate(w0(size(y)),qw(size(y)));w0=1.0_dp;if(present(weights))w0=weights;qw=w0
    do it=0,niter
      call locfit_fit(x,y,x,at_data,opts,prior_weights=qw,base=base,censored=censored,scale=scale,style=style)
      if(it<niter)then
        do i=1,size(y)
          vv=variance(at_data%fit(i))
          if(vv>sqrt(tiny(1.0_dp)))then
            qw(i)=w0(i)/vv
          else
            qw(i)=0.0_dp
          end if
        end do
      end if
    end do
    call locfit_fit(x,y,x_eval,result,opts,prior_weights=qw,base=base,censored=censored,scale=scale,style=style)
  end subroutine locfit_quasi_fit
end module locfit_robust

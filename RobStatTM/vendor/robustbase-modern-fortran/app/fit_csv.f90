! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
program fit_csv
   use robustbase
   implicit none
   character(len=1024)::filename,method
   real(dp),allocatable::xv(:),y(:),x(:,:)
   type(robust_regression_result)::fit
   type(fast_lts_result)::fast_fit
   type(partitioned_lts_result)::part_lts
   type(by_logistic_result)::by_fit
   type(lmrob_result)::lm_fit
   type(glmrob_result)::glm_fit
   integer::n
   if(command_argument_count()<1)then
      write(*,'(a)')'usage: fit_csv FILE [lts|fastlts|partlts|mm|lar|lmrob|smdm|by|mqle-binomial|mqle-poisson|mt]'
      error stop 2
   end if
   call get_command_argument(1,filename)
   method='mm';if(command_argument_count()>=2)call get_command_argument(2,method)
   call read_two_column_csv(trim(filename),xv,y)
   n=size(y);allocate(x(n,2));x(:,1)=1.0_dp;x(:,2)=xv
   write(*,'(a,1x,a)')'method',trim(method)
   write(*,'(a,1x,i0)')'nobs',n
   select case(trim(method))
   case('lts')
      call lts_regression(x,y,fit,n_starts=250)
      call print_standard_fit(fit)
   case('fastlts')
      call fast_lts_regression(x,y,fast_fit,alpha=0.5_dp,sampling='best',n_starts=500)
      call print_fast_fit(fast_fit)
   case('partlts')
      call fast_lts_partitioned(x,y,part_lts,alpha=0.5_dp,n_partitions=5,n_starts=100)
      call print_fast_fit(part_lts%estimate)
      write(*,'(a,1x,i0)')'partitions',part_lts%partitions
   case('mm')
      call mm_regression(x,y,fit,n_starts=200)
      call print_standard_fit(fit)
   case('lar')
      call lmrob_lar_fit(x,y,lm_fit)
      call print_lmrob_fit(lm_fit)
   case('lmrob')
      call lmrob_fit(x,y,lm_fit,method='MM',n_resample=500,sampling='nonsingular')
      call print_lmrob_fit(lm_fit)
   case('smdm')
      call lmrob_fit(x,y,lm_fit,method='SMDM',n_resample=500,sampling='nonsingular')
      call print_lmrob_fit(lm_fit)
   case('by')
      call require_binary(y,'BY')
      call by_logistic_fit(x,y,by_fit)
      write(*,'(a,*(1x,es18.10))')'coefficients',by_fit%coefficients
      write(*,'(a,*(1x,es18.10))')'standard_errors',by_fit%standard_errors
      write(*,'(a,1x,es18.10)')'objective',by_fit%objective
      write(*,'(a,1x,l1)')'converged',by_fit%converged
   case('mqle-binomial')
      call require_binary(y,'Mqle binomial')
      call glmrob_mqle_fit(x,y,'binomial',glm_fit)
      call print_glm_fit(glm_fit)
   case('mqle-poisson')
      if(any(y<0.0_dp))error stop 'fit_csv: Poisson response must be nonnegative'
      call glmrob_mqle_fit(x,y,'poisson',glm_fit)
      call print_glm_fit(glm_fit)
   case('mt')
      call require_binary(y,'MT')
      call glmrob_mt_fit(x,y,glm_fit)
      call print_glm_fit(glm_fit)
   case default
      error stop 'fit_csv: unsupported method'
   end select
contains
   subroutine require_binary(values,name)
      real(dp),intent(in)::values(:)
      character(len=*),intent(in)::name
      if(any(values<0.0_dp) .or. any(values>1.0_dp) .or. any(abs(values-real(nint(values),dp))>10.0_dp*epsilon(1.0_dp))) &
         error stop 'fit_csv: '//trim(name)//' response must contain only zero and one'
   end subroutine require_binary
   subroutine print_standard_fit(model)
      type(robust_regression_result),intent(in)::model
      write(*,'(a,*(1x,es18.10))')'coefficients',model%coefficients
      write(*,'(a,1x,es18.10)')'scale',model%scale
      write(*,'(a,1x,l1)')'converged',model%converged
   end subroutine print_standard_fit
   subroutine print_fast_fit(model)
      type(fast_lts_result),intent(in)::model
      write(*,'(a,*(1x,es18.10))')'coefficients',model%coefficients
      write(*,'(a,1x,es18.10)')'scale',model%scale
      write(*,'(a,1x,es18.10)')'raw_objective',model%objective
      write(*,'(a,1x,i0)')'h',model%h
      write(*,'(a,1x,l1)')'converged',model%converged
   end subroutine print_fast_fit
   subroutine print_lmrob_fit(model)
      type(lmrob_result),intent(in)::model
      write(*,'(a,*(1x,es18.10))')'coefficients',model%coefficients
      write(*,'(a,*(1x,es18.10))')'standard_errors',model%standard_errors
      write(*,'(a,1x,es18.10)')'scale',model%scale
      write(*,'(a,1x,a)')'chain',trim(model%method)
      write(*,'(a,1x,l1)')'converged',model%converged
   end subroutine print_lmrob_fit
   subroutine print_glm_fit(model)
      type(glmrob_result),intent(in)::model
      write(*,'(a,*(1x,es18.10))')'coefficients',model%coefficients
      write(*,'(a,*(1x,es18.10))')'standard_errors',model%standard_errors
      write(*,'(a,1x,es18.10)')'objective',model%objective
      write(*,'(a,1x,l1)')'converged',model%converged
   end subroutine print_glm_fit
end program fit_csv

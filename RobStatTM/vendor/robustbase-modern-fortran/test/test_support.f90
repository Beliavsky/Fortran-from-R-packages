! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module test_support
   use robustbase_kinds, only: dp
   implicit none
   private
   public :: assert_close, assert_true, seed_rng, exp_model
contains
   subroutine assert_close(actual,expected,tol,message)
      real(dp),intent(in)::actual,expected,tol
      character(len=*),intent(in)::message
      if(abs(actual-expected)>tol) then
         write(*,'(a,2(1x,es24.15),1x,a,1x,es12.4)') trim(message)//':',actual,expected,'tol',tol
         error stop 1
      end if
   end subroutine
   subroutine assert_true(value,message)
      logical,intent(in)::value
      character(len=*),intent(in)::message
      if(.not.value) then;write(*,'(a)')trim(message);error stop 1;end if
   end subroutine
   subroutine seed_rng(seed)
      integer,intent(in)::seed
      integer::n,i
      integer,allocatable::put(:)
      call random_seed(size=n);allocate(put(n));put=[(seed+104729*i,i=1,n)];call random_seed(put=put)
   end subroutine
   subroutine exp_model(theta,xx,yhat)
      real(dp),intent(in)::theta(:),xx(:,:)
      real(dp),intent(out)::yhat(:)
      yhat=theta(1)*exp(theta(2)*xx(:,1))
   end subroutine exp_model
end module test_support

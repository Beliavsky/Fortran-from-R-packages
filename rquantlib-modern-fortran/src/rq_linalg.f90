! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
module rq_linalg
  use rq_kinds, only: dp
  implicit none
  private
  public :: least_squares_solve
  interface
    subroutine dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
      import dp
      character(len=1),intent(in)::trans
      integer,intent(in)::m,n,nrhs,lda,ldb,lwork
      real(dp),intent(inout)::a(lda,*),b(ldb,*),work(*)
      integer,intent(out)::info
    end subroutine dgels
  end interface
contains
  subroutine least_squares_solve(x,y,beta,status)
    real(dp),intent(in)::x(:,:),y(:)
    real(dp),intent(out)::beta(:)
    integer,intent(out),optional::status
    real(dp),allocatable::a(:,:),b(:,:),work(:)
    integer::m,n,nrhs,lda,ldb,lwork,info
    m=size(x,1); n=size(x,2); nrhs=1; lda=max(1,m); ldb=max(m,n)
    allocate(a(lda,n),b(ldb,1)); a=0.0_dp; a(1:m,:)=x; b=0.0_dp; b(1:m,1)=y
    lwork=max(1,4*max(m,n)); allocate(work(lwork))
    call dgels('N',m,n,nrhs,a,lda,b,ldb,work,lwork,info)
    if(info==0) then; beta=b(1:n,1); else; beta=0.0_dp; end if
    if(present(status)) status=info
  end subroutine least_squares_solve
end module rq_linalg

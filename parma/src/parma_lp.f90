! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! Native simplex support replacing the upstream Rglpk binding for continuous LPs.
module parma_lp
   use parma_kinds, only: dp
   use parma_types, only: lp_result
   implicit none
   private
   public :: lp_simplex

contains

   subroutine lp_simplex(c,a,b,result,maximize,max_iter,tol)
      real(dp), intent(in) :: c(:),a(:,:),b(:)
      type(lp_result), intent(out) :: result
      logical, intent(in), optional :: maximize
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: tableau(:,:),cost(:)
      real(dp) :: epsilonx,ratio,best_ratio,pivot
      integer :: m,n,i,j,row,col,iter,maxit
      logical :: do_max

      m = size(b)
      n = size(c)
      maxit = 10000
      if (present(max_iter)) maxit = max_iter
      epsilonx = 1.0e-10_dp
      if (present(tol)) epsilonx = tol
      do_max = .true.
      if (present(maximize)) do_max = maximize
      allocate(result%x(n),tableau(m+1,n+m+1),cost(n))
      result%x = 0.0_dp
      if (size(a,1) /= m .or. size(a,2) /= n .or. any(b < -epsilonx)) then
         result%status = 2
         result%message = 'simplex requires A*x <= b, b >= 0, x >= 0'
         return
      end if
      cost = c
      if (.not. do_max) cost = -cost
      tableau = 0.0_dp
      tableau(1:m,1:n) = a
      do i = 1,m
         tableau(i,n+i) = 1.0_dp
      end do
      tableau(1:m,n+m+1) = b
      tableau(m+1,1:n) = -cost
      do iter = 1,maxit
         col = 0
         do j = 1,n+m
            if (tableau(m+1,j) < -epsilonx) then
               if (col == 0) then
                  col = j
               else if (tableau(m+1,j) < tableau(m+1,col)) then
                  col = j
               end if
            end if
         end do
         if (col == 0) exit
         row = 0
         best_ratio = huge(1.0_dp)
         do i = 1,m
            if (tableau(i,col) > epsilonx) then
               ratio = tableau(i,n+m+1)/tableau(i,col)
               if (ratio < best_ratio) then
                  best_ratio = ratio
                  row = i
               end if
            end if
         end do
         if (row == 0) then
            result%status = 3
            result%message = 'unbounded linear program'
            return
         end if
         pivot = tableau(row,col)
         tableau(row,:) = tableau(row,:)/pivot
         do i = 1,m+1
            if (i /= row) tableau(i,:) = tableau(i,:)-tableau(i,col)*tableau(row,:)
         end do
      end do
      do j = 1,n
         row = 0
         do i = 1,m
            if (abs(tableau(i,j)-1.0_dp) <= epsilonx .and. &
                count(abs(tableau(1:m,j)) > epsilonx) == 1) then
               row = i
               exit
            end if
         end do
         if (row > 0) result%x(j) = tableau(row,n+m+1)
      end do
      result%objective = dot_product(c,result%x)
      result%iterations = min(iter,maxit)
      if (iter > maxit) then
         result%status = 1
         result%message = 'iteration limit'
      else
         result%status = 0
         result%message = 'optimal'
      end if
   end subroutine lp_simplex

end module parma_lp

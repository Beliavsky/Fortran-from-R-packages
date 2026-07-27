! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! Matrix, simulation, and compatibility helpers derived from parma 1.7.
module parma_helpers
   use parma_kinds, only: dp
   use parma_rng, only: seed_rng
   use parma_linalg, only: project_box_budget, jacobi_eigen
   use parma_risk, only: empirical_quantile
   implicit none
   private
   public :: eye, ones_vector, zeros_vector, repmat, tril, triu
   public :: condition_number, is_even, log_transform, simweights
   public :: percentile_values, x_into_bounds

contains

   function eye(n,m) result(a)
      integer, intent(in) :: n
      integer, intent(in), optional :: m
      integer :: ncol,i
      real(dp), allocatable :: a(:,:)
      ncol = n
      if (present(m)) ncol = m
      allocate(a(n,ncol))
      a = 0.0_dp
      do i = 1,min(n,ncol)
         a(i,i) = 1.0_dp
      end do
   end function eye

   function ones_vector(n) result(x)
      integer, intent(in) :: n
      real(dp) :: x(n)
      x = 1.0_dp
   end function ones_vector

   function zeros_vector(n) result(x)
      integer, intent(in) :: n
      real(dp) :: x(n)
      x = 0.0_dp
   end function zeros_vector

   function repmat(a,nrow,ncol) result(out)
      real(dp), intent(in) :: a(:,:)
      integer, intent(in) :: nrow,ncol
      real(dp) :: out(size(a,1)*nrow,size(a,2)*ncol)
      integer :: i,j,r1,r2,c1,c2
      do i = 1,nrow
         r1 = (i-1)*size(a,1)+1
         r2 = i*size(a,1)
         do j = 1,ncol
            c1 = (j-1)*size(a,2)+1
            c2 = j*size(a,2)
            out(r1:r2,c1:c2) = a
         end do
      end do
   end function repmat

   function tril(a,k) result(out)
      real(dp), intent(in) :: a(:,:)
      integer, intent(in), optional :: k
      real(dp) :: out(size(a,1),size(a,2))
      integer :: offset,i,j
      offset = 0
      if (present(k)) offset = k
      out = 0.0_dp
      do i = 1,size(a,1)
         do j = 1,size(a,2)
            if (j-i <= offset) out(i,j) = a(i,j)
         end do
      end do
   end function tril

   function triu(a,k) result(out)
      real(dp), intent(in) :: a(:,:)
      integer, intent(in), optional :: k
      real(dp) :: out(size(a,1),size(a,2))
      integer :: offset,i,j
      offset = 0
      if (present(k)) offset = k
      out = 0.0_dp
      do i = 1,size(a,1)
         do j = 1,size(a,2)
            if (j-i >= offset) out(i,j) = a(i,j)
         end do
      end do
   end function triu

   function condition_number(a) result(value)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: value
      real(dp), allocatable :: ata(:,:),eig(:),vectors(:,:)
      integer :: n,info
      n = size(a,2)
      allocate(ata(n,n),eig(n),vectors(n,n))
      ata = matmul(transpose(a),a)
      call jacobi_eigen(ata,eig,vectors,info)
      if (info < 0 .or. minval(eig) <= tiny(1.0_dp)) then
         value = huge(1.0_dp)
      else
         value = sqrt(maxval(eig)/minval(eig))
      end if
   end function condition_number

   elemental function is_even(x) result(answer)
      real(dp), intent(in) :: x
      logical :: answer
      answer = modulo(nint(x),2) == 0
   end function is_even

   elemental function log_transform(x,lb,ub) result(value)
      real(dp), intent(in) :: x,lb,ub
      real(dp) :: value
      value = lb+(ub-lb)/(1.0_dp+exp(-x))
   end function log_transform

   function x_into_bounds(x,lb,ub) result(y)
      real(dp), intent(in) :: x(:),lb(:),ub(:)
      real(dp) :: y(size(x)),width
      integer :: i,count
      do i = 1,size(x)
         width = ub(i)-lb(i)
         if (width <= 0.0_dp) then
            y(i) = lb(i)
         else
            y(i) = x(i)
            count = 0
            do while ((y(i) < lb(i) .or. y(i) > ub(i)) .and. count < 50)
               if (y(i) < lb(i)) y(i) = 2.0_dp*lb(i)-y(i)
               if (y(i) > ub(i)) y(i) = 2.0_dp*ub(i)-y(i)
               count = count+1
            end do
            y(i) = min(max(y(i),lb(i)),ub(i))
         end if
      end do
   end function x_into_bounds

   subroutine simweights(weights,lb,ub,budget,forecast,constrain_positive,seed,info)
      real(dp), intent(out) :: weights(:)
      real(dp), intent(in) :: lb(:),ub(:),budget
      real(dp), intent(in), optional :: forecast(:)
      logical, intent(in), optional :: constrain_positive
      integer, intent(in), optional :: seed
      integer, intent(out), optional :: info
      real(dp), allocatable :: trial(:),projected(:)
      logical :: require_positive
      integer :: iter,ierr

      allocate(trial(size(weights)),projected(size(weights)))
      require_positive = .true.
      if (present(constrain_positive)) require_positive = constrain_positive
      if (present(seed)) then
         call seed_rng(seed)
      else
         call seed_rng(12345)
      end if
      ierr = 1
      do iter = 1,10000
         call random_number(trial)
         trial = lb+(ub-lb)*trial
         call project_box_budget(trial,lb,ub,budget,projected,ierr)
         if (ierr /= 0) exit
         if (require_positive .and. present(forecast)) then
            if (dot_product(projected,forecast) < 0.0_dp) cycle
         end if
         weights = projected
         ierr = 0
         exit
      end do
      if (present(info)) info = ierr
   end subroutine simweights

   function percentile_values(x,percentages) result(values)
      real(dp), intent(in) :: x(:),percentages(:)
      real(dp) :: values(size(percentages))
      integer :: i
      do i = 1,size(percentages)
         values(i) = empirical_quantile(x,percentages(i)/100.0_dp)
      end do
   end function percentile_values

end module parma_helpers

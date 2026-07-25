! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module gogarch_ica
   use gogarch_kinds, only : dp
   use gogarch_linalg, only : covariance_matrix, symmetric_sqrt, symmetric_invsqrt, identity_matrix
   use gogarch_orthogonal, only : uprod_r
   implicit none
   private

   type, public :: ica_result
      real(dp), allocatable :: rotation(:,:)
      real(dp), allocatable :: mixing(:,:)
      real(dp), allocatable :: sources(:,:)
      real(dp), allocatable :: whitening(:,:)
      real(dp), allocatable :: covariance_sqrt(:,:)
      integer :: iterations = 0
      integer :: status = 1
   end type ica_result

   public :: fastica

contains

   function fastica(data, max_iterations, tolerance) result(result)
      real(dp), intent(in) :: data(:,:)
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(ica_result) :: result
      real(dp) :: covariance(size(data,2),size(data,2)), invroot(size(data,2),size(data,2))
      real(dp) :: covroot(size(data,2),size(data,2)), whitened(size(data,1),size(data,2))
      real(dp) :: w(size(data,2),size(data,2)), wnew(size(data,2),size(data,2))
      real(dp) :: y(size(data,1),size(data,2)), g(size(data,1),size(data,2))
      real(dp) :: gp(size(data,1),size(data,2)), decor(size(data,2),size(data,2))
      real(dp) :: check(size(data,2),size(data,2)), theta(max(1,size(data,2)*(size(data,2)-1)/2))
      real(dp) :: tol, convergence
      logical :: ok1, ok2, ok3
      integer :: n, m, maxit, iter, j
      n = size(data,1)
      m = size(data,2)
      maxit = 500
      if (present(max_iterations)) maxit = max(1,max_iterations)
      tol = 1.0e-7_dp
      if (present(tolerance)) tol = max(tolerance,epsilon(1.0_dp))
      allocate(result%rotation(m,m),result%mixing(m,m),result%sources(n,m),result%whitening(m,m),result%covariance_sqrt(m,m))
      covariance = covariance_matrix(data,center=.false.)
      invroot = symmetric_invsqrt(covariance,1.0e-12_dp,ok1)
      covroot = symmetric_sqrt(covariance,1.0e-12_dp,ok2)
      if (.not. (ok1 .and. ok2)) then
         result%rotation = 0.0_dp
         result%mixing = 0.0_dp
         result%sources = 0.0_dp
         result%whitening = 0.0_dp
         result%covariance_sqrt = 0.0_dp
         result%status = 2
         return
      end if
      whitened = matmul(data,invroot)
      if (m == 1) then
         w = identity_matrix(1)
      else
         do j = 1, m*(m-1)/2
            theta(j) = 0.13_dp+0.07_dp*real(j,dp)
         end do
         w = transpose(uprod_r(theta(1:m*(m-1)/2),ok3))
         if (.not. ok3) w = identity_matrix(m)
      end if
      do iter = 1, maxit
         y = matmul(whitened,transpose(w))
         g = tanh(y)
         gp = 1.0_dp-g*g
         wnew = matmul(transpose(g),whitened)/real(n,dp)
         do j = 1, m
            wnew(j,:) = wnew(j,:)-sum(gp(:,j))/real(n,dp)*w(j,:)
         end do
         decor = symmetric_invsqrt(matmul(wnew,transpose(wnew)),1.0e-12_dp,ok3)
         if (.not. ok3) exit
         wnew = matmul(decor,wnew)
         check = matmul(wnew,transpose(w))
         convergence = maxval(abs(abs([(check(j,j),j=1,m)])-1.0_dp))
         w = wnew
         if (convergence < tol) then
            result%status = 0
            exit
         end if
      end do
      result%iterations = min(iter,maxit)
      result%rotation = transpose(w)
      result%whitening = invroot
      result%covariance_sqrt = covroot
      result%sources = matmul(whitened,result%rotation)
      result%mixing = matmul(covroot,result%rotation)
      if (result%status /= 0 .and. iter > maxit) result%status = 1
   end function fastica

end module gogarch_ica

! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch_ica
   use rmgarch_kinds, only : dp
   use rmgarch_math, only : pi, covariance_matrix, inverse_sqrt_symmetric, &
      sqrt_symmetric, inverse_matrix, identity_matrix
   use rmgarch_rng, only : random_normal
   use rmgarch_types, only : ica_result
   implicit none
   private

   public :: fastica, radical

contains

   function fastica(x, max_iterations, tolerance, nonlinearity) result(result)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      character(len=*), intent(in), optional :: nonlinearity
      type(ica_result) :: result
      real(dp), allocatable :: centered(:,:), whitened(:,:), w(:,:), wnew(:,:), y(:,:), g(:,:), gp(:,:)
      real(dp), allocatable :: decor(:,:), dewhite(:,:), cross(:,:), covariance(:,:), temp(:,:)
      real(dp) :: tol, convergence
      integer :: n, m, maxit, i, j, iter
      logical :: ok
      character(len=16) :: gfun

      n = size(x,1)
      m = size(x,2)
      maxit = 500
      if (present(max_iterations)) maxit = max(1,max_iterations)
      tol = 1.0e-7_dp
      if (present(tolerance)) tol = max(tolerance,epsilon(1.0_dp))
      gfun = 'tanh'
      if (present(nonlinearity)) gfun = lowercase(adjustl(nonlinearity))

      if (n <= 1 .or. m < 1) then
         call allocate_empty_result(result,n,m,2)
         return
      end if

      allocate(result%center(m),centered(n,m),whitened(n,m),w(m,m),wnew(m,m),y(n,m),g(n,m),gp(n,m))
      allocate(decor(m,m),dewhite(m,m),cross(m,m),covariance(m,m),temp(m,m))
      result%center = sum(x,dim=1)/real(n,dp)
      centered = x-spread(result%center,1,n)
      covariance = covariance_matrix(centered)
      decor = inverse_sqrt_symmetric(covariance,1.0e-10_dp,ok)
      if (.not. ok) then
         call allocate_empty_arrays(result,n,m,2)
         return
      end if
      dewhite = sqrt_symmetric(covariance,0.0_dp,ok)
      if (.not. ok) then
         call allocate_empty_arrays(result,n,m,2)
         return
      end if
      whitened = matmul(centered,decor)

      do j = 1, m
         do i = 1, m
            w(i,j) = random_normal()
         end do
      end do
      w = symmetric_decorrelation(w,ok)
      if (.not. ok) then
         call allocate_empty_arrays(result,n,m,3)
         return
      end if

      convergence = huge(1.0_dp)
      do iter = 1, maxit
         y = matmul(whitened,transpose(w))
         call nonlinear_transform(y,gfun,g,gp)
         wnew = matmul(transpose(g),whitened)/real(n,dp)
         do i = 1, m
            wnew(i,:) = wnew(i,:)-sum(gp(:,i))/real(n,dp)*w(i,:)
         end do
         wnew = symmetric_decorrelation(wnew,ok)
         if (.not. ok) exit
         cross = matmul(wnew,transpose(w))
         convergence = 0.0_dp
         do i = 1, m
            convergence = max(convergence,abs(abs(cross(i,i))-1.0_dp))
         end do
         w = wnew
         if (convergence < tol) exit
      end do

      call finalize_ica_result(result,centered,covariance,decor,dewhite,w,iter,maxit, &
         ok .and. convergence < max(tol,1.0e-5_dp))
   end function fastica

   function radical(x, max_sweeps, angle_points, spacing, tolerance, demean) result(result)
      !! Pairwise-rotation RADICAL ICA using a Vasicek spacing-entropy objective.
      !! The implementation follows the computational idea of RADICAL: whiten,
      !! repeatedly optimize every two-dimensional rotation, and minimize the
      !! sum of marginal entropy estimates.
      real(dp), intent(in) :: x(:,:)
      integer, intent(in), optional :: max_sweeps, angle_points, spacing
      real(dp), intent(in), optional :: tolerance
      logical, intent(in), optional :: demean
      type(ica_result) :: result
      real(dp), allocatable :: centered(:,:), whitened(:,:), covariance(:,:), decor(:,:), dewhite(:,:)
      real(dp), allocatable :: w(:,:), col_i(:), col_j(:), trial_i(:), trial_j(:)
      real(dp) :: theta, best_theta, objective, best_objective, c, s, tol, max_change
      integer :: n, m, nsweep, nangle, mspace, sweep, i, j, k
      logical :: ok, remove_mean

      n = size(x,1)
      m = size(x,2)
      nsweep = 12
      if (present(max_sweeps)) nsweep = max(1,max_sweeps)
      nangle = 61
      if (present(angle_points)) nangle = max(5,angle_points)
      if (mod(nangle,2) == 0) nangle = nangle+1
      mspace = max(1,nint(sqrt(real(max(n,1),dp))))
      if (present(spacing)) mspace = max(1,spacing)
      tol = 2.0e-3_dp
      if (present(tolerance)) tol = max(tolerance,epsilon(1.0_dp))
      remove_mean = .true.
      if (present(demean)) remove_mean = demean

      if (n <= 2 .or. m < 1) then
         call allocate_empty_result(result,n,m,2)
         return
      end if
      mspace = min(mspace,max(1,(n-1)/2))

      allocate(result%center(m),centered(n,m),whitened(n,m),covariance(m,m),decor(m,m),dewhite(m,m))
      allocate(w(m,m),col_i(n),col_j(n),trial_i(n),trial_j(n))
      if (remove_mean) then
         result%center = sum(x,dim=1)/real(n,dp)
      else
         result%center = 0.0_dp
      end if
      centered = x-spread(result%center,1,n)
      covariance = covariance_matrix(centered)
      decor = inverse_sqrt_symmetric(covariance,1.0e-10_dp,ok)
      if (.not. ok) then
         call allocate_empty_arrays(result,n,m,2)
         return
      end if
      dewhite = sqrt_symmetric(covariance,0.0_dp,ok)
      if (.not. ok) then
         call allocate_empty_arrays(result,n,m,2)
         return
      end if
      whitened = matmul(centered,decor)
      w = identity_matrix(m)

      do sweep = 1, nsweep
         max_change = 0.0_dp
         do i = 1, m-1
            do j = i+1, m
               col_i = whitened(:,i)
               col_j = whitened(:,j)
               best_theta = 0.0_dp
               best_objective = entropy_spacing(col_i,mspace)+entropy_spacing(col_j,mspace)
               do k = 1, nangle
                  theta = -0.25_dp*pi+0.5_dp*pi*real(k-1,dp)/real(nangle-1,dp)
                  c = cos(theta)
                  s = sin(theta)
                  trial_i = c*col_i+s*col_j
                  trial_j = -s*col_i+c*col_j
                  objective = entropy_spacing(trial_i,mspace)+entropy_spacing(trial_j,mspace)
                  if (objective < best_objective) then
                     best_objective = objective
                     best_theta = theta
                  end if
               end do
               c = cos(best_theta)
               s = sin(best_theta)
               whitened(:,i) = c*col_i+s*col_j
               whitened(:,j) = -s*col_i+c*col_j
               call rotate_rows(w,i,j,c,s)
               max_change = max(max_change,abs(best_theta))
            end do
         end do
         if (max_change < tol) exit
      end do

      call finalize_ica_result(result,centered,covariance,decor,dewhite,w,sweep,nsweep,.true.)
   end function radical

   subroutine nonlinear_transform(y, name, g, gp)
      real(dp), intent(in) :: y(:,:)
      character(len=*), intent(in) :: name
      real(dp), intent(out) :: g(size(y,1),size(y,2)), gp(size(y,1),size(y,2))
      select case (trim(name))
      case ('pow3','cube','cubic')
         g = y**3
         gp = 3.0_dp*y*y
      case ('gauss','gaussian')
         g = y*exp(-0.5_dp*y*y)
         gp = (1.0_dp-y*y)*exp(-0.5_dp*y*y)
      case ('skew')
         g = y*y
         gp = 2.0_dp*y
      case default
         g = tanh(y)
         gp = 1.0_dp-g*g
      end select
   end subroutine nonlinear_transform

   subroutine finalize_ica_result(result, centered, covariance, decor, dewhite, w, iterations, maxit, converged)
      type(ica_result), intent(inout) :: result
      real(dp), intent(in) :: centered(:,:), covariance(:,:), decor(:,:), dewhite(:,:), w(:,:)
      integer, intent(in) :: iterations, maxit
      logical, intent(in) :: converged
      real(dp) :: temp(size(w,1),size(w,2))
      logical :: ok
      integer :: n, m

      n = size(centered,1)
      m = size(centered,2)
      allocate(result%unmixing(m,m),result%mixing(m,m),result%rotation(m,m),result%sources(n,m))
      allocate(result%covariance(m,m),result%whitening(m,m),result%dewhitening(m,m))
      result%rotation = w
      result%covariance = covariance
      result%whitening = decor
      result%dewhitening = dewhite
      result%unmixing = matmul(w,decor)
      temp = inverse_matrix(result%unmixing,ok)
      if (ok) then
         result%mixing = temp
      else
         result%mixing = 0.0_dp
      end if
      result%sources = matmul(centered,transpose(result%unmixing))
      result%iterations = min(iterations,maxit)
      result%status = merge(0,1,ok .and. converged)
   end subroutine finalize_ica_result

   subroutine allocate_empty_result(result,n,m,status)
      type(ica_result), intent(inout) :: result
      integer, intent(in) :: n, m, status
      allocate(result%center(m))
      result%center = 0.0_dp
      call allocate_empty_arrays(result,n,m,status)
   end subroutine allocate_empty_result

   subroutine allocate_empty_arrays(result,n,m,status)
      type(ica_result), intent(inout) :: result
      integer, intent(in) :: n, m, status
      if (.not. allocated(result%mixing)) allocate(result%mixing(m,m))
      if (.not. allocated(result%unmixing)) allocate(result%unmixing(m,m))
      if (.not. allocated(result%rotation)) allocate(result%rotation(m,m))
      if (.not. allocated(result%sources)) allocate(result%sources(n,m))
      if (.not. allocated(result%covariance)) allocate(result%covariance(m,m))
      if (.not. allocated(result%whitening)) allocate(result%whitening(m,m))
      if (.not. allocated(result%dewhitening)) allocate(result%dewhitening(m,m))
      result%mixing = 0.0_dp
      result%unmixing = 0.0_dp
      result%rotation = 0.0_dp
      result%sources = 0.0_dp
      result%covariance = 0.0_dp
      result%whitening = 0.0_dp
      result%dewhitening = 0.0_dp
      result%status = status
   end subroutine allocate_empty_arrays

   function symmetric_decorrelation(w, ok) result(out)
      real(dp), intent(in) :: w(:,:)
      logical, intent(out) :: ok
      real(dp) :: out(size(w,1),size(w,2)), invsqrt(size(w,1),size(w,1))
      invsqrt = inverse_sqrt_symmetric(matmul(w,transpose(w)),1.0e-12_dp,ok)
      if (ok) then
         out = matmul(invsqrt,w)
      else
         out = 0.0_dp
      end if
   end function symmetric_decorrelation

   subroutine rotate_rows(w,i,j,c,s)
      real(dp), intent(inout) :: w(:,:)
      integer, intent(in) :: i, j
      real(dp), intent(in) :: c, s
      real(dp) :: row_i(size(w,2)), row_j(size(w,2))
      row_i = w(i,:)
      row_j = w(j,:)
      w(i,:) = c*row_i+s*row_j
      w(j,:) = -s*row_i+c*row_j
   end subroutine rotate_rows

   function entropy_spacing(x, spacing) result(value)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: spacing
      real(dp) :: value, sorted(size(x)), width
      integer :: n, i, lo, hi, count

      n = size(x)
      if (n <= 2) then
         value = huge(1.0_dp)
         return
      end if
      sorted = x
      call sort_real(sorted)
      value = 0.0_dp
      count = 0
      do i = 1, n
         lo = max(1,i-spacing)
         hi = min(n,i+spacing)
         width = max(sorted(hi)-sorted(lo),1.0e-12_dp)
         value = value+log(real(n,dp)*width/real(max(hi-lo,1),dp))
         count = count+1
      end do
      value = value/real(count,dp)
   end function entropy_spacing

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      integer :: i, j
      real(dp) :: key
      do i = 2, size(x)
         key = x(i)
         j = i-1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j+1) = x(j)
            j = j-1
         end do
         x(j+1) = key
      end do
   end subroutine sort_real

   pure function lowercase(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i, code
      out = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) out(i:i) = achar(code+32)
      end do
   end function lowercase

end module rmgarch_ica

! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from ghyp 1.6.5 by Marc Weibel, David Luethi, and Henriette-Elise Breymann.
module ghyp_utilities
   use ghyp_kinds, only : dp
   use ghyp_model, only : ghyp_model_type, make_ghyp, transform_ghyp, &
      ghyp_moments, moments_result, model_gaussian
   use ghyp_distribution, only : qghyp
   use ghyp_linalg, only : inverse_spd, logdet_spd
   implicit none
   private

   type, public :: alpha_delta_result
      real(dp) :: lambda = 0.0_dp
      real(dp) :: alpha = 0.0_dp
      real(dp) :: delta = 0.0_dp
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: delta_matrix(:,:)
      logical :: ok = .false.
      character(len=160) :: message = ''
   end type alpha_delta_result

   type, public :: qq_result
      real(dp), allocatable :: theoretical(:)
      real(dp), allocatable :: sample(:)
      logical :: ok = .false.
      character(len=160) :: message = ''
   end type qq_result

   public :: subset_ghyp, standardize_ghyp, ghyp_alpha_delta, qqghyp_data

contains

   function subset_ghyp(model, indices) result(output)
      type(ghyp_model_type), intent(in) :: model
      integer, intent(in) :: indices(:)
      type(ghyp_model_type) :: output
      real(dp), allocatable :: mu(:), gamma(:), scatter(:,:)
      integer :: i, j, m
      m = size(indices)
      if (.not. model%ok .or. m < 1 .or. any(indices < 1) .or. &
          any(indices > model%dimension())) then
         output%message = 'invalid marginal indices'
         return
      end if
      allocate(mu(m),gamma(m),scatter(m,m))
      do i = 1, m
         mu(i) = model%mu(indices(i))
         gamma(i) = model%gamma(indices(i))
         do j = 1, m
            scatter(i,j) = model%scatter(indices(i),indices(j))
         end do
      end do
      output = make_ghyp(model%lambda,model%chi,model%psi,mu,scatter,gamma)
      if (model%family == model_gaussian .and. output%ok) output%family = model_gaussian
   end function subset_ghyp

   function standardize_ghyp(model, center, scale) result(output)
      type(ghyp_model_type), intent(in) :: model
      logical, intent(in), optional :: center, scale
      type(ghyp_model_type) :: output
      type(moments_result) :: moments
      real(dp), allocatable :: a(:,:), shift(:)
      logical :: do_center, do_scale
      integer :: i, d
      do_center = .true.; do_scale = .true.
      if (present(center)) do_center = center
      if (present(scale)) do_scale = scale
      moments = ghyp_moments(model)
      if (.not. moments%ok) then
         output%message = 'model moments are unavailable'
         return
      end if
      d = model%dimension()
      allocate(a(d,d),shift(d))
      a = 0.0_dp; shift = 0.0_dp
      do i = 1, d
         if (do_scale) then
            if (moments%covariance(i,i) <= 0.0_dp) then
               output%message = 'nonpositive marginal variance'
               return
            end if
            a(i,i) = 1.0_dp/sqrt(moments%covariance(i,i))
         else
            a(i,i) = 1.0_dp
         end if
      end do
      if (do_center) shift = -matmul(a,moments%mean)
      output = transform_ghyp(model,a,shift)
   end function standardize_ghyp

   function ghyp_alpha_delta(model) result(result)
      type(ghyp_model_type), intent(in) :: model
      type(alpha_delta_result) :: result
      real(dp), allocatable :: inverse(:,:)
      real(dp) :: logdet, detroot, psi_term
      logical :: ok
      integer :: d
      if (.not. model%ok .or. model%family == model_gaussian) then
         result%message = 'alpha-delta parameters require a non-Gaussian model'
         return
      end if
      d = model%dimension()
      allocate(result%beta(d),result%mu(d),result%delta_matrix(d,d))
      call inverse_spd(model%scatter,inverse,ok)
      if (.not. ok) then
         result%message = 'singular scatter matrix'
         return
      end if
      logdet = logdet_spd(model%scatter,ok)
      if (.not. ok) then
         result%message = 'invalid scatter determinant'
         return
      end if
      detroot = exp(logdet/real(d,dp))
      result%lambda = model%lambda
      result%mu = model%mu
      result%beta = matmul(inverse,model%gamma)
      result%delta = sqrt(max(0.0_dp,model%chi*detroot))
      result%delta_matrix = model%scatter/detroot
      psi_term = model%psi+dot_product(model%gamma,matmul(inverse,model%gamma))
      result%alpha = sqrt(max(0.0_dp,psi_term/detroot))
      result%ok = .true.
   end function ghyp_alpha_delta

   function qqghyp_data(data, model) result(result)
      real(dp), intent(in) :: data(:)
      type(ghyp_model_type), intent(in) :: model
      type(qq_result) :: result
      integer :: i, n
      if (model%dimension() /= 1 .or. size(data) < 1) then
         result%message = 'univariate data and model are required'
         return
      end if
      n = size(data)
      allocate(result%sample(n),result%theoretical(n))
      result%sample = data
      call insertion_sort(result%sample)
      do i = 1, n
         result%theoretical(i) = qghyp((real(i,dp)-0.5_dp)/real(n,dp),model)
      end do
      result%ok = .true.
   end function qqghyp_data

   subroutine insertion_sort(x)
      real(dp), intent(inout) :: x(:)
      real(dp) :: key
      integer :: i, j
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
   end subroutine insertion_sort

end module ghyp_utilities

! SPDX-License-Identifier: GPL-3.0-or-later
module cla_core
   use kind_mod, only: dp
   use cla_types, only: cla_result_t, cla_success, cla_invalid_input, &
      cla_infeasible_bounds, cla_singular_system, cla_no_improvement
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   public :: cla_solve
   public :: cla_mean_sigma

   interface
      subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         import dp
         integer, intent(in) :: n, nrhs, lda, ldb
         integer, intent(out) :: ipiv(*), info
         real(dp), intent(inout) :: a(lda,*), b(ldb,*)
      end subroutine dgesv

      subroutine dpotrf(uplo, n, a, lda, info)
         import dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, lda
         integer, intent(out) :: info
         real(dp), intent(inout) :: a(lda,*)
      end subroutine dpotrf
   end interface

contains

   function cla_solve(mu, covar, lower, upper, tol_lambda, check_covariance) result(out)
      !! Markowitz critical-line algorithm with individual box constraints.
      real(dp), intent(in) :: mu(:)
      real(dp), intent(in) :: covar(:,:)
      real(dp), intent(in) :: lower(:)
      real(dp), intent(in) :: upper(:)
      real(dp), intent(in), optional :: tol_lambda
      logical, intent(in), optional :: check_covariance
      type(cla_result_t) :: out

      real(dp), allocatable :: w(:), inv(:,:), inv_out(:,:), inv_candidates(:,:,:), candidate_inv(:,:)
      real(dp), allocatable :: work_weights(:,:), work_lambda(:), work_gamma(:)
      logical, allocatable :: work_mask(:,:), valid(:)
      integer, allocatable :: free(:), bounded(:), candidate_free(:)
      real(dp), allocatable :: lambda_out(:), boundary_values(:)
      real(dp) :: lambda, lambda_in, lambda_leave, gamma, tolerance
      real(dp) :: tol_work
      real(dp), allocatable :: free_weights(:)
      integer :: n, nf, nb, k, candidate, entering, leaving, info
      integer :: capacity, nturn, max_steps, step
      logical :: use_check, has_leave, remove_free

      n = size(mu)
      out%n_assets = n
      tolerance = 1.0e-7_dp
      if (present(tol_lambda)) tolerance = tol_lambda
      use_check = .true.
      if (present(check_covariance)) use_check = check_covariance

      if (n < 1 .or. size(covar,1) /= n .or. size(covar,2) /= n .or. &
          size(lower) /= n .or. size(upper) /= n .or. tolerance <= 0.0_dp .or. &
          .not. all(ieee_is_finite(mu)) .or. .not. all(ieee_is_finite(covar)) .or. &
          .not. all(ieee_is_finite(lower)) .or. .not. all(ieee_is_finite(upper)) .or. &
          any(lower > upper)) then
         out%info = cla_invalid_input
         return
      end if
      if (sum(lower) > 1.0_dp + 100.0_dp*epsilon(1.0_dp) .or. &
          sum(upper) < 1.0_dp - 100.0_dp*epsilon(1.0_dp)) then
         out%info = cla_infeasible_bounds
         return
      end if
      if (use_check) then
         call check_positive_definite(covar, info)
         if (info /= 0) out%warnings = out%warnings + 1
      end if

      call initialize_weights(mu, lower, upper, w, entering, info)
      if (info /= 0) then
         out%info = info
         return
      end if
      allocate(free(1))
      free(1) = entering

      capacity = max(8, 2*n + 4)
      allocate(work_weights(n,capacity), work_mask(n,capacity), &
         work_lambda(capacity), work_gamma(capacity))
      work_weights = 0.0_dp
      work_mask = .false.
      work_lambda = 0.0_dp
      work_gamma = 0.0_dp
      nturn = 0
      lambda = 1.0_dp
      max_steps = max(100, 20*n*n)

      if (n == 1) then
         call append_turning_point(w, free, 0.0_dp, 0.0_dp)
      else
         do step = 1, max_steps
            nf = size(free)
            if (nf < 1 .or. nf > n .or. lambda <= 0.0_dp) exit

            lambda_in = 0.0_dp
            entering = 0
            if (nf > 1) then
               if (.not. allocated(inv)) then
                  call compute_inverse_columns(mu, covar, w, free, inv, info)
                  if (info /= 0) then
                     out%info = cla_singular_system
                     exit
                  end if
               end if
               call free_to_bound_candidates(inv, w, free, lower, upper, &
                  lambda_in, entering, boundary_values)
            else
               allocate(boundary_values(0))
            end if

            bounded = complement_indices(n, free)
            nb = size(bounded)
            lambda_leave = -huge(1.0_dp)
            leaving = 0
            has_leave = .false.
            if (nb > 0) then
               allocate(lambda_out(nb), inv_candidates(nf+1,3,nb), valid(nb))
               lambda_out = -huge(1.0_dp)
               inv_candidates = 0.0_dp
               do candidate = 1, nb
                  allocate(candidate_free(nf+1))
                  candidate_free(1:nf) = free
                  candidate_free(nf+1) = bounded(candidate)
                  call compute_inverse_columns(mu, covar, w, candidate_free, &
                     candidate_inv, info)
                  if (info == 0) then
                     inv_candidates(:,:,candidate) = candidate_inv
                     lambda_out(candidate) = bound_to_free_lambda( &
                        inv_candidates(:,:,candidate), nf+1, w, bounded, candidate)
                  end if
                  deallocate(candidate_free)
                  if (allocated(candidate_inv)) deallocate(candidate_inv)
               end do

               valid = ieee_is_finite(lambda_out) .and. &
                  lambda_out > -0.5_dp*huge(1.0_dp)
               if (nturn > 0 .and. any(valid .and. &
                   lambda_out >= lambda*(1.0_dp - tolerance))) then
                  tol_work = tolerance
                  do
                     valid = ieee_is_finite(lambda_out) .and. &
                        lambda_out > -0.5_dp*huge(1.0_dp) .and. &
                        lambda_out < lambda*(1.0_dp - tol_work)
                     if (any(valid)) exit
                     tol_work = 0.5_dp*tol_work
                     if (tol_work <= 100.0_dp*epsilon(1.0_dp)) exit
                  end do
               end if
               if (any(valid)) then
                  k = maxloc(merge(lambda_out, -huge(1.0_dp), valid), dim=1)
                  lambda_leave = lambda_out(k)
                  leaving = bounded(k)
                  allocate(inv_out(nf+1,3))
                  inv_out = inv_candidates(:,:,k)
                  has_leave = .true.
               end if
               deallocate(lambda_out, inv_candidates, valid)
            end if

            remove_free = lambda_in > lambda_leave
            lambda = max(lambda_in, lambda_leave)
            if (.not. ieee_is_finite(lambda) .or. lambda <= 0.0_dp) lambda = 0.0_dp

            if (lambda > 0.0_dp) then
               if (remove_free) then
                  if (entering <= 0) then
                     out%info = cla_no_improvement
                     exit
                  end if
                  k = find_position(free, entering)
                  w(entering) = boundary_values(k)
                  call remove_integer(free, k)
                  if (allocated(inv)) deallocate(inv)
                  call compute_inverse_columns(mu, covar, w, free, inv, info)
                  if (info /= 0) then
                     out%info = cla_singular_system
                     exit
                  end if
               else if (has_leave) then
                  call append_integer(free, leaving)
                  if (allocated(inv)) deallocate(inv)
                  call move_alloc(inv_out, inv)
               else
                  out%warnings = out%warnings + 1
                  lambda = 0.0_dp
               end if
            end if

            if (.not. allocated(inv)) then
               call compute_inverse_columns(mu, covar, w, free, inv, info)
               if (info /= 0) then
                  out%info = cla_singular_system
                  exit
               end if
            end if
            bounded = complement_indices(n, free)
            call compute_free_weights(lambda, inv, w(bounded), free_weights, gamma)
            w(free) = free_weights
            where (abs(w - lower) <= 100.0_dp*epsilon(1.0_dp)) w = lower
            where (abs(w - upper) <= 100.0_dp*epsilon(1.0_dp)) w = upper
            call append_turning_point(w, free, lambda, gamma)

            if (allocated(boundary_values)) deallocate(boundary_values)
            if (allocated(inv_out)) deallocate(inv_out)
            if (lambda <= 0.0_dp) exit
         end do
      end if

      if (out%info == cla_success .and. nturn == 0) out%info = cla_no_improvement
      if (nturn > 0) then
         out%n_turning = nturn
         allocate(out%weights(n,nturn), out%free_mask(n,nturn), &
            out%lambdas(nturn), out%gammas(nturn), out%sigma(nturn), out%mu(nturn))
         out%weights = work_weights(:,1:nturn)
         out%free_mask = work_mask(:,1:nturn)
         out%lambdas = work_lambda(1:nturn)
         out%gammas = work_gamma(1:nturn)
         call cla_mean_sigma(out%weights, mu, covar, out%sigma, out%mu)
      end if

   contains

      subroutine append_turning_point(weights, indices, lam, gam)
         real(dp), intent(in) :: weights(:), lam, gam
         integer, intent(in) :: indices(:)
         if (nturn == capacity) call grow_storage()
         nturn = nturn + 1
         work_weights(:,nturn) = weights
         work_mask(:,nturn) = .false.
         work_mask(indices,nturn) = .true.
         work_lambda(nturn) = lam
         work_gamma(nturn) = gam
      end subroutine append_turning_point

      subroutine grow_storage()
         real(dp), allocatable :: w2(:,:), l2(:), g2(:)
         logical, allocatable :: m2(:,:)
         integer :: new_capacity
         new_capacity = 2*capacity
         allocate(w2(n,new_capacity), m2(n,new_capacity), l2(new_capacity), g2(new_capacity))
         w2 = 0.0_dp; m2 = .false.; l2 = 0.0_dp; g2 = 0.0_dp
         w2(:,1:capacity) = work_weights
         m2(:,1:capacity) = work_mask
         l2(1:capacity) = work_lambda
         g2(1:capacity) = work_gamma
         call move_alloc(w2, work_weights)
         call move_alloc(m2, work_mask)
         call move_alloc(l2, work_lambda)
         call move_alloc(g2, work_gamma)
         capacity = new_capacity
      end subroutine grow_storage

   end function cla_solve

   subroutine initialize_weights(mu, lower, upper, weights, free_index, info)
      real(dp), intent(in) :: mu(:), lower(:), upper(:)
      real(dp), allocatable, intent(out) :: weights(:)
      integer, intent(out) :: free_index, info
      integer, allocatable :: order(:)
      real(dp), allocatable :: sorted_weights(:)
      integer :: n, i, j, tmp

      n = size(mu)
      allocate(order(n), sorted_weights(n), weights(n))
      order = [(i, i=1,n)]
      do i = 1, n-1
         j = i - 1 + maxloc(mu(order(i:n)), dim=1)
         if (j /= i) then
            tmp = order(i); order(i) = order(j); order(j) = tmp
         end if
      end do
      sorted_weights = lower(order)
      i = 0
      do while (sum(sorted_weights) < 1.0_dp - 100.0_dp*epsilon(1.0_dp))
         i = i + 1
         if (i > n) then
            info = cla_infeasible_bounds
            return
         end if
         sorted_weights(i) = upper(order(i))
      end do
      if (i == 0) i = 1
      sorted_weights(i) = 1.0_dp - sum(sorted_weights, mask=[(j /= i, j=1,n)])
      weights(order) = sorted_weights
      free_index = order(i)
      info = cla_success
   end subroutine initialize_weights

   subroutine compute_inverse_columns(mu, covar, weights, free, inverse, info)
      real(dp), intent(in) :: mu(:), covar(:,:), weights(:)
      integer, intent(in) :: free(:)
      real(dp), allocatable, intent(out) :: inverse(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: a(:,:), rhs(:,:)
      integer, allocatable :: bounded(:), piv(:)
      integer :: nf, i, j, lapack_info

      nf = size(free)
      if (nf < 1) then
         info = 1
         allocate(inverse(0,0))
         return
      end if
      bounded = complement_indices(size(mu), free)
      allocate(a(nf,nf), rhs(nf,3), piv(nf))
      do j = 1, nf
         do i = 1, nf
            a(i,j) = covar(free(i),free(j))
         end do
      end do
      rhs(:,1) = 1.0_dp
      rhs(:,2) = mu(free)
      rhs(:,3) = 0.0_dp
      do j = 1, size(bounded)
         rhs(:,3) = rhs(:,3) + covar(free,bounded(j))*weights(bounded(j))
      end do
      call dgesv(nf, 3, a, nf, piv, rhs, nf, lapack_info)
      if (lapack_info /= 0 .or. .not. all(ieee_is_finite(rhs))) then
         info = lapack_info
         allocate(inverse(0,0))
      else
         allocate(inverse(nf,3))
         inverse = rhs
         info = 0
      end if
   end subroutine compute_inverse_columns

   subroutine compute_free_weights(lambda, inverse, bounded_weights, weights, gamma)
      real(dp), intent(in) :: lambda
      real(dp), intent(in) :: inverse(:,:)
      real(dp), intent(in) :: bounded_weights(:)
      real(dp), allocatable, intent(out) :: weights(:)
      real(dp), intent(out) :: gamma
      real(dp) :: sums(3)
      sums = sum(inverse, dim=1)
      gamma = (-lambda*sums(2) + (1.0_dp - sum(bounded_weights) + sums(3)))/sums(1)
      allocate(weights(size(inverse,1)))
      weights = -inverse(:,3) + gamma*inverse(:,1) + lambda*inverse(:,2)
   end subroutine compute_free_weights

   subroutine free_to_bound_candidates(inverse, weights, free, lower, upper, &
      best_lambda, best_asset, boundaries)
      real(dp), intent(in) :: inverse(:,:), weights(:), lower(:), upper(:)
      integer, intent(in) :: free(:)
      real(dp), intent(out) :: best_lambda
      integer, intent(out) :: best_asset
      real(dp), allocatable, intent(out) :: boundaries(:)
      integer, allocatable :: bounded(:)
      real(dp), allocatable :: candidates(:)
      real(dp) :: sums(3), c1, c4, ci, sum_bounded
      integer :: i, k

      bounded = complement_indices(size(weights), free)
      sum_bounded = sum(weights(bounded))
      sums = sum(inverse, dim=1)
      c1 = sums(1)
      allocate(boundaries(size(free)), candidates(size(free)))
      candidates = -huge(1.0_dp)
      do i = 1, size(free)
         c4 = inverse(i,1)
         ci = -c1*inverse(i,2) + sums(2)*c4
         if (ci > 0.0_dp) then
            boundaries(i) = upper(free(i))
         else if (ci < 0.0_dp) then
            boundaries(i) = lower(free(i))
         else
            boundaries(i) = 0.0_dp
         end if
         if (abs(ci) > 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(c1))) then
            candidates(i) = ((1.0_dp - sum_bounded + sums(3))*c4 - &
               c1*(boundaries(i) + inverse(i,3)))/ci
         end if
      end do
      k = maxloc(candidates, dim=1)
      best_lambda = candidates(k)
      best_asset = free(k)
   end subroutine free_to_bound_candidates

   real(dp) function bound_to_free_lambda(inverse, row, weights, bounded, position) result(lambda)
      real(dp), intent(in) :: inverse(:,:), weights(:)
      integer, intent(in) :: row, bounded(:), position
      integer, allocatable :: other(:)
      real(dp) :: sums(3), c1, c4, ci, sum_other
      integer :: j, count

      allocate(other(max(0,size(bounded)-1)))
      count = 0
      do j = 1, size(bounded)
         if (j /= position) then
            count = count + 1
            other(count) = bounded(j)
         end if
      end do
      sum_other = sum(weights(other))
      sums = sum(inverse, dim=1)
      c1 = sums(1)
      c4 = inverse(row,1)
      ci = -c1*inverse(row,2) + sums(2)*c4
      if (abs(ci) <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(c1))) then
         lambda = 0.0_dp
      else
         lambda = ((1.0_dp - sum_other + sums(3))*c4 - &
            c1*(weights(bounded(position)) + inverse(row,3)))/ci
      end if
   end function bound_to_free_lambda

   pure function complement_indices(n, selected) result(other)
      integer, intent(in) :: n, selected(:)
      integer, allocatable :: other(:)
      logical :: used(n)
      integer :: i, n_other
      used = .false.
      used(selected) = .true.
      allocate(other(count(.not. used)))
      n_other = 0
      do i = 1, n
         if (.not. used(i)) then
            n_other = n_other + 1
            other(n_other) = i
         end if
      end do
   end function complement_indices

   pure integer function find_position(values, target) result(position)
      integer, intent(in) :: values(:), target
      integer :: i
      position = 0
      do i = 1, size(values)
         if (values(i) == target) then
            position = i
            return
         end if
      end do
   end function find_position

   subroutine append_integer(values, value)
      integer, allocatable, intent(inout) :: values(:)
      integer, intent(in) :: value
      integer, allocatable :: tmp(:)
      allocate(tmp(size(values)+1))
      tmp(1:size(values)) = values
      tmp(size(tmp)) = value
      call move_alloc(tmp, values)
   end subroutine append_integer

   subroutine remove_integer(values, position)
      integer, allocatable, intent(inout) :: values(:)
      integer, intent(in) :: position
      integer, allocatable :: tmp(:)
      integer :: n
      n = size(values)
      allocate(tmp(n-1))
      if (position > 1) tmp(1:position-1) = values(1:position-1)
      if (position < n) tmp(position:n-1) = values(position+1:n)
      call move_alloc(tmp, values)
   end subroutine remove_integer

   subroutine check_positive_definite(matrix, info)
      real(dp), intent(in) :: matrix(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: a(:,:)
      integer :: n
      n = size(matrix,1)
      allocate(a(n,n))
      a = 0.5_dp*(matrix + transpose(matrix))
      call dpotrf('L', n, a, n, info)
   end subroutine check_positive_definite

   subroutine cla_mean_sigma(weights, asset_mu, covar, sigma, mean_return)
      real(dp), intent(in) :: weights(:,:), asset_mu(:), covar(:,:)
      real(dp), intent(out) :: sigma(:), mean_return(:)
      integer :: j
      do j = 1, size(weights,2)
         sigma(j) = sqrt(max(dot_product(weights(:,j), matmul(covar,weights(:,j))),0.0_dp))
         mean_return(j) = dot_product(weights(:,j), asset_mu)
      end do
   end subroutine cla_mean_sigma

end module cla_core

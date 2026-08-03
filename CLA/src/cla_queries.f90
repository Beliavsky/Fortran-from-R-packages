! SPDX-License-Identifier: GPL-3.0-or-later
module cla_queries
   use kind_mod, only: dp
   use cla_types, only: cla_result_t, cla_path_query_t, cla_success, &
      cla_invalid_input, cla_out_of_range
   implicit none
   private
   public :: cla_find_sigma, cla_find_mu

contains

   function cla_find_sigma(target_mu, result, covar, equal_tolerance) result(out)
      !! Interpolate portfolio standard deviation and weights at target means.
      real(dp), intent(in) :: target_mu(:)
      type(cla_result_t), intent(in) :: result
      real(dp), intent(in) :: covar(:,:)
      real(dp), intent(in), optional :: equal_tolerance
      type(cla_path_query_t) :: out
      integer, allocatable :: order(:)
      real(dp) :: tol, a
      integer :: n, m, j, i

      tol = 1.0e-6_dp
      if (present(equal_tolerance)) tol = equal_tolerance
      n = result%n_assets
      m = result%n_turning
      if (m < 1 .or. size(covar,1) /= n .or. size(covar,2) /= n .or. &
          .not. allocated(result%weights) .or. .not. allocated(result%mu)) then
         out%info = cla_invalid_input
         return
      end if
      order = sort_indices(result%mu)
      allocate(out%value(size(target_mu)), out%weights(n,size(target_mu)))
      do j = 1, size(target_mu)
         if (target_mu(j) < result%mu(order(1)) - tol .or. &
             target_mu(j) > result%mu(order(m)) + tol) then
            out%info = cla_out_of_range
            return
         end if
         i = lower_interval(target_mu(j), result%mu(order))
         if (i >= m .or. abs(result%mu(order(i)) - result%mu(order(i+1))) <= tol) then
            out%weights(:,j) = result%weights(:,order(i))
         else
            a = (target_mu(j) - result%mu(order(i+1)))/ &
               (result%mu(order(i)) - result%mu(order(i+1)))
            out%weights(:,j) = a*result%weights(:,order(i)) + &
               (1.0_dp-a)*result%weights(:,order(i+1))
         end if
         out%value(j) = sqrt(max(dot_product(out%weights(:,j), &
            matmul(covar,out%weights(:,j))),0.0_dp))
      end do
      out%info = cla_success
   end function cla_find_sigma

   function cla_find_mu(target_sigma, result, covar, tolerance, equal_tolerance) result(out)
      !! Invert the piecewise-hyperbolic frontier to obtain means and weights.
      real(dp), intent(in) :: target_sigma(:)
      type(cla_result_t), intent(in) :: result
      real(dp), intent(in) :: covar(:,:)
      real(dp), intent(in), optional :: tolerance, equal_tolerance
      type(cla_path_query_t) :: out
      integer, allocatable :: order(:)
      real(dp) :: root_tol, equal_tol, a, lo, hi, mid, fmid
      integer :: n, m, j, i, iteration

      root_tol = 1.0e-6_dp
      equal_tol = 1.0e-6_dp
      if (present(tolerance)) root_tol = tolerance
      if (present(equal_tolerance)) equal_tol = equal_tolerance
      n = result%n_assets
      m = result%n_turning
      if (m < 1 .or. size(covar,1) /= n .or. size(covar,2) /= n .or. &
          .not. allocated(result%weights) .or. .not. allocated(result%sigma)) then
         out%info = cla_invalid_input
         return
      end if
      order = sort_indices(result%sigma)
      allocate(out%value(size(target_sigma)), out%weights(n,size(target_sigma)))
      do j = 1, size(target_sigma)
         if (target_sigma(j) < result%sigma(order(1)) - equal_tol .or. &
             target_sigma(j) > result%sigma(order(m)) + equal_tol) then
            out%info = cla_out_of_range
            return
         end if
         i = lower_interval(target_sigma(j), result%sigma(order))
         if (i >= m .or. abs(result%sigma(order(i))-result%sigma(order(i+1))) <= equal_tol) then
            out%value(j) = result%mu(order(i))
            out%weights(:,j) = result%weights(:,order(i))
         else
            lo = min(result%mu(order(i)),result%mu(order(i+1)))
            hi = max(result%mu(order(i)),result%mu(order(i+1)))
            do iteration = 1, 100
               mid = 0.5_dp*(lo+hi)
               fmid = segment_sigma(mid, order(i), order(i+1)) - target_sigma(j)
               if (abs(fmid) <= root_tol .or. abs(hi-lo) <= root_tol) exit
               if ((segment_sigma(lo,order(i),order(i+1))-target_sigma(j))*fmid <= 0.0_dp) then
                  hi = mid
               else
                  lo = mid
               end if
            end do
            out%value(j) = mid
            a = (mid - result%mu(order(i+1)))/ &
               (result%mu(order(i)) - result%mu(order(i+1)))
            out%weights(:,j) = a*result%weights(:,order(i)) + &
               (1.0_dp-a)*result%weights(:,order(i+1))
         end if
      end do
      out%info = cla_success

   contains

      real(dp) function segment_sigma(target_mu, left_index, right_index) result(value)
         real(dp), intent(in) :: target_mu
         integer, intent(in) :: left_index, right_index
         real(dp) :: alpha
         real(dp) :: weights(n)
         if (abs(result%mu(left_index)-result%mu(right_index)) <= equal_tol) then
            weights = result%weights(:,left_index)
         else
            alpha = (target_mu-result%mu(right_index))/ &
               (result%mu(left_index)-result%mu(right_index))
            weights = alpha*result%weights(:,left_index) + &
               (1.0_dp-alpha)*result%weights(:,right_index)
         end if
         value = sqrt(max(dot_product(weights,matmul(covar,weights)),0.0_dp))
      end function segment_sigma

   end function cla_find_mu

   pure function sort_indices(values) result(order)
      real(dp), intent(in) :: values(:)
      integer, allocatable :: order(:)
      integer :: i, j, tmp
      allocate(order(size(values)))
      order = [(i,i=1,size(values))]
      do i = 2, size(values)
         tmp = order(i)
         j = i-1
         do while (j >= 1)
            if (values(order(j)) <= values(tmp)) exit
            order(j+1) = order(j)
            j = j-1
         end do
         order(j+1) = tmp
      end do
   end function sort_indices

   pure integer function lower_interval(value, sorted_values) result(index)
      real(dp), intent(in) :: value, sorted_values(:)
      integer :: lo, hi, mid
      if (value <= sorted_values(1)) then
         index = 1
         return
      end if
      if (value >= sorted_values(size(sorted_values))) then
         index = size(sorted_values)
         return
      end if
      lo = 1
      hi = size(sorted_values)
      do while (hi-lo > 1)
         mid = (lo+hi)/2
         if (sorted_values(mid) <= value) then
            lo = mid
         else
            hi = mid
         end if
      end do
      index = lo
   end function lower_interval

end module cla_queries

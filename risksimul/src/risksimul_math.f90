! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from riskSimul 0.1.2 by Wolfgang Hormann and Ismail Basoglu.
module risksimul_math
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use ghyp_kinds, only : dp
   use ghyp_special, only : normal_quantile, student_cdf
   use risksimul_types, only : allocation_result
   implicit none
   private

   public :: student_quantile, gamma_cdf, gamma_quantile
   public :: orthogonal_completion, integer_allocation
   public :: optimal_allocation_heuristic, inverse_table_quantile
   public :: sample_variance_from_sums, sample_covariance_from_sums

contains

   function student_quantile(p, nu) result(x)
      real(dp), intent(in) :: p, nu
      real(dp) :: x, lo, hi, mid, value
      integer :: iter

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
         return
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
         return
      else if (nu <= 0.0_dp) then
         x = ieee_value(1.0_dp,ieee_quiet_nan)
         return
      else if (abs(p-0.5_dp) <= epsilon(1.0_dp)) then
         x = 0.0_dp
         return
      end if

      x = normal_quantile(p)
      if (nu > 2.0_dp) x = x*sqrt(nu/max(nu-2.0_dp,tiny(1.0_dp)))
      lo = min(-1.0_dp, x-1.0_dp)
      hi = max( 1.0_dp, x+1.0_dp)
      do while (student_cdf(lo,nu) > p)
         hi = lo
         lo = 2.0_dp*lo
      end do
      do while (student_cdf(hi,nu) < p)
         lo = hi
         hi = 2.0_dp*hi
      end do
      do iter = 1, 100
         mid = 0.5_dp*(lo+hi)
         value = student_cdf(mid,nu)
         if (value < p) then
            lo = mid
         else
            hi = mid
         end if
         if (abs(hi-lo) <= 4.0e-13_dp*max(1.0_dp,abs(mid))) exit
      end do
      x = 0.5_dp*(lo+hi)
   end function student_quantile

   pure function regularized_gamma_p(a, x) result(value)
      real(dp), intent(in) :: a, x
      real(dp) :: value
      real(dp) :: sum_value, del, ap, b, c, d, h, an
      integer :: n
      real(dp), parameter :: eps = 4.0e-15_dp
      real(dp), parameter :: fpmin = 1.0e-300_dp

      if (a <= 0.0_dp .or. x < 0.0_dp) then
         value = ieee_value(1.0_dp,ieee_quiet_nan)
         return
      else if (x <= 0.0_dp) then
         value = 0.0_dp
         return
      end if

      if (x < a+1.0_dp) then
         ap = a
         sum_value = 1.0_dp/a
         del = sum_value
         do n = 1, 10000
            ap = ap+1.0_dp
            del = del*x/ap
            sum_value = sum_value+del
            if (abs(del) <= abs(sum_value)*eps) exit
         end do
         value = sum_value*exp(-x+a*log(x)-log_gamma(a))
      else
         b = x+1.0_dp-a
         c = 1.0_dp/fpmin
         d = 1.0_dp/max(abs(b),fpmin)
         if (b < 0.0_dp) d = -d
         h = d
         do n = 1, 10000
            an = -real(n,dp)*(real(n,dp)-a)
            b = b+2.0_dp
            d = an*d+b
            if (abs(d) < fpmin) d = sign(fpmin,d)
            c = b+an/c
            if (abs(c) < fpmin) c = sign(fpmin,c)
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del-1.0_dp) <= eps) exit
         end do
         value = 1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
      end if
      value = min(1.0_dp,max(0.0_dp,value))
   end function regularized_gamma_p

   pure function gamma_cdf(x, shape, scale) result(value)
      real(dp), intent(in) :: x, shape, scale
      real(dp) :: value
      if (scale <= 0.0_dp) then
         value = ieee_value(1.0_dp,ieee_quiet_nan)
      else
         value = regularized_gamma_p(shape,x/scale)
      end if
   end function gamma_cdf

   function gamma_quantile(p, shape, scale) result(x)
      real(dp), intent(in) :: p, shape, scale
      real(dp) :: x, lo, hi, mid, cdf, density, step
      integer :: iter

      if (p <= 0.0_dp) then
         x = 0.0_dp
         return
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
         return
      else if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         x = ieee_value(1.0_dp,ieee_quiet_nan)
         return
      end if

      lo = 0.0_dp
      hi = max(scale,shape*scale)
      do while (gamma_cdf(hi,shape,scale) < p)
         hi = 2.0_dp*hi
         if (hi >= huge(1.0_dp)/4.0_dp) exit
      end do
      x = min(hi,max(lo,shape*scale))
      do iter = 1, 100
         cdf = gamma_cdf(x,shape,scale)
         if (cdf < p) then
            lo = x
         else
            hi = x
         end if
         if (x > 0.0_dp) then
            density = exp((shape-1.0_dp)*log(x)-x/scale- &
               log_gamma(shape)-shape*log(scale))
         else
            density = 0.0_dp
         end if
         if (density > tiny(1.0_dp)) then
            step = (cdf-p)/density
            mid = x-step
         else
            mid = 0.5_dp*(lo+hi)
         end if
         if (mid <= lo .or. mid >= hi .or. .not. ieee_is_finite(mid)) &
            mid = 0.5_dp*(lo+hi)
         x = mid
         if (abs(hi-lo) <= 8.0e-13_dp*max(1.0_dp,x)) exit
      end do
   end function gamma_quantile

   subroutine orthogonal_completion(direction, basis, ok)
      real(dp), intent(in) :: direction(:)
      real(dp), allocatable, intent(out) :: basis(:,:)
      logical, intent(out) :: ok
      real(dp), allocatable :: candidate(:)
      real(dp) :: norm_value
      integer :: d, j, k, selected

      d = size(direction)
      allocate(basis(d,d),candidate(d))
      basis = 0.0_dp
      norm_value = sqrt(dot_product(direction,direction))
      if (d < 1 .or. norm_value <= tiny(1.0_dp)) then
         ok = .false.
         return
      end if
      basis(:,1) = direction/norm_value
      do j = 2, d
         selected = 0
         do k = 1, d
            candidate = 0.0_dp
            candidate(k) = 1.0_dp
            candidate = candidate-matmul(basis(:,1:j-1), &
               matmul(transpose(basis(:,1:j-1)),candidate))
            norm_value = sqrt(dot_product(candidate,candidate))
            if (norm_value > 1.0e-10_dp) then
               selected = k
               exit
            end if
         end do
         if (selected == 0) then
            ok = .false.
            return
         end if
         basis(:,j) = candidate/norm_value
      end do
      ok = .true.
   end subroutine orthogonal_completion

   subroutine integer_allocation(fractions, total, minimum_count, counts)
      real(dp), intent(in) :: fractions(:)
      integer, intent(in) :: total, minimum_count
      integer, allocatable, intent(out) :: counts(:)
      real(dp), allocatable :: f(:), cumulative(:)
      real(dp) :: denom
      integer :: i, previous

      allocate(counts(size(fractions)),f(size(fractions)),cumulative(size(fractions)))
      f = max(fractions,0.0_dp)
      denom = sum(f)
      if (denom <= tiny(1.0_dp)) then
         f = 1.0_dp/real(size(f),dp)
      else
         f = f/denom
      end if
      cumulative = real(total,dp)*[ (sum(f(1:i)), i=1,size(f)) ]
      previous = 0
      do i = 1, size(f)
         counts(i) = floor(cumulative(i))-previous
         previous = floor(cumulative(i))
         counts(i) = max(counts(i),minimum_count)
      end do
   end subroutine integer_allocation

   function optimal_allocation_heuristic(a, tolerance, upstream_compatibility) result(result)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tolerance
      logical, intent(in), optional :: upstream_compatibility
      type(allocation_result) :: result
      real(dp), allocatable :: extreme(:,:), current(:), best(:), objectives(:)
      real(dp) :: eps, delta, eta, best_obj, current_obj, candidate_obj
      real(dp) :: denom
      integer :: j, current_index, candidate_index, iterations
      logical :: compat

      eps = 1.0e-6_dp
      if (present(tolerance)) eps = tolerance
      compat = .false.
      if (present(upstream_compatibility)) compat = upstream_compatibility
      if (size(a,1) < 1 .or. size(a,2) < 1 .or. any(a < 0.0_dp)) then
         result%message = 'allocation matrix must be nonnegative and nonempty'
         return
      end if

      allocate(extreme(size(a,1),size(a,2)),current(size(a,1)),best(size(a,1)))
      allocate(objectives(size(a,2)))
      do j = 1, size(a,2)
         extreme(:,j) = sqrt(a(:,j))
         denom = sum(extreme(:,j))
         if (denom <= tiny(1.0_dp)) then
            extreme(:,j) = 1.0_dp/real(size(a,1),dp)
         else
            extreme(:,j) = extreme(:,j)/denom
         end if
      end do
      current = sum(extreme,dim=2)/real(size(a,2),dp)
      best = current
      objectives = [ (sum(a(:,j)/(current+1.0e-16_dp)), j=1,size(a,2)) ]
      current_index = maxloc(objectives,dim=1)
      current_obj = objectives(current_index)
      best_obj = current_obj
      delta = huge(1.0_dp)
      eta = 1.0_dp
      iterations = 0

      do while (abs(delta)/max(best_obj,tiny(1.0_dp)) > eps .and. iterations < 100000)
         iterations = iterations+1
         current = (eta*current+extreme(:,current_index))/(eta+1.0_dp)
         objectives = [ (sum(a(:,j)/(current+1.0e-16_dp)), j=1,size(a,2)) ]
         candidate_index = maxloc(objectives,dim=1)
         candidate_obj = objectives(candidate_index)
         if (candidate_obj <= best_obj) then
            best = current
            if (compat) then
               best_obj = current_obj
            else
               best_obj = candidate_obj
            end if
         end if
         if (candidate_index /= current_index) then
            eta = eta+1.0_dp
            delta = current_obj-candidate_obj
            current_obj = candidate_obj
            current_index = candidate_index
         else if (iterations > 10) then
            delta = 0.0_dp
         end if
      end do
      result%fractions = best/sum(best)
      result%objectives = [ (sum(a(:,j)/(result%fractions+1.0e-16_dp)), &
         j=1,size(a,2)) ]
      result%worst_objective = maxval(result%objectives)
      result%ok = .true.
   end function optimal_allocation_heuristic

   function inverse_table_quantile(p, x, cdf) result(value)
      real(dp), intent(in) :: p
      real(dp), intent(in) :: x(:), cdf(:)
      real(dp) :: value, fraction
      integer :: lo, hi, mid

      if (size(x) /= size(cdf) .or. size(x) < 2) then
         value = ieee_value(1.0_dp,ieee_quiet_nan)
         return
      end if
      if (p <= cdf(1)) then
         value = x(1)
         return
      else if (p >= cdf(size(cdf))) then
         value = x(size(x))
         return
      end if
      lo = 1
      hi = size(cdf)
      do while (hi-lo > 1)
         mid = (lo+hi)/2
         if (cdf(mid) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      if (cdf(hi) <= cdf(lo)+tiny(1.0_dp)) then
         value = 0.5_dp*(x(lo)+x(hi))
      else
         fraction = (p-cdf(lo))/(cdf(hi)-cdf(lo))
         value = x(lo)+fraction*(x(hi)-x(lo))
      end if
   end function inverse_table_quantile

   pure function sample_variance_from_sums(sum_x, sum_x2, n) result(value)
      real(dp), intent(in) :: sum_x, sum_x2
      integer, intent(in) :: n
      real(dp) :: value
      if (n <= 1) then
         value = 0.0_dp
      else
         value = max(0.0_dp,(sum_x2-sum_x*sum_x/real(n,dp))/real(n-1,dp))
      end if
   end function sample_variance_from_sums

   pure function sample_covariance_from_sums(sum_x, sum_y, sum_xy, n) result(value)
      real(dp), intent(in) :: sum_x, sum_y, sum_xy
      integer, intent(in) :: n
      real(dp) :: value
      if (n <= 1) then
         value = 0.0_dp
      else
         value = (sum_xy-sum_x*sum_y/real(n,dp))/real(n-1,dp)
      end if
   end function sample_covariance_from_sums

end module risksimul_math

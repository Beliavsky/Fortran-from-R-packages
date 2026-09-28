module psych_stats
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   implicit none
   private
   public :: finite_count, finite_mean, finite_variance, finite_sd, finite_covariance, finite_correlation
   public :: finite_quantile, finite_median, rank_average, covariance_matrix, correlation_matrix
   public :: center_columns, scale_columns, trace_matrix, identity_matrix

contains

   pure integer function finite_count(x) result(n)
      real(dp), intent(in) :: x(:) !! Numeric vector whose non-NaN elements are counted.
      integer :: i
      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) n = n + 1
      end do
   end function finite_count

   pure real(dp) function finite_mean(x) result(m)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN entries are ignored.
      real(dp) :: s
      integer :: i, n
      s = 0.0_dp
      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) then
            s = s + x(i)
            n = n + 1
         end if
      end do
      if (n > 0) then
         m = s/real(n,dp)
      else
         m = 0.0_dp
      end if
   end function finite_mean

   pure real(dp) function finite_variance(x) result(v)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN entries are ignored and denominator is `n-1`.
      real(dp) :: m, s
      integer :: i, n
      m = finite_mean(x)
      s = 0.0_dp
      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) then
            s = s + (x(i)-m)**2
            n = n + 1
         end if
      end do
      if (n > 1) then
         v = s/real(n-1,dp)
      else
         v = 0.0_dp
      end if
   end function finite_variance

   pure real(dp) function finite_sd(x) result(s)
      real(dp), intent(in) :: x(:) !! Numeric vector; NaN entries are ignored.
      s = sqrt(max(0.0_dp, finite_variance(x)))
   end function finite_sd

   pure real(dp) function finite_covariance(x, y) result(c)
      real(dp), intent(in) :: x(:) !! First vector; paired NaNs are omitted.
      real(dp), intent(in) :: y(:) !! Second vector with the same length as `x`.
      real(dp) :: mx, my, sx, sy, sxy
      integer :: i, n
      sx = 0.0_dp
      sy = 0.0_dp
      n = 0
      do i = 1, min(size(x),size(y))
         if (.not. ieee_is_nan(x(i)) .and. .not. ieee_is_nan(y(i))) then
            sx = sx + x(i)
            sy = sy + y(i)
            n = n + 1
         end if
      end do
      if (n <= 1) then
         c = 0.0_dp
         return
      end if
      mx = sx/real(n,dp)
      my = sy/real(n,dp)
      sxy = 0.0_dp
      do i = 1, min(size(x),size(y))
         if (.not. ieee_is_nan(x(i)) .and. .not. ieee_is_nan(y(i))) then
            sxy = sxy + (x(i)-mx)*(y(i)-my)
         end if
      end do
      c = sxy/real(n-1,dp)
   end function finite_covariance

   pure real(dp) function finite_correlation(x, y) result(r)
      real(dp), intent(in) :: x(:) !! First vector; paired NaNs are omitted.
      real(dp), intent(in) :: y(:) !! Second vector with the same length as `x`.
      real(dp) :: mx, my, sx, sy, sxx, syy, sxy
      integer :: i, n
      sx = 0.0_dp
      sy = 0.0_dp
      n = 0
      do i = 1, min(size(x),size(y))
         if (.not. ieee_is_nan(x(i)) .and. .not. ieee_is_nan(y(i))) then
            sx = sx + x(i)
            sy = sy + y(i)
            n = n + 1
         end if
      end do
      if (n <= 1) then
         r = 0.0_dp
         return
      end if
      mx = sx/real(n,dp)
      my = sy/real(n,dp)
      sxx = 0.0_dp
      syy = 0.0_dp
      sxy = 0.0_dp
      do i = 1, min(size(x),size(y))
         if (.not. ieee_is_nan(x(i)) .and. .not. ieee_is_nan(y(i))) then
            sxx = sxx + (x(i)-mx)**2
            syy = syy + (y(i)-my)**2
            sxy = sxy + (x(i)-mx)*(y(i)-my)
         end if
      end do
      if (sxx <= 0.0_dp .or. syy <= 0.0_dp) then
         r = 0.0_dp
      else
         r = sxy/sqrt(sxx*syy)
      end if
   end function finite_correlation

   pure subroutine sort_real(x)
      real(dp), intent(inout) :: x(:) !! Vector sorted into ascending order in place.
      real(dp) :: key
      integer :: i, j
      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j+1) = x(j)
            j = j - 1
         end do
         x(j+1) = key
      end do
   end subroutine sort_real

   pure real(dp) function finite_quantile(x, probability) result(q)
      real(dp), intent(in) :: x(:) !! Numeric sample; NaN entries are omitted.
      real(dp), intent(in) :: probability !! Quantile probability in `[0,1]`; R type-7 interpolation is used.
      real(dp), allocatable :: work(:)
      real(dp) :: h, fraction
      integer :: i, j, n
      n = finite_count(x)
      if (n == 0) then
         q = 0.0_dp
         return
      end if
      allocate(work(n))
      j = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) then
            j = j + 1
            work(j) = x(i)
         end if
      end do
      call sort_real(work)
      if (probability <= 0.0_dp) then
         q = work(1)
      else if (probability >= 1.0_dp) then
         q = work(n)
      else
         h = 1.0_dp + real(n-1,dp)*probability
         i = floor(h)
         fraction = h-real(i,dp)
         if (i >= n) then
            q = work(n)
         else
            q = work(i) + fraction*(work(i+1)-work(i))
         end if
      end if
   end function finite_quantile

   pure real(dp) function finite_median(x) result(m)
      real(dp), intent(in) :: x(:) !! Numeric sample; NaN entries are omitted.
      m = finite_quantile(x,0.5_dp)
   end function finite_median

   pure subroutine rank_average(x, ranks)
      real(dp), intent(in) :: x(:) !! Values to rank in ascending order; NaNs receive NaN-like zero rank.
      real(dp), intent(out) :: ranks(:) !! Average ranks for ties, with length equal to `x`.
      integer, allocatable :: order(:)
      integer :: i, j, k, n, tmp
      real(dp) :: av
      n = size(x)
      ranks = 0.0_dp
      allocate(order(n))
      do i = 1, n
         order(i) = i
      end do
      do i = 2, n
         tmp = order(i)
         j = i - 1
         do while (j >= 1)
            if (x(order(j)) <= x(tmp)) exit
            order(j+1) = order(j)
            j = j - 1
         end do
         order(j+1) = tmp
      end do
      i = 1
      do while (i <= n)
         if (ieee_is_nan(x(order(i)))) then
            ranks(order(i)) = 0.0_dp
            i = i + 1
         else
            j = i
            do while (j < n)
               if (abs(x(order(j+1)) - x(order(i))) > epsilon(1.0_dp) * &
                  max(1.0_dp, abs(x(order(j+1))), abs(x(order(i))))) exit
               j = j + 1
            end do
            av = 0.5_dp*real(i+j,dp)
            do k = i, j
               ranks(order(k)) = av
            end do
            i = j + 1
         end if
      end do
   end subroutine rank_average

   pure subroutine covariance_matrix(x, cov)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix; pairwise NaNs are omitted.
      real(dp), intent(out) :: cov(:, :) !! Pairwise covariance matrix with shape `(p,p)`.
      integer :: i, j
      do j = 1, size(x,2)
         do i = 1, size(x,2)
            cov(i,j) = finite_covariance(x(:,i),x(:,j))
         end do
      end do
   end subroutine covariance_matrix

   pure subroutine correlation_matrix(x, cor)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix; pairwise NaNs are omitted.
      real(dp), intent(out) :: cor(:, :) !! Pairwise Pearson-correlation matrix with shape `(p,p)`.
      integer :: i, j
      do j = 1, size(x,2)
         do i = 1, size(x,2)
            cor(i,j) = finite_correlation(x(:,i),x(:,j))
         end do
      end do
   end subroutine correlation_matrix

   pure subroutine center_columns(x, centered, means)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix.
      real(dp), intent(out) :: centered(:, :) !! Column-centered matrix, preserving NaN entries.
      real(dp), intent(out) :: means(:) !! Column means after omitting NaNs.
      integer :: i, j
      do j = 1, size(x,2)
         means(j) = finite_mean(x(:,j))
         do i = 1, size(x,1)
            if (ieee_is_nan(x(i,j))) then
               centered(i,j) = x(i,j)
            else
               centered(i,j) = x(i,j)-means(j)
            end if
         end do
      end do
   end subroutine center_columns

   pure subroutine scale_columns(x, scaled, means, sds)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix.
      real(dp), intent(out) :: scaled(:, :) !! Standardized matrix, preserving NaN entries.
      real(dp), intent(out) :: means(:) !! Column means after omitting NaNs.
      real(dp), intent(out) :: sds(:) !! Column sample standard deviations after omitting NaNs.
      integer :: i, j
      do j = 1, size(x,2)
         means(j) = finite_mean(x(:,j))
         sds(j) = finite_sd(x(:,j))
         do i = 1, size(x,1)
            if (ieee_is_nan(x(i,j))) then
               scaled(i,j) = x(i,j)
            else if (sds(j) > 0.0_dp) then
               scaled(i,j) = (x(i,j)-means(j))/sds(j)
            else
               scaled(i,j) = 0.0_dp
            end if
         end do
      end do
   end subroutine scale_columns

   pure real(dp) function trace_matrix(a) result(value)
      real(dp), intent(in) :: a(:, :) !! Square matrix whose trace is required.
      integer :: i
      value = 0.0_dp
      do i = 1, min(size(a,1),size(a,2))
         value = value + a(i,i)
      end do
   end function trace_matrix

   pure subroutine identity_matrix(a)
      real(dp), intent(out) :: a(:, :) !! Square output matrix set to the identity.
      integer :: i
      a = 0.0_dp
      do i = 1, min(size(a,1),size(a,2))
         a(i,i) = 1.0_dp
      end do
   end subroutine identity_matrix

end module psych_stats

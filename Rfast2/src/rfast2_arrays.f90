module rfast2_arrays
   use rfast_special, only : dp
   use rfast_arrays, only : sort_real, order_real, median_r, colmeans, rowsums, colsums, colvars
   implicit none
   private

   public :: quantile_rfast2, row_quantile, col_quantile
   public :: trim_mean, row_trim_mean, col_trim_mean
   public :: intersect_real, merge_sorted
   public :: is_lower_tri, is_upper_tri, is_skew_symmetric, lud_parts
   public :: col_group_sum, col_group_mean, col_group_min, col_group_max, col_group_median, col_group_var
   public :: jack_mean, coljack_means, rowjack_means, colmeansvars
   public :: col_mse, col_mae, col_pkl, col_ukl
   public :: col_accuracy, col_sensitivity, col_specificity, col_precision, col_fscore, col_fmi, col_fbscore

contains

   real(dp) function quantile_rfast2(x, prob) result(q)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: prob
      real(dp), allocatable :: y(:)
      real(dp) :: h, frac
      integer :: n, lo, hi

      n = size(x)
      if (n == 0) then
         q = huge(1.0_dp)
         return
      end if
      if (prob <= 0.0_dp) then
         q = minval(x)
         return
      end if
      if (prob >= 1.0_dp) then
         q = maxval(x)
         return
      end if
      allocate(y(n))
      y = x
      call sort_real(y)
      h = 1.0_dp + real(n - 1, dp) * prob
      lo = floor(h)
      hi = ceiling(h)
      lo = max(1, min(n, lo))
      hi = max(1, min(n, hi))
      frac = h - real(lo, dp)
      q = (1.0_dp - frac) * y(lo) + frac * y(hi)
   end function quantile_rfast2

   function col_quantile(x, probs) result(out)
      real(dp), intent(in) :: x(:,:), probs(:)
      real(dp) :: out(size(probs), size(x,2))
      integer :: i, j

      do j = 1, size(x,2)
         do i = 1, size(probs)
            out(i,j) = quantile_rfast2(x(:,j), probs(i))
         end do
      end do
   end function col_quantile

   function row_quantile(x, probs) result(out)
      real(dp), intent(in) :: x(:,:), probs(:)
      real(dp) :: out(size(x,1), size(probs))
      integer :: i, j

      do i = 1, size(x,1)
         do j = 1, size(probs)
            out(i,j) = quantile_rfast2(x(i,:), probs(j))
         end do
      end do
   end function row_quantile

   real(dp) function trim_mean(x, a) result(m)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: a
      real(dp), allocatable :: y(:)
      real(dp) :: trim
      integer :: n, k, lo, hi

      n = size(x)
      if (n == 0) then
         m = huge(1.0_dp)
         return
      end if
      trim = 0.05_dp
      if (present(a)) trim = a
      trim = max(0.0_dp, min(0.499999999_dp, trim))
      allocate(y(n))
      y = x
      call sort_real(y)
      k = floor(trim * real(n,dp))
      lo = k + 1
      hi = n - k
      m = sum(y(lo:hi)) / real(hi - lo + 1, dp)
   end function trim_mean

   function col_trim_mean(x, a) result(out)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: a
      real(dp) :: out(size(x,2))
      integer :: j

      do j = 1, size(x,2)
         if (present(a)) then
            out(j) = trim_mean(x(:,j), a)
         else
            out(j) = trim_mean(x(:,j))
         end if
      end do
   end function col_trim_mean

   function row_trim_mean(x, a) result(out)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: a
      real(dp) :: out(size(x,1))
      integer :: i

      do i = 1, size(x,1)
         if (present(a)) then
            out(i) = trim_mean(x(i,:), a)
         else
            out(i) = trim_mean(x(i,:))
         end if
      end do
   end function row_trim_mean

   function intersect_real(x, y) result(z)
      real(dp), intent(in) :: x(:), y(:)
      real(dp), allocatable :: z(:), tmp(:)
      integer :: i, n

      allocate(tmp(min(size(x), size(y))))
      n = 0
      do i = 1, size(x)
         if (any(real_equal(y,x(i)))) then
            if (n == 0 .or. .not. any(real_equal(tmp(1:n),x(i)))) then
               n = n + 1
               tmp(n) = x(i)
            end if
         end if
      end do
      allocate(z(n))
      if (n > 0) z = tmp(1:n)
   end function intersect_real

   function merge_sorted(x, y) result(z)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: z(size(x) + size(y))
      integer :: i, j, k

      i = 1
      j = 1
      k = 1
      do while (i <= size(x) .and. j <= size(y))
         if (x(i) <= y(j)) then
            z(k) = x(i)
            i = i + 1
         else
            z(k) = y(j)
            j = j + 1
         end if
         k = k + 1
      end do
      if (i <= size(x)) z(k:) = x(i:)
      if (j <= size(y)) z(k:) = y(j:)
   end function merge_sorted

   logical function is_lower_tri(x, diag) result(ok)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: diag
      logical :: d
      integer :: i, j

      d = .false.
      if (present(diag)) d = diag
      ok = .true.
      do j = 1, size(x,2)
         do i = 1, min(size(x,1), j - merge(0,1,d))
            if (d) then
               if (i <= j .and. abs(x(i,j)) > 0.0_dp) then
                  ok = .false.
                  return
               end if
            else
               if (i < j .and. abs(x(i,j)) > 0.0_dp) then
                  ok = .false.
                  return
               end if
            end if
         end do
      end do
   end function is_lower_tri

   logical function is_upper_tri(x, diag) result(ok)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: diag
      logical :: d
      integer :: i, j

      d = .false.
      if (present(diag)) d = diag
      ok = .true.
      do j = 1, size(x,2)
         do i = 1, size(x,1)
            if (d) then
               if (i >= j .and. abs(x(i,j)) > 0.0_dp) then
                  ok = .false.
                  return
               end if
            else
               if (i > j .and. abs(x(i,j)) > 0.0_dp) then
                  ok = .false.
                  return
               end if
            end if
         end do
      end do
   end function is_upper_tri

   logical function is_skew_symmetric(x, tol) result(ok)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: tol
      real(dp) :: eps

      eps = 0.0_dp
      if (present(tol)) eps = tol
      if (size(x,1) /= size(x,2)) then
         ok = .false.
      else
         ok = maxval(abs(x + transpose(x))) <= eps
      end if
   end function is_skew_symmetric

   subroutine lud_parts(x, lower, upper, diagv)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: lower(:), upper(:), diagv(:)
      integer :: nl, nu, nd, i, j, k

      nd = min(size(x,1), size(x,2))
      nl = 0
      nu = 0
      do j = 1, size(x,2)
         do i = j + 1, size(x,1)
            nl = nl + 1
         end do
         do i = 1, min(j - 1, size(x,1))
            nu = nu + 1
         end do
      end do
      allocate(lower(nl), upper(nu), diagv(nd))
      k = 0
      do j = 1, size(x,2)
         do i = j + 1, size(x,1)
            k = k + 1
            lower(k) = x(i,j)
         end do
      end do
      k = 0
      do j = 1, size(x,2)
         do i = 1, min(j - 1, size(x,1))
            k = k + 1
            upper(k) = x(i,j)
         end do
      end do
      do i = 1, nd
         diagv(i) = x(i,i)
      end do
   end subroutine lud_parts

   function col_group_sum(x, group, ngroup) result(out)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: group(:), ngroup
      real(dp) :: out(ngroup, size(x,2))
      integer :: i

      out = 0.0_dp
      do i = 1, size(x,1)
         if (group(i) >= 1 .and. group(i) <= ngroup) out(group(i),:) = out(group(i),:) + x(i,:)
      end do
   end function col_group_sum

   function col_group_mean(x, group, ngroup) result(out)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: group(:), ngroup
      real(dp) :: out(ngroup, size(x,2))
      integer :: cnt(ngroup), i, g

      out = col_group_sum(x, group, ngroup)
      cnt = 0
      do i = 1, size(group)
         g = group(i)
         if (g >= 1 .and. g <= ngroup) cnt(g) = cnt(g) + 1
      end do
      do g = 1, ngroup
         if (cnt(g) > 0) out(g,:) = out(g,:) / real(cnt(g),dp)
      end do
   end function col_group_mean

   function col_group_min(x, group, ngroup) result(out)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: group(:), ngroup
      real(dp) :: out(ngroup, size(x,2))
      logical :: seen(ngroup)
      integer :: i, g

      out = huge(1.0_dp)
      seen = .false.
      do i = 1, size(x,1)
         g = group(i)
         if (g >= 1 .and. g <= ngroup) then
            out(g,:) = min(out(g,:), x(i,:))
            seen(g) = .true.
         end if
      end do
      do g = 1, ngroup
         if (.not. seen(g)) out(g,:) = 0.0_dp
      end do
   end function col_group_min

   function col_group_max(x, group, ngroup) result(out)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: group(:), ngroup
      real(dp) :: out(ngroup, size(x,2))
      logical :: seen(ngroup)
      integer :: i, g

      out = -huge(1.0_dp)
      seen = .false.
      do i = 1, size(x,1)
         g = group(i)
         if (g >= 1 .and. g <= ngroup) then
            out(g,:) = max(out(g,:), x(i,:))
            seen(g) = .true.
         end if
      end do
      do g = 1, ngroup
         if (.not. seen(g)) out(g,:) = 0.0_dp
      end do
   end function col_group_max

   function col_group_median(x, group, ngroup) result(out)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: group(:), ngroup
      real(dp) :: out(ngroup, size(x,2))
      integer :: g, j

      out = 0.0_dp
      do g = 1, ngroup
         do j = 1, size(x,2)
            if (count(group == g) > 0) out(g,j) = median_r(pack(x(:,j), group == g))
         end do
      end do
   end function col_group_median

   function col_group_var(x, group, ngroup, std) result(out)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: group(:), ngroup
      logical, intent(in), optional :: std
      real(dp) :: out(ngroup, size(x,2))
      real(dp) :: m
      integer :: g, j, n
      logical :: sd

      sd = .false.
      if (present(std)) sd = std
      out = 0.0_dp
      do g = 1, ngroup
         n = count(group == g)
         if (n <= 1) cycle
         do j = 1, size(x,2)
            m = sum(pack(x(:,j), group == g)) / real(n,dp)
            out(g,j) = sum((pack(x(:,j), group == g) - m)**2) / real(n - 1,dp)
            if (sd) out(g,j) = sqrt(out(g,j))
         end do
      end do
   end function col_group_var

   pure real(dp) function jack_mean(x) result(v)
      real(dp), intent(in) :: x(:)
      real(dp) :: n

      n = real(size(x),dp)
      v = sum(x) * (n - 1.0_dp) / (n*n)
   end function jack_mean

   pure function coljack_means(x) result(v)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: v(size(x,2)), n

      n = real(size(x,1),dp)
      v = sum(x,dim=1) * (n - 1.0_dp) / (n*n)
   end function coljack_means

   pure function rowjack_means(x) result(v)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: v(size(x,1)), n

      n = real(size(x,2),dp)
      v = sum(x,dim=2) * (n - 1.0_dp) / (n*n)
   end function rowjack_means

   function colmeansvars(x, std) result(out)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: std
      real(dp) :: out(2,size(x,2))
      logical :: sd

      sd = .false.
      if (present(std)) sd = std
      out(1,:) = colmeans(x)
      out(2,:) = colvars(x)
      if (sd) out(2,:) = sqrt(out(2,:))
   end function colmeansvars

   pure function col_mse(y, yhat) result(v)
      real(dp), intent(in) :: y(:), yhat(:,:)
      real(dp) :: v(size(yhat,2))
      integer :: j

      do j = 1, size(yhat,2)
         v(j) = sum((y - yhat(:,j))**2) / real(size(y),dp)
      end do
   end function col_mse

   pure function col_mae(y, yhat) result(v)
      real(dp), intent(in) :: y(:), yhat(:,:)
      real(dp) :: v(size(yhat,2))
      integer :: j

      do j = 1, size(yhat,2)
         v(j) = sum(abs(y - yhat(:,j))) / real(size(y),dp)
      end do
   end function col_mae

   function col_pkl(y, yhat) result(v)
      real(dp), intent(in) :: y(:), yhat(:,:)
      real(dp) :: v(size(yhat,2)), yh
      integer :: i, j

      v = 0.0_dp
      do j = 1, size(yhat,2)
         do i = 1, size(y)
            yh = min(1.0_dp - epsilon(1.0_dp), max(tiny(1.0_dp), yhat(i,j)))
            if (y(i) > 0.0_dp) v(j) = v(j) + y(i) * log(y(i)/yh)
            if (y(i) < 1.0_dp) v(j) = v(j) + (1.0_dp-y(i))*log((1.0_dp-y(i))/(1.0_dp-yh))
         end do
      end do
   end function col_pkl

   function col_ukl(y, yhat) result(v)
      real(dp), intent(in) :: y(:), yhat(:,:)
      real(dp) :: v(size(yhat,2)), yh
      integer :: i, j

      v = 0.0_dp
      do j = 1, size(yhat,2)
         do i = 1, size(y)
            yh = max(tiny(1.0_dp), yhat(i,j))
            if (y(i) > 0.0_dp) v(j) = v(j) + y(i) * log(y(i)/yh)
         end do
      end do
   end function col_ukl

   subroutine confusion_counts(group, preds, tp, fp, tn, fn)
      integer, intent(in) :: group(:), preds(:,:)
      integer, intent(out) :: tp(size(preds,2)), fp(size(preds,2)), tn(size(preds,2)), fn(size(preds,2))
      integer :: j

      do j = 1, size(preds,2)
         tp(j) = count(group == 1 .and. preds(:,j) == 1)
         fp(j) = count(group == 0 .and. preds(:,j) == 1)
         tn(j) = count(group == 0 .and. preds(:,j) == 0)
         fn(j) = count(group == 1 .and. preds(:,j) == 0)
      end do
   end subroutine confusion_counts

   function col_accuracy(group, preds) result(v)
      integer, intent(in) :: group(:), preds(:,:)
      real(dp) :: v(size(preds,2))
      integer :: j

      do j = 1, size(preds,2)
         v(j) = real(count(group == preds(:,j)),dp) / real(size(group),dp)
      end do
   end function col_accuracy

   function col_sensitivity(group, preds) result(v)
      integer, intent(in) :: group(:), preds(:,:)
      real(dp) :: v(size(preds,2))
      integer :: tp(size(preds,2)), fp(size(preds,2)), tn(size(preds,2)), fn(size(preds,2)), j

      call confusion_counts(group,preds,tp,fp,tn,fn)
      do j = 1, size(v)
         v(j) = real(tp(j),dp) / real(max(1,tp(j)+fn(j)),dp)
      end do
   end function col_sensitivity

   function col_specificity(group, preds) result(v)
      integer, intent(in) :: group(:), preds(:,:)
      real(dp) :: v(size(preds,2))
      integer :: tp(size(preds,2)), fp(size(preds,2)), tn(size(preds,2)), fn(size(preds,2)), j

      call confusion_counts(group,preds,tp,fp,tn,fn)
      do j = 1, size(v)
         v(j) = real(tn(j),dp) / real(max(1,tn(j)+fp(j)),dp)
      end do
   end function col_specificity

   function col_precision(group, preds) result(v)
      integer, intent(in) :: group(:), preds(:,:)
      real(dp) :: v(size(preds,2))
      integer :: tp(size(preds,2)), fp(size(preds,2)), tn(size(preds,2)), fn(size(preds,2)), j

      call confusion_counts(group,preds,tp,fp,tn,fn)
      do j = 1, size(v)
         v(j) = real(tp(j),dp) / real(max(1,tp(j)+fp(j)),dp)
      end do
   end function col_precision

   function col_fscore(group, preds) result(v)
      integer, intent(in) :: group(:), preds(:,:)
      real(dp) :: v(size(preds,2)), p(size(preds,2)), r(size(preds,2))

      p = col_precision(group,preds)
      r = col_sensitivity(group,preds)
      where (p + r > 0.0_dp)
         v = 2.0_dp * p * r / (p + r)
      elsewhere
         v = 0.0_dp
      end where
   end function col_fscore

   function col_fmi(group, preds) result(v)
      integer, intent(in) :: group(:), preds(:,:)
      real(dp) :: v(size(preds,2)), p(size(preds,2)), r(size(preds,2))

      p = col_precision(group,preds)
      r = col_sensitivity(group,preds)
      v = sqrt(max(0.0_dp,p*r))
   end function col_fmi

   function col_fbscore(group, preds, b) result(v)
      integer, intent(in) :: group(:), preds(:,:)
      real(dp), intent(in) :: b
      real(dp) :: v(size(preds,2)), p(size(preds,2)), r(size(preds,2)), den(size(preds,2))

      p = col_precision(group,preds)
      r = col_sensitivity(group,preds)
      den = b*b*p + r
      where (den > 0.0_dp)
         v = (1.0_dp + b)**2 * p * r / den
      elsewhere
         v = 0.0_dp
      end where
   end function col_fbscore

   elemental pure logical function real_equal(a,b) result(eq)
      real(dp), intent(in) :: a,b
      eq = (a <= b .and. a >= b)
   end function real_equal

end module rfast2_arrays

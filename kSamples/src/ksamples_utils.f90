module ksamples_utils
   use r_kinds, only : dp
   use r_distributions, only : r_pchisq, r_pnorm, r_qnorm, r_dnorm
   use r_sorting, only : r_average_ranks
   use r_descriptive, only : r_variance
   implicit none
   private

   public :: combination_count, log_choose, upper_chisq, upper_normal
   public :: average_ranks, score_average_ties, sample_variance
   public :: shuffle_real, shuffle_integer
   public :: initial_group_labels, next_group_labels, grouped_values
   public :: convolve_distribution, convolve_values, sort_real_in_place
   public :: normal_cdf, normal_pdf, normal_quantile
   public :: unique_sorted, relative_frequencies, multinomial_count


contains

   pure real(dp) function log_choose(n, k) result(value)
      integer, intent(in) :: n !! Population size; must be nonnegative.
      integer, intent(in) :: k !! Number selected; values outside 0:n yield negative infinity.
      integer :: kk, j
      if (n < 0 .or. k < 0 .or. k > n) then
         value = -huge(1.0_dp)
         return
      end if
      kk = min(k, n - k)
      value = 0.0_dp
      do j = 1, kk
         value = value + log(real(n - kk + j, dp)) - log(real(j, dp))
      end do
   end function log_choose

   pure real(dp) function combination_count(n, k) result(value)
      integer, intent(in) :: n !! Population size.
      integer, intent(in) :: k !! Number selected.
      real(dp) :: lc
      lc = log_choose(n, k)
      if (lc <= -0.5_dp*huge(1.0_dp)) then
         value = 0.0_dp
      else if (lc > log(huge(1.0_dp))) then
         value = huge(1.0_dp)
      else
         value = exp(lc)
      end if
   end function combination_count

   pure real(dp) function multinomial_count(ns) result(value)
      integer, intent(in) :: ns(:) !! Positive group sizes whose sum is the pooled size.
      integer :: i, nleft
      value = 1.0_dp
      nleft = sum(ns)
      do i = 1, max(0, size(ns) - 1)
         value = value*combination_count(nleft, ns(i))
         if (value >= huge(1.0_dp)) then
            value = huge(1.0_dp)
            return
         end if
         nleft = nleft - ns(i)
      end do
   end function multinomial_count

   pure real(dp) function upper_chisq(x, df) result(p)
      real(dp), intent(in) :: x !! Chi-square statistic.
      real(dp), intent(in) :: df !! Positive chi-square degrees of freedom.
      p = r_pchisq(x, df, lower_tail=.false.)
   end function upper_chisq

   pure real(dp) function upper_normal(x) result(p)
      real(dp), intent(in) :: x !! Standard-normal variate.
      p = r_pnorm(x, lower_tail=.false.)
   end function upper_normal

   pure elemental real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x !! Standard-normal variate.
      p = r_pnorm(x)
   end function normal_cdf

   pure elemental real(dp) function normal_pdf(x) result(p)
      real(dp), intent(in) :: x !! Standard-normal variate.
      p = r_dnorm(x)
   end function normal_pdf

   pure elemental real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p !! Probability in [0,1].
      x = r_qnorm(p)
   end function normal_quantile

   pure subroutine average_ranks(x, ranks)
      real(dp), intent(in) :: x(:) !! Values to rank; exact equal values are tied.
      real(dp), allocatable, intent(out) :: ranks(:) !! Average 1-based ranks in original order.
      call r_average_ranks(x, ranks)
   end subroutine average_ranks

   pure real(dp) function sample_variance(x) result(value)
      real(dp), intent(in) :: x(:) !! Finite sample values.
      value = r_variance(x)
   end function sample_variance

   pure subroutine score_average_ties(x, scores, averaged)
      real(dp), intent(in) :: x(:) !! Observations whose equal values define tie groups.
      real(dp), intent(in) :: scores(:) !! Ordered scores for ranks 1:size(x).
      real(dp), allocatable, intent(out) :: averaged(:) !! Score assigned to each observation after tie averaging.
      integer, allocatable :: idx(:)
      integer :: i, j, n, left, right, tmp
      real(dp) :: av
      n = size(x)
      if (size(scores) /= n) error stop 'score_average_ties: nonconforming scores'
      allocate(idx(n), averaged(n))
      idx = [(i, i = 1, n)]
      do i = 2, n
         tmp = idx(i)
         j = i - 1
         do while (j >= 1)
            if (x(idx(j)) <= x(tmp)) exit
            idx(j + 1) = idx(j)
            j = j - 1
         end do
         idx(j + 1) = tmp
      end do
      left = 1
      do while (left <= n)
         right = left
         do while (right < n)
            if (x(idx(right + 1)) /= x(idx(left))) exit
            right = right + 1
         end do
         av = sum(scores(left:right))/real(right - left + 1, dp)
         do i = left, right
            averaged(idx(i)) = av
         end do
         left = right + 1
      end do
   end subroutine score_average_ties

   subroutine shuffle_real(x)
      real(dp), intent(inout) :: x(:) !! Values shuffled in place using Fisher-Yates and random_number.
      integer :: i, j
      real(dp) :: u, tmp
      do i = size(x), 2, -1
         call random_number(u)
         j = 1 + int(u*real(i, dp))
         if (j > i) j = i
         tmp = x(i)
         x(i) = x(j)
         x(j) = tmp
      end do
   end subroutine shuffle_real

   subroutine shuffle_integer(x)
      integer, intent(inout) :: x(:) !! Integer values shuffled in place using Fisher-Yates and random_number.
      integer :: i, j, tmp
      real(dp) :: u
      do i = size(x), 2, -1
         call random_number(u)
         j = 1 + int(u*real(i, dp))
         if (j > i) j = i
         tmp = x(i)
         x(i) = x(j)
         x(j) = tmp
      end do
   end subroutine shuffle_integer

   pure subroutine initial_group_labels(ns, labels)
      integer, intent(in) :: ns(:) !! Nonnegative labeled-group sizes whose sum defines the pooled length.
      integer, allocatable, intent(out) :: labels(:) !! Initial nondecreasing multiset of group labels 1:size(ns).
      integer :: g, first, last
      if (any(ns < 0)) error stop 'initial_group_labels: negative group size'
      allocate(labels(sum(ns)))
      first = 1
      do g = 1, size(ns)
         last = first + ns(g) - 1
         if (last >= first) labels(first:last) = g
         first = last + 1
      end do
   end subroutine initial_group_labels

   pure subroutine next_group_labels(labels, has_next)
      integer, intent(inout) :: labels(:) !! Current nondecreasing/lexicographic multiset permutation, updated in place.
      logical, intent(out) :: has_next !! True when labels was advanced to another distinct permutation.
      integer :: i, j, left, right, tmp
      i = size(labels) - 1
      do while (i >= 1)
         if (labels(i) < labels(i + 1)) exit
         i = i - 1
      end do
      if (i < 1) then
         has_next = .false.
         return
      end if
      j = size(labels)
      do while (labels(j) <= labels(i))
         j = j - 1
      end do
      tmp = labels(i)
      labels(i) = labels(j)
      labels(j) = tmp
      left = i + 1
      right = size(labels)
      do while (left < right)
         tmp = labels(left)
         labels(left) = labels(right)
         labels(right) = tmp
         left = left + 1
         right = right - 1
      end do
      has_next = .true.
   end subroutine next_group_labels

   pure subroutine grouped_values(x, labels, ngroups, grouped)
      real(dp), intent(in) :: x(:) !! Pooled values in their fixed original order.
      integer, intent(in) :: labels(:) !! One group label for each pooled value.
      integer, intent(in) :: ngroups !! Number of labeled groups; valid labels are 1:ngroups.
      real(dp), allocatable, intent(out) :: grouped(:) !! Values concatenated by increasing group label.
      integer :: g, i, out_pos
      if (size(labels) /= size(x)) error stop 'grouped_values: shape mismatch'
      if (any(labels < 1) .or. any(labels > ngroups)) error stop 'grouped_values: invalid label'
      allocate(grouped(size(x)))
      out_pos = 0
      do g = 1, ngroups
         do i = 1, size(x)
            if (labels(i) == g) then
               out_pos = out_pos + 1
               grouped(out_pos) = x(i)
            end if
         end do
      end do
   end subroutine grouped_values

   pure subroutine convolve_values(x1, x2, x)
      real(dp), intent(in) :: x1(:) !! First set of values.
      real(dp), intent(in) :: x2(:) !! Second set of values.
      real(dp), allocatable, intent(out) :: x(:) !! All pairwise sums in first-major order.
      integer :: i, j, k
      allocate(x(size(x1)*size(x2)))
      k = 0
      do i = 1, size(x1)
         do j = 1, size(x2)
            k = k + 1
            x(k) = x1(i) + x2(j)
         end do
      end do
   end subroutine convolve_values

   pure subroutine convolve_distribution(x1, p1, x2, p2, x, p)
      real(dp), intent(in) :: x1(:) !! Strictly increasing support of the first distribution.
      real(dp), intent(in) :: p1(:) !! Probabilities corresponding to x1.
      real(dp), intent(in) :: x2(:) !! Strictly increasing support of the second distribution.
      real(dp), intent(in) :: p2(:) !! Probabilities corresponding to x2.
      real(dp), allocatable, intent(out) :: x(:) !! Sorted unique support of the convolution, rounded to 1e-8.
      real(dp), allocatable, intent(out) :: p(:) !! Probabilities corresponding to x.
      real(dp), allocatable :: tx(:), tp(:)
      real(dp) :: value, prob
      integer :: i, j, n, pos
      if (size(x1) /= size(p1) .or. size(x2) /= size(p2)) error stop 'convolve_distribution: shape mismatch'
      allocate(tx(size(x1)*size(x2)), tp(size(x1)*size(x2)))
      n = 0
      do i = 1, size(x1)
         do j = 1, size(x2)
            value = anint(1.0e8_dp*(x1(i) + x2(j)))/1.0e8_dp
            prob = p1(i)*p2(j)
            pos = 1
            do while (pos <= n)
               if (tx(pos) >= value) exit
               pos = pos + 1
            end do
            if (pos <= n) then
               if (tx(pos) == value) then
                  tp(pos) = tp(pos) + prob
                  cycle
               end if
            end if
               if (pos <= n) then
                  tx(pos + 1:n + 1) = tx(pos:n)
                  tp(pos + 1:n + 1) = tp(pos:n)
               end if
               tx(pos) = value
               tp(pos) = prob
               n = n + 1
         end do
      end do
      allocate(x(n), p(n))
      x = tx(:n)
      p = tp(:n)
   end subroutine convolve_distribution

   pure subroutine sort_real_in_place(x)
      real(dp), intent(inout) :: x(:) !! Real vector sorted into ascending order in place.
      integer :: i, j
      real(dp) :: key
      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine sort_real_in_place

   pure subroutine unique_sorted(x, unique_values)
      real(dp), intent(in) :: x(:) !! Input values; order is not significant.
      real(dp), allocatable, intent(out) :: unique_values(:) !! Sorted distinct exact values.
      real(dp), allocatable :: tmp(:)
      integer :: i, n
      tmp = x
      call sort_real_in_place(tmp)
      if (size(tmp) == 0) then
         allocate(unique_values(0))
         return
      end if
      n = 1
      do i = 2, size(tmp)
         if (tmp(i) /= tmp(n)) then
            n = n + 1
            tmp(n) = tmp(i)
         end if
      end do
      allocate(unique_values(n))
      unique_values = tmp(:n)
   end subroutine unique_sorted

   pure subroutine relative_frequencies(values, support, probability, digits)
      real(dp), intent(in) :: values(:) !! Values whose empirical distribution is tabulated.
      real(dp), allocatable, intent(out) :: support(:) !! Sorted unique rounded values.
      real(dp), allocatable, intent(out) :: probability(:) !! Relative frequencies at support points.
      integer, intent(in), optional :: digits !! Decimal digits used before tabulation; default 8.
      real(dp), allocatable :: tmp(:), uniq(:)
      real(dp) :: scale
      integer :: i, j, d
      d = 8
      if (present(digits)) d = digits
      scale = 10.0_dp**d
      allocate(tmp(size(values)))
      tmp = anint(values*scale)/scale
      call unique_sorted(tmp, uniq)
      allocate(support(size(uniq)), probability(size(uniq)))
      support = uniq
      probability = 0.0_dp
      do i = 1, size(tmp)
         do j = 1, size(uniq)
            if (tmp(i) == uniq(j)) then
               probability(j) = probability(j) + 1.0_dp
               exit
            end if
         end do
      end do
      if (size(tmp) > 0) probability = probability/real(size(tmp), dp)
   end subroutine relative_frequencies


end module ksamples_utils

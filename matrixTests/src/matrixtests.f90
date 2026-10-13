! SPDX-License-Identifier: GPL-2.0-only
! Row/column numerical adapters for matrixTests 0.2.3.1.
module matrixtests
   use matrixtests_types
   use matrixtests_kernels
   use matrixtests_ks, only : ks_test
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none
   private
   public :: dp, test_result, mt_ok, mt_invalid, mt_insufficient, mt_degenerate, mt_limit, mt_nonfinite
   public :: row_t_equalvar, col_t_equalvar
   public :: row_t_paired, col_t_paired
   public :: row_t_onesample, col_t_onesample
   public :: row_t_welch, col_t_welch
   public :: row_cor_pearson, col_cor_pearson
   public :: row_f_var, col_f_var
   public :: row_oneway_equalvar, col_oneway_equalvar
   public :: row_oneway_welch, col_oneway_welch
   public :: row_kruskalwallis, col_kruskalwallis
   public :: row_waerden, col_waerden
   public :: row_bartlett, col_bartlett
   public :: row_jarquebera, col_jarquebera
   public :: row_flignerkilleen, col_flignerkilleen
   public :: row_levene, col_levene
   public :: row_brownforsythe, col_brownforsythe
   public :: row_wilcoxon_onesample, col_wilcoxon_onesample
   public :: row_wilcoxon_paired, col_wilcoxon_paired
   public :: row_wilcoxon_twosample, col_wilcoxon_twosample
   public :: row_cosinor, col_cosinor
   public :: row_andersondarling, col_andersondarling
   public :: row_kolmogorovsmirnov_twosample, col_kolmogorovsmirnov_twosample
contains

   pure function row_t_equalvar(x, y, null, alternative, conf_level) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      real(dp), intent(in) :: y(:, :) !! Second matrix; one row broadcasts over x, otherwise test counts must agree.
      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      real(dp), intent(in), optional :: conf_level !! Common confidence level in [0,1], default 0.95; NaN disables intervals.
      type(test_result), allocatable :: r(:)
      integer :: i
      integer :: j
      real(dp) :: mu
      character(9) :: alt
      real(dp) :: conf
      allocate(r(size(x, 1)))
      r = empty_result()
      if (size(y, 1) /= size(x, 1) .and. size(y, 1) /= 1) return
      mu = 0.0_dp
      if (present(null)) mu = null
      alt = 'two.sided'
      if (present(alternative)) then
         if (.not. valid_alternative(alternative)) return
         alt = alternative
      end if
      conf = 0.95_dp
      if (present(conf_level)) conf = conf_level
      do i = 1, size(x, 1)
         j = i
         if (size(y, 1) == 1) j = 1
         r(i) = t_kernel(x(i, :), y(j, :), 2, mu, alt, conf)
      end do
   end function row_t_equalvar

   pure function col_t_equalvar(x, y, null, alternative, conf_level) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      real(dp), intent(in) :: y(:, :) !! Second matrix; one col broadcasts over x, otherwise test counts must agree.
      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      real(dp), intent(in), optional :: conf_level !! Common confidence level in [0,1], default 0.95; NaN disables intervals.
      type(test_result), allocatable :: r(:)
      r = row_t_equalvar(transpose(x), transpose(y), null, alternative, conf_level)
   end function col_t_equalvar

   pure function row_t_paired(x, y, null, alternative, conf_level) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      real(dp), intent(in) :: y(:, :) !! Second matrix; one row broadcasts over x, otherwise test counts must agree.
      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      real(dp), intent(in), optional :: conf_level !! Common confidence level in [0,1], default 0.95; NaN disables intervals.
      type(test_result), allocatable :: r(:)
      integer :: i
      integer :: j
      real(dp) :: mu
      character(9) :: alt
      real(dp) :: conf
      allocate(r(size(x, 1)))
      r = empty_result()
      if (size(y, 1) /= size(x, 1) .and. size(y, 1) /= 1) return
      if (size(x, 2) /= size(y, 2)) return
      mu = 0.0_dp
      if (present(null)) mu = null
      alt = 'two.sided'
      if (present(alternative)) then
         if (.not. valid_alternative(alternative)) return
         alt = alternative
      end if
      conf = 0.95_dp
      if (present(conf_level)) conf = conf_level
      do i = 1, size(x, 1)
         j = i
         if (size(y, 1) == 1) j = 1
         r(i) = t_kernel(x(i, :), y(j, :), 4, mu, alt, conf)
      end do
   end function row_t_paired

   pure function col_t_paired(x, y, null, alternative, conf_level) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      real(dp), intent(in) :: y(:, :) !! Second matrix; one col broadcasts over x, otherwise test counts must agree.
      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      real(dp), intent(in), optional :: conf_level !! Common confidence level in [0,1], default 0.95; NaN disables intervals.
      type(test_result), allocatable :: r(:)
      r = row_t_paired(transpose(x), transpose(y), null, alternative, conf_level)
   end function col_t_paired

   pure function row_t_onesample(x, null, alternative, conf_level) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      real(dp), intent(in), optional :: conf_level !! Common confidence level in [0,1], default 0.95; NaN disables intervals.
      type(test_result), allocatable :: r(:)
      integer :: i
      real(dp) :: mu
      character(9) :: alt
      real(dp) :: conf
      allocate(r(size(x, 1)))
      r = empty_result()
      mu = 0.0_dp
      if (present(null)) mu = null
      alt = 'two.sided'
      if (present(alternative)) then
         if (.not. valid_alternative(alternative)) return
         alt = alternative
      end if
      conf = 0.95_dp
      if (present(conf_level)) conf = conf_level
      do i = 1, size(x, 1)
         r(i) = t_kernel(x(i, :), [real(dp) ::], 1, mu, alt, conf)
      end do
   end function row_t_onesample

   pure function col_t_onesample(x, null, alternative, conf_level) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      real(dp), intent(in), optional :: conf_level !! Common confidence level in [0,1], default 0.95; NaN disables intervals.
      type(test_result), allocatable :: r(:)
      r = row_t_onesample(transpose(x), null, alternative, conf_level)
   end function col_t_onesample

   pure function row_t_welch(x, y, null, alternative, conf_level) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      real(dp), intent(in) :: y(:, :) !! Second matrix; one row broadcasts over x, otherwise test counts must agree.
      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      real(dp), intent(in), optional :: conf_level !! Common confidence level in [0,1], default 0.95; NaN disables intervals.
      type(test_result), allocatable :: r(:)
      integer :: i
      integer :: j
      real(dp) :: mu
      character(9) :: alt
      real(dp) :: conf
      allocate(r(size(x, 1)))
      r = empty_result()
      if (size(y, 1) /= size(x, 1) .and. size(y, 1) /= 1) return
      mu = 0.0_dp
      if (present(null)) mu = null
      alt = 'two.sided'
      if (present(alternative)) then
         if (.not. valid_alternative(alternative)) return
         alt = alternative
      end if
      conf = 0.95_dp
      if (present(conf_level)) conf = conf_level
      do i = 1, size(x, 1)
         j = i
         if (size(y, 1) == 1) j = 1
         r(i) = t_kernel(x(i, :), y(j, :), 3, mu, alt, conf)
      end do
   end function row_t_welch

   pure function col_t_welch(x, y, null, alternative, conf_level) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      real(dp), intent(in) :: y(:, :) !! Second matrix; one col broadcasts over x, otherwise test counts must agree.
      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      real(dp), intent(in), optional :: conf_level !! Common confidence level in [0,1], default 0.95; NaN disables intervals.
      type(test_result), allocatable :: r(:)
      r = row_t_welch(transpose(x), transpose(y), null, alternative, conf_level)
   end function col_t_welch

   pure function row_cor_pearson(x, y, alternative, conf_level) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      real(dp), intent(in) :: y(:, :) !! Second matrix; one row broadcasts over x, otherwise test counts must agree.
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      real(dp), intent(in), optional :: conf_level !! Common confidence level in [0,1], default 0.95; NaN disables intervals.
      type(test_result), allocatable :: r(:)
      integer :: i
      integer :: j
      character(9) :: alt
      real(dp) :: conf
      allocate(r(size(x, 1)))
      r = empty_result()
      if (size(y, 1) /= size(x, 1) .and. size(y, 1) /= 1) return
      if (size(x, 2) /= size(y, 2)) return
      alt = 'two.sided'
      if (present(alternative)) then
         if (.not. valid_alternative(alternative)) return
         alt = alternative
      end if
      conf = 0.95_dp
      if (present(conf_level)) conf = conf_level
      do i = 1, size(x, 1)
         j = i
         if (size(y, 1) == 1) j = 1
         r(i) = correlation_kernel(x(i, :), y(j, :), alt, conf)
      end do
   end function row_cor_pearson

   pure function col_cor_pearson(x, y, alternative, conf_level) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      real(dp), intent(in) :: y(:, :) !! Second matrix; one col broadcasts over x, otherwise test counts must agree.
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      real(dp), intent(in), optional :: conf_level !! Common confidence level in [0,1], default 0.95; NaN disables intervals.
      type(test_result), allocatable :: r(:)
      r = row_cor_pearson(transpose(x), transpose(y), alternative, conf_level)
   end function col_cor_pearson

   pure function row_f_var(x, y, null, alternative, conf_level) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      real(dp), intent(in) :: y(:, :) !! Second matrix; one row broadcasts over x, otherwise test counts must agree.
      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      real(dp), intent(in), optional :: conf_level !! Common confidence level in [0,1], default 0.95; NaN disables intervals.
      type(test_result), allocatable :: r(:)
      integer :: i
      integer :: j
      real(dp) :: mu
      character(9) :: alt
      real(dp) :: conf
      allocate(r(size(x, 1)))
      r = empty_result()
      if (size(y, 1) /= size(x, 1) .and. size(y, 1) /= 1) return
      mu = 1.0_dp
      if (present(null)) mu = null
      alt = 'two.sided'
      if (present(alternative)) then
         if (.not. valid_alternative(alternative)) return
         alt = alternative
      end if
      conf = 0.95_dp
      if (present(conf_level)) conf = conf_level
      do i = 1, size(x, 1)
         j = i
         if (size(y, 1) == 1) j = 1
         r(i) = f_kernel(x(i, :), y(j, :), mu, alt, conf)
      end do
   end function row_f_var

   pure function col_f_var(x, y, null, alternative, conf_level) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      real(dp), intent(in) :: y(:, :) !! Second matrix; one col broadcasts over x, otherwise test counts must agree.
      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      real(dp), intent(in), optional :: conf_level !! Common confidence level in [0,1], default 0.95; NaN disables intervals.
      type(test_result), allocatable :: r(:)
      r = row_f_var(transpose(x), transpose(y), null, alternative, conf_level)
   end function col_f_var

   pure function row_oneway_equalvar(x, g) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.
      type(test_result), allocatable :: r(:)
      integer :: i
      allocate(r(size(x, 1)))
      r = empty_result()
      do i = 1, size(x, 1)
         r(i) = group_kernel(x(i, :), g, 1)
      end do
   end function row_oneway_equalvar

   pure function col_oneway_equalvar(x, g) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.
      type(test_result), allocatable :: r(:)
      r = row_oneway_equalvar(transpose(x), g)
   end function col_oneway_equalvar

   pure function row_oneway_welch(x, g) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.
      type(test_result), allocatable :: r(:)
      integer :: i
      allocate(r(size(x, 1)))
      r = empty_result()
      do i = 1, size(x, 1)
         r(i) = group_kernel(x(i, :), g, 2)
      end do
   end function row_oneway_welch

   pure function col_oneway_welch(x, g) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.
      type(test_result), allocatable :: r(:)
      r = row_oneway_welch(transpose(x), g)
   end function col_oneway_welch

   pure function row_kruskalwallis(x, g) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.
      type(test_result), allocatable :: r(:)
      integer :: i
      allocate(r(size(x, 1)))
      r = empty_result()
      do i = 1, size(x, 1)
         r(i) = group_kernel(x(i, :), g, 3)
      end do
   end function row_kruskalwallis

   pure function col_kruskalwallis(x, g) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.
      type(test_result), allocatable :: r(:)
      r = row_kruskalwallis(transpose(x), g)
   end function col_kruskalwallis

   pure function row_waerden(x, g) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.
      type(test_result), allocatable :: r(:)
      integer :: i
      allocate(r(size(x, 1)))
      r = empty_result()
      do i = 1, size(x, 1)
         r(i) = group_kernel(x(i, :), g, 4)
      end do
   end function row_waerden

   pure function col_waerden(x, g) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.
      type(test_result), allocatable :: r(:)
      r = row_waerden(transpose(x), g)
   end function col_waerden

   pure function row_bartlett(x, g) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.
      type(test_result), allocatable :: r(:)
      integer :: i
      allocate(r(size(x, 1)))
      r = empty_result()
      do i = 1, size(x, 1)
         r(i) = group_kernel(x(i, :), g, 5)
      end do
   end function row_bartlett

   pure function col_bartlett(x, g) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.
      type(test_result), allocatable :: r(:)
      r = row_bartlett(transpose(x), g)
   end function col_bartlett

   function row_jarquebera(x) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      type(test_result), allocatable :: r(:)
      integer :: i
      allocate(r(size(x, 1)))
      r = empty_result()
      do i = 1, size(x, 1)
         r(i) = jarque_kernel(x(i, :))
      end do
   end function row_jarquebera

   function col_jarquebera(x) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      type(test_result), allocatable :: r(:)
      r = row_jarquebera(transpose(x))
   end function col_jarquebera

   pure function row_flignerkilleen(x, g) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.
      type(test_result), allocatable :: r(:)
      integer :: i
      allocate(r(size(x, 1)))
      r = empty_result()
      do i = 1, size(x, 1)
         r(i) = group_kernel(x(i, :), g, 6)
      end do
   end function row_flignerkilleen

   pure function col_flignerkilleen(x, g) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.
      type(test_result), allocatable :: r(:)
      r = row_flignerkilleen(transpose(x), g)
   end function col_flignerkilleen

   pure function row_levene(x, g) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.
      type(test_result), allocatable :: r(:)
      integer :: i
      allocate(r(size(x, 1)))
      r = empty_result()
      do i = 1, size(x, 1)
         r(i) = group_kernel(x(i, :), g, 7)
      end do
   end function row_levene

   pure function col_levene(x, g) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.
      type(test_result), allocatable :: r(:)
      r = row_levene(transpose(x), g)
   end function col_levene

   pure function row_brownforsythe(x, g) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.
      type(test_result), allocatable :: r(:)
      integer :: i
      allocate(r(size(x, 1)))
      r = empty_result()
      do i = 1, size(x, 1)
         r(i) = group_kernel(x(i, :), g, 8)
      end do
   end function row_brownforsythe

   pure function col_brownforsythe(x, g) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      integer, intent(in) :: g(:) !! Group labels along the observation dimension; nonpositive labels omit observations.
      type(test_result), allocatable :: r(:)
      r = row_brownforsythe(transpose(x), g)
   end function col_brownforsythe

   pure function row_wilcoxon_onesample(x, null, alternative, exact, correct) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      logical, intent(in), optional :: exact !! Exact inference request; absent chooses the upstream sample-size rule.
      logical, intent(in), optional :: correct !! Apply continuity correction to approximate Wilcoxon p-values, default true.
      type(test_result), allocatable :: r(:)
      integer :: i
      real(dp) :: mu
      character(9) :: alt
      allocate(r(size(x, 1)))
      r = empty_result()
      mu = 0.0_dp
      if (present(null)) mu = null
      alt = 'two.sided'
      if (present(alternative)) then
         if (.not. valid_alternative(alternative)) return
         alt = alternative
      end if
      do i = 1, size(x, 1)
         r(i) = wilcoxon_kernel(x(i, :), [real(dp) ::], 1, &
            mu, alt, exact, correct)
      end do
   end function row_wilcoxon_onesample

   pure function col_wilcoxon_onesample(x, null, alternative, exact, correct) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      logical, intent(in), optional :: exact !! Exact inference request; absent chooses the upstream sample-size rule.
      logical, intent(in), optional :: correct !! Apply continuity correction to approximate Wilcoxon p-values, default true.
      type(test_result), allocatable :: r(:)
      r = row_wilcoxon_onesample(transpose(x), null, alternative, exact, correct)
   end function col_wilcoxon_onesample

   pure function row_wilcoxon_paired(x, y, null, alternative, exact, correct) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      real(dp), intent(in) :: y(:, :) !! Second matrix; one row broadcasts over x, otherwise test counts must agree.
      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      logical, intent(in), optional :: exact !! Exact inference request; absent chooses the upstream sample-size rule.
      logical, intent(in), optional :: correct !! Apply continuity correction to approximate Wilcoxon p-values, default true.
      type(test_result), allocatable :: r(:)
      integer :: i
      integer :: j
      real(dp) :: mu
      character(9) :: alt
      allocate(r(size(x, 1)))
      r = empty_result()
      if (size(y, 1) /= size(x, 1) .and. size(y, 1) /= 1) return
      if (size(x, 2) /= size(y, 2)) return
      mu = 0.0_dp
      if (present(null)) mu = null
      alt = 'two.sided'
      if (present(alternative)) then
         if (.not. valid_alternative(alternative)) return
         alt = alternative
      end if
      do i = 1, size(x, 1)
         j = i
         if (size(y, 1) == 1) j = 1
         r(i) = wilcoxon_kernel(x(i, :), y(j, :), 3, &
            mu, alt, exact, correct)
      end do
   end function row_wilcoxon_paired

   pure function col_wilcoxon_paired(x, y, null, alternative, exact, correct) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      real(dp), intent(in) :: y(:, :) !! Second matrix; one col broadcasts over x, otherwise test counts must agree.
      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      logical, intent(in), optional :: exact !! Exact inference request; absent chooses the upstream sample-size rule.
      logical, intent(in), optional :: correct !! Apply continuity correction to approximate Wilcoxon p-values, default true.
      type(test_result), allocatable :: r(:)
      r = row_wilcoxon_paired(transpose(x), transpose(y), null, alternative, exact, correct)
   end function col_wilcoxon_paired

   pure function row_wilcoxon_twosample(x, y, null, alternative, exact, correct) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      real(dp), intent(in) :: y(:, :) !! Second matrix; one row broadcasts over x, otherwise test counts must agree.
      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      logical, intent(in), optional :: exact !! Exact inference request; absent chooses the upstream sample-size rule.
      logical, intent(in), optional :: correct !! Apply continuity correction to approximate Wilcoxon p-values, default true.
      type(test_result), allocatable :: r(:)
      integer :: i
      integer :: j
      real(dp) :: mu
      character(9) :: alt
      allocate(r(size(x, 1)))
      r = empty_result()
      if (size(y, 1) /= size(x, 1) .and. size(y, 1) /= 1) return
      mu = 0.0_dp
      if (present(null)) mu = null
      alt = 'two.sided'
      if (present(alternative)) then
         if (.not. valid_alternative(alternative)) return
         alt = alternative
      end if
      do i = 1, size(x, 1)
         j = i
         if (size(y, 1) == 1) j = 1
         r(i) = wilcoxon_kernel(x(i, :), y(j, :), 2, &
            mu, alt, exact, correct)
      end do
   end function row_wilcoxon_twosample

   pure function col_wilcoxon_twosample(x, y, null, alternative, exact, correct) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      real(dp), intent(in) :: y(:, :) !! Second matrix; one col broadcasts over x, otherwise test counts must agree.
      real(dp), intent(in), optional :: null !! Common null mean/location shift (default 0), or variance ratio (default 1).
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      logical, intent(in), optional :: exact !! Exact inference request; absent chooses the upstream sample-size rule.
      logical, intent(in), optional :: correct !! Apply continuity correction to approximate Wilcoxon p-values, default true.
      type(test_result), allocatable :: r(:)
      r = row_wilcoxon_twosample(transpose(x), transpose(y), null, alternative, exact, correct)
   end function col_wilcoxon_twosample

   pure function row_cosinor(x, t, period) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      real(dp), intent(in) :: t(:) !! Times along the observation dimension; NaN times are omitted.
      real(dp), intent(in), optional :: period !! Positive period in time units, default 24.
      type(test_result), allocatable :: r(:)
      integer :: i
      real(dp) :: cycle
      allocate(r(size(x, 1)))
      r = empty_result()
      cycle = 24.0_dp
      if (present(period)) cycle = period
      do i = 1, size(x, 1)
         r(i) = cosinor_kernel(x(i, :), t, cycle)
      end do
   end function row_cosinor

   pure function col_cosinor(x, t, period) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      real(dp), intent(in) :: t(:) !! Times along the observation dimension; NaN times are omitted.
      real(dp), intent(in), optional :: period !! Positive period in time units, default 24.
      type(test_result), allocatable :: r(:)
      r = row_cosinor(transpose(x), t, period)
   end function col_cosinor

   pure function row_andersondarling(x) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      type(test_result), allocatable :: r(:)
      integer :: i
      allocate(r(size(x, 1)))
      r = empty_result()
      do i = 1, size(x, 1)
         r(i) = ad_kernel(x(i, :))
      end do
   end function row_andersondarling

   pure function col_andersondarling(x) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      type(test_result), allocatable :: r(:)
      r = row_andersondarling(transpose(x))
   end function col_andersondarling

   pure function row_kolmogorovsmirnov_twosample(x, y, alternative, exact) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by row; IEEE NaNs are missing.
      real(dp), intent(in) :: y(:, :) !! Second matrix; one row broadcasts over x, otherwise test counts must agree.
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      logical, intent(in), optional :: exact !! Exact inference request; absent chooses the upstream sample-size rule.
      type(test_result), allocatable :: r(:)
      integer :: i
      integer :: j
      real(dp), allocatable :: a(:), b(:)
      integer :: status
      character(9) :: alt
      allocate(r(size(x, 1)))
      r = empty_result()
      if (size(y, 1) /= size(x, 1) .and. size(y, 1) /= 1) return
      alt = 'two.sided'
      if (present(alternative)) then
         if (.not. valid_alternative(alternative)) return
         alt = alternative
      end if
      do i = 1, size(x, 1)
         j = i
         if (size(y, 1) == 1) j = 1
         a = pack(x(i, :), .not. ieee_is_nan(x(i, :)))
         b = pack(y(j, :), .not. ieee_is_nan(y(j, :)))
         r(i)%obs_x = size(a)
         r(i)%obs_y = size(b)
         r(i)%obs_tot = size(a) + size(b)
         r(i)%alternative = alt
         call ks_test(a, b, r(i)%statistic, r(i)%pvalue, status, alt, exact, r(i)%exact)
         r(i)%status = mt_ok
         if (status /= 0) r(i)%status = mt_invalid
         if (status == 3) r(i)%status = mt_limit
      end do
   end function row_kolmogorovsmirnov_twosample

   pure function col_kolmogorovsmirnov_twosample(x, y, alternative, exact) result(r)
      real(dp), intent(in) :: x(:, :) !! Data matrix, tested independently by col; IEEE NaNs are missing.
      real(dp), intent(in) :: y(:, :) !! Second matrix; one col broadcasts over x, otherwise test counts must agree.
      character(*), intent(in), optional :: alternative !! Common alternative: two.sided (default), less, or greater.
      logical, intent(in), optional :: exact !! Exact inference request; absent chooses the upstream sample-size rule.
      type(test_result), allocatable :: r(:)
      r = row_kolmogorovsmirnov_twosample(transpose(x), transpose(y), alternative, exact)
   end function col_kolmogorovsmirnov_twosample
end module matrixtests

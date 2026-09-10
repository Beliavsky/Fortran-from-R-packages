! SPDX-License-Identifier: GPL-2.0-or-later
! Wavelet cross-validation and denoising wrappers translated from wavethresh 4.7.3.
module wavethresh_cv
   use wavethresh_types, only : dp, wd_t, rss_result_t, wavelet_cv_result_t
   use wavethresh_transform_1d, only : wd, wr_wd, wst, wr_wst, is_power_of_two
   use wavethresh_threshold, only : threshold_wd, threshold_wd_levels
   use wavethresh_stats, only : ssq, dof, l2norm, madmad
   implicit none
   private

   public :: denwr, rsswav, crsswav, wavelet_cv, cwcv, full_wavelet_cv, get_rss_wst, wvcvlrss, wst_cv
   public :: get_rss_wst_levels, wvcvlrss_levels, wst_cvl

contains

   function denwr(object, start_level) result(data)
      type(wd_t), intent(in) :: object !! Decimated wavelet object reconstructed by the legacy denwr operation.
      integer, intent(in), optional :: start_level !! Coarsest zero-based level included in reconstruction; default is zero.
      real(dp), allocatable :: data(:)
      data = wr_wd(object, start_level)
   end function denwr

   function rsswav(noisy, value, filter_number, family, threshold_type, ll) result(result)
      real(dp), intent(in) :: noisy(:) !! Noisy dyadic series used by odd/even wavelet cross-validation.
      real(dp), intent(in), optional :: value !! Manual hard/soft threshold; default is 1.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      character(len=*), intent(in), optional :: threshold_type !! hard or soft threshold rule; default is hard.
      integer, intent(in), optional :: ll !! First zero-based detail level thresholded; default is 3.
      type(rss_result_t) :: result
      type(wd_t) :: odd_object
      type(wd_t) :: even_object
      type(wd_t) :: full_object
      type(wd_t) :: thresholded
      real(dp), allocatable :: odd(:)
      real(dp), allocatable :: even(:)
      real(dp), allocatable :: odd_fit(:)
      real(dp), allocatable :: even_fit(:)
      real(dp), allocatable :: interp(:)
      integer, allocatable :: levels(:)
      real(dp) :: threshold
      real(dp) :: fnum
      integer :: first_level
      integer :: m
      integer :: i
      integer :: nlev
      character(len=24) :: fam
      character(len=16) :: kind

      result%threshold = 1.0_dp
      if (present(value)) result%threshold = value
      threshold = result%threshold
      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      kind = "hard"
      if (present(threshold_type)) kind = threshold_type
      first_level = 3
      if (present(ll)) first_level = ll

      if (size(noisy) < 4 .or. mod(size(noisy), 2) /= 0) then
         result%message = "rsswav requires an even data length of at least four"
         return
      end if
      m = size(noisy) / 2
      if (.not. is_power_of_two(m)) then
         result%message = "rsswav requires half the data length to be a power of two"
         return
      end if
      allocate(odd(m), even(m), interp(m))
      odd = noisy(1:size(noisy):2)
      even = noisy(2:size(noisy):2)
      odd_object = wd(odd, filter_number=fnum, family=trim(fam))
      even_object = wd(even, filter_number=fnum, family=trim(fam))
      if (.not. odd_object%ok .or. .not. even_object%ok) then
         result%message = "wavelet decomposition failed in rsswav"
         return
      end if
      nlev = odd_object%nlevels
      call make_levels(first_level, nlev, levels)
      thresholded = threshold_wd(odd_object, levels=levels, threshold_type=trim(kind), &
         policy="manual", value=threshold)
      odd_fit = wr_wd(thresholded)
      thresholded = threshold_wd(even_object, levels=levels, threshold_type=trim(kind), &
         policy="manual", value=threshold)
      even_fit = wr_wd(thresholded)

      interp(1) = even(1)
      do i = 2, m
         interp(i) = 0.5_dp * (even(i - 1) + even(i))
      end do
      result%ssq = 0.5_dp * ssq(interp, odd_fit)
      interp(1) = odd(1)
      do i = 2, m
         interp(i) = 0.5_dp * (odd(i - 1) + odd(i))
      end do
      result%ssq = result%ssq + 0.5_dp * ssq(interp, even_fit)

      full_object = wd(noisy, filter_number=fnum, family=trim(fam))
      if (.not. full_object%ok) then
         result%message = "full wavelet decomposition failed in rsswav"
         return
      end if
      call make_levels(first_level, full_object%nlevels, levels)
      thresholded = threshold_wd(full_object, levels=levels, threshold_type=trim(kind), &
         policy="manual", value=threshold)
      result%df = dof(thresholded)
      result%ok = .true.
      result%message = "ok"
   end function rsswav

   function crsswav(noisy, value, filter_number, family, threshold_type, ll) result(result)
      real(dp), intent(in) :: noisy(:) !! Noisy dyadic series used by the compiled-style cross-validation objective.
      real(dp), intent(in), optional :: value !! Manual hard/soft threshold; default is 1.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      character(len=*), intent(in), optional :: threshold_type !! hard or soft threshold rule; default is hard.
      integer, intent(in), optional :: ll !! First zero-based detail level thresholded; default is 3.
      type(rss_result_t) :: result
      result = rsswav(noisy, value, filter_number, family, threshold_type, ll)
   end function crsswav

   function get_rss_wst(ndata, threshold, levels, filter_number, family, threshold_type) result(value)
      real(dp), intent(in) :: ndata(:) !! Even-length dyadic series whose odd/even stationary-wavelet prediction error is computed.
      real(dp), intent(in) :: threshold !! Manual threshold applied uniformly to the selected stationary detail levels.
      integer, intent(in) :: levels(:) !! Zero-based stationary detail levels thresholded in both odd and even fits.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      character(len=*), intent(in), optional :: threshold_type !! hard or soft threshold rule; default is soft.
      real(dp) :: value
      type(wd_t) :: odd_object
      type(wd_t) :: even_object
      type(wd_t) :: thresholded
      real(dp), allocatable :: odd(:)
      real(dp), allocatable :: even(:)
      real(dp), allocatable :: odd_fit(:)
      real(dp), allocatable :: even_fit(:)
      real(dp), allocatable :: prediction(:)
      real(dp) :: fnum
      character(len=24) :: fam
      character(len=16) :: kind
      integer :: i
      integer :: m

      value = huge(1.0_dp)
      if (size(ndata) < 4 .or. mod(size(ndata), 2) /= 0) return
      m = size(ndata) / 2
      if (.not. is_power_of_two(m)) return
      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      kind = "soft"
      if (present(threshold_type)) kind = threshold_type

      odd = ndata(1:size(ndata):2)
      even = ndata(2:size(ndata):2)
      odd_object = wst(odd, filter_number=fnum, family=trim(fam))
      even_object = wst(even, filter_number=fnum, family=trim(fam))
      if (.not. odd_object%ok .or. .not. even_object%ok) return
      do i = 1, size(levels)
         if (levels(i) < 0 .or. levels(i) >= odd_object%nlevels) return
      end do

      thresholded = threshold_wd(odd_object, levels=levels, threshold_type=trim(kind), &
         policy="manual", value=threshold)
      odd_fit = wr_wst(thresholded)
      thresholded = threshold_wd(even_object, levels=levels, threshold_type=trim(kind), &
         policy="manual", value=threshold)
      even_fit = wr_wst(thresholded)
      if (size(odd_fit) /= m .or. size(even_fit) /= m) return

      allocate(prediction(m))
      do i = 1, m - 1
         prediction(i) = 0.5_dp * (odd_fit(i) + odd_fit(i + 1))
      end do
      prediction(m) = 0.5_dp * (odd_fit(m) + odd_fit(1))
      value = l2norm(prediction, even)
      do i = 1, m - 1
         prediction(i) = 0.5_dp * (even_fit(i) + even_fit(i + 1))
      end do
      prediction(m) = 0.5_dp * (even_fit(m) + even_fit(1))
      value = 0.5_dp * (value + l2norm(prediction, odd))
   end function get_rss_wst

   function wvcvlrss(threshold, ndata, levels, filter_number, family, threshold_type) result(value)
      real(dp), intent(in) :: threshold !! Manual scalar threshold forwarded to the stationary-wavelet CV objective.
      real(dp), intent(in) :: ndata(:) !! Even-length dyadic series used to evaluate the CV prediction error.
      integer, intent(in) :: levels(:) !! Zero-based stationary detail levels included in the CV objective.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      character(len=*), intent(in), optional :: threshold_type !! hard or soft threshold rule; default is soft.
      real(dp) :: value

      value = get_rss_wst(ndata, threshold, levels, filter_number, family, threshold_type)
   end function wvcvlrss

   function get_rss_wst_levels(ndata, thresholds, levels, filter_number, family, threshold_type) result(value)
      real(dp), intent(in) :: ndata(:) !! Even-length dyadic series whose odd/even stationary-wavelet CV error is evaluated.
      real(dp), intent(in) :: thresholds(:) !! Manual threshold values paired elementwise with the selected levels.
      integer, intent(in) :: levels(:) !! Zero-based stationary detail levels paired elementwise with thresholds.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      character(len=*), intent(in), optional :: threshold_type !! hard or soft threshold rule; default is soft.
      real(dp) :: value
      type(wd_t) :: odd_object
      type(wd_t) :: even_object
      type(wd_t) :: thresholded
      real(dp), allocatable :: odd(:)
      real(dp), allocatable :: even(:)
      real(dp), allocatable :: odd_fit(:)
      real(dp), allocatable :: even_fit(:)
      real(dp), allocatable :: prediction(:)
      real(dp) :: fnum
      character(len=24) :: fam
      character(len=16) :: kind
      integer :: i
      integer :: m

      value = huge(1.0_dp)
      if (size(thresholds) /= size(levels) .or. any(thresholds < 0.0_dp)) return
      if (size(ndata) < 4 .or. mod(size(ndata), 2) /= 0) return
      m = size(ndata) / 2
      if (.not. is_power_of_two(m)) return
      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      kind = "soft"
      if (present(threshold_type)) kind = threshold_type

      odd = ndata(1:size(ndata):2)
      even = ndata(2:size(ndata):2)
      odd_object = wst(odd, filter_number=fnum, family=trim(fam))
      even_object = wst(even, filter_number=fnum, family=trim(fam))
      if (.not. odd_object%ok .or. .not. even_object%ok) return
      do i = 1, size(levels)
         if (levels(i) < 0 .or. levels(i) >= odd_object%nlevels) return
      end do

      thresholded = threshold_wd_levels(odd_object, levels, thresholds, trim(kind))
      if (.not. thresholded%ok) return
      odd_fit = wr_wst(thresholded)
      thresholded = threshold_wd_levels(even_object, levels, thresholds, trim(kind))
      if (.not. thresholded%ok) return
      even_fit = wr_wst(thresholded)
      if (size(odd_fit) /= m .or. size(even_fit) /= m) return

      allocate(prediction(m))
      do i = 1, m - 1
         prediction(i) = 0.5_dp * (odd_fit(i) + odd_fit(i + 1))
      end do
      prediction(m) = 0.5_dp * (odd_fit(m) + odd_fit(1))
      value = l2norm(prediction, even)
      do i = 1, m - 1
         prediction(i) = 0.5_dp * (even_fit(i) + even_fit(i + 1))
      end do
      prediction(m) = 0.5_dp * (even_fit(m) + even_fit(1))
      value = 0.5_dp * (value + l2norm(prediction, odd))
   end function get_rss_wst_levels

   function wvcvlrss_levels(thresholds, ndata, levels, filter_number, family, threshold_type) result(value)
      real(dp), intent(in) :: thresholds(:) !! Manual level-specific thresholds forwarded to the stationary CV objective.
      real(dp), intent(in) :: ndata(:) !! Even-length dyadic series used to evaluate the CV prediction error.
      integer, intent(in) :: levels(:) !! Zero-based stationary detail levels paired with thresholds.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      character(len=*), intent(in), optional :: threshold_type !! hard or soft threshold rule; default is soft.
      real(dp) :: value

      value = get_rss_wst_levels(ndata, thresholds, levels, filter_number, family, threshold_type)
   end function wvcvlrss_levels

   function wavelet_cv(noisy, filter_number, family, threshold_type, tolerance, ll) result(result)
      real(dp), intent(in) :: noisy(:) !! Noisy dyadic series for golden-section wavelet cross-validation.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      character(len=*), intent(in), optional :: threshold_type !! hard or soft threshold rule; default is soft.
      real(dp), intent(in), optional :: tolerance !! Relative golden-section convergence tolerance; default is 0.01.
      integer, intent(in), optional :: ll !! First zero-based detail level thresholded; default is 3.
      type(wavelet_cv_result_t) :: result
      type(wd_t) :: object
      type(wd_t) :: thresholded
      type(rss_result_t) :: trial
      integer, allocatable :: levels(:)
      real(dp), allocatable :: thresholds(:)
      real(dp), allocatable :: errors(:)
      real(dp) :: fnum
      real(dp) :: tol
      real(dp) :: ax
      real(dp) :: bx
      real(dp) :: cx
      real(dp) :: x0
      real(dp) :: x1
      real(dp) :: x2
      real(dp) :: x3
      real(dp) :: f1
      real(dp) :: f2
      real(dp) :: ratio
      real(dp) :: complement
      real(dp) :: correction
      integer :: first_level
      integer :: iterations
      character(len=24) :: fam
      character(len=16) :: kind

      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      kind = "soft"
      if (present(threshold_type)) kind = threshold_type
      tol = 0.01_dp
      if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
      first_level = 3
      if (present(ll)) first_level = ll

      object = wd(noisy, filter_number=fnum, family=trim(fam))
      if (.not. object%ok) then
         result%message = "wavelet decomposition failed in WaveletCV"
         return
      end if
      call make_levels(first_level, object%nlevels, levels)
      thresholded = threshold_wd(object, levels=levels, threshold_type=trim(kind), policy="universal")
      result%universal_df = dof(thresholded)
      result%universal_reconstruction = wr_wd(thresholded)
      result%universal_threshold = pooled_universal_threshold(object, levels)
      if (result%universal_threshold <= 0.0_dp) then
         result%cv_reconstruction = result%universal_reconstruction
         result%cv_threshold = 0.0_dp
         result%cv_df = result%universal_df
         result%trial_thresholds = [0.0_dp]
         result%trial_errors = [0.0_dp]
         result%ok = .true.
         result%message = "ok"
         return
      end if

      ratio = 0.61803399_dp
      complement = 1.0_dp - ratio
      ax = 0.0_dp
      bx = 0.5_dp * result%universal_threshold
      cx = result%universal_threshold
      x0 = ax
      x3 = cx
      if (abs(cx - bx) > abs(bx - ax)) then
         x1 = bx
         x2 = bx + complement * (cx - bx)
      else
         x2 = bx
         x1 = bx - complement * (bx - ax)
      end if
      trial = rsswav(noisy, x1, fnum, trim(fam), trim(kind), first_level)
      if (.not. trial%ok) then
         result%message = trial%message
         return
      end if
      f1 = trial%ssq
      trial = rsswav(noisy, x2, fnum, trim(fam), trim(kind), first_level)
      if (.not. trial%ok) then
         result%message = trial%message
         return
      end if
      f2 = trial%ssq
      thresholds = [ax, cx, x1, x2]
      errors = [rss_value(noisy, ax, fnum, fam, kind, first_level), &
         rss_value(noisy, cx, fnum, fam, kind, first_level), f1, f2]
      iterations = 0
      do while (abs(x3 - x0) > tol * max(abs(x1) + abs(x2), 1.0_dp))
         iterations = iterations + 1
         if (iterations > 200) exit
         if (f2 < f1) then
            x0 = x1
            x1 = x2
            x2 = ratio * x1 + complement * x3
            f1 = f2
            f2 = rss_value(noisy, x2, fnum, fam, kind, first_level)
            thresholds = [thresholds, x2]
            errors = [errors, f2]
         else
            x3 = x2
            x2 = x1
            x1 = ratio * x2 + complement * x0
            f2 = f1
            f1 = rss_value(noisy, x1, fnum, fam, kind, first_level)
            thresholds = [thresholds, x1]
            errors = [errors, f1]
         end if
      end do
      if (f1 < f2) then
         result%cv_threshold = x1
      else
         result%cv_threshold = x2
      end if
      if (size(noisy) > 2) then
         correction = sqrt(max(1.0_dp - log(2.0_dp) / log(real(size(noisy), dp)), epsilon(1.0_dp)))
         result%cv_threshold = result%cv_threshold / correction
      end if
      thresholded = threshold_wd(object, levels=levels, threshold_type=trim(kind), &
         policy="manual", value=result%cv_threshold)
      result%cv_df = dof(thresholded)
      result%cv_reconstruction = wr_wd(thresholded)
      result%trial_thresholds = thresholds
      result%trial_errors = errors
      result%ok = .true.
      result%message = "ok"
   end function wavelet_cv

   function cwcv(noisy, ll, filter_number, family, threshold_type, tolerance) result(result)
      real(dp), intent(in) :: noisy(:) !! Noisy dyadic series for the compiled-style CWCV cross-validation workflow.
      integer, intent(in), optional :: ll !! First zero-based detail level thresholded; default is three.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      character(len=*), intent(in), optional :: threshold_type !! Hard or soft thresholding rule; default is soft.
      real(dp), intent(in), optional :: tolerance !! Relative golden-section convergence tolerance; default is 0.01.
      type(wavelet_cv_result_t) :: result

      result = wavelet_cv(noisy, filter_number, family, threshold_type, tolerance, ll)
   end function cwcv

   function full_wavelet_cv(noisy, ll, threshold_type, filter_number, family, tolerance) result(value)
      real(dp), intent(in) :: noisy(:) !! Noisy dyadic series whose cross-validated wavelet threshold is requested.
      integer, intent(in), optional :: ll !! First zero-based detail level included in cross-validation; default is three.
      character(len=*), intent(in), optional :: threshold_type !! hard or soft thresholding rule; default is soft.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      real(dp), intent(in), optional :: tolerance !! Relative golden-section convergence tolerance; default is 0.01.
      real(dp) :: value
      type(wavelet_cv_result_t) :: result

      result = wavelet_cv(noisy, filter_number, family, threshold_type, tolerance, ll)
      if (result%ok) then
         value = result%cv_threshold
      else
         value = huge(1.0_dp)
      end if
   end function full_wavelet_cv

   function wst_cv(ndata, ll, threshold_type, filter_number, family, tolerance) result(result)
      real(dp), intent(in) :: ndata(:) !! Dyadic series for stationary-wavelet golden-section cross-validation.
      integer, intent(in), optional :: ll !! First zero-based detail level included in cross-validation; default is three.
      character(len=*), intent(in), optional :: threshold_type !! hard or soft thresholding rule; default is soft.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      real(dp), intent(in), optional :: tolerance !! Relative golden-section convergence tolerance; default is 0.01.
      type(wavelet_cv_result_t) :: result
      type(wd_t) :: object
      type(wd_t) :: thresholded
      integer, allocatable :: all_levels(:)
      integer, allocatable :: cv_levels(:)
      real(dp), allocatable :: thresholds(:)
      real(dp), allocatable :: errors(:)
      real(dp) :: fnum
      real(dp) :: tol
      real(dp) :: ax
      real(dp) :: bx
      real(dp) :: cx
      real(dp) :: x0
      real(dp) :: x1
      real(dp) :: x2
      real(dp) :: x3
      real(dp) :: f1
      real(dp) :: f2
      real(dp) :: ratio
      real(dp) :: complement
      real(dp) :: correction
      integer :: first_level
      integer :: iterations
      character(len=24) :: fam
      character(len=16) :: kind

      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      kind = "soft"
      if (present(threshold_type)) kind = threshold_type
      tol = 0.01_dp
      if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
      first_level = 3
      if (present(ll)) first_level = ll

      object = wst(ndata, filter_number=fnum, family=trim(fam))
      if (.not. object%ok) then
         result%message = "stationary wavelet decomposition failed in wstCV"
         return
      end if
      call make_levels(first_level, object%nlevels, all_levels)
      if (size(all_levels) == 0) then
         result%message = "wstCV requires at least one thresholdable detail level"
         return
      end if
      thresholded = threshold_wd(object, levels=all_levels, threshold_type=trim(kind), policy="universal")
      result%universal_threshold = pooled_universal_threshold(object, all_levels)
      result%universal_df = dof(thresholded)
      result%universal_reconstruction = wr_wst(thresholded)

      if (first_level >= object%nlevels - 1) then
         allocate(cv_levels(0))
      else
         allocate(cv_levels(object%nlevels - first_level - 1))
         cv_levels = [(iterations, iterations = first_level, object%nlevels - 2)]
      end if
      if (size(cv_levels) == 0 .or. result%universal_threshold <= 0.0_dp) then
         result%cv_threshold = 0.0_dp
         result%cv_df = result%universal_df
         result%cv_reconstruction = result%universal_reconstruction
         result%trial_thresholds = [0.0_dp]
         result%trial_errors = [0.0_dp]
         result%ok = .true.
         result%message = "ok"
         return
      end if

      ratio = 0.61803399_dp
      complement = 1.0_dp - ratio
      ax = 0.0_dp
      bx = 0.5_dp * result%universal_threshold
      cx = result%universal_threshold
      x0 = ax
      x3 = cx
      if (abs(cx - bx) > abs(bx - ax)) then
         x1 = bx
         x2 = bx + complement * (cx - bx)
      else
         x2 = bx
         x1 = bx - complement * (bx - ax)
      end if
      f1 = get_rss_wst(ndata, x1, cv_levels, fnum, trim(fam), trim(kind))
      f2 = get_rss_wst(ndata, x2, cv_levels, fnum, trim(fam), trim(kind))
      if (f1 >= huge(1.0_dp) .or. f2 >= huge(1.0_dp)) then
         result%message = "stationary-wavelet CV objective failed"
         return
      end if
      thresholds = [ax, cx, x1, x2]
      errors = [get_rss_wst(ndata, ax, cv_levels, fnum, trim(fam), trim(kind)), &
         get_rss_wst(ndata, cx, cv_levels, fnum, trim(fam), trim(kind)), f1, f2]
      iterations = 0
      do while (abs(x3 - x0) > tol * max(abs(x1) + abs(x2), epsilon(1.0_dp)))
         iterations = iterations + 1
         if (iterations > 200) exit
         if (f2 < f1) then
            x0 = x1
            x1 = x2
            x2 = ratio * x1 + complement * x3
            f1 = f2
            f2 = get_rss_wst(ndata, x2, cv_levels, fnum, trim(fam), trim(kind))
            thresholds = [thresholds, x2]
            errors = [errors, f2]
         else
            x3 = x2
            x2 = x1
            x1 = ratio * x2 + complement * x0
            f2 = f1
            f1 = get_rss_wst(ndata, x1, cv_levels, fnum, trim(fam), trim(kind))
            thresholds = [thresholds, x1]
            errors = [errors, f1]
         end if
      end do
      if (f1 < f2) then
         result%cv_threshold = x1
      else
         result%cv_threshold = x2
      end if
      correction = sqrt(max(1.0_dp - log(2.0_dp) / log(real(size(ndata), dp)), epsilon(1.0_dp)))
      result%cv_threshold = result%cv_threshold / correction

      thresholded = threshold_wd(object, levels=cv_levels, threshold_type=trim(kind), &
         policy="manual", value=result%cv_threshold)
      thresholded = threshold_wd(thresholded, levels=[object%nlevels - 1], &
         threshold_type=trim(kind), policy="universal", by_level=.true.)
      result%cv_df = dof(thresholded)
      result%cv_reconstruction = wr_wst(thresholded)
      result%trial_thresholds = thresholds
      result%trial_errors = errors
      result%ok = .true.
      result%message = "ok"
   end function wst_cv

   function wst_cvl(ndata, ll, threshold_type, filter_number, family, tolerance) result(result)
      real(dp), intent(in) :: ndata(:) !! Dyadic series for level-specific stationary-wavelet cross-validation.
      integer, intent(in), optional :: ll !! First zero-based detail level included in cross-validation; default is three.
      character(len=*), intent(in), optional :: threshold_type !! hard or soft thresholding rule; default is soft.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      real(dp), intent(in), optional :: tolerance !! Relative coordinate-search convergence tolerance; default is 0.01.
      type(wavelet_cv_result_t) :: result
      type(wd_t) :: object
      type(wd_t) :: thresholded
      integer, allocatable :: levels(:)
      real(dp), allocatable :: thresholds(:)
      real(dp), allocatable :: previous(:)
      real(dp), allocatable :: trial(:)
      real(dp) :: fnum
      real(dp) :: tol
      real(dp) :: universal
      real(dp) :: top_threshold
      real(dp) :: left
      real(dp) :: right
      real(dp) :: x1
      real(dp) :: x2
      real(dp) :: f1
      real(dp) :: f2
      real(dp) :: ratio
      real(dp) :: complement
      real(dp) :: relative_change
      integer :: first_level
      integer :: coordinate
      integer :: iteration
      integer :: sweep
      character(len=24) :: fam
      character(len=16) :: kind

      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      kind = "soft"
      if (present(threshold_type)) kind = threshold_type
      tol = 0.01_dp
      if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
      first_level = 3
      if (present(ll)) first_level = ll

      object = wst(ndata, filter_number=fnum, family=trim(fam))
      if (.not. object%ok) then
         result%message = "wstCVl stationary wavelet decomposition failed"
         return
      end if
      if (first_level >= object%nlevels - 1) then
         result%message = "wstCVl requires at least one level below the finest detail level"
         return
      end if
      allocate(levels(object%nlevels - first_level - 1))
      levels = [(iteration, iteration = first_level, object%nlevels - 2)]
      universal = pooled_universal_threshold(object, levels)
      result%universal_threshold = universal
      if (universal <= 0.0_dp) then
         allocate(thresholds(size(levels)), source=0.0_dp)
      else
         allocate(thresholds(size(levels)), source=0.5_dp * universal)
      end if
      ratio = 0.61803399_dp
      complement = 1.0_dp - ratio
      allocate(trial(size(thresholds)))

      do sweep = 1, 20
         previous = thresholds
         do coordinate = 1, size(thresholds)
            left = 0.0_dp
            right = universal
            if (right <= 0.0_dp) cycle
            x1 = right * complement
            x2 = right * ratio
            trial = thresholds
            trial(coordinate) = x1
            f1 = get_rss_wst_levels(ndata, trial, levels, fnum, trim(fam), trim(kind))
            trial(coordinate) = x2
            f2 = get_rss_wst_levels(ndata, trial, levels, fnum, trim(fam), trim(kind))
            do iteration = 1, 80
               if (abs(right - left) <= tol * max(abs(x1) + abs(x2), epsilon(1.0_dp))) exit
               if (f2 < f1) then
                  left = x1
                  x1 = x2
                  f1 = f2
                  x2 = ratio * x1 + complement * right
                  trial = thresholds
                  trial(coordinate) = x2
                  f2 = get_rss_wst_levels(ndata, trial, levels, fnum, trim(fam), trim(kind))
               else
                  right = x2
                  x2 = x1
                  f2 = f1
                  x1 = ratio * x2 + complement * left
                  trial = thresholds
                  trial(coordinate) = x1
                  f1 = get_rss_wst_levels(ndata, trial, levels, fnum, trim(fam), trim(kind))
               end if
            end do
            if (f1 < f2) then
               thresholds(coordinate) = x1
            else
               thresholds(coordinate) = x2
            end if
         end do
         relative_change = maxval(abs(thresholds - previous)) / max(maxval(abs(previous)), 1.0_dp)
         if (relative_change <= tol) exit
      end do

      thresholded = threshold_wd_levels(object, levels, thresholds, trim(kind))
      if (.not. thresholded%ok) then
         result%message = thresholded%message
         return
      end if
      top_threshold = pooled_universal_threshold(thresholded, [object%nlevels - 1])
      thresholded = threshold_wd(thresholded, levels=[object%nlevels - 1], threshold_type=trim(kind), &
         policy="manual", value=top_threshold)
      result%cv_thresholds = [thresholds, top_threshold]
      result%cv_threshold = thresholds(1)
      result%cv_df = dof(thresholded)
      result%cv_reconstruction = wr_wst(thresholded)
      result%universal_df = dof(threshold_wd(object, levels=levels, threshold_type=trim(kind), policy="universal"))
      result%universal_reconstruction = wr_wst(threshold_wd(object, levels=levels, threshold_type=trim(kind), &
         policy="universal"))
      result%trial_errors = [get_rss_wst_levels(ndata, previous, levels, fnum, trim(fam), trim(kind)), &
         get_rss_wst_levels(ndata, thresholds, levels, fnum, trim(fam), trim(kind))]
      result%ok = .true.
      result%message = "ok"
   end function wst_cvl

   subroutine make_levels(first_level, nlevels, levels)
      integer, intent(in) :: first_level !! First zero-based detail level selected.
      integer, intent(in) :: nlevels !! Number of available detail levels.
      integer, allocatable, intent(out) :: levels(:) !! Consecutive selected levels from first_level through nlevels-1.
      integer :: first
      integer :: i
      first = max(first_level, 0)
      if (first >= nlevels) then
         allocate(levels(0))
      else
         allocate(levels(nlevels - first))
         levels = [(i, i = first, nlevels - 1)]
      end if
   end subroutine make_levels

   function pooled_universal_threshold(object, levels) result(value)
      type(wd_t), intent(in) :: object !! Wavelet decomposition providing the selected detail coefficients.
      integer, intent(in) :: levels(:) !! Zero-based detail levels pooled for the universal threshold.
      real(dp) :: value
      real(dp), allocatable :: pooled(:)
      real(dp) :: noise
      integer :: total
      integer :: position
      integer :: i
      integer :: level
      total = 0
      do i = 1, size(levels)
         level = levels(i)
         if (level >= 0 .and. level < object%nlevels) total = total + size(object%detail(level)%values)
      end do
      if (total <= 1) then
         value = 0.0_dp
         return
      end if
      allocate(pooled(total))
      position = 1
      do i = 1, size(levels)
         level = levels(i)
         if (level < 0 .or. level >= object%nlevels) cycle
         pooled(position:position + size(object%detail(level)%values) - 1) = object%detail(level)%values
         position = position + size(object%detail(level)%values)
      end do
      noise = sqrt(madmad(pooled))
      if (noise <= tiny(1.0_dp)) then
         value = 0.0_dp
      else
         value = sqrt(2.0_dp * log(real(total, dp))) * noise
      end if
   end function pooled_universal_threshold

   function rss_value(noisy, threshold, filter_number, family, threshold_type, first_level) result(value)
      real(dp), intent(in) :: noisy(:) !! Noisy series passed through rsswav.
      real(dp), intent(in) :: threshold !! Trial manual threshold.
      real(dp), intent(in) :: filter_number !! Wavethresh filter number.
      character(len=*), intent(in) :: family !! Wavethresh filter family.
      character(len=*), intent(in) :: threshold_type !! hard or soft threshold rule.
      integer, intent(in) :: first_level !! First zero-based detail level thresholded.
      real(dp) :: value
      type(rss_result_t) :: trial
      trial = rsswav(noisy, threshold, filter_number, family, threshold_type, first_level)
      if (trial%ok) then
         value = trial%ssq
      else
         value = huge(1.0_dp)
      end if
   end function rss_value

end module wavethresh_cv

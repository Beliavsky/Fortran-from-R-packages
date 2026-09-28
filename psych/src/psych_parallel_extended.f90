module psych_parallel_extended
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use r_linalg, only : symmetric_eigen
   use psych_correlations, only : squared_multiple_correlations
   use psych_factor_extended, only : minimum_residual_factor
   use psych_stats, only : correlation_matrix, finite_correlation, finite_mean, finite_quantile
   implicit none
   private
   public :: parallel_analysis_normal, random_effect_factor

contains

   subroutine parallel_analysis_normal(r, n_observations, n_iter, quantile, observed_pc, reference_pc_mean, &
      reference_pc_quantile, observed_fa, reference_fa_mean, reference_fa_quantile, n_components, n_factors, &
      info, smc)
      real(dp), intent(in) :: r(:, :) !! Observed item correlation matrix used for Horn-style parallel analysis.
      integer, intent(in) :: n_observations !! Number of simulated observations per normal-reference replicate; must exceed one.
      integer, intent(in) :: n_iter !! Number of independent normal-reference replicates; must be positive.
      real(dp), intent(in) :: quantile !! Upper reference quantile in `[0,1]`, commonly 0.95.
      real(dp), allocatable, intent(out) :: observed_pc(:) !! Eigenvalues of the observed correlation matrix.
      real(dp), allocatable, intent(out) :: reference_pc_mean(:) !! Mean simulated principal-component eigenvalues.
      real(dp), allocatable, intent(out) :: reference_pc_quantile(:) !! Requested simulated PC eigenvalue quantiles.
      real(dp), allocatable, intent(out) :: observed_fa(:) !! Observed factor-reference eigenvalues, optionally with SMC diagonal.
      real(dp), allocatable, intent(out) :: reference_fa_mean(:) !! Mean simulated factor-reference eigenvalues.
      real(dp), allocatable, intent(out) :: reference_fa_quantile(:) !! Requested simulated factor-reference quantiles.
      integer, intent(out) :: n_components !! Leading observed PC eigenvalues exceeding their reference quantiles.
      integer, intent(out) :: n_factors !! Leading observed factor-reference eigenvalues exceeding their reference quantiles.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions, SMC, or eigensolution failure.
      logical, intent(in), optional :: smc !! Replace factor-analysis diagonals by squared multiple correlations; default false.
      real(dp), allocatable :: values(:), vectors(:, :), simulated(:, :), simulated_r(:, :), work(:, :), smc_values(:)
      real(dp), allocatable :: pc_replicates(:, :), fa_replicates(:, :), column(:)
      integer :: i, iter, local_info, p
      logical :: use_smc

      p = size(r, 1)
      use_smc = .false.
      if (present(smc)) use_smc = smc
      if (p < 1 .or. size(r, 2) /= p .or. n_observations < 2 .or. n_iter < 1) then
         info = 1
         n_components = 0
         n_factors = 0
         allocate(observed_pc(0), reference_pc_mean(0), reference_pc_quantile(0))
         allocate(observed_fa(0), reference_fa_mean(0), reference_fa_quantile(0))
         return
      end if

      call symmetric_eigen(r, values, vectors, local_info, descending=.true.)
      if (local_info /= 0) then
         info = 10 + local_info
         n_components = 0
         n_factors = 0
         allocate(observed_pc(0), reference_pc_mean(0), reference_pc_quantile(0))
         allocate(observed_fa(0), reference_fa_mean(0), reference_fa_quantile(0))
         return
      end if
      allocate(observed_pc(p))
      observed_pc = values
      allocate(work(p, p), smc_values(p))
      work = r
      if (use_smc) then
         call squared_multiple_correlations(r, smc_values, local_info)
         if (local_info /= 0) then
            info = 20 + local_info
            n_components = 0
            n_factors = 0
            allocate(observed_fa(0), reference_pc_mean(0), reference_pc_quantile(0))
            allocate(reference_fa_mean(0), reference_fa_quantile(0))
            return
         end if
         do i = 1, p
            work(i, i) = smc_values(i)
         end do
      end if
      call symmetric_eigen(work, values, vectors, local_info, descending=.true.)
      if (local_info /= 0) then
         info = 30 + local_info
         n_components = 0
         n_factors = 0
         allocate(observed_fa(0), reference_pc_mean(0), reference_pc_quantile(0))
         allocate(reference_fa_mean(0), reference_fa_quantile(0))
         return
      end if
      allocate(observed_fa(p))
      observed_fa = values

      allocate(simulated(n_observations, p), simulated_r(p, p), pc_replicates(n_iter, p), fa_replicates(n_iter, p))
      do iter = 1, n_iter
         call fill_standard_normal(simulated)
         call correlation_matrix(simulated, simulated_r)
         call symmetric_eigen(simulated_r, values, vectors, local_info, descending=.true.)
         if (local_info /= 0) then
            info = 40 + local_info
            return
         end if
         pc_replicates(iter, :) = values
         work = simulated_r
         if (use_smc) then
            call squared_multiple_correlations(simulated_r, smc_values, local_info)
            if (local_info /= 0) then
               info = 50 + local_info
               return
            end if
            do i = 1, p
               work(i, i) = smc_values(i)
            end do
         end if
         call symmetric_eigen(work, values, vectors, local_info, descending=.true.)
         if (local_info /= 0) then
            info = 60 + local_info
            return
         end if
         fa_replicates(iter, :) = values
      end do

      allocate(reference_pc_mean(p), reference_pc_quantile(p), reference_fa_mean(p), reference_fa_quantile(p))
      allocate(column(n_iter))
      do i = 1, p
         reference_pc_mean(i) = sum(pc_replicates(:, i)) / real(n_iter, dp)
         reference_fa_mean(i) = sum(fa_replicates(:, i)) / real(n_iter, dp)
         column = pc_replicates(:, i)
         reference_pc_quantile(i) = finite_quantile(column, min(1.0_dp, max(0.0_dp, quantile)))
         column = fa_replicates(:, i)
         reference_fa_quantile(i) = finite_quantile(column, min(1.0_dp, max(0.0_dp, quantile)))
      end do
      n_components = leading_above(observed_pc, reference_pc_quantile)
      n_factors = leading_above(observed_fa, reference_fa_quantile)
      info = 0
   end subroutine parallel_analysis_normal

   subroutine random_effect_factor(data, n_factors, loadings, uniqueness, subject_means, within_correlations, &
      objective, info, add_fix_noise)
      real(dp), intent(in) :: data(:, :) !! Complete subject-by-item data matrix for within-subject factor analysis.
      integer, intent(in) :: n_factors !! Number of MINRES factors fitted after removing subject mean effects.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! MINRES loadings of subject-centered item data.
      real(dp), allocatable, intent(out) :: uniqueness(:) !! MINRES uniquenesses of the subject-centered item data.
      real(dp), allocatable, intent(out) :: subject_means(:) !! Mean item score for each subject before centering.
      real(dp), allocatable, intent(out) :: within_correlations(:) !! Correlation of each item with the subject mean score.
      real(dp), intent(out) :: objective !! MINRES off-diagonal residual objective for the centered correlation matrix.
      integer, intent(out) :: info !! Zero on success; nonzero for missing values or factor-extraction failure.
      logical, intent(in), optional :: add_fix_noise !! Add Gaussian noise with SD 0.03 after centering, matching R `fix=TRUE`.
      real(dp), allocatable :: adjusted(:, :), r(:, :), noise(:, :)
      integer :: i, j, k
      logical :: add_noise

      if (size(data, 1) < 2 .or. size(data, 2) < 2) then
         info = 1
         objective = 0.0_dp
         allocate(loadings(0, 0), uniqueness(0), subject_means(0), within_correlations(0))
         return
      end if
      do j = 1, size(data, 2)
         do i = 1, size(data, 1)
            if (ieee_is_nan(data(i, j))) then
               info = 2
               objective = 0.0_dp
               allocate(loadings(0, 0), uniqueness(0), subject_means(0), within_correlations(0))
               return
            end if
         end do
      end do
      allocate(subject_means(size(data, 1)), within_correlations(size(data, 2)))
      do i = 1, size(data, 1)
         subject_means(i) = finite_mean(data(i, :))
      end do
      do j = 1, size(data, 2)
         within_correlations(j) = finite_correlation(data(:, j), subject_means)
      end do
      allocate(adjusted(size(data, 1), size(data, 2)))
      do j = 1, size(data, 2)
         adjusted(:, j) = data(:, j) - subject_means
      end do
      add_noise = .true.
      if (present(add_fix_noise)) add_noise = add_fix_noise
      if (add_noise) then
         allocate(noise(size(data, 1), size(data, 2)))
         call fill_standard_normal(noise)
         adjusted = adjusted + 0.03_dp * noise
      end if
      allocate(r(size(data, 2), size(data, 2)))
      call correlation_matrix(adjusted, r)
      k = min(max(1, n_factors), size(data, 2))
      allocate(loadings(size(data, 2), k), uniqueness(size(data, 2)))
      call minimum_residual_factor(r, k, loadings, uniqueness, objective, info)
   end subroutine random_effect_factor

   pure integer function leading_above(observed, reference) result(count_above)
      real(dp), intent(in) :: observed(:) !! Ordered observed eigenvalues.
      real(dp), intent(in) :: reference(:) !! Corresponding ordered reference thresholds.
      integer :: i

      count_above = 0
      do i = 1, min(size(observed), size(reference))
         if (observed(i) > reference(i)) then
            count_above = count_above + 1
         else
            exit
         end if
      end do
   end function leading_above

   subroutine fill_standard_normal(x)
      real(dp), intent(out) :: x(:, :) !! Matrix filled with independent standard-normal pseudorandom values.
      real(dp) :: angle, radius, u1, u2
      integer :: i, j

      do j = 1, size(x, 2)
         i = 1
         do while (i <= size(x, 1))
            call random_number(u1)
            call random_number(u2)
            u1 = max(u1, tiny(1.0_dp))
            radius = sqrt(-2.0_dp * log(u1))
            angle = 2.0_dp * acos(-1.0_dp) * u2
            x(i, j) = radius * cos(angle)
            if (i + 1 <= size(x, 1)) x(i + 1, j) = radius * sin(angle)
            i = i + 2
         end do
      end do
   end subroutine fill_standard_normal

end module psych_parallel_extended

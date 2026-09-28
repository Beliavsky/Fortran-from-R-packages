module psych_simulation
   use r_kinds, only : dp, r_pi
   use r_linalg, only : symmetric_eigen
   use psych_irt, only : irt_probability
   use psych_distributions, only : normal_cdf
   use psych_stats, only : correlation_matrix
   implicit none
   private
   public :: simulate_correlated_normal, simplex_correlation, congeneric_correlation
   public :: structure_correlation, simulate_rasch, simulate_irt_2pl, simulate_dichotomous
   public :: make_congeneric_model, make_hierarchical_model, simulate_item_dichot
   public :: simulate_npl_items, simulate_npn_items, simulate_polychoric_marginals
contains

   subroutine make_congeneric_model(loadings, model, error_sd, pattern, n, observed, latent)
      real(dp), intent(in), optional :: loadings(:) !! Common-factor loadings; defaults to `(0.8, 0.7, 0.6, 0.5)`.
      real(dp), allocatable, intent(out) :: model(:, :) !! Population covariance/correlation matrix implied by the pattern.
      real(dp), intent(in), optional :: error_sd(:) !! Item-specific residual standard deviations; defaults to `sqrt(1-load^2)`.
      real(dp), allocatable, intent(out), optional :: pattern(:, :) !! Item-by-latent pattern matrix of common and error factors.
      integer, intent(in), optional :: n !! Number of simulated cases when `observed` or `latent` is requested; default 1000.
      real(dp), allocatable, intent(out), optional :: observed(:, :) !! Optional simulated case-by-item observations.
      real(dp), allocatable, intent(out), optional :: latent(:, :) !! Optional simulated common-factor and residual scores.
      real(dp), allocatable :: load(:), local_pattern(:, :), z(:, :)
      integer :: i, p, sample_size

      if (present(loadings)) then
         allocate(load(size(loadings)))
         load = loadings
      else
         allocate(load(4))
         load = [0.8_dp, 0.7_dp, 0.6_dp, 0.5_dp]
      end if
      p = size(load)
      allocate(local_pattern(p, p + 1))
      local_pattern = 0.0_dp
      local_pattern(:, 1) = load
      do i = 1, p
         if (present(error_sd)) then
            local_pattern(i, i + 1) = error_sd(i)
         else
            local_pattern(i, i + 1) = sqrt(max(0.0_dp, 1.0_dp - load(i) * load(i)))
         end if
      end do
      model = matmul(local_pattern, transpose(local_pattern))
      if (present(pattern)) pattern = local_pattern

      if (present(observed) .or. present(latent)) then
         sample_size = 1000
         if (present(n)) sample_size = max(0, n)
         allocate(z(sample_size, p + 1))
         call fill_standard_normal(z)
         if (present(latent)) latent = z
         if (present(observed)) observed = matmul(z, transpose(local_pattern))
      end if
   end subroutine make_congeneric_model

   subroutine make_hierarchical_model(result, gload, fload, n, raw, info)
      real(dp), allocatable, intent(out) :: result(:, :) !! Model matrix, simulated raw data, or simulated sample correlation.
      real(dp), intent(in), optional :: gload(:, :) !! Higher-order loadings; defaults to the upstream 3-by-1 vector.
      real(dp), intent(in), optional :: fload(:, :) !! Item-by-group-factor loadings; defaults to the upstream 9-by-3 matrix.
      integer, intent(in), optional :: n !! Simulated case count; zero or absent returns the population model.
      logical, intent(in), optional :: raw !! With positive `n`, return raw cases when true; otherwise return correlations.
      integer, intent(out), optional :: info !! Zero on success or the eigensolver status from simulation.
      real(dp), allocatable :: g(:, :), f(:, :), phi(:, :), model(:, :), sample(:, :), sample_cor(:, :)
      integer :: i, local_info, n_cases, p
      logical :: return_raw

      if (present(gload)) then
         g = gload
      else
         allocate(g(3, 1))
         g(:, 1) = [0.9_dp, 0.8_dp, 0.7_dp]
      end if
      if (present(fload)) then
         f = fload
      else
         allocate(f(9, 3))
         f = reshape([0.8_dp, 0.7_dp, 0.6_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.7_dp, 0.6_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.6_dp, 0.5_dp, 0.4_dp], [9, 3])
      end if

      phi = matmul(g, transpose(g))
      do i = 1, size(phi, 1)
         phi(i, i) = 1.0_dp
      end do
      model = matmul(matmul(f, phi), transpose(f))
      do i = 1, size(model, 1)
         model(i, i) = 1.0_dp
      end do

      n_cases = 0
      if (present(n)) n_cases = max(0, n)
      if (n_cases == 0) then
         result = model
         if (present(info)) info = 0
         return
      end if

      call simulate_correlated_normal(n_cases, model, sample, local_info)
      if (present(info)) info = local_info
      if (local_info /= 0) then
         allocate(result(0, 0))
         return
      end if
      return_raw = .false.
      if (present(raw)) return_raw = raw
      if (return_raw) then
         result = sample
      else
         p = size(sample, 2)
         allocate(sample_cor(p, p))
         call correlation_matrix(sample, sample_cor)
         result = sample_cor
      end if
   end subroutine make_hierarchical_model

   subroutine simulate_item_dichot(nvar, nsub, responses, circum, xloading, yloading, gloading, xbias, ybias, low, high)
      integer, intent(in) :: nvar !! Number of dichotomous items; the non-circumplex layout is intended for multiples of four.
      integer, intent(in) :: nsub !! Number of simulated subjects.
      real(dp), allocatable, intent(out) :: responses(:, :) !! Simulated zero/one subject-by-item response matrix.
      logical, intent(in), optional :: circum !! Use evenly spaced circumplex angles instead of four quarter-circle directions.
      real(dp), intent(in), optional :: xloading !! Loading multiplier for the first latent coordinate; default `0.6`.
      real(dp), intent(in), optional :: yloading !! Loading multiplier for the second latent coordinate; default `0.6`.
      real(dp), intent(in), optional :: gloading !! General-factor loading common to every item; default zero.
      real(dp), intent(in), optional :: xbias !! Additive bias for the first latent coordinate; default zero.
      real(dp), intent(in), optional :: ybias !! Additive bias for the second latent coordinate; default zero.
      real(dp), intent(in), optional :: low !! Threshold applied to odd source items before column regrouping; default zero.
      real(dp), intent(in), optional :: high !! Threshold applied to even source items before column regrouping; default zero.
      real(dp), allocatable :: error(:, :), g(:), item(:, :), rad(:), truex(:), truey(:)
      real(dp) :: avloading, errorweight, gl, hi, lo, xl, xb, yl, yb
      integer :: i, j, nvar2, quarter
      logical :: use_circum

      xl = 0.6_dp
      yl = 0.6_dp
      gl = 0.0_dp
      xb = 0.0_dp
      yb = 0.0_dp
      lo = 0.0_dp
      hi = 0.0_dp
      use_circum = .false.
      if (present(xloading)) xl = xloading
      if (present(yloading)) yl = yloading
      if (present(gloading)) gl = gloading
      if (present(xbias)) xb = xbias
      if (present(ybias)) yb = ybias
      if (present(low)) lo = low
      if (present(high)) hi = high
      if (present(circum)) use_circum = circum

      allocate(g(nsub), truex(nsub), truey(nsub), error(nsub, nvar), item(nsub, nvar), rad(nvar))
      call fill_standard_normal_vector(g)
      call fill_standard_normal_vector(truex)
      call fill_standard_normal_vector(truey)
      call fill_standard_normal(error)
      truex = truex * xl + xb
      truey = truey * yl + yb
      avloading = 0.5_dp * (xl + yl)
      errorweight = sqrt(max(0.0_dp, 1.0_dp - avloading * avloading - gl * gl))

      if (use_circum) then
         do j = 1, nvar
            rad(j) = 2.0_dp * r_pi * real(j - 1, dp) / real(max(1, nvar), dp)
         end do
      else
         quarter = max(1, nvar / 4)
         do j = 1, nvar
            if (j <= quarter) then
               rad(j) = 0.0_dp
            else if (j <= 2 * quarter) then
               rad(j) = 0.5_dp * r_pi
            else if (j <= 3 * quarter) then
               rad(j) = r_pi
            else
               rad(j) = 1.5_dp * r_pi
            end if
         end do
      end if

      do j = 1, nvar
         do i = 1, nsub
            item(i, j) = gl * g(i) + truex(i) * cos(rad(j)) + truey(i) * sin(rad(j)) + errorweight * error(i, j)
         end do
      end do

      nvar2 = nvar / 2
      allocate(responses(nsub, 2 * nvar2))
      do j = 1, nvar2
         responses(:, j) = merge(1.0_dp, 0.0_dp, item(:, 2 * j - 1) >= lo)
         responses(:, nvar2 + j) = merge(1.0_dp, 0.0_dp, item(:, 2 * j) >= hi)
      end do
   end subroutine simulate_item_dichot

   subroutine simulate_correlated_normal(n, correlation, x, info)
      integer, intent(in) :: n !! Number of multivariate observations to generate.
      real(dp), intent(in) :: correlation(:, :) !! Target symmetric positive-semidefinite correlation/covariance matrix.
      real(dp), allocatable, intent(out) :: x(:, :) !! Generated observation-by-variable normal matrix.
      integer, intent(out) :: info !! Zero on successful eigendecomposition.
      real(dp), allocatable :: values(:), vectors(:, :), root(:, :), z(:, :)
      integer :: j, p
      p = size(correlation, 1)
      call symmetric_eigen(correlation, values, vectors, info, descending=.true.)
      if (info /= 0) then
         allocate(x(0, 0))
         return
      end if
      allocate(root(p, p), z(n, p), x(n, p))
      root = vectors
      do j = 1, p
         root(:, j) = root(:, j) * sqrt(max(values(j), 0.0_dp))
      end do
      call fill_standard_normal(z)
      x = matmul(z, transpose(root))
   end subroutine simulate_correlated_normal

   pure subroutine simplex_correlation(p, rho, correlation)
      integer, intent(in) :: p !! Number of variables in the simplex/AR(1)-style correlation matrix.
      real(dp), intent(in) :: rho !! Adjacent-variable correlation parameter.
      real(dp), allocatable, intent(out) :: correlation(:, :) !! Matrix with entry rho**abs(i-j).
      integer :: i, j
      allocate(correlation(p, p))
      do j = 1, p
         do i = 1, p
            correlation(i, j) = rho ** abs(i - j)
         end do
      end do
   end subroutine simplex_correlation

   pure subroutine congeneric_correlation(loadings, correlation)
      real(dp), intent(in) :: loadings(:) !! Standardized loadings of variables on one common factor.
      real(dp), allocatable, intent(out) :: correlation(:, :) !! Congeneric population correlation matrix with unit diagonal.
      integer :: i, j, p
      p = size(loadings)
      allocate(correlation(p, p))
      do j = 1, p
         do i = 1, p
            correlation(i, j) = loadings(i) * loadings(j)
         end do
         correlation(j, j) = 1.0_dp
      end do
   end subroutine congeneric_correlation

   pure subroutine structure_correlation(loadings, phi, correlation)
      real(dp), intent(in) :: loadings(:, :) !! Variable-by-factor loading matrix.
      real(dp), intent(in) :: phi(:, :) !! Factor correlation matrix.
      real(dp), allocatable, intent(out) :: correlation(:, :) !! Model-implied standardized variable correlation matrix.
      integer :: j
      correlation = matmul(matmul(loadings, phi), transpose(loadings))
      do j = 1, size(correlation, 1)
         correlation(j, j) = 1.0_dp
      end do
   end subroutine structure_correlation

   subroutine simulate_rasch(n, difficulty, responses)
      integer, intent(in) :: n !! Number of subjects to simulate.
      real(dp), intent(in) :: difficulty(:) !! Item difficulty parameters.
      real(dp), allocatable, intent(out) :: responses(:, :) !! Simulated zero/one subject-by-item responses.
      real(dp), allocatable :: theta(:)
      real(dp) :: u
      integer :: i, j
      allocate(theta(n), responses(n, size(difficulty)))
      call fill_standard_normal_vector(theta)
      do j = 1, size(difficulty)
         do i = 1, n
            call random_number(u)
            responses(i, j) = merge(1.0_dp, 0.0_dp, u < irt_probability(theta(i), difficulty(j), 1.0_dp))
         end do
      end do
   end subroutine simulate_rasch

   subroutine simulate_irt_2pl(n, difficulty, discrimination, responses)
      integer, intent(in) :: n !! Number of subjects to simulate.
      real(dp), intent(in) :: difficulty(:) !! Item difficulty parameters.
      real(dp), intent(in) :: discrimination(:) !! Item discrimination parameters.
      real(dp), allocatable, intent(out) :: responses(:, :) !! Simulated zero/one subject-by-item responses.
      real(dp), allocatable :: theta(:)
      real(dp) :: u
      integer :: i, j
      allocate(theta(n), responses(n, size(difficulty)))
      call fill_standard_normal_vector(theta)
      do j = 1, size(difficulty)
         do i = 1, n
            call random_number(u)
            responses(i, j) = merge(1.0_dp, 0.0_dp, &
               u < irt_probability(theta(i), difficulty(j), discrimination(j)))
         end do
      end do
   end subroutine simulate_irt_2pl

   subroutine simulate_dichotomous(probability, n, responses)
      real(dp), intent(in) :: probability(:) !! Marginal success probability for each simulated item.
      integer, intent(in) :: n !! Number of independent response rows to generate.
      real(dp), allocatable, intent(out) :: responses(:, :) !! Independent Bernoulli response matrix.
      real(dp) :: u
      integer :: i, j
      allocate(responses(n, size(probability)))
      do j = 1, size(probability)
         do i = 1, n
            call random_number(u)
            responses(i, j) = merge(1.0_dp, 0.0_dp, u < probability(j))
         end do
      end do
   end subroutine simulate_dichotomous

   subroutine simulate_npl_items(n, difficulty, discrimination, items, theta, gamma, zeta, mu, sd, theta_input, trials, info)
      integer, intent(in) :: n !! Number of subjects to simulate when no explicit `theta_input` vector is supplied.
      real(dp), intent(in) :: difficulty(:) !! Per-item logistic difficulty parameters.
      real(dp), intent(in) :: discrimination(:) !! Per-item logistic discrimination parameters; length must match difficulty.
      real(dp), allocatable, intent(out) :: items(:, :) !! Simulated subject-by-item binomial response counts.
      real(dp), allocatable, intent(out) :: theta(:) !! Latent trait values used to generate the responses.
      real(dp), intent(in), optional :: gamma(:) !! Lower asymptotes, either length one or one value per item; default zero.
      real(dp), intent(in), optional :: zeta(:) !! Upper asymptotes, either length one or one value per item; default one.
      real(dp), intent(in), optional :: mu !! Mean of generated normal traits when `theta_input` is absent; default zero.
      real(dp), intent(in), optional :: sd !! Standard deviation of generated traits when `theta_input` is absent; default one.
      real(dp), intent(in), optional :: theta_input(:) !! Optional explicit latent trait vector replacing random trait generation.
      integer, intent(in), optional :: trials !! Binomial trials per item; default one for dichotomous responses.
      integer, intent(out), optional :: info !! Zero on success; nonzero for incompatible parameter dimensions.

      call simulate_irt_link(.true., n, difficulty, discrimination, items, theta, gamma, zeta, mu, sd, &
         theta_input, trials, info)
   end subroutine simulate_npl_items

   subroutine simulate_npn_items(n, difficulty, discrimination, items, theta, gamma, zeta, mu, sd, theta_input, trials, info)
      integer, intent(in) :: n !! Number of subjects to simulate when no explicit `theta_input` vector is supplied.
      real(dp), intent(in) :: difficulty(:) !! Per-item normal-ogive difficulty parameters.
      real(dp), intent(in) :: discrimination(:) !! Per-item normal-ogive discrimination parameters; length must match difficulty.
      real(dp), allocatable, intent(out) :: items(:, :) !! Simulated subject-by-item binomial response counts.
      real(dp), allocatable, intent(out) :: theta(:) !! Latent trait values used to generate the responses.
      real(dp), intent(in), optional :: gamma(:) !! Lower asymptotes, either length one or one value per item; default zero.
      real(dp), intent(in), optional :: zeta(:) !! Upper asymptotes, either length one or one value per item; default one.
      real(dp), intent(in), optional :: mu !! Mean of generated normal traits when `theta_input` is absent; default zero.
      real(dp), intent(in), optional :: sd !! Standard deviation of generated traits when `theta_input` is absent; default one.
      real(dp), intent(in), optional :: theta_input(:) !! Optional explicit latent trait vector replacing random trait generation.
      integer, intent(in), optional :: trials !! Binomial trials per item; default one for dichotomous responses.
      integer, intent(out), optional :: info !! Zero on success; nonzero for incompatible parameter dimensions.

      call simulate_irt_link(.false., n, difficulty, discrimination, items, theta, gamma, zeta, mu, sd, &
         theta_input, trials, info)
   end subroutine simulate_npn_items

   subroutine simulate_polychoric_marginals(correlation, thresholds, n, responses, info)
      real(dp), intent(in) :: correlation(:, :) !! Target latent correlation matrix for the simulated normal variables.
      real(dp), intent(in) :: thresholds(:, :) !! Item-by-threshold latent cut points defining ordinal categories.
      integer, intent(in) :: n !! Number of response rows to simulate.
      real(dp), allocatable, intent(out) :: responses(:, :) !! Simulated zero-based ordinal subject-by-item responses.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid dimensions or eigensolver failure.
      real(dp), allocatable :: latent(:, :)
      integer :: i, j, k, local_info

      local_info = 0
      if (n < 0 .or. size(correlation, 1) /= size(correlation, 2) .or. &
         size(thresholds, 1) /= size(correlation, 1)) then
         allocate(responses(0, 0))
         local_info = 1
         if (present(info)) info = local_info
         return
      end if
      call simulate_correlated_normal(n, correlation, latent, local_info)
      if (local_info /= 0) then
         allocate(responses(0, 0))
         if (present(info)) info = local_info
         return
      end if
      allocate(responses(n, size(correlation, 1)))
      responses = 0.0_dp
      do j = 1, size(correlation, 1)
         do i = 1, n
            do k = 1, size(thresholds, 2)
               if (latent(i, j) > thresholds(j, k)) responses(i, j) = real(k, dp)
            end do
         end do
      end do
      if (present(info)) info = local_info
   end subroutine simulate_polychoric_marginals

   subroutine simulate_irt_link(logistic_link, n, difficulty, discrimination, items, theta, gamma, zeta, mu, sd, &
      theta_input, trials, info)
      logical, intent(in) :: logistic_link !! Use logistic probabilities when true and normal-ogive probabilities when false.
      integer, intent(in) :: n !! Number of generated subjects when explicit traits are not supplied.
      real(dp), intent(in) :: difficulty(:) !! Per-item difficulty parameters.
      real(dp), intent(in) :: discrimination(:) !! Per-item discrimination parameters.
      real(dp), allocatable, intent(out) :: items(:, :) !! Simulated subject-by-item binomial response counts.
      real(dp), allocatable, intent(out) :: theta(:) !! Latent traits used for response generation.
      real(dp), intent(in), optional :: gamma(:) !! Lower asymptotes with length one or number of items; default zero.
      real(dp), intent(in), optional :: zeta(:) !! Upper asymptotes with length one or number of items; default one.
      real(dp), intent(in), optional :: mu !! Mean for generated normal traits; default zero.
      real(dp), intent(in), optional :: sd !! Standard deviation for generated normal traits; default one.
      real(dp), intent(in), optional :: theta_input(:) !! Explicit latent traits; when present their length defines subjects.
      integer, intent(in), optional :: trials !! Number of independent Bernoulli trials accumulated per cell; default one.
      integer, intent(out), optional :: info !! Zero on success; nonzero for parameter dimension errors.
      real(dp) :: g, location, p, s, u, upper_value
      integer :: i, j, local_info, n_subjects, n_trials, t

      local_info = 0
      if (size(discrimination) /= size(difficulty) .or. size(difficulty) < 1) local_info = 1
      if (present(gamma)) then
         if (size(gamma) /= 1 .and. size(gamma) /= size(difficulty)) local_info = 2
      end if
      if (present(zeta)) then
         if (size(zeta) /= 1 .and. size(zeta) /= size(difficulty)) local_info = 3
      end if
      if (present(theta_input)) then
         n_subjects = size(theta_input)
      else
         n_subjects = n
      end if
      if (n_subjects < 0) local_info = 4
      n_trials = 1
      if (present(trials)) n_trials = trials
      if (n_trials < 1) local_info = 5
      if (local_info /= 0) then
         allocate(items(0, 0), theta(0))
         if (present(info)) info = local_info
         return
      end if

      allocate(theta(n_subjects), items(n_subjects, size(difficulty)))
      if (present(theta_input)) then
         theta = theta_input
      else
         call fill_standard_normal_vector(theta)
         location = 0.0_dp
         if (present(mu)) location = mu
         s = 1.0_dp
         if (present(sd)) s = sd
         theta = location + s * theta
      end if

      do j = 1, size(difficulty)
         g = 0.0_dp
         if (present(gamma)) g = gamma(min(j, size(gamma)))
         upper_value = 1.0_dp
         if (present(zeta)) upper_value = zeta(min(j, size(zeta)))
         do i = 1, n_subjects
            if (logistic_link) then
               p = irt_probability(theta(i), difficulty(j), discrimination(j))
            else
               p = normal_cdf(discrimination(j) * (theta(i) - difficulty(j)))
            end if
            p = max(0.0_dp, min(1.0_dp, g + (upper_value - g) * p))
            items(i, j) = 0.0_dp
            do t = 1, n_trials
               call random_number(u)
               if (u < p) items(i, j) = items(i, j) + 1.0_dp
            end do
         end do
      end do
      if (present(info)) info = local_info
   end subroutine simulate_irt_link

   subroutine fill_standard_normal(x)
      real(dp), intent(out) :: x(:, :) !! Matrix filled with independent standard-normal draws.
      real(dp) :: u1, u2, radius, angle
      integer :: i, j, k, total
      total = size(x)
      k = 0
      do while (k < total)
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1, tiny(1.0_dp))
         radius = sqrt(-2.0_dp * log(u1))
         angle = 2.0_dp * r_pi * u2
         k = k + 1
         call flat_assign(x, k, radius * cos(angle))
         if (k < total) then
            k = k + 1
            call flat_assign(x, k, radius * sin(angle))
         end if
      end do
      i = size(x, 1)
      j = size(x, 2)
      if (i + j < 0) x = 0.0_dp
   end subroutine fill_standard_normal

   subroutine fill_standard_normal_vector(x)
      real(dp), intent(out) :: x(:) !! Vector filled with independent standard-normal draws.
      real(dp) :: u1, u2, radius, angle
      integer :: k
      k = 1
      do while (k <= size(x))
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1, tiny(1.0_dp))
         radius = sqrt(-2.0_dp * log(u1))
         angle = 2.0_dp * r_pi * u2
         x(k) = radius * cos(angle)
         if (k + 1 <= size(x)) x(k + 1) = radius * sin(angle)
         k = k + 2
      end do
   end subroutine fill_standard_normal_vector

   subroutine flat_assign(x, index, value)
      real(dp), intent(inout) :: x(:, :) !! Matrix receiving one value in Fortran column-major order.
      integer, intent(in) :: index !! One-based flattened column-major index.
      real(dp), intent(in) :: value !! Value assigned at the flattened position.
      integer :: row, column
      row = mod(index - 1, size(x, 1)) + 1
      column = (index - 1) / size(x, 1) + 1
      x(row, column) = value
   end subroutine flat_assign
end module psych_simulation

module funData_simulation
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use funData_kinds, only : dp
   use funData_types, only : basis_spec, fun_data, int_vector, multi_fun_data, real_vector, sim_fun_result, sim_multi_result
   use funData_core, only : create_fun_data_1d, fun_value, nobs_fun_data, support_dim, support_size
   use funData_numeric, only : tensor_product2, tensor_product3
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)

   public :: add_error_fun_data
   public :: add_error_multi_fun_data
   public :: eigenfunctions
   public :: eigenvalues
   public :: sim_fun_data
   public :: sim_multi_fun_data
   public :: sparsify_fun_data
   public :: sparsify_multi_fun_data

contains

   pure subroutine eigenvalues(m, value_type, values, ok)
      integer, intent(in) :: m !! Number of decreasing eigenvalues requested; must be positive.
      character(len=*), intent(in) :: value_type !! Eigenvalue family: linear, exponential, or wiener.
      real(dp), allocatable, intent(out) :: values(:) !! Generated eigenvalues in decreasing order.
      logical, intent(out) :: ok !! True when m is positive and the eigenvalue family is recognized.
      integer :: j
      character(len=:), allocatable :: kind

      ok = .false.
      if (m <= 0) return
      allocate(values(m))
      kind = trim(to_lower(value_type))
      select case (kind)
      case ("linear")
         do j = 1, m
            values(j) = real(m + 1 - j, dp) / real(m, dp)
         end do
      case ("exponential")
         do j = 1, m
            values(j) = exp(-real(j - 1, dp) / 2.0_dp)
         end do
      case ("wiener")
         do j = 1, m
            values(j) = 1.0_dp / (0.5_dp * pi * real(2 * j - 1, dp))**2
         end do
      case default
         return
      end select
      ok = .true.
   end subroutine eigenvalues

   pure recursive subroutine eigenfunctions(argvals, m, basis_type, result, ignore_deg, ok)
      real(dp), intent(in) :: argvals(:) !! Fine one-dimensional grid on the interval supporting the basis functions.
      integer, intent(in) :: m !! Number of basis functions returned after any PolyHigh exclusions.
      character(len=*), intent(in) :: basis_type !! Basis family: Poly, PolyHigh, Fourier, FourierLin, or Wiener.
      type(fun_data), intent(out) :: result !! One-dimensional functional-data object whose observations are basis functions.
      integer, intent(in), optional :: ignore_deg(:) !! One-based polynomial-function indices excluded for PolyHigh.
      logical, intent(out) :: ok !! True when the requested basis family and parameters are valid.
      type(fun_data) :: full_poly
      real(dp), allocatable :: phi(:, :)
      integer, allocatable :: keep(:)
      integer :: j
      integer :: k
      integer :: total_m
      real(dp) :: xmin
      real(dp) :: span
      real(dp) :: z
      real(dp) :: norm_const
      character(len=:), allocatable :: kind
      logical :: base_ok

      ok = .false.
      if (size(argvals) == 0 .or. m <= 0) return
      xmin = minval(argvals)
      span = maxval(argvals) - xmin
      if (span <= 0.0_dp) return
      kind = trim(to_lower(basis_type))

      if (kind == "polyhigh") then
         if (.not. present(ignore_deg)) return
         if (size(ignore_deg) == 0 .or. any(ignore_deg <= 0)) return
         total_m = m + size(ignore_deg)
         call eigenfunctions(argvals, total_m, "poly", full_poly, ok = base_ok)
         if (.not. base_ok) return
         allocate(keep(m))
         k = 0
         do j = 1, total_m
            if (.not. any(ignore_deg == j)) then
               k = k + 1
               if (k <= m) keep(k) = j
            end if
         end do
         if (k /= m) return
         allocate(phi(m, size(argvals)))
         do j = 1, m
            do k = 1, size(argvals)
               phi(j, k) = fun_value(full_poly, keep(j), [k])
            end do
         end do
         call create_fun_data_1d(argvals, phi, result, ok)
         return
      end if

      allocate(phi(m, size(argvals)))
      select case (kind)
      case ("poly")
         phi(1, :) = 1.0_dp
         if (m >= 2) then
            do k = 1, size(argvals)
               phi(2, k) = 2.0_dp * (argvals(k) - xmin) / span - 1.0_dp
            end do
         end if
         do j = 3, m
            do k = 1, size(argvals)
               z = 2.0_dp * (argvals(k) - xmin) / span - 1.0_dp
               phi(j, k) = real(2 * j - 3, dp) / real(j - 1, dp) * z * phi(j - 1, k) - &
                  real(j - 2, dp) / real(j - 1, dp) * phi(j - 2, k)
            end do
         end do
         do j = 1, m
            phi(j, :) = sqrt(real(2 * j - 1, dp) / span) * phi(j, :)
         end do
      case ("fourier", "fourierlin")
         phi(1, :) = sqrt(1.0_dp / span)
         do j = 2, m
            if (mod(j, 2) == 0) then
               do k = 1, size(argvals)
                  phi(j, k) = sqrt(2.0_dp / span) * cos(real(j / 2, dp) * &
                     (2.0_dp * pi * (argvals(k) - xmin) / span - pi))
               end do
            else
               do k = 1, size(argvals)
                  phi(j, k) = sqrt(2.0_dp / span) * sin(real(j / 2, dp) * &
                     (2.0_dp * pi * (argvals(k) - xmin) / span - pi))
               end do
            end if
         end do
         if (kind == "fourierlin") then
            if (abs(xmin) > epsilon(1.0_dp) .or. abs(maxval(argvals) - 1.0_dp) > epsilon(1.0_dp)) return
            do k = 1, size(argvals)
               phi(m, k) = argvals(k) - 0.5_dp
               do j = 1, (m - 1) / 2
                  phi(m, k) = phi(m, k) + (-1.0_dp)**j / (pi * real(j, dp)) * &
                     sin(real(j, dp) * (2.0_dp * pi * argvals(k) - pi))
               end do
            end do
            norm_const = 1.0_dp / 3.0_dp - 1.0_dp / 4.0_dp
            do j = 1, (m - 1) / 2
               norm_const = norm_const - 1.0_dp / (2.0_dp * pi**2 * real(j * j, dp))
            end do
            if (norm_const <= 0.0_dp) return
            phi(m, :) = phi(m, :) / sqrt(norm_const)
         end if
      case ("wiener")
         do j = 1, m
            do k = 1, size(argvals)
               phi(j, k) = sqrt(2.0_dp / span) * sin(0.5_dp * pi * real(2 * j - 1, dp) * &
                  (argvals(k) - xmin) / span)
            end do
         end do
      case default
         return
      end select
      call create_fun_data_1d(argvals, phi, result, ok)
   end subroutine eigenfunctions

   subroutine sim_fun_data(argvals, m, basis_type, value_type, n, result, ignore_deg, ok)
      type(real_vector), intent(in) :: argvals(:) !! Sampling grid vectors, one per support dimension.
      integer, intent(in) :: m(:) !! Number of marginal basis functions; length one may be replicated to all dimensions.
      character(len=*), intent(in) :: basis_type(:) !! Marginal basis families; length one may be replicated to all dimensions.
      character(len=*), intent(in) :: value_type !! Eigenvalue family used as score variances.
      integer, intent(in) :: n !! Number of functional observations to simulate.
      type(sim_fun_result), intent(out) :: result !! Simulated data, true tensor-product basis, and true eigenvalues.
      type(int_vector), intent(in), optional :: ignore_deg(:) !! Optional PolyHigh exclusions by support dimension.
      logical, intent(out) :: ok !! True when basis construction and simulation succeed for one to three support dimensions.
      type(fun_data), allocatable :: marginal(:)
      type(fun_data) :: temp2
      integer, allocatable :: use_m(:)
      character(len=16), allocatable :: use_type(:)
      real(dp), allocatable :: scores(:, :)
      real(dp), allocatable :: matrix(:, :)
      integer :: d
      integer :: i
      integer :: j
      integer :: p
      integer :: m_total
      logical :: basis_ok
      logical :: tensor_ok

      ok = .false.
      d = size(argvals)
      if (d < 1 .or. d > 3 .or. n <= 0) return
      allocate(use_m(d))
      if (size(m) == 1) then
         use_m = m(1)
      else if (size(m) == d) then
         use_m = m
      else
         return
      end if
      if (any(use_m <= 0)) return
      allocate(use_type(d))
      if (size(basis_type) == 1) then
         use_type = basis_type(1)
      else if (size(basis_type) == d) then
         use_type = basis_type
      else
         return
      end if
      allocate(marginal(d))
      do j = 1, d
         if (present(ignore_deg)) then
            if (size(ignore_deg) /= d) return
            if (allocated(ignore_deg(j)%values)) then
               call eigenfunctions(argvals(j)%values, use_m(j), use_type(j), marginal(j), &
                  ignore_deg(j)%values, basis_ok)
            else
               call eigenfunctions(argvals(j)%values, use_m(j), use_type(j), marginal(j), ok = basis_ok)
            end if
         else
            call eigenfunctions(argvals(j)%values, use_m(j), use_type(j), marginal(j), ok = basis_ok)
         end if
         if (.not. basis_ok) return
      end do

      if (d == 1) then
         result%true_funs = marginal(1)
      else if (d == 2) then
         call tensor_product2(marginal(1), marginal(2), result%true_funs, tensor_ok)
         if (.not. tensor_ok) return
      else
         call tensor_product2(marginal(1), marginal(2), temp2, tensor_ok)
         if (.not. tensor_ok) return
         call tensor_basis_with_third(temp2, marginal(3), result%true_funs, tensor_ok)
         if (.not. tensor_ok) return
      end if

      m_total = nobs_fun_data(result%true_funs)
      call eigenvalues(m_total, value_type, result%true_vals, basis_ok)
      if (.not. basis_ok) return
      allocate(scores(n, m_total))
      do i = 1, n
         do j = 1, m_total
            scores(i, j) = normal_random() * sqrt(result%true_vals(j))
         end do
      end do
      allocate(matrix(n, support_size(result%true_funs)))
      matrix = 0.0_dp
      do p = 1, support_size(result%true_funs)
         do i = 1, n
            do j = 1, m_total
               matrix(i, p) = matrix(i, p) + scores(i, j) * result%true_funs%x(j + (p - 1) * m_total)
            end do
         end do
      end do
      call simulation_matrix_to_fun(result%true_funs, matrix, result%sim_data, basis_ok)
      if (.not. basis_ok) return
      ok = .true.
   end subroutine sim_fun_data

   subroutine sim_multi_fun_data(mode, specs, value_type, n, result, split_m, split_basis_type, split_ignore_deg, ok)
      character(len=*), intent(in) :: mode !! Multivariate basis construction mode: split or weighted.
      type(basis_spec), intent(in) :: specs(:) !! Component grids and weighted-mode marginal basis specifications.
      character(len=*), intent(in) :: value_type !! Eigenvalue family used for multivariate score variances.
      integer, intent(in) :: n !! Number of multivariate observations to simulate.
      type(sim_multi_result), intent(out) :: result !! Simulated multivariate data, true basis, and score variances.
      integer, intent(in), optional :: split_m !! Number of split-basis functions for split mode.
      character(len=*), intent(in), optional :: split_basis_type !! Basis family used on the concatenated split support.
      integer, intent(in), optional :: split_ignore_deg(:) !! Optional PolyHigh exclusions for split mode.
      logical, intent(out) :: ok !! True when the selected construction mode and all component specifications are valid.
      character(len=:), allocatable :: kind
      type(multi_fun_data) :: basis

      ok = .false.
      if (size(specs) < 1 .or. n <= 0) return
      kind = trim(to_lower(mode))
      select case (kind)
      case ("split")
         if (.not. present(split_m) .or. .not. present(split_basis_type)) return
         if (present(split_ignore_deg)) then
            call build_split_basis(specs, split_m, split_basis_type, basis, split_ignore_deg, ok)
         else
            call build_split_basis(specs, split_m, split_basis_type, basis, ok = ok)
         end if
      case ("weighted")
         call build_weighted_basis(specs, basis, ok)
      case default
         return
      end select
      if (.not. ok) return
      call simulate_multi_from_basis(basis, value_type, n, result, ok)
   end subroutine sim_multi_fun_data

   subroutine sparsify_fun_data(object, min_obs, max_obs, result, ok)
      type(fun_data), intent(in) :: object !! One-dimensional regular data to sparsify using NaNs for omitted points.
      integer, intent(in) :: min_obs !! Minimum number of retained grid points for each observation.
      integer, intent(in) :: max_obs !! Maximum number of retained grid points for each observation.
      type(fun_data), intent(out) :: result !! Sparse copy of object using IEEE NaNs for omitted values.
      logical, intent(out) :: ok !! True when one-dimensional support and observation-count limits are valid.
      integer, allocatable :: order(:)
      logical, allocatable :: keep(:)
      real(dp) :: qnan
      real(dp) :: u
      integer :: i
      integer :: j
      integer :: k
      integer :: m
      integer :: ni
      integer :: n

      ok = .false.
      if (support_dim(object) /= 1) return
      m = object%dims(2)
      if (min_obs < 1 .or. max_obs < min_obs .or. max_obs > m) return
      result = object
      qnan = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(order(m))
      allocate(keep(m))
      n = nobs_fun_data(object)
      do i = 1, n
         do j = 1, m
            order(j) = j
         end do
         if (min_obs == max_obs) then
            ni = min_obs
         else
            call random_number(u)
            ni = min_obs + int(u * real(max_obs - min_obs + 1, dp))
            ni = min(ni, max_obs)
         end if
         call shuffle_int(order)
         keep = .false.
         do k = 1, ni
            keep(order(k)) = .true.
         end do
         do j = 1, m
            if (.not. keep(j)) result%x(i + (j - 1) * n) = qnan
         end do
      end do
      ok = .true.
   end subroutine sparsify_fun_data

   subroutine sparsify_multi_fun_data(object, min_obs, max_obs, result, ok)
      type(multi_fun_data), intent(in) :: object !! Multivariate regular functional data to sparsify component by component.
      integer, intent(in) :: min_obs(:) !! Minimum retained grid-point counts, one per component.
      integer, intent(in) :: max_obs(:) !! Maximum retained grid-point counts, one per component.
      type(multi_fun_data), intent(out) :: result !! Sparse multivariate copy.
      logical, intent(out) :: ok !! True when limits are supplied for every component and all component calls succeed.
      logical :: component_ok
      integer :: j

      ok = .false.
      if (.not. allocated(object%components)) return
      if (size(min_obs) /= size(object%components) .or. size(max_obs) /= size(object%components)) return
      allocate(result%components(size(object%components)))
      do j = 1, size(object%components)
         call sparsify_fun_data(object%components(j), min_obs(j), max_obs(j), result%components(j), component_ok)
         if (.not. component_ok) return
      end do
      ok = .true.
   end subroutine sparsify_multi_fun_data

   subroutine add_error_fun_data(object, sd, result, ok)
      type(fun_data), intent(in) :: object !! Regular functional data to perturb by independent Gaussian white noise.
      real(dp), intent(in) :: sd !! Positive standard deviation of the additive Gaussian noise.
      type(fun_data), intent(out) :: result !! Noisy copy of object with original support metadata preserved.
      logical, intent(out) :: ok !! True when sd is positive.
      integer :: i

      ok = .false.
      if (sd <= 0.0_dp) return
      result = object
      do i = 1, size(result%x)
         result%x(i) = result%x(i) + sd * normal_random()
      end do
      ok = .true.
   end subroutine add_error_fun_data

   subroutine add_error_multi_fun_data(object, sd, result, ok)
      type(multi_fun_data), intent(in) :: object !! Multivariate regular functional data to perturb independently by component.
      real(dp), intent(in) :: sd(:) !! Positive noise standard deviations, either length one or one per component.
      type(multi_fun_data), intent(out) :: result !! Noisy multivariate copy.
      logical, intent(out) :: ok !! True when standard deviations are positive and dimensionally compatible.
      logical :: component_ok
      real(dp) :: use_sd
      integer :: j

      ok = .false.
      if (.not. allocated(object%components)) return
      if (size(sd) /= 1 .and. size(sd) /= size(object%components)) return
      if (any(sd <= 0.0_dp)) return
      allocate(result%components(size(object%components)))
      do j = 1, size(object%components)
         if (size(sd) == 1) then
            use_sd = sd(1)
         else
            use_sd = sd(j)
         end if
         call add_error_fun_data(object%components(j), use_sd, result%components(j), component_ok)
         if (.not. component_ok) return
      end do
      ok = .true.
   end subroutine add_error_multi_fun_data

   pure subroutine tensor_basis_with_third(first_two, third, result, ok)
      type(fun_data), intent(in) :: first_two !! Two-dimensional tensor basis whose observations index the first two marginal bases.
      type(fun_data), intent(in) :: third !! Third one-dimensional marginal basis.
      type(fun_data), intent(out) :: result !! Three-dimensional tensor basis with all observation combinations.
      logical, intent(out) :: ok !! True when dimensions have the expected 2D and 1D form.
      type(real_vector) :: grids(3)
      integer :: dims(4)
      real(dp), allocatable :: out_x(:)
      integer :: a
      integer :: b
      integer :: c
      integer :: i
      integer :: j
      integer :: obs
      integer :: n12
      integer :: n3
      integer :: nout

      ok = .false.
      if (support_dim(first_two) /= 2 .or. support_dim(third) /= 1) return
      n12 = nobs_fun_data(first_two)
      n3 = nobs_fun_data(third)
      nout = n12 * n3
      dims = [nout, first_two%dims(2), first_two%dims(3), third%dims(2)]
      grids(1)%values = first_two%argvals(1)%values
      grids(2)%values = first_two%argvals(2)%values
      grids(3)%values = third%argvals(1)%values
      allocate(out_x(product(dims)))
      do i = 1, n12
         do j = 1, n3
            obs = (i - 1) * n3 + j
            do c = 1, third%dims(2)
               do b = 1, first_two%dims(3)
                  do a = 1, first_two%dims(2)
                     out_x(obs + (a - 1) * nout + (b - 1) * nout * first_two%dims(2) + &
                        (c - 1) * nout * first_two%dims(2) * first_two%dims(3)) = &
                        fun_value(first_two, i, [a, b]) * fun_value(third, j, [c])
                  end do
               end do
            end do
         end do
      end do
      call create_from_grids(grids, dims, out_x, result, ok)
   end subroutine tensor_basis_with_third

   pure subroutine simulation_matrix_to_fun(template, matrix, result, ok)
      type(fun_data), intent(in) :: template !! Functional-data object supplying support grids and dimensions.
      real(dp), intent(in) :: matrix(:, :) !! Observation-by-flattened-support matrix to reshape into template support dimensions.
      type(fun_data), intent(out) :: result !! Functional-data object with matrix rows as observations on template support.
      logical, intent(out) :: ok !! True when matrix width equals the template support size.
      type(real_vector), allocatable :: grids(:)
      integer, allocatable :: dims(:)
      integer :: j

      ok = .false.
      if (size(matrix, 2) /= support_size(template)) return
      allocate(grids(size(template%argvals)))
      do j = 1, size(grids)
         grids(j)%values = template%argvals(j)%values
      end do
      allocate(dims(size(template%dims)))
      dims = template%dims
      dims(1) = size(matrix, 1)
      call create_from_grids(grids, dims, reshape(matrix, [size(matrix)]), result, ok)
   end subroutine simulation_matrix_to_fun

   pure subroutine create_from_grids(grids, dims, x, result, ok)
      type(real_vector), intent(in) :: grids(:) !! Support grids copied into the constructed functional-data object.
      integer, intent(in) :: dims(:) !! Observation and support dimensions for flattened values.
      real(dp), intent(in) :: x(:) !! Flattened column-major functional values.
      type(fun_data), intent(out) :: result !! Constructed regular functional-data object.
      logical, intent(out) :: ok !! True when dimensions, grids, and data size agree.
      integer :: j

      ok = .false.
      if (size(dims) /= size(grids) + 1) return
      if (product(dims) /= size(x)) return
      allocate(result%argvals(size(grids)))
      do j = 1, size(grids)
         if (size(grids(j)%values) /= dims(j + 1)) return
         result%argvals(j)%values = grids(j)%values
      end do
      result%dims = dims
      result%x = x
      ok = .true.
   end subroutine create_from_grids

   subroutine build_split_basis(specs, m, basis_type, result, ignore_deg, ok)
      type(basis_spec), intent(in) :: specs(:) !! One-dimensional component grids that form consecutive pieces of the split support.
      integer, intent(in) :: m !! Number of multivariate split-basis functions.
      character(len=*), intent(in) :: basis_type !! Basis family used on the concatenated translated grid.
      type(multi_fun_data), intent(out) :: result !! Multivariate basis obtained by splitting and randomly signing a large basis.
      integer, intent(in), optional :: ignore_deg(:) !! Optional PolyHigh exclusions for the large basis.
      logical, intent(out) :: ok !! True when every component is one-dimensional and basis construction succeeds.
      type(fun_data) :: big_basis
      real(dp), allocatable :: big_grid(:)
      real(dp), allocatable :: piece(:, :)
      integer :: i
      integer :: j
      integer :: offset
      integer :: total_points
      real(dp) :: shift
      real(dp) :: u
      real(dp) :: sign_value
      logical :: basis_ok

      ok = .false.
      total_points = 0
      do j = 1, size(specs)
         if (.not. allocated(specs(j)%argvals)) return
         if (size(specs(j)%argvals) /= 1) return
         total_points = total_points + size(specs(j)%argvals(1)%values)
      end do
      allocate(big_grid(total_points))
      offset = 0
      do j = 1, size(specs)
         if (j == 1) then
            big_grid(1:size(specs(j)%argvals(1)%values)) = specs(j)%argvals(1)%values
         else
            shift = big_grid(offset) - minval(specs(j)%argvals(1)%values)
            big_grid(offset + 1:offset + size(specs(j)%argvals(1)%values)) = &
               specs(j)%argvals(1)%values + shift
         end if
         offset = offset + size(specs(j)%argvals(1)%values)
      end do
      if (present(ignore_deg)) then
         call eigenfunctions(big_grid, m, basis_type, big_basis, ignore_deg, basis_ok)
      else
         call eigenfunctions(big_grid, m, basis_type, big_basis, ok = basis_ok)
      end if
      if (.not. basis_ok) return
      allocate(result%components(size(specs)))
      offset = 0
      do j = 1, size(specs)
         allocate(piece(m, size(specs(j)%argvals(1)%values)))
         do i = 1, m
            piece(i, :) = big_basis%x(i + offset * m:i + (offset + size(piece, 2) - 1) * m:m)
         end do
         call random_number(u)
         if (u < 0.5_dp) then
            sign_value = -1.0_dp
         else
            sign_value = 1.0_dp
         end if
         piece = sign_value * piece
         call create_fun_data_1d(specs(j)%argvals(1)%values, piece, result%components(j), basis_ok)
         if (.not. basis_ok) return
         offset = offset + size(piece, 2)
      end do
      ok = .true.
   end subroutine build_split_basis

   subroutine build_weighted_basis(specs, result, ok)
      type(basis_spec), intent(in) :: specs(:) !! Per-component one- or two-dimensional grids and marginal basis specifications.
      type(multi_fun_data), intent(out) :: result !! Weighted multivariate orthonormal basis with equal basis counts per component.
      logical, intent(out) :: ok !! True when each component has one or two support dimensions and equal total basis size.
      type(fun_data) :: first_basis
      type(fun_data) :: second_basis
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: weights(:)
      integer :: j
      integer :: m_total
      logical :: basis_ok

      ok = .false.
      allocate(alpha(size(specs)))
      call random_number(alpha)
      alpha = 0.2_dp + 0.6_dp * alpha
      allocate(weights(size(specs)))
      weights = sqrt(alpha / sum(alpha))
      allocate(result%components(size(specs)))
      m_total = -1
      do j = 1, size(specs)
         if (.not. allocated(specs(j)%argvals) .or. .not. allocated(specs(j)%m) .or. &
            .not. allocated(specs(j)%basis_type)) return
         if (size(specs(j)%argvals) < 1 .or. size(specs(j)%argvals) > 2) return
         if (size(specs(j)%m) /= size(specs(j)%argvals)) return
         if (size(specs(j)%basis_type) /= size(specs(j)%argvals)) return
         call build_spec_basis(specs(j), first_basis, second_basis, result%components(j), basis_ok)
         if (.not. basis_ok) return
         if (m_total < 0) then
            m_total = nobs_fun_data(result%components(j))
         else if (nobs_fun_data(result%components(j)) /= m_total) then
            return
         end if
         result%components(j)%x = weights(j) * result%components(j)%x
      end do
      ok = .true.
   end subroutine build_weighted_basis

   pure subroutine build_spec_basis(spec, first_basis, second_basis, result, ok)
      type(basis_spec), intent(in) :: spec !! One component's grid, marginal basis counts, families, and optional exclusions.
      type(fun_data), intent(out) :: first_basis !! Workspace/result for the first marginal basis.
      type(fun_data), intent(out) :: second_basis !! Workspace for an optional second marginal basis.
      type(fun_data), intent(out) :: result !! One- or two-dimensional component basis.
      logical, intent(out) :: ok !! True when all marginal bases and an optional tensor product are constructed.
      logical :: basis_ok
      logical :: tensor_ok

      ok = .false.
      if (allocated(spec%ignore_deg)) then
         if (size(spec%ignore_deg) >= 1) then
            if (allocated(spec%ignore_deg(1)%values)) then
               call eigenfunctions(spec%argvals(1)%values, spec%m(1), spec%basis_type(1), first_basis, &
                  spec%ignore_deg(1)%values, basis_ok)
            else
               call eigenfunctions(spec%argvals(1)%values, spec%m(1), spec%basis_type(1), first_basis, ok = basis_ok)
            end if
         else
            call eigenfunctions(spec%argvals(1)%values, spec%m(1), spec%basis_type(1), first_basis, ok = basis_ok)
         end if
      else
         call eigenfunctions(spec%argvals(1)%values, spec%m(1), spec%basis_type(1), first_basis, ok = basis_ok)
      end if
      if (.not. basis_ok) return
      if (size(spec%argvals) == 1) then
         result = first_basis
         ok = .true.
         return
      end if
      if (allocated(spec%ignore_deg)) then
         if (size(spec%ignore_deg) >= 2) then
            if (allocated(spec%ignore_deg(2)%values)) then
               call eigenfunctions(spec%argvals(2)%values, spec%m(2), spec%basis_type(2), second_basis, &
                  spec%ignore_deg(2)%values, basis_ok)
            else
               call eigenfunctions(spec%argvals(2)%values, spec%m(2), spec%basis_type(2), second_basis, ok = basis_ok)
            end if
         else
            call eigenfunctions(spec%argvals(2)%values, spec%m(2), spec%basis_type(2), second_basis, ok = basis_ok)
         end if
      else
         call eigenfunctions(spec%argvals(2)%values, spec%m(2), spec%basis_type(2), second_basis, ok = basis_ok)
      end if
      if (.not. basis_ok) return
      call tensor_product2(first_basis, second_basis, result, tensor_ok)
      if (.not. tensor_ok) return
      ok = .true.
   end subroutine build_spec_basis

   subroutine simulate_multi_from_basis(true_funs, value_type, n, result, ok)
      type(multi_fun_data), intent(in) :: true_funs !! Multivariate basis with equal observation count in every component.
      character(len=*), intent(in) :: value_type !! Eigenvalue family used as variances for common scores.
      integer, intent(in) :: n !! Number of multivariate observations to simulate.
      type(sim_multi_result), intent(inout) :: result !! Result receiving simulated data and eigenvalues while preserving true_funs.
      logical, intent(out) :: ok !! True when components have a common basis size and simulation succeeds.
      real(dp), allocatable :: scores(:, :)
      real(dp), allocatable :: matrix(:, :)
      integer :: i
      integer :: j
      integer :: p
      integer :: m_total
      logical :: build_ok

      ok = .false.
      if (.not. allocated(true_funs%components)) return
      if (size(true_funs%components) == 0) return
      m_total = nobs_fun_data(true_funs%components(1))
      do j = 2, size(true_funs%components)
         if (nobs_fun_data(true_funs%components(j)) /= m_total) return
      end do
      result%true_funs = true_funs
      call eigenvalues(m_total, value_type, result%true_vals, build_ok)
      if (.not. build_ok) return
      allocate(scores(n, m_total))
      do i = 1, n
         do j = 1, m_total
            scores(i, j) = normal_random() * sqrt(result%true_vals(j))
         end do
      end do
      allocate(result%sim_data%components(size(true_funs%components)))
      do j = 1, size(true_funs%components)
         allocate(matrix(n, support_size(true_funs%components(j))))
         matrix = 0.0_dp
         do p = 1, support_size(true_funs%components(j))
            do i = 1, n
               matrix(i, p) = dot_product(scores(i, :), &
                  true_funs%components(j)%x(1 + (p - 1) * m_total:p * m_total))
            end do
         end do
         call simulation_matrix_to_fun(true_funs%components(j), matrix, result%sim_data%components(j), build_ok)
         if (.not. build_ok) return
         deallocate(matrix)
      end do
      ok = .true.
   end subroutine simulate_multi_from_basis

   subroutine shuffle_int(x)
      integer, intent(inout) :: x(:) !! Integer vector shuffled in place using Fisher-Yates and Fortran random_number.
      real(dp) :: u
      integer :: i
      integer :: j
      integer :: tmp

      do i = size(x), 2, -1
         call random_number(u)
         j = 1 + int(u * real(i, dp))
         j = min(j, i)
         tmp = x(i)
         x(i) = x(j)
         x(j) = tmp
      end do
   end subroutine shuffle_int

   real(dp) function normal_random() result(value)
      real(dp) :: u1
      real(dp) :: u2

      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      value = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
   end function normal_random

   pure function to_lower(text) result(lower)
      character(len=*), intent(in) :: text !! ASCII option text normalized to lowercase.
      character(len=len(text)) :: lower
      integer :: code
      integer :: i

      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function to_lower

end module funData_simulation

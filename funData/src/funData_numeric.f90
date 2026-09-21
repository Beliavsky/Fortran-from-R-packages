module funData_numeric
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use funData_kinds, only : dp
   use funData_types, only : fun_data, int_vector, irreg_curve, irreg_fun_data, multi_fun_data, real_vector
   use funData_core, only : create_fun_data, create_fun_data_1d, extract_fun_data, fun_value, nobs_fun_data, &
      nobs_irreg_fun_data, nobs_multi_fun_data, support_dim, support_size
   implicit none
   private

   public :: approx_na
   public :: flip_fun_data
   public :: flip_fun_irreg_data
   public :: flip_irreg_fun_data
   public :: flip_multi_fun_data
   public :: int_weights
   public :: integrate_fun_data
   public :: integrate_irreg_fun_data
   public :: integrate_multi_fun_data
   public :: mean_fun_data
   public :: mean_irreg_fun_data
   public :: mean_multi_fun_data
   public :: norm_fun_data
   public :: norm_irreg_fun_data
   public :: norm_multi_fun_data
   public :: scalar_product_fun_data
   public :: scalar_product_fun_irreg
   public :: scalar_product_irreg_fun
   public :: scalar_product_irreg_fun_data
   public :: scalar_product_multi_fun_data
   public :: tensor_product2
   public :: tensor_product3

contains

   pure subroutine int_weights(argvals, method, weights, ok)
      real(dp), intent(in) :: argvals(:) !! One-dimensional observation grid used by the quadrature rule.
      character(len=*), intent(in) :: method !! Quadrature rule, either "trapezoidal" or "midpoint".
      real(dp), allocatable, intent(out) :: weights(:) !! Integration weights corresponding to each grid point.
      logical, intent(out) :: ok !! True when the method is recognized and weights are returned.
      character(len=:), allocatable :: actual_method
      integer :: d
      integer :: i

      ok = .false.
      d = size(argvals)
      allocate(weights(d))
      if (d == 0) then
         weights = 0.0_dp
         ok = .true.
         return
      end if

      actual_method = trim(to_lower(method))
      if (actual_method == "trapezoidal" .and. d < 3) actual_method = "midpoint"

      select case (actual_method)
      case ("midpoint")
         weights = 0.0_dp
         do i = 2, d
            weights(i) = argvals(i) - argvals(i - 1)
         end do
      case ("trapezoidal")
         weights(1) = 0.5_dp * (argvals(2) - argvals(1))
         do i = 2, d - 1
            weights(i) = 0.5_dp * (argvals(i + 1) - argvals(i - 1))
         end do
         weights(d) = 0.5_dp * (argvals(d) - argvals(d - 1))
      case default
         return
      end select
      ok = .true.
   end subroutine int_weights

   pure subroutine integrate_fun_data(object, values, method, ok)
      type(fun_data), intent(in) :: object !! Regular functional-data object with one to three support dimensions.
      real(dp), allocatable, intent(out) :: values(:) !! Numerical integral of each observation over the full rectangular support.
      character(len=*), intent(in), optional :: method !! Quadrature rule passed to int_weights; default is trapezoidal.
      logical, intent(out) :: ok !! True when the support dimension and quadrature method are supported.
      type(real_vector), allocatable :: weights(:)
      integer, allocatable :: sub(:)
      integer :: d
      integer :: i
      integer :: j
      integer :: p
      integer :: q
      integer :: n_support
      real(dp) :: w
      character(len=16) :: use_method
      logical :: weight_ok

      ok = .false.
      use_method = "trapezoidal"
      if (present(method)) use_method = method
      d = support_dim(object)
      if (d < 1 .or. d > 3) return
      allocate(weights(d))
      do j = 1, d
         call int_weights(object%argvals(j)%values, trim(use_method), weights(j)%values, weight_ok)
         if (.not. weight_ok) return
      end do

      allocate(values(nobs_fun_data(object)))
      values = 0.0_dp
      n_support = support_size(object)
      allocate(sub(d))
      do p = 0, n_support - 1
         q = p
         w = 1.0_dp
         do j = 1, d
            sub(j) = modulo(q, object%dims(j + 1)) + 1
            q = q / object%dims(j + 1)
            w = w * weights(j)%values(sub(j))
         end do
         do i = 1, nobs_fun_data(object)
            values(i) = values(i) + w * fun_value(object, i, sub)
         end do
      end do
      ok = .true.
   end subroutine integrate_fun_data

   pure subroutine integrate_irreg_fun_data(object, values, method, full_domain, ok)
      type(irreg_fun_data), intent(in) :: object !! Irregular one-dimensional functional data to integrate curve by curve.
      real(dp), allocatable, intent(out) :: values(:) !! Integral for every irregular observation.
      character(len=*), intent(in), optional :: method !! Quadrature rule; default is trapezoidal with short-grid fallback.
      logical, intent(in), optional :: full_domain !! If true, extrapolate every curve linearly to the global observed range.
      logical, intent(out) :: ok !! True when quadrature succeeds for every curve.
      type(irreg_curve) :: curve
      real(dp), allocatable :: weights(:)
      real(dp) :: domain(2)
      character(len=16) :: use_method
      logical :: do_full
      logical :: weight_ok
      integer :: i

      ok = .false.
      use_method = "trapezoidal"
      if (present(method)) use_method = method
      do_full = .false.
      if (present(full_domain)) do_full = full_domain
      allocate(values(nobs_irreg_fun_data(object)))
      values = 0.0_dp
      if (size(values) == 0) then
         ok = .true.
         return
      end if
      if (do_full) domain = global_irreg_range(object)

      do i = 1, size(values)
         curve = object%curves(i)
         if (do_full) call extrapolate_curve(curve, domain)
         call int_weights(curve%argvals, trim(use_method), weights, weight_ok)
         if (.not. weight_ok) return
         values(i) = sum(weights * curve%x)
      end do
      ok = .true.
   end subroutine integrate_irreg_fun_data

   pure subroutine integrate_multi_fun_data(object, values, method, ok)
      type(multi_fun_data), intent(in) :: object !! Multivariate functional-data object to integrate componentwise and sum.
      real(dp), allocatable, intent(out) :: values(:) !! Sum of component integrals for each observation.
      character(len=*), intent(in), optional :: method !! Quadrature rule used for every component.
      logical, intent(out) :: ok !! True when every component has a common observation count and integrates successfully.
      real(dp), allocatable :: component_values(:)
      logical :: component_ok
      integer :: j

      ok = .false.
      if (.not. allocated(object%components)) return
      allocate(values(nobs_multi_fun_data(object)))
      values = 0.0_dp
      do j = 1, size(object%components)
         if (nobs_fun_data(object%components(j)) /= size(values)) return
         if (present(method)) then
            call integrate_fun_data(object%components(j), component_values, method, component_ok)
         else
            call integrate_fun_data(object%components(j), component_values, ok = component_ok)
         end if
         if (.not. component_ok) return
         values = values + component_values
      end do
      ok = .true.
   end subroutine integrate_multi_fun_data

   pure subroutine norm_fun_data(object, values, squared, weight, method, obs, ok)
      type(fun_data), intent(in) :: object !! Regular functional-data object whose observation norms are requested.
      real(dp), allocatable, intent(out) :: values(:) !! Weighted squared norms or norms for selected observations.
      logical, intent(in), optional :: squared !! Return squared norms when true; default follows R and is true.
      real(dp), intent(in), optional :: weight !! Positive scalar multiplier applied to each squared norm; default is one.
      character(len=*), intent(in), optional :: method !! Integration method, defaulting to trapezoidal.
      integer, intent(in), optional :: obs(:) !! Optional one-based observation indices to evaluate.
      logical, intent(out) :: ok !! True when selection, integration, and weight validation succeed.
      type(fun_data) :: selected
      integer, allocatable :: selected_obs(:)
      real(dp), allocatable :: integrals(:)
      real(dp) :: use_weight
      logical :: use_squared
      logical :: select_ok
      integer :: i

      ok = .false.
      use_squared = .true.
      if (present(squared)) use_squared = squared
      use_weight = 1.0_dp
      if (present(weight)) use_weight = weight
      if (use_weight <= 0.0_dp) return
      if (present(obs)) then
         selected_obs = obs
      else
         allocate(selected_obs(nobs_fun_data(object)))
         do i = 1, size(selected_obs)
            selected_obs(i) = i
         end do
      end if
      call extract_fun_data(object, selected_obs, result = selected, ok = select_ok)
      if (.not. select_ok) return
      selected%x = selected%x * selected%x
      if (present(method)) then
         call integrate_fun_data(selected, integrals, method, ok)
      else
         call integrate_fun_data(selected, integrals, ok = ok)
      end if
      if (.not. ok) return
      values = use_weight * integrals
      if (.not. use_squared) values = sqrt(values)
   end subroutine norm_fun_data

   pure subroutine norm_irreg_fun_data(object, values, squared, weight, method, full_domain, obs, ok)
      type(irreg_fun_data), intent(in) :: object !! Irregular functional-data object whose observation norms are requested.
      real(dp), allocatable, intent(out) :: values(:) !! Weighted squared norms or norms for selected irregular observations.
      logical, intent(in), optional :: squared !! Return squared norms when true; default is true.
      real(dp), intent(in), optional :: weight !! Positive scalar multiplier applied to each squared norm; default is one.
      character(len=*), intent(in), optional :: method !! Integration method, defaulting to trapezoidal.
      logical, intent(in), optional :: full_domain !! Extrapolate selected curves to their global domain before integration.
      integer, intent(in), optional :: obs(:) !! Optional one-based observation indices to evaluate.
      logical, intent(out) :: ok !! True when indices, weight, and numerical integration are valid.
      type(irreg_fun_data) :: selected
      real(dp), allocatable :: integrals(:)
      integer, allocatable :: selected_obs(:)
      real(dp) :: use_weight
      real(dp) :: domain(2)
      logical :: use_squared
      logical :: do_full
      integer :: i
      integer :: j

      ok = .false.
      use_squared = .true.
      if (present(squared)) use_squared = squared
      use_weight = 1.0_dp
      if (present(weight)) use_weight = weight
      if (use_weight <= 0.0_dp) return
      do_full = .false.
      if (present(full_domain)) do_full = full_domain
      if (present(obs)) then
         selected_obs = obs
      else
         allocate(selected_obs(nobs_irreg_fun_data(object)))
         do i = 1, size(selected_obs)
            selected_obs(i) = i
         end do
      end if
      if (any(selected_obs < 1) .or. any(selected_obs > nobs_irreg_fun_data(object))) return
      allocate(selected%curves(size(selected_obs)))
      do i = 1, size(selected_obs)
         selected%curves(i) = object%curves(selected_obs(i))
      end do
      if (do_full) then
         domain = global_irreg_range(selected)
         do i = 1, size(selected%curves)
            call extrapolate_curve(selected%curves(i), domain)
         end do
      end if
      do i = 1, size(selected%curves)
         do j = 1, size(selected%curves(i)%x)
            selected%curves(i)%x(j) = selected%curves(i)%x(j) * selected%curves(i)%x(j)
         end do
      end do
      if (present(method)) then
         call integrate_irreg_fun_data(selected, integrals, method, .false., ok)
      else
         call integrate_irreg_fun_data(selected, integrals, full_domain = .false., ok = ok)
      end if
      if (.not. ok) return
      values = use_weight * integrals
      if (.not. use_squared) values = sqrt(values)
   end subroutine norm_irreg_fun_data

   pure subroutine norm_multi_fun_data(object, values, squared, weights, method, obs, ok)
      type(multi_fun_data), intent(in) :: object !! Multivariate functional-data object whose weighted norms are requested.
      real(dp), allocatable, intent(out) :: values(:) !! Weighted squared norms or norms across all components.
      logical, intent(in), optional :: squared !! Return squared multivariate norms when true; default is true.
      real(dp), intent(in), optional :: weights(:) !! Positive component weights; defaults to one for every component.
      character(len=*), intent(in), optional :: method !! Common integration method for all components.
      integer, intent(in), optional :: obs(:) !! Optional observation indices evaluated consistently in every component.
      logical, intent(out) :: ok !! True when weights and all component norm calculations succeed.
      real(dp), allocatable :: component_norms(:)
      real(dp), allocatable :: use_weights(:)
      logical :: component_ok
      logical :: use_squared
      integer :: j

      ok = .false.
      if (.not. allocated(object%components)) return
      allocate(use_weights(size(object%components)))
      use_weights = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= size(use_weights)) return
         use_weights = weights
      end if
      if (any(use_weights <= 0.0_dp)) return
      use_squared = .true.
      if (present(squared)) use_squared = squared
      if (present(obs)) then
         allocate(values(size(obs)))
      else
         allocate(values(nobs_multi_fun_data(object)))
      end if
      values = 0.0_dp
      do j = 1, size(object%components)
         if (present(method) .and. present(obs)) then
            call norm_fun_data(object%components(j), component_norms, squared = .true., weight = use_weights(j), &
               method = method, obs = obs, ok = component_ok)
         else if (present(method)) then
            call norm_fun_data(object%components(j), component_norms, squared = .true., weight = use_weights(j), &
               method = method, ok = component_ok)
         else if (present(obs)) then
            call norm_fun_data(object%components(j), component_norms, squared = .true., weight = use_weights(j), &
               obs = obs, ok = component_ok)
         else
            call norm_fun_data(object%components(j), component_norms, squared = .true., weight = use_weights(j), &
               ok = component_ok)
         end if
         if (.not. component_ok) return
         values = values + component_norms
      end do
      if (.not. use_squared) values = sqrt(values)
      ok = .true.
   end subroutine norm_multi_fun_data

   pure subroutine scalar_product_fun_data(object1, object2, values, method, ok)
      type(fun_data), intent(in) :: object1 !! First regular functional-data object.
      type(fun_data), intent(in) :: object2 !! Second regular object on the same support, allowing one-observation broadcasting.
      real(dp), allocatable, intent(out) :: values(:) !! Pairwise scalar products for the broadcasted observations.
      character(len=*), intent(in), optional :: method !! Integration method for the product functions.
      logical, intent(out) :: ok !! True when domains and observation counts are compatible.
      type(fun_data) :: product_data
      integer :: n
      integer :: n_support
      integer :: i
      integer :: p
      integer :: i1
      integer :: i2
      logical :: same_domain

      ok = .false.
      same_domain = regular_domains_equal(object1, object2)
      if (.not. same_domain) return
      if (nobs_fun_data(object1) /= nobs_fun_data(object2) .and. nobs_fun_data(object1) /= 1 .and. &
         nobs_fun_data(object2) /= 1) return
      n = max(nobs_fun_data(object1), nobs_fun_data(object2))
      product_data = object1
      product_data%dims(1) = n
      n_support = support_size(object1)
      if (allocated(product_data%x)) deallocate(product_data%x)
      allocate(product_data%x(n * n_support))
      do p = 0, n_support - 1
         do i = 1, n
            i1 = min(i, nobs_fun_data(object1))
            i2 = min(i, nobs_fun_data(object2))
            product_data%x(i + p * n) = object1%x(i1 + p * nobs_fun_data(object1)) * &
               object2%x(i2 + p * nobs_fun_data(object2))
         end do
      end do
      if (present(method)) then
         call integrate_fun_data(product_data, values, method, ok)
      else
         call integrate_fun_data(product_data, values, ok = ok)
      end if
   end subroutine scalar_product_fun_data

   pure subroutine scalar_product_irreg_fun_data(object1, object2, values, method, ok)
      type(irreg_fun_data), intent(in) :: object1 !! First irregular functional-data object.
      type(irreg_fun_data), intent(in) :: object2 !! Second irregular object with matching grids or one covering broadcast grid.
      real(dp), allocatable, intent(out) :: values(:) !! Pairwise scalar products for compatible irregular observations.
      character(len=*), intent(in), optional :: method !! Integration method used on each resulting irregular grid.
      logical, intent(out) :: ok !! True when grids and observation counts permit pointwise products.
      type(irreg_fun_data) :: product_data
      integer :: n
      integer :: i
      integer :: i1
      integer :: i2
      integer :: j
      integer :: k
      logical :: curve_ok

      ok = .false.
      if (nobs_irreg_fun_data(object1) /= nobs_irreg_fun_data(object2) .and. &
         nobs_irreg_fun_data(object1) /= 1 .and. nobs_irreg_fun_data(object2) /= 1) return
      n = max(nobs_irreg_fun_data(object1), nobs_irreg_fun_data(object2))
      allocate(product_data%curves(n))
      do i = 1, n
         i1 = min(i, nobs_irreg_fun_data(object1))
         i2 = min(i, nobs_irreg_fun_data(object2))
         call product_irreg_curves(object1%curves(i1), object2%curves(i2), product_data%curves(i), curve_ok)
         if (.not. curve_ok) return
         do j = 1, size(product_data%curves(i)%x)
            k = find_grid_value(object2%curves(i2)%argvals, product_data%curves(i)%argvals(j))
            if (k == 0) return
         end do
      end do
      if (present(method)) then
         call integrate_irreg_fun_data(product_data, values, method, .false., ok)
      else
         call integrate_irreg_fun_data(product_data, values, ok = ok)
      end if
   end subroutine scalar_product_irreg_fun_data

   pure subroutine scalar_product_fun_irreg(object1, object2, values, method, ok)
      type(fun_data), intent(in) :: object1 !! Regular one-dimensional reference object containing all irregular grid points.
      type(irreg_fun_data), intent(in) :: object2 !! Irregular object whose grids define the product domain.
      real(dp), allocatable, intent(out) :: values(:) !! Pairwise scalar products evaluated on each irregular observation grid.
      character(len=*), intent(in), optional :: method !! Integration method on irregular grids.
      logical, intent(out) :: ok !! True when regular and irregular observation counts and domains are compatible.
      type(irreg_fun_data) :: product_data
      integer :: n
      integer :: i
      integer :: i1
      integer :: i2
      integer :: j
      integer :: k

      ok = .false.
      if (support_dim(object1) /= 1) return
      if (nobs_fun_data(object1) /= nobs_irreg_fun_data(object2) .and. nobs_fun_data(object1) /= 1) return
      n = nobs_irreg_fun_data(object2)
      allocate(product_data%curves(n))
      do i = 1, n
         i1 = min(i, nobs_fun_data(object1))
         i2 = i
         product_data%curves(i)%argvals = object2%curves(i2)%argvals
         allocate(product_data%curves(i)%x(size(object2%curves(i2)%x)))
         do j = 1, size(object2%curves(i2)%x)
            k = find_grid_value(object1%argvals(1)%values, object2%curves(i2)%argvals(j))
            if (k == 0) return
            product_data%curves(i)%x(j) = fun_value(object1, i1, [k]) * object2%curves(i2)%x(j)
         end do
      end do
      if (present(method)) then
         call integrate_irreg_fun_data(product_data, values, method, .false., ok)
      else
         call integrate_irreg_fun_data(product_data, values, ok = ok)
      end if
   end subroutine scalar_product_fun_irreg

   pure subroutine scalar_product_irreg_fun(object1, object2, values, method, ok)
      type(irreg_fun_data), intent(in) :: object1 !! Irregular first object whose grids define the product domain.
      type(fun_data), intent(in) :: object2 !! Regular one-dimensional object containing all irregular grid points.
      real(dp), allocatable, intent(out) :: values(:) !! Pairwise scalar products evaluated on irregular grids.
      character(len=*), intent(in), optional :: method !! Integration method on irregular grids.
      logical, intent(out) :: ok !! True when observation counts and domains are compatible.

      if (present(method)) then
         call scalar_product_fun_irreg(object2, object1, values, method, ok)
      else
         call scalar_product_fun_irreg(object2, object1, values, ok = ok)
      end if
   end subroutine scalar_product_irreg_fun

   pure subroutine scalar_product_multi_fun_data(object1, object2, values, weights, method, ok)
      type(multi_fun_data), intent(in) :: object1 !! First multivariate functional-data object.
      type(multi_fun_data), intent(in) :: object2 !! Second multivariate object with matching component structure.
      real(dp), allocatable, intent(out) :: values(:) !! Weighted sum of component scalar products by observation.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative component weights with at least one positive entry.
      character(len=*), intent(in), optional :: method !! Common numerical integration method.
      logical, intent(out) :: ok !! True when component structures, weights, and domains are compatible.
      real(dp), allocatable :: component_values(:)
      real(dp), allocatable :: use_weights(:)
      logical :: component_ok
      integer :: j
      integer :: n

      ok = .false.
      if (.not. allocated(object1%components) .or. .not. allocated(object2%components)) return
      if (size(object1%components) /= size(object2%components)) return
      allocate(use_weights(size(object1%components)))
      use_weights = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= size(use_weights)) return
         use_weights = weights
      end if
      if (any(use_weights < 0.0_dp) .or. all(abs(use_weights) <= 0.0_dp)) return
      n = max(nobs_multi_fun_data(object1), nobs_multi_fun_data(object2))
      allocate(values(n))
      values = 0.0_dp
      do j = 1, size(object1%components)
         if (present(method)) then
            call scalar_product_fun_data(object1%components(j), object2%components(j), component_values, &
               method, component_ok)
         else
            call scalar_product_fun_data(object1%components(j), object2%components(j), component_values, ok = component_ok)
         end if
         if (.not. component_ok) return
         values = values + use_weights(j) * component_values
      end do
      ok = .true.
   end subroutine scalar_product_multi_fun_data

   pure subroutine mean_fun_data(object, result, remove_nan, ok)
      type(fun_data), intent(in) :: object !! Regular functional-data object whose pointwise mean is requested.
      type(fun_data), intent(out) :: result !! One-observation regular object containing the pointwise mean function.
      logical, intent(in), optional :: remove_nan !! Ignore NaNs when true; otherwise any NaN in a grid cell yields NaN.
      logical, intent(out) :: ok !! True when the input has at least one observation and mean construction succeeds.
      type(real_vector), allocatable :: grids(:)
      integer, allocatable :: dims(:)
      real(dp), allocatable :: means(:)
      real(dp) :: qnan
      real(dp) :: total
      integer :: count
      integer :: i
      integer :: j
      integer :: n
      integer :: p
      logical :: na_rm
      logical :: found_nan

      ok = .false.
      n = nobs_fun_data(object)
      if (n <= 0) return
      na_rm = .false.
      if (present(remove_nan)) na_rm = remove_nan
      qnan = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(means(support_size(object)))
      do p = 0, support_size(object) - 1
         total = 0.0_dp
         count = 0
         found_nan = .false.
         do i = 1, n
            if (ieee_is_nan(object%x(i + p * n))) then
               found_nan = .true.
            else
               total = total + object%x(i + p * n)
               count = count + 1
            end if
         end do
         if ((found_nan .and. .not. na_rm) .or. count == 0) then
            means(p + 1) = qnan
         else
            means(p + 1) = total / real(count, dp)
         end if
      end do
      allocate(grids(size(object%argvals)))
      do j = 1, size(grids)
         grids(j)%values = object%argvals(j)%values
      end do
      allocate(dims(size(object%dims)))
      dims = object%dims
      dims(1) = 1
      call create_fun_data(grids, dims, means, result, ok)
   end subroutine mean_fun_data

   pure subroutine mean_irreg_fun_data(object, result, ok)
      type(irreg_fun_data), intent(in) :: object !! Irregular data whose curves must share an identical grid.
      type(irreg_fun_data), intent(out) :: result !! One-curve irregular object containing the pointwise mean.
      logical, intent(out) :: ok !! True when all observation grids are identical and at least one observation exists.
      integer :: i

      ok = .false.
      if (nobs_irreg_fun_data(object) < 1) return
      do i = 2, nobs_irreg_fun_data(object)
         if (size(object%curves(i)%argvals) /= size(object%curves(1)%argvals)) return
         if (any(abs(object%curves(i)%argvals - object%curves(1)%argvals) > 0.0_dp)) return
      end do
      allocate(result%curves(1))
      result%curves(1)%argvals = object%curves(1)%argvals
      allocate(result%curves(1)%x(size(object%curves(1)%x)))
      result%curves(1)%x = 0.0_dp
      do i = 1, nobs_irreg_fun_data(object)
         result%curves(1)%x = result%curves(1)%x + object%curves(i)%x
      end do
      result%curves(1)%x = result%curves(1)%x / real(nobs_irreg_fun_data(object), dp)
      ok = .true.
   end subroutine mean_irreg_fun_data

   pure subroutine mean_multi_fun_data(object, result, remove_nan, ok)
      type(multi_fun_data), intent(in) :: object !! Multivariate functional data whose componentwise means are requested.
      type(multi_fun_data), intent(out) :: result !! Multivariate one-observation mean function.
      logical, intent(in), optional :: remove_nan !! Ignore NaNs within regular component means when true.
      logical, intent(out) :: ok !! True when every component mean succeeds.
      logical :: component_ok
      integer :: j

      ok = .false.
      if (.not. allocated(object%components)) return
      allocate(result%components(size(object%components)))
      do j = 1, size(object%components)
         if (present(remove_nan)) then
            call mean_fun_data(object%components(j), result%components(j), remove_nan, component_ok)
         else
            call mean_fun_data(object%components(j), result%components(j), ok = component_ok)
         end if
         if (.not. component_ok) return
      end do
      ok = .true.
   end subroutine mean_multi_fun_data

   pure subroutine approx_na(object, result, ok)
      type(fun_data), intent(in) :: object !! One-dimensional regular functional data containing missing values represented by NaNs.
      type(fun_data), intent(out) :: result !! Copy with interior missing runs linearly interpolated and edge NaNs retained.
      logical, intent(out) :: ok !! True when the input has one-dimensional support.
      integer :: i
      integer :: j
      integer :: left
      integer :: right
      integer :: m
      integer :: n
      real(dp) :: t

      ok = .false.
      if (support_dim(object) /= 1) return
      result = object
      n = nobs_fun_data(object)
      m = object%dims(2)
      do i = 1, n
         j = 1
         do while (j <= m)
            if (.not. ieee_is_nan(result%x(i + (j - 1) * n))) then
               j = j + 1
               cycle
            end if
            left = j - 1
            right = j
            do while (right <= m)
               if (.not. ieee_is_nan(result%x(i + (right - 1) * n))) exit
               right = right + 1
            end do
            if (left >= 1 .and. right <= m) then
               do j = left + 1, right - 1
                  t = (object%argvals(1)%values(j) - object%argvals(1)%values(left)) / &
                     (object%argvals(1)%values(right) - object%argvals(1)%values(left))
                  result%x(i + (j - 1) * n) = (1.0_dp - t) * result%x(i + (left - 1) * n) + &
                     t * result%x(i + (right - 1) * n)
               end do
            end if
            j = max(right, j + 1)
         end do
      end do
      ok = .true.
   end subroutine approx_na

   pure subroutine flip_fun_data(ref_object, new_object, result, ok)
      type(fun_data), intent(in) :: ref_object !! Regular reference functions defining the desired orientation.
      type(fun_data), intent(in) :: new_object !! Regular functions to flip when their scalar product with reference is negative.
      type(fun_data), intent(out) :: result !! Reoriented copy of new_object.
      logical, intent(out) :: ok !! True when domains and observation counts permit pairwise comparison.
      real(dp), allocatable :: products(:)
      integer :: i
      integer :: p
      integer :: n

      if (nobs_fun_data(ref_object) /= nobs_fun_data(new_object) .and. nobs_fun_data(ref_object) /= 1) then
         ok = .false.
         return
      end if
      call scalar_product_fun_data(ref_object, new_object, products, ok = ok)
      if (.not. ok) return
      result = new_object
      n = nobs_fun_data(new_object)
      do i = 1, n
         if (products(i) < 0.0_dp) then
            do p = 0, support_size(new_object) - 1
               result%x(i + p * n) = -result%x(i + p * n)
            end do
         end if
      end do
   end subroutine flip_fun_data

   pure subroutine flip_fun_irreg_data(ref_object, new_object, result, ok)
      type(fun_data), intent(in) :: ref_object !! Regular reference functions whose one-dimensional grid contains irregular grids.
      type(irreg_fun_data), intent(in) :: new_object !! Irregular functions to orient to the regular reference.
      type(irreg_fun_data), intent(out) :: result !! Reoriented irregular copy of new_object.
      logical, intent(out) :: ok !! True when the reference grid contains all irregular observation points.
      real(dp), allocatable :: products(:)
      integer :: i

      ok = .false.
      if (support_dim(ref_object) /= 1) return
      if (nobs_fun_data(ref_object) /= nobs_irreg_fun_data(new_object) .and. nobs_fun_data(ref_object) /= 1) return
      call scalar_product_fun_irreg(ref_object, new_object, products, ok = ok)
      if (.not. ok) return
      result = new_object
      do i = 1, size(products)
         if (products(i) < 0.0_dp) result%curves(i)%x = -result%curves(i)%x
      end do
   end subroutine flip_fun_irreg_data

   pure subroutine flip_irreg_fun_data(ref_object, new_object, result, ok)
      type(irreg_fun_data), intent(in) :: ref_object !! Irregular reference object, possibly containing one curve for broadcasting.
      type(irreg_fun_data), intent(in) :: new_object !! Irregular functions to orient to the reference.
      type(irreg_fun_data), intent(out) :: result !! Reoriented copy of new_object.
      logical, intent(out) :: ok !! True when grids and observation counts permit scalar-product comparison.
      real(dp), allocatable :: products(:)
      integer :: i

      if (nobs_irreg_fun_data(ref_object) /= nobs_irreg_fun_data(new_object) .and. &
         nobs_irreg_fun_data(ref_object) /= 1) then
         ok = .false.
         return
      end if
      call scalar_product_irreg_fun_data(ref_object, new_object, products, ok = ok)
      if (.not. ok) return
      result = new_object
      do i = 1, size(products)
         if (products(i) < 0.0_dp) result%curves(i)%x = -result%curves(i)%x
      end do
   end subroutine flip_irreg_fun_data

   pure subroutine flip_multi_fun_data(ref_object, new_object, result, ok)
      type(multi_fun_data), intent(in) :: ref_object !! Multivariate reference functions defining orientation.
      type(multi_fun_data), intent(in) :: new_object !! Multivariate functions to flip jointly across components.
      type(multi_fun_data), intent(out) :: result !! Reoriented multivariate copy of new_object.
      logical, intent(out) :: ok !! True when the two multivariate objects are compatible.
      real(dp), allocatable :: products(:)
      integer :: i
      integer :: j
      integer :: p
      integer :: n

      call scalar_product_multi_fun_data(ref_object, new_object, products, ok = ok)
      if (.not. ok) return
      result = new_object
      do j = 1, size(result%components)
         n = nobs_fun_data(result%components(j))
         do i = 1, n
            if (products(i) < 0.0_dp) then
               do p = 0, support_size(result%components(j)) - 1
                  result%components(j)%x(i + p * n) = -result%components(j)%x(i + p * n)
               end do
            end if
         end do
      end do
   end subroutine flip_multi_fun_data

   pure subroutine tensor_product2(object1, object2, result, ok)
      type(fun_data), intent(in) :: object1 !! First one-dimensional regular functional-data object.
      type(fun_data), intent(in) :: object2 !! Second one-dimensional regular functional-data object.
      type(fun_data), intent(out) :: result !! Two-dimensional tensor-product object for all observation combinations.
      logical, intent(out) :: ok !! True when both inputs have one-dimensional support.
      type(real_vector) :: grids(2)
      integer :: dims(3)
      real(dp), allocatable :: out_x(:)
      integer :: a
      integer :: b
      integer :: i
      integer :: j
      integer :: obs
      integer :: n1
      integer :: n2
      integer :: nout

      ok = .false.
      if (support_dim(object1) /= 1 .or. support_dim(object2) /= 1) return
      n1 = nobs_fun_data(object1)
      n2 = nobs_fun_data(object2)
      nout = n1 * n2
      dims = [nout, object1%dims(2), object2%dims(2)]
      grids(1)%values = object1%argvals(1)%values
      grids(2)%values = object2%argvals(1)%values
      allocate(out_x(product(dims)))
      do i = 1, n1
         do j = 1, n2
            obs = (i - 1) * n2 + j
            do b = 1, object2%dims(2)
               do a = 1, object1%dims(2)
                  out_x(obs + (a - 1) * nout + (b - 1) * nout * object1%dims(2)) = &
                     fun_value(object1, i, [a]) * fun_value(object2, j, [b])
               end do
            end do
         end do
      end do
      call create_fun_data(grids, dims, out_x, result, ok)
   end subroutine tensor_product2

   pure subroutine tensor_product3(object1, object2, object3, result, ok)
      type(fun_data), intent(in) :: object1 !! First one-dimensional regular functional-data object.
      type(fun_data), intent(in) :: object2 !! Second one-dimensional regular functional-data object.
      type(fun_data), intent(in) :: object3 !! Third one-dimensional regular functional-data object.
      type(fun_data), intent(out) :: result !! Three-dimensional tensor-product object for every observation combination.
      logical, intent(out) :: ok !! True when all three inputs have one-dimensional support.
      type(real_vector) :: grids(3)
      integer :: dims(4)
      real(dp), allocatable :: out_x(:)
      integer :: a
      integer :: b
      integer :: c
      integer :: i
      integer :: j
      integer :: k
      integer :: obs
      integer :: n1
      integer :: n2
      integer :: n3
      integer :: nout

      ok = .false.
      if (support_dim(object1) /= 1 .or. support_dim(object2) /= 1 .or. support_dim(object3) /= 1) return
      n1 = nobs_fun_data(object1)
      n2 = nobs_fun_data(object2)
      n3 = nobs_fun_data(object3)
      nout = n1 * n2 * n3
      dims = [nout, object1%dims(2), object2%dims(2), object3%dims(2)]
      grids(1)%values = object1%argvals(1)%values
      grids(2)%values = object2%argvals(1)%values
      grids(3)%values = object3%argvals(1)%values
      allocate(out_x(product(dims)))
      do i = 1, n1
         do j = 1, n2
            do k = 1, n3
               obs = (i - 1) * n2 * n3 + (j - 1) * n3 + k
               do c = 1, object3%dims(2)
                  do b = 1, object2%dims(2)
                     do a = 1, object1%dims(2)
                        out_x(obs + (a - 1) * nout + (b - 1) * nout * object1%dims(2) + &
                           (c - 1) * nout * object1%dims(2) * object2%dims(2)) = &
                           fun_value(object1, i, [a]) * fun_value(object2, j, [b]) * fun_value(object3, k, [c])
                     end do
                  end do
               end do
            end do
         end do
      end do
      call create_fun_data(grids, dims, out_x, result, ok)
   end subroutine tensor_product3

   pure logical function regular_domains_equal(object1, object2) result(equal)
      type(fun_data), intent(in) :: object1 !! First regular functional-data object in a domain comparison.
      type(fun_data), intent(in) :: object2 !! Second regular functional-data object in a domain comparison.
      integer :: j

      equal = .false.
      if (support_dim(object1) /= support_dim(object2)) return
      if (size(object1%dims) /= size(object2%dims)) return
      if (any(object1%dims(2:) /= object2%dims(2:))) return
      do j = 1, size(object1%argvals)
         if (size(object1%argvals(j)%values) /= size(object2%argvals(j)%values)) return
         if (any(abs(object1%argvals(j)%values - object2%argvals(j)%values) > 0.0_dp)) return
      end do
      equal = .true.
   end function regular_domains_equal

   pure integer function find_grid_value(grid, value) result(pos)
      real(dp), intent(in) :: grid(:) !! Grid searched by exact numerical equality.
      real(dp), intent(in) :: value !! Grid value whose position is requested.
      integer :: i

      pos = 0
      do i = 1, size(grid)
         if (abs(grid(i) - value) <= 0.0_dp) then
            pos = i
            return
         end if
      end do
   end function find_grid_value

   pure subroutine product_irreg_curves(curve1, curve2, product_curve, ok)
      type(irreg_curve), intent(in) :: curve1 !! First irregular curve, possibly defined on a superset of curve2's grid.
      type(irreg_curve), intent(in) :: curve2 !! Second irregular curve, possibly defined on a superset of curve1's grid.
      type(irreg_curve), intent(out) :: product_curve !! Pointwise product on the smaller compatible grid.
      logical, intent(out) :: ok !! True when one grid is identical to or contains the other grid.
      integer :: j
      integer :: k
      logical :: first_contains_second
      logical :: second_contains_first

      ok = .false.
      first_contains_second = .true.
      do j = 1, size(curve2%argvals)
         if (find_grid_value(curve1%argvals, curve2%argvals(j)) == 0) first_contains_second = .false.
      end do
      second_contains_first = .true.
      do j = 1, size(curve1%argvals)
         if (find_grid_value(curve2%argvals, curve1%argvals(j)) == 0) second_contains_first = .false.
      end do
      if (first_contains_second) then
         product_curve%argvals = curve2%argvals
         allocate(product_curve%x(size(curve2%x)))
         do j = 1, size(curve2%x)
            k = find_grid_value(curve1%argvals, curve2%argvals(j))
            product_curve%x(j) = curve1%x(k) * curve2%x(j)
         end do
      else if (second_contains_first) then
         product_curve%argvals = curve1%argvals
         allocate(product_curve%x(size(curve1%x)))
         do j = 1, size(curve1%x)
            k = find_grid_value(curve2%argvals, curve1%argvals(j))
            product_curve%x(j) = curve1%x(j) * curve2%x(k)
         end do
      else
         return
      end if
      ok = .true.
   end subroutine product_irreg_curves

   pure function global_irreg_range(object) result(rng)
      type(irreg_fun_data), intent(in) :: object !! Irregular object whose global minimum and maximum grid values are requested.
      real(dp) :: rng(2)
      integer :: i

      rng = [huge(0.0_dp), -huge(0.0_dp)]
      do i = 1, nobs_irreg_fun_data(object)
         if (size(object%curves(i)%argvals) > 0) then
            rng(1) = min(rng(1), minval(object%curves(i)%argvals))
            rng(2) = max(rng(2), maxval(object%curves(i)%argvals))
         end if
      end do
   end function global_irreg_range

   pure subroutine extrapolate_curve(curve, domain)
      type(irreg_curve), intent(inout) :: curve !! Irregular curve extended in place to the supplied domain endpoints.
      real(dp), intent(in) :: domain(2) !! Global lower and upper support boundaries used for extrapolation.
      real(dp), allocatable :: new_t(:)
      real(dp), allocatable :: new_x(:)
      integer :: n

      n = size(curve%argvals)
      if (n == 0) return
      allocate(new_t(n + 2))
      allocate(new_x(n + 2))
      new_t(1) = domain(1)
      new_t(2:n + 1) = curve%argvals
      new_t(n + 2) = domain(2)
      if (n == 1) then
         new_x = curve%x(1)
      else
         new_x(1) = curve%x(1) + (curve%x(2) - curve%x(1)) / &
            (curve%argvals(2) - curve%argvals(1)) * (domain(1) - curve%argvals(1))
         new_x(2:n + 1) = curve%x
         new_x(n + 2) = curve%x(n - 1) + (curve%x(n) - curve%x(n - 1)) / &
            (curve%argvals(n) - curve%argvals(n - 1)) * (domain(2) - curve%argvals(n - 1))
      end if
      curve%argvals = new_t
      curve%x = new_x
   end subroutine extrapolate_curve

   pure function to_lower(text) result(lower)
      character(len=*), intent(in) :: text !! ASCII text normalized to lowercase for option comparisons.
      character(len=len(text)) :: lower
      integer :: i
      integer :: code

      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function to_lower

end module funData_numeric

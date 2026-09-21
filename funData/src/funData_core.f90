module funData_core
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use funData_kinds, only : dp
   use funData_types, only : fun_data, int_vector, irreg_curve, irreg_fun_data, multi_fun_data, real_vector
   implicit none
   private

   public :: as_fun_data
   public :: as_irreg_fun_data
   public :: as_multi_fun_data
   public :: create_fun_data
   public :: create_fun_data_1d
   public :: create_irreg_fun_data
   public :: extract_fun_data
   public :: extract_irreg_fun_data
   public :: extract_multi_fun_data
   public :: fun_value
   public :: nobs_fun_data
   public :: nobs_irreg_fun_data
   public :: nobs_multi_fun_data
   public :: support_dim
   public :: support_size

contains

   pure subroutine create_fun_data(argvals, dims, x, result, ok)
      type(real_vector), intent(in) :: argvals(:) !! Sampling grids, one numeric vector per support dimension.
      integer, intent(in) :: dims(:) !! Array dimensions in R order: number of observations followed by support sizes.
      real(dp), intent(in) :: x(:) !! Flattened functional values in column-major order matching dims.
      type(fun_data), intent(out) :: result !! Constructed regular functional-data object.
      logical, intent(out) :: ok !! True when dimensions and sampling grids are mutually consistent.
      integer :: j

      ok = .false.
      if (size(dims) < 2) return
      if (dims(1) < 0) return
      if (size(argvals) /= size(dims) - 1) return
      if (product(dims) /= size(x)) return
      do j = 1, size(argvals)
         if (.not. allocated(argvals(j)%values)) return
         if (size(argvals(j)%values) /= dims(j + 1)) return
      end do

      allocate(result%argvals(size(argvals)))
      do j = 1, size(argvals)
         result%argvals(j)%values = argvals(j)%values
      end do
      result%dims = dims
      result%x = x
      ok = .true.
   end subroutine create_fun_data

   pure subroutine create_fun_data_1d(argvals, x, result, ok)
      real(dp), intent(in) :: argvals(:) !! One-dimensional sampling grid in increasing or otherwise user-supplied order.
      real(dp), intent(in) :: x(:, :) !! Observation-by-grid matrix of functional values.
      type(fun_data), intent(out) :: result !! Constructed one-dimensional regular functional-data object.
      logical, intent(out) :: ok !! True when the matrix width matches the grid length.
      type(real_vector) :: grids(1)
      integer :: dims(2)

      grids(1)%values = argvals
      dims = [size(x, 1), size(x, 2)]
      call create_fun_data(grids, dims, reshape(x, [size(x)]), result, ok)
   end subroutine create_fun_data_1d

   pure subroutine create_irreg_fun_data(argvals, values, result, ok)
      type(real_vector), intent(in) :: argvals(:) !! Per-observation sampling grids for irregular one-dimensional functions.
      type(real_vector), intent(in) :: values(:) !! Per-observation values corresponding one-for-one to argvals entries.
      type(irreg_fun_data), intent(out) :: result !! Constructed irregular functional-data object.
      logical, intent(out) :: ok !! True when every observation has matching grid and value lengths.
      integer :: i

      ok = .false.
      if (size(argvals) /= size(values)) return
      allocate(result%curves(size(argvals)))
      do i = 1, size(argvals)
         if (.not. allocated(argvals(i)%values)) return
         if (.not. allocated(values(i)%values)) return
         if (size(argvals(i)%values) /= size(values(i)%values)) return
         result%curves(i)%argvals = argvals(i)%values
         result%curves(i)%x = values(i)%values
      end do
      ok = .true.
   end subroutine create_irreg_fun_data

   pure integer function nobs_fun_data(object) result(n)
      type(fun_data), intent(in) :: object !! Regular functional-data object whose observation count is requested.
      if (allocated(object%dims)) then
         n = object%dims(1)
      else
         n = 0
      end if
   end function nobs_fun_data

   pure integer function nobs_irreg_fun_data(object) result(n)
      type(irreg_fun_data), intent(in) :: object !! Irregular functional-data object whose observation count is requested.
      if (allocated(object%curves)) then
         n = size(object%curves)
      else
         n = 0
      end if
   end function nobs_irreg_fun_data

   pure integer function nobs_multi_fun_data(object) result(n)
      type(multi_fun_data), intent(in) :: object !! Multivariate functional-data object whose common observation count is requested.
      if (.not. allocated(object%components)) then
         n = 0
      else if (size(object%components) == 0) then
         n = 0
      else
         n = nobs_fun_data(object%components(1))
      end if
   end function nobs_multi_fun_data

   pure integer function support_dim(object) result(d)
      type(fun_data), intent(in) :: object !! Regular functional-data object whose support dimension is requested.
      if (allocated(object%dims)) then
         d = max(0, size(object%dims) - 1)
      else
         d = 0
      end if
   end function support_dim

   pure integer function support_size(object) result(n)
      type(fun_data), intent(in) :: object !! Regular functional-data object whose total number of grid points is requested.
      if (.not. allocated(object%dims)) then
         n = 0
      else if (size(object%dims) <= 1) then
         n = 0
      else
         n = product(object%dims(2:))
      end if
   end function support_size

   pure integer function flat_index(dims, indices) result(idx)
      integer, intent(in) :: dims(:) !! Dimensions of a column-major array.
      integer, intent(in) :: indices(:) !! One-based subscript in every dimension.
      integer :: j
      integer :: stride

      idx = indices(1)
      stride = dims(1)
      do j = 2, size(dims)
         idx = idx + (indices(j) - 1) * stride
         stride = stride * dims(j)
      end do
   end function flat_index

   pure real(dp) function fun_value(object, obs, support_indices) result(value)
      type(fun_data), intent(in) :: object !! Regular functional-data object to inspect.
      integer, intent(in) :: obs !! One-based observation index.
      integer, intent(in) :: support_indices(:) !! One-based grid index in each support dimension.
      integer, allocatable :: indices(:)

      allocate(indices(size(object%dims)))
      indices(1) = obs
      indices(2:) = support_indices
      value = object%x(flat_index(object%dims, indices))
   end function fun_value

   pure subroutine as_multi_fun_data(object, result)
      type(fun_data), intent(in) :: object !! Regular functional-data object to wrap as one multivariate component.
      type(multi_fun_data), intent(out) :: result !! One-component multivariate functional-data object.

      allocate(result%components(1))
      result%components(1) = object
   end subroutine as_multi_fun_data

   pure subroutine as_irreg_fun_data(object, result, ok)
      type(fun_data), intent(in) :: object !! One-dimensional regular functional data, possibly containing NaNs.
      type(irreg_fun_data), intent(out) :: result !! Irregular object formed by dropping NaN values independently by observation.
      logical, intent(out) :: ok !! True when the input has one-dimensional support.
      integer :: i
      integer :: j
      integer :: k
      integer :: m
      integer :: n
      real(dp) :: value

      ok = .false.
      if (support_dim(object) /= 1) return
      n = nobs_fun_data(object)
      m = object%dims(2)
      allocate(result%curves(n))
      do i = 1, n
         k = 0
         do j = 1, m
            value = fun_value(object, i, [j])
            if (.not. ieee_is_nan(value)) k = k + 1
         end do
         allocate(result%curves(i)%argvals(k))
         allocate(result%curves(i)%x(k))
         k = 0
         do j = 1, m
            value = fun_value(object, i, [j])
            if (.not. ieee_is_nan(value)) then
               k = k + 1
               result%curves(i)%argvals(k) = object%argvals(1)%values(j)
               result%curves(i)%x(k) = value
            end if
         end do
      end do
      ok = .true.
   end subroutine as_irreg_fun_data

   pure subroutine as_fun_data(object, result, ok)
      type(irreg_fun_data), intent(in) :: object !! Irregular one-dimensional functional data to place on the union grid.
      type(fun_data), intent(out) :: result !! Regular object on the sorted union grid with missing entries represented by NaN.
      logical, intent(out) :: ok !! True when conversion succeeds, including empty inputs.
      real(dp), allocatable :: grid(:)
      real(dp), allocatable :: temp(:)
      real(dp), allocatable :: mat(:, :)
      integer :: i
      integer :: j
      integer :: k
      integer :: n_total
      integer :: n_unique
      real(dp) :: qnan

      ok = .false.
      qnan = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. allocated(object%curves)) then
         allocate(grid(0))
         allocate(mat(0, 0))
         call create_fun_data_1d(grid, mat, result, ok)
         return
      end if

      n_total = 0
      do i = 1, size(object%curves)
         n_total = n_total + size(object%curves(i)%argvals)
      end do
      allocate(temp(n_total))
      k = 0
      do i = 1, size(object%curves)
         do j = 1, size(object%curves(i)%argvals)
            k = k + 1
            temp(k) = object%curves(i)%argvals(j)
         end do
      end do
      call sort_real(temp)
      n_unique = 0
      do i = 1, size(temp)
         if (i == 1) then
            n_unique = 1
         else if (abs(temp(i) - temp(i - 1)) > 0.0_dp) then
            n_unique = n_unique + 1
         end if
      end do
      allocate(grid(n_unique))
      k = 0
      do i = 1, size(temp)
         if (i == 1 .or. abs(temp(i) - temp(max(1, i - 1))) > 0.0_dp) then
            k = k + 1
            grid(k) = temp(i)
         end if
      end do

      allocate(mat(size(object%curves), n_unique))
      mat = qnan
      do i = 1, size(object%curves)
         do j = 1, size(object%curves(i)%argvals)
            k = find_exact(grid, object%curves(i)%argvals(j))
            if (k > 0) mat(i, k) = object%curves(i)%x(j)
         end do
      end do
      call create_fun_data_1d(grid, mat, result, ok)
   end subroutine as_fun_data

   pure integer function find_exact(x, value) result(pos)
      real(dp), intent(in) :: x(:) !! Numeric vector to search using exact equality, matching R's grid membership behavior.
      real(dp), intent(in) :: value !! Value whose first matching position is requested.
      integer :: i

      pos = 0
      do i = 1, size(x)
         if (abs(x(i) - value) <= 0.0_dp) then
            pos = i
            return
         end if
      end do
   end function find_exact

   pure subroutine sort_real(x)
      real(dp), intent(inout) :: x(:) !! Numeric vector sorted into nondecreasing order in place.
      integer :: i
      integer :: j
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
   end subroutine sort_real

   pure subroutine extract_fun_data(object, obs, support_indices, result, ok)
      type(fun_data), intent(in) :: object !! Regular functional-data object to subset.
      integer, intent(in) :: obs(:) !! One-based observation indices to retain, in requested order.
      type(int_vector), intent(in), optional :: support_indices(:) !! Optional one-based indices retained in each support dimension.
      type(fun_data), intent(out) :: result !! Extracted regular functional-data object.
      logical, intent(out) :: ok !! True when all requested indices are valid and the extraction succeeds.
      type(real_vector), allocatable :: grids(:)
      type(int_vector), allocatable :: sel(:)
      integer, allocatable :: out_dims(:)
      integer, allocatable :: out_sub(:)
      integer, allocatable :: in_sub(:)
      real(dp), allocatable :: out_x(:)
      integer :: d
      integer :: i
      integer :: j
      integer :: p
      integer :: q
      integer :: total_support

      ok = .false.
      d = support_dim(object)
      if (any(obs < 1) .or. any(obs > nobs_fun_data(object))) return
      allocate(sel(d))
      do j = 1, d
         if (present(support_indices)) then
            if (size(support_indices) /= d) return
            if (.not. allocated(support_indices(j)%values)) return
            sel(j)%values = support_indices(j)%values
         else
            allocate(sel(j)%values(object%dims(j + 1)))
            do i = 1, object%dims(j + 1)
               sel(j)%values(i) = i
            end do
         end if
         if (any(sel(j)%values < 1) .or. any(sel(j)%values > object%dims(j + 1))) return
      end do

      allocate(grids(d))
      allocate(out_dims(d + 1))
      out_dims(1) = size(obs)
      do j = 1, d
         grids(j)%values = object%argvals(j)%values(sel(j)%values)
         out_dims(j + 1) = size(sel(j)%values)
      end do
      total_support = product(out_dims(2:))
      allocate(out_x(size(obs) * total_support))
      allocate(out_sub(d))
      allocate(in_sub(d))

      do p = 0, total_support - 1
         q = p
         do j = 1, d
            out_sub(j) = modulo(q, out_dims(j + 1)) + 1
            q = q / out_dims(j + 1)
            in_sub(j) = sel(j)%values(out_sub(j))
         end do
         do i = 1, size(obs)
            out_x(i + p * size(obs)) = fun_value(object, obs(i), in_sub)
         end do
      end do
      call create_fun_data(grids, out_dims, out_x, result, ok)
   end subroutine extract_fun_data

   pure subroutine extract_irreg_fun_data(object, obs, requested_argvals, result, ok)
      type(irreg_fun_data), intent(in) :: object !! Irregular functional-data object to subset.
      integer, intent(in) :: obs(:) !! One-based observation indices to retain.
      real(dp), intent(in), optional :: requested_argvals(:) !! Optional exact grid values retained within each curve.
      type(irreg_fun_data), intent(out) :: result !! Extracted irregular object; empty curves are omitted after domain filtering.
      logical, intent(out) :: ok !! True when observation indices are valid and extraction succeeds.
      integer :: i
      integer :: j
      integer :: k
      integer :: keep_count
      integer :: curve_count
      logical :: keep

      ok = .false.
      if (any(obs < 1) .or. any(obs > nobs_irreg_fun_data(object))) return
      curve_count = 0
      do i = 1, size(obs)
         if (.not. present(requested_argvals)) then
            curve_count = curve_count + 1
         else
            keep_count = 0
            do j = 1, size(object%curves(obs(i))%argvals)
               keep = any(abs(object%curves(obs(i))%argvals(j) - requested_argvals) <= 0.0_dp)
               if (keep) keep_count = keep_count + 1
            end do
            if (keep_count > 0) curve_count = curve_count + 1
         end if
      end do
      allocate(result%curves(curve_count))
      k = 0
      do i = 1, size(obs)
         if (.not. present(requested_argvals)) then
            k = k + 1
            result%curves(k) = object%curves(obs(i))
         else
            keep_count = 0
            do j = 1, size(object%curves(obs(i))%argvals)
               if (any(abs(object%curves(obs(i))%argvals(j) - requested_argvals) <= 0.0_dp)) keep_count = keep_count + 1
            end do
            if (keep_count > 0) then
               k = k + 1
               allocate(result%curves(k)%argvals(keep_count))
               allocate(result%curves(k)%x(keep_count))
               keep_count = 0
               do j = 1, size(object%curves(obs(i))%argvals)
                  if (any(abs(object%curves(obs(i))%argvals(j) - requested_argvals) <= 0.0_dp)) then
                     keep_count = keep_count + 1
                     result%curves(k)%argvals(keep_count) = object%curves(obs(i))%argvals(j)
                     result%curves(k)%x(keep_count) = object%curves(obs(i))%x(j)
                  end if
               end do
            end if
         end if
      end do
      ok = .true.
   end subroutine extract_irreg_fun_data

   pure subroutine extract_multi_fun_data(object, obs, result, ok)
      type(multi_fun_data), intent(in) :: object !! Multivariate functional-data object whose observations are to be subset.
      integer, intent(in) :: obs(:) !! One-based observation indices retained in every component.
      type(multi_fun_data), intent(out) :: result !! Multivariate object with the selected observations.
      logical, intent(out) :: ok !! True when each component can be subset with the requested observation indices.
      integer :: j
      logical :: component_ok

      ok = .false.
      if (.not. allocated(object%components)) return
      allocate(result%components(size(object%components)))
      do j = 1, size(object%components)
         call extract_fun_data(object%components(j), obs, result = result%components(j), ok = component_ok)
         if (.not. component_ok) return
      end do
      ok = .true.
   end subroutine extract_multi_fun_data

end module funData_core

! desirability-fortran
!
! Modern Fortran translation of the computational code in the R package
! "desirability" by Max Kuhn, version 2.1.
!
! The original package is distributed under the GNU General Public License,
! version 2. This translation is distributed under the same license.
! See COPYING.

module desirability
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    implicit none
    private

    integer, parameter, public :: dp = real64
    integer, parameter :: input_missing = 0
    integer, parameter :: input_numeric = 1
    integer, parameter :: input_categorical = 2
    integer, parameter :: n_noninformative = 100

    type, public :: desirability_input
        private
        integer :: kind = input_missing
        real(dp) :: number = 0.0_dp
        character(len=:), allocatable :: text
    end type desirability_input

    type, abstract, public :: desirability_function
        real(dp) :: missing = 0.0_dp
        logical :: has_tol = .false.
        real(dp) :: tol = 0.0_dp
    contains
        procedure(raw_value_interface), deferred :: raw_value
        procedure :: value => evaluate_input
    end type desirability_function

    type, extends(desirability_function), public :: d_max_type
        real(dp) :: low = 0.0_dp
        real(dp) :: high = 1.0_dp
        real(dp) :: scale = 1.0_dp
    contains
        procedure :: raw_value => d_max_raw
    end type d_max_type

    type, extends(desirability_function), public :: d_min_type
        real(dp) :: low = 0.0_dp
        real(dp) :: high = 1.0_dp
        real(dp) :: scale = 1.0_dp
    contains
        procedure :: raw_value => d_min_raw
    end type d_min_type

    type, extends(desirability_function), public :: d_target_type
        real(dp) :: low = 0.0_dp
        real(dp) :: target = 0.5_dp
        real(dp) :: high = 1.0_dp
        real(dp) :: low_scale = 1.0_dp
        real(dp) :: high_scale = 1.0_dp
    contains
        procedure :: raw_value => d_target_raw
    end type d_target_type

    type, extends(desirability_function), public :: d_box_type
        real(dp) :: low = 0.0_dp
        real(dp) :: high = 1.0_dp
    contains
        procedure :: raw_value => d_box_raw
    end type d_box_type

    type, extends(desirability_function), public :: d_arb_type
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: d(:)
        real(dp), allocatable :: interp_x(:)
        real(dp), allocatable :: interp_d(:)
    contains
        procedure :: raw_value => d_arb_raw
    end type d_arb_type

    type, extends(desirability_function), public :: d_categorical_type
        character(len=:), allocatable :: labels(:)
        real(dp), allocatable :: values(:)
    contains
        procedure :: raw_value => d_categorical_raw
    end type d_categorical_type

    type, public :: desirability_holder
        class(desirability_function), allocatable :: fun
    end type desirability_holder

    type, public :: d_overall_type
        type(desirability_holder), allocatable :: d(:)
    contains
        procedure :: add => d_overall_add
        procedure :: size => d_overall_size
    end type d_overall_type

    abstract interface
        function raw_value_interface(self, input) result(out)
            import :: dp, desirability_function, desirability_input
            class(desirability_function), intent(in) :: self
            type(desirability_input), intent(in) :: input
            real(dp) :: out
        end function raw_value_interface
    end interface

    interface predict
        module procedure predict_numeric_scalar
        module procedure predict_numeric_vector
        module procedure predict_category_scalar
        module procedure predict_category_vector
        module procedure predict_overall_numeric
        module procedure predict_overall_mixed
    end interface predict

    interface predict_all
        module procedure predict_all_numeric
        module procedure predict_all_mixed
    end interface predict_all

    public :: d_max, d_min, d_target, d_box, d_arb, d_categorical
    public :: d_overall, hold, predict, predict_all
    public :: numeric_input, categorical_input, missing_input, quiet_nan

contains

    function quiet_nan() result(x)
        real(dp) :: x

        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function quiet_nan

    function numeric_input(x) result(input)
        real(dp), intent(in) :: x
        type(desirability_input) :: input

        if (ieee_is_nan(x)) then
            input%kind = input_missing
        else
            input%kind = input_numeric
            input%number = x
        end if
    end function numeric_input

    function categorical_input(x) result(input)
        character(len=*), intent(in) :: x
        type(desirability_input) :: input

        if (len_trim(x) == 0) then
            input%kind = input_missing
        else
            input%kind = input_categorical
            input%text = trim(x)
        end if
    end function categorical_input

    function missing_input() result(input)
        type(desirability_input) :: input

        input%kind = input_missing
    end function missing_input

    function evaluate_input(self, input, missing) result(out)
        class(desirability_function), intent(in) :: self
        type(desirability_input), intent(in) :: input
        real(dp), intent(in), optional :: missing
        real(dp) :: out

        if (input%kind == input_missing) then
            if (present(missing)) then
                out = missing
            else
                out = self%missing
            end if
        else
            out = self%raw_value(input)
        end if

        if (self%has_tol) then
            if (out == 0.0_dp) out = self%tol
        end if
    end function evaluate_input

    function d_max(low, high, scale, tol) result(obj)
        real(dp), intent(in) :: low, high
        real(dp), intent(in), optional :: scale, tol
        type(d_max_type) :: obj

        if (low >= high) error stop "d_max: low must be less than high"

        obj%low = low
        obj%high = high
        if (present(scale)) obj%scale = scale
        if (obj%scale <= 0.0_dp) error stop "d_max: scale must be positive"
        call set_tolerance(obj, tol)
        obj%missing = mean_on_grid(obj, low, high)
    end function d_max

    function d_min(low, high, scale, tol) result(obj)
        real(dp), intent(in) :: low, high
        real(dp), intent(in), optional :: scale, tol
        type(d_min_type) :: obj

        if (low >= high) error stop "d_min: low must be less than high"

        obj%low = low
        obj%high = high
        if (present(scale)) obj%scale = scale
        if (obj%scale <= 0.0_dp) error stop "d_min: scale must be positive"
        call set_tolerance(obj, tol)
        obj%missing = mean_on_grid(obj, low, high)
    end function d_min

    function d_target(low, target, high, low_scale, high_scale, tol) result(obj)
        real(dp), intent(in) :: low, target, high
        real(dp), intent(in), optional :: low_scale, high_scale, tol
        type(d_target_type) :: obj

        if (low >= high) error stop "d_target: low must be less than high"
        if (low >= target) error stop "d_target: low must be less than target"
        if (target >= high) error stop "d_target: target must be less than high"

        obj%low = low
        obj%target = target
        obj%high = high
        if (present(low_scale)) obj%low_scale = low_scale
        if (present(high_scale)) obj%high_scale = high_scale
        if (obj%low_scale <= 0.0_dp .or. obj%high_scale <= 0.0_dp) then
            error stop "d_target: scale parameters must be positive"
        end if
        call set_tolerance(obj, tol)
        obj%missing = mean_on_grid(obj, low, high)
    end function d_target

    function d_box(low, high, tol) result(obj)
        real(dp), intent(in) :: low, high
        real(dp), intent(in), optional :: tol
        type(d_box_type) :: obj

        if (low >= high) error stop "d_box: low must be less than high"

        obj%low = low
        obj%high = high
        call set_tolerance(obj, tol)
        obj%missing = mean_on_grid(obj, low, high)
    end function d_box

    function d_arb(x, d, tol) result(obj)
        real(dp), intent(in) :: x(:), d(:)
        real(dp), intent(in), optional :: tol
        type(d_arb_type) :: obj
        integer :: i

        if (size(x) /= size(d)) error stop "d_arb: x and d must have the same length"
        if (size(x) < 2) error stop "d_arb: x and d must contain at least two values"
        if (any(ieee_is_nan(x)) .or. any(ieee_is_nan(d))) then
            error stop "d_arb: x and d may not contain NaN"
        end if
        if (any(d < 0.0_dp) .or. any(d > 1.0_dp)) then
            error stop "d_arb: desirabilities must satisfy 0 <= d <= 1"
        end if

        allocate(obj%x(size(x)), obj%d(size(d)))
        obj%x = x
        obj%d = d
        call stable_sort_pairs(obj%x, obj%d)
        call build_interpolation_knots(obj%x, obj%d, obj%interp_x, obj%interp_d)
        if (size(obj%interp_x) < 2) then
            error stop "d_arb: at least two distinct x values are required"
        end if

        call set_tolerance(obj, tol)
        obj%missing = mean_on_grid(obj, obj%x(1), obj%x(size(obj%x)))

        do i = 2, size(obj%x)
            if (obj%x(i) < obj%x(i - 1)) error stop "d_arb: internal sorting failure"
        end do
    end function d_arb

    function d_categorical(labels, values, tol) result(obj)
        character(len=*), intent(in) :: labels(:)
        real(dp), intent(in) :: values(:)
        real(dp), intent(in), optional :: tol
        type(d_categorical_type) :: obj
        integer :: i, j, max_len

        if (size(values) < 2) then
            error stop "d_categorical: at least two values are required"
        end if
        if (size(labels) /= size(values)) then
            error stop "d_categorical: labels and values must have the same length"
        end if
        if (any([(len_trim(labels(i)) == 0, i = 1, size(labels))])) then
            error stop "d_categorical: labels must be nonempty"
        end if
        do i = 1, size(labels) - 1
            do j = i + 1, size(labels)
                if (trim(labels(i)) == trim(labels(j))) then
                    error stop "d_categorical: labels must be unique"
                end if
            end do
        end do

        max_len = 1
        do i = 1, size(labels)
            max_len = max(max_len, len_trim(labels(i)))
        end do
        allocate(character(len=max_len) :: obj%labels(size(labels)))
        allocate(obj%values(size(values)))
        do i = 1, size(labels)
            obj%labels(i) = trim(labels(i))
        end do
        obj%values = values
        call set_tolerance(obj, tol)
        obj%missing = sum(values) / real(size(values), dp)
    end function d_categorical

    function d_overall(functions) result(obj)
        type(desirability_holder), intent(in), optional :: functions(:)
        type(d_overall_type) :: obj
        integer :: i

        if (.not. present(functions)) then
            allocate(obj%d(0))
            return
        end if

        allocate(obj%d(size(functions)))
        do i = 1, size(functions)
            if (.not. allocated(functions(i)%fun)) then
                error stop "d_overall: every holder must contain a function"
            end if
            allocate(obj%d(i)%fun, source=functions(i)%fun)
        end do
    end function d_overall

    function hold(fun) result(holder)
        class(desirability_function), intent(in) :: fun
        type(desirability_holder) :: holder

        allocate(holder%fun, source=fun)
    end function hold

    subroutine d_overall_add(self, fun)
        class(d_overall_type), intent(inout) :: self
        class(desirability_function), intent(in) :: fun
        type(desirability_holder), allocatable :: tmp(:)
        integer :: i, n

        if (allocated(self%d)) then
            n = size(self%d)
        else
            n = 0
        end if

        allocate(tmp(n + 1))
        do i = 1, n
            allocate(tmp(i)%fun, source=self%d(i)%fun)
        end do
        allocate(tmp(n + 1)%fun, source=fun)
        call move_alloc(tmp, self%d)
    end subroutine d_overall_add

    integer function d_overall_size(self) result(n)
        class(d_overall_type), intent(in) :: self

        if (allocated(self%d)) then
            n = size(self%d)
        else
            n = 0
        end if
    end function d_overall_size

    function predict_numeric_scalar(object, newdata, missing) result(out)
        class(desirability_function), intent(in) :: object
        real(dp), intent(in) :: newdata
        real(dp), intent(in), optional :: missing
        real(dp) :: out

        select type (object)
        type is (d_categorical_type)
            error stop "predict: categorical desirability requires character input"
        class default
            out = object%value(numeric_input(newdata), missing)
        end select
    end function predict_numeric_scalar

    function predict_numeric_vector(object, newdata, missing) result(out)
        class(desirability_function), intent(in) :: object
        real(dp), intent(in) :: newdata(:)
        real(dp), intent(in), optional :: missing
        real(dp), allocatable :: out(:)
        integer :: i

        allocate(out(size(newdata)))
        do i = 1, size(newdata)
            out(i) = predict_numeric_scalar(object, newdata(i), missing)
        end do
    end function predict_numeric_vector

    function predict_category_scalar(object, newdata, missing) result(out)
        class(desirability_function), intent(in) :: object
        character(len=*), intent(in) :: newdata
        real(dp), intent(in), optional :: missing
        real(dp) :: out

        select type (object)
        type is (d_categorical_type)
            out = object%value(categorical_input(newdata), missing)
        class default
            error stop "predict: numeric desirability requires real input"
        end select
    end function predict_category_scalar

    function predict_category_vector(object, newdata, missing) result(out)
        class(desirability_function), intent(in) :: object
        character(len=*), intent(in) :: newdata(:)
        real(dp), intent(in), optional :: missing
        real(dp), allocatable :: out(:)
        integer :: i

        allocate(out(size(newdata)))
        do i = 1, size(newdata)
            out(i) = predict_category_scalar(object, newdata(i), missing)
        end do
    end function predict_category_vector

    function predict_overall_numeric(object, newdata) result(out)
        type(d_overall_type), intent(in) :: object
        real(dp), intent(in) :: newdata(:, :)
        real(dp), allocatable :: out(:)
        real(dp), allocatable :: all_values(:, :)
        integer :: nfun

        nfun = object%size()
        if (nfun == 0) error stop "predict: d_overall contains no functions"
        if (size(newdata, 2) /= nfun) then
            error stop "predict: newdata columns must match the number of functions"
        end if

        all_values = predict_all_numeric(object, newdata)
        allocate(out(size(newdata, 1)))
        out = all_values(:, nfun + 1)
    end function predict_overall_numeric

    function predict_all_numeric(object, newdata) result(out)
        type(d_overall_type), intent(in) :: object
        real(dp), intent(in) :: newdata(:, :)
        real(dp), allocatable :: out(:, :)
        integer :: i, j, nfun

        nfun = object%size()
        if (nfun == 0) error stop "predict_all: d_overall contains no functions"
        if (size(newdata, 2) /= nfun) then
            error stop "predict_all: newdata columns must match the number of functions"
        end if

        allocate(out(size(newdata, 1), nfun + 1))
        do j = 1, nfun
            select type (f => object%d(j)%fun)
            type is (d_categorical_type)
                error stop "predict_all: categorical function requires mixed input"
            class default
                do i = 1, size(newdata, 1)
                    out(i, j) = f%value(numeric_input(newdata(i, j)))
                end do
            end select
        end do
        call fill_overall_column(out, nfun)
    end function predict_all_numeric

    function predict_overall_mixed(object, newdata) result(out)
        type(d_overall_type), intent(in) :: object
        type(desirability_input), intent(in) :: newdata(:, :)
        real(dp), allocatable :: out(:)
        real(dp), allocatable :: all_values(:, :)
        integer :: nfun

        nfun = object%size()
        if (nfun == 0) error stop "predict: d_overall contains no functions"
        if (size(newdata, 2) /= nfun) then
            error stop "predict: newdata columns must match the number of functions"
        end if

        all_values = predict_all_mixed(object, newdata)
        allocate(out(size(newdata, 1)))
        out = all_values(:, nfun + 1)
    end function predict_overall_mixed

    function predict_all_mixed(object, newdata) result(out)
        type(d_overall_type), intent(in) :: object
        type(desirability_input), intent(in) :: newdata(:, :)
        real(dp), allocatable :: out(:, :)
        integer :: i, j, nfun

        nfun = object%size()
        if (nfun == 0) error stop "predict_all: d_overall contains no functions"
        if (size(newdata, 2) /= nfun) then
            error stop "predict_all: newdata columns must match the number of functions"
        end if

        allocate(out(size(newdata, 1), nfun + 1))
        do j = 1, nfun
            do i = 1, size(newdata, 1)
                out(i, j) = object%d(j)%fun%value(newdata(i, j))
            end do
        end do
        call fill_overall_column(out, nfun)
    end function predict_all_mixed

    subroutine fill_overall_column(values, nfun)
        real(dp), intent(inout) :: values(:, :)
        integer, intent(in) :: nfun
        integer :: i

        do i = 1, size(values, 1)
            values(i, nfun + 1) = product(values(i, 1:nfun)) ** &
                (1.0_dp / real(nfun, dp))
        end do
    end subroutine fill_overall_column

    function d_max_raw(self, input) result(out)
        class(d_max_type), intent(in) :: self
        type(desirability_input), intent(in) :: input
        real(dp) :: out, x

        call require_numeric(input, "d_max")
        x = input%number
        if (x < self%low) then
            out = 0.0_dp
        else if (x > self%high) then
            out = 1.0_dp
        else
            out = ((x - self%low) / (self%high - self%low)) ** self%scale
        end if
    end function d_max_raw

    function d_min_raw(self, input) result(out)
        class(d_min_type), intent(in) :: self
        type(desirability_input), intent(in) :: input
        real(dp) :: out, x

        call require_numeric(input, "d_min")
        x = input%number
        if (x < self%low) then
            out = 1.0_dp
        else if (x > self%high) then
            out = 0.0_dp
        else
            out = ((x - self%high) / (self%low - self%high)) ** self%scale
        end if
    end function d_min_raw

    function d_target_raw(self, input) result(out)
        class(d_target_type), intent(in) :: self
        type(desirability_input), intent(in) :: input
        real(dp) :: out, x

        call require_numeric(input, "d_target")
        x = input%number
        if (x < self%low .or. x > self%high) then
            out = 0.0_dp
        else if (x <= self%target) then
            out = ((x - self%low) / (self%target - self%low)) ** self%low_scale
        else
            out = ((x - self%high) / (self%target - self%high)) ** self%high_scale
        end if
    end function d_target_raw

    function d_box_raw(self, input) result(out)
        class(d_box_type), intent(in) :: self
        type(desirability_input), intent(in) :: input
        real(dp) :: out, x

        call require_numeric(input, "d_box")
        x = input%number
        if (x < self%low .or. x > self%high) then
            out = 0.0_dp
        else
            out = 1.0_dp
        end if
    end function d_box_raw

    function d_arb_raw(self, input) result(out)
        class(d_arb_type), intent(in) :: self
        type(desirability_input), intent(in) :: input
        real(dp) :: out, x, weight
        integer :: lo, hi, mid, n

        call require_numeric(input, "d_arb")
        x = input%number
        n = size(self%interp_x)

        if (x < self%x(1)) then
            out = self%d(1)
            return
        else if (x > self%x(size(self%x))) then
            out = self%d(size(self%d))
            return
        end if

        if (x <= self%interp_x(1)) then
            out = self%interp_d(1)
            return
        else if (x >= self%interp_x(n)) then
            out = self%interp_d(n)
            return
        end if

        lo = 1
        hi = n
        do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (x >= self%interp_x(mid)) then
                lo = mid
            else
                hi = mid
            end if
        end do

        weight = (x - self%interp_x(lo)) / &
            (self%interp_x(hi) - self%interp_x(lo))
        out = self%interp_d(lo) + weight * &
            (self%interp_d(hi) - self%interp_d(lo))
    end function d_arb_raw

    function d_categorical_raw(self, input) result(out)
        class(d_categorical_type), intent(in) :: self
        type(desirability_input), intent(in) :: input
        real(dp) :: out
        integer :: i

        if (input%kind /= input_categorical) then
            error stop "d_categorical: categorical input required"
        end if

        do i = 1, size(self%labels)
            if (trim(input%text) == trim(self%labels(i))) then
                out = self%values(i)
                return
            end if
        end do
        error stop "d_categorical: input label is not defined"
    end function d_categorical_raw

    subroutine require_numeric(input, caller)
        type(desirability_input), intent(in) :: input
        character(len=*), intent(in) :: caller

        if (input%kind /= input_numeric) then
            error stop trim(caller) // ": numeric input required"
        end if
    end subroutine require_numeric

    subroutine set_tolerance(obj, tol)
        class(desirability_function), intent(inout) :: obj
        real(dp), intent(in), optional :: tol

        if (present(tol)) then
            obj%has_tol = .true.
            obj%tol = tol
        else
            obj%has_tol = .false.
            obj%tol = 0.0_dp
        end if
    end subroutine set_tolerance

    function mean_on_grid(obj, low, high) result(mean_value)
        class(desirability_function), intent(in) :: obj
        real(dp), intent(in) :: low, high
        real(dp) :: mean_value, x
        integer :: i

        mean_value = 0.0_dp
        do i = 1, n_noninformative
            x = low + real(i - 1, dp) * (high - low) / &
                real(n_noninformative - 1, dp)
            mean_value = mean_value + obj%raw_value(numeric_input(x))
        end do
        mean_value = mean_value / real(n_noninformative, dp)
    end function mean_on_grid

    subroutine stable_sort_pairs(x, y)
        real(dp), intent(inout) :: x(:), y(:)
        integer :: i, j
        real(dp) :: key_x, key_y

        do i = 2, size(x)
            key_x = x(i)
            key_y = y(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= key_x) exit
                x(j + 1) = x(j)
                y(j + 1) = y(j)
                j = j - 1
            end do
            x(j + 1) = key_x
            y(j + 1) = key_y
        end do
    end subroutine stable_sort_pairs

    subroutine build_interpolation_knots(x, d, knot_x, knot_d)
        real(dp), intent(in) :: x(:), d(:)
        real(dp), allocatable, intent(out) :: knot_x(:), knot_d(:)
        real(dp), allocatable :: temp_x(:), temp_d(:)
        real(dp) :: sum_d
        integer :: i, j, n_unique, count

        allocate(temp_x(size(x)), temp_d(size(d)))
        n_unique = 0
        i = 1
        do while (i <= size(x))
            j = i
            sum_d = 0.0_dp
            count = 0
            do while (j <= size(x))
                if (x(j) /= x(i)) exit
                sum_d = sum_d + d(j)
                count = count + 1
                j = j + 1
            end do
            n_unique = n_unique + 1
            temp_x(n_unique) = x(i)
            temp_d(n_unique) = sum_d / real(count, dp)
            i = j
        end do

        allocate(knot_x(n_unique), knot_d(n_unique))
        knot_x = temp_x(1:n_unique)
        knot_d = temp_d(1:n_unique)
    end subroutine build_interpolation_knots

end module desirability

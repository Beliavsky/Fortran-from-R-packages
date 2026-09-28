module car_kernels
    use car_kinds, only: dp
    use car_utils, only: mean_valid, quantile_type7, sample_variance
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_quiet_nan, ieee_value
    implicit none
    private

    public :: adaptive_kernel
    public :: biweight_kernel
    public :: bw_nrd0
    public :: epanechnikov_kernel
    public :: normal_kernel

contains

    pure elemental function normal_kernel(x) result(value)
        real(dp), intent(in) :: x !! Standardized kernel coordinate.
        real(dp) :: value
        real(dp), parameter :: inv_sqrt_two_pi = 0.39894228040143267794_dp

        value = inv_sqrt_two_pi * exp(-0.5_dp * x * x)
    end function normal_kernel

    pure elemental function epanechnikov_kernel(x) result(value)
        real(dp), intent(in) :: x !! Standardized kernel coordinate; support is [-1, 1].
        real(dp) :: value

        if (abs(x) > 1.0_dp) then
            value = 0.0_dp
        else
            value = 0.75_dp * (1.0_dp - x * x)
        end if
    end function epanechnikov_kernel

    pure elemental function biweight_kernel(x) result(value)
        real(dp), intent(in) :: x !! Standardized kernel coordinate; support is [-1, 1].
        real(dp) :: value

        if (abs(x) > 1.0_dp) then
            value = 0.0_dp
        else
            value = (15.0_dp / 16.0_dp) * (1.0_dp - x * x)**2
        end if
    end function biweight_kernel

    pure subroutine bw_nrd0(x, bandwidth, status)
        real(dp), intent(in) :: x(:) !! Numeric sample; NaNs are omitted, and at least two valid observations are required.
        real(dp), intent(out) :: bandwidth !! R-compatible normal-reference bandwidth from stats::bw.nrd0.
        integer, intent(out) :: status !! Zero on success; nonzero when fewer than two valid observations exist.
        real(dp) :: hi
        real(dp) :: iqr
        real(dp) :: lo
        real(dp) :: q25
        real(dp) :: q75
        real(dp) :: variance
        integer :: n
        integer :: quantile_status
        integer :: variance_status

        call sample_variance(x, variance, n, variance_status)
        if (variance_status /= 0 .or. n < 2) then
            bandwidth = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 1
            return
        end if

        hi = sqrt(max(0.0_dp, variance))
        call quantile_type7(x, 0.25_dp, q25, quantile_status)
        if (quantile_status /= 0) then
            bandwidth = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 1
            return
        end if
        call quantile_type7(x, 0.75_dp, q75, quantile_status)
        if (quantile_status /= 0) then
            bandwidth = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 1
            return
        end if
        iqr = q75 - q25
        lo = min(hi, iqr / 1.34_dp)
        if (lo <= 0.0_dp) lo = hi
        if (lo <= 0.0_dp) then
            lo = maxval(abs(x), mask=.not. ieee_is_nan(x))
        end if
        if (lo <= 0.0_dp) lo = 1.0_dp

        bandwidth = 0.9_dp * lo * real(n, dp)**(-0.2_dp)
        status = 0
    end subroutine bw_nrd0

    pure function kernel_from_id(kernel_id, x) result(value)
        integer, intent(in) :: kernel_id !! Kernel selector: 1 normal, 2 Epanechnikov, 3 biweight.
        real(dp), intent(in) :: x !! Standardized kernel coordinate.
        real(dp) :: value

        select case (kernel_id)
        case (1)
            value = normal_kernel(x)
        case (2)
            value = epanechnikov_kernel(x)
        case (3)
            value = biweight_kernel(x)
        case default
            value = ieee_value(0.0_dp, ieee_quiet_nan)
        end select
    end function kernel_from_id

    subroutine adaptive_kernel(x, grid, density, initial_density, bandwidth_used, status, kernel, bandwidth, &
            adjust, n_points, from, to, cut, na_rm)
        real(dp), intent(in) :: x(:) !! Input sample for adaptive density estimation.
        real(dp), allocatable, intent(out) :: grid(:) !! Equally spaced evaluation grid.
        real(dp), allocatable, intent(out) :: density(:) !! Adaptive-kernel density estimates on grid.
        real(dp), allocatable, intent(out) :: initial_density(:) !! Pilot fixed-bandwidth density estimates on grid.
        real(dp), intent(out) :: bandwidth_used !! Actual fixed pilot bandwidth after kernel scaling and adjustment.
        integer, intent(out) :: status !! Zero on success; nonzero for invalid input, kernel, bandwidth, or pilot density.
        character(len=*), intent(in), optional :: kernel !! Kernel name: "normal", "epanechnikov"/"depan", or "biweight"/"dbiwt"; default normal.
        real(dp), intent(in), optional :: bandwidth !! Numeric pilot bandwidth; omitted uses bw.nrd0 and car kernel-scale correction.
        real(dp), intent(in), optional :: adjust !! Positive multiplicative bandwidth adjustment; default one.
        integer, intent(in), optional :: n_points !! Number of grid points; default 500 and must be at least two.
        real(dp), intent(in), optional :: from !! Lower grid endpoint; default min(x) - cut*bandwidth.
        real(dp), intent(in), optional :: to !! Upper grid endpoint; default max(x) + cut*bandwidth.
        real(dp), intent(in), optional :: cut !! Number of bandwidths added beyond the sample range; default three.
        logical, intent(in), optional :: na_rm !! If true omit NaNs; default true. If false any NaN is an error.
        real(dp), allocatable :: clean(:)
        real(dp), allocatable :: local_density(:)
        real(dp), allocatable :: pilot_at_data(:)
        real(dp), allocatable :: factors(:)
        real(dp), allocatable :: log_pilot(:)
        real(dp) :: adjust_value
        real(dp) :: bw
        real(dp) :: cut_value
        real(dp) :: fac
        real(dp) :: from_value
        real(dp) :: kernel_scale
        real(dp) :: pbar
        real(dp) :: to_value
        real(dp) :: distance
        real(dp) :: best_distance
        logical :: remove_na
        integer :: bandwidth_status
        integer :: i
        integer :: j
        integer :: kernel_id
        integer :: mean_status
        integer :: n
        integer :: n_grid
        integer :: n_valid
        integer :: nearest

        remove_na = .true.
        if (present(na_rm)) remove_na = na_rm
        if (.not. remove_na .and. any(ieee_is_nan(x))) then
            allocate(grid(0), density(0), initial_density(0))
            bandwidth_used = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 1
            return
        end if

        n = count(.not. ieee_is_nan(x))
        if (n < 2) then
            allocate(grid(0), density(0), initial_density(0))
            bandwidth_used = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 1
            return
        end if
        allocate(clean(n))
        clean = pack(x, .not. ieee_is_nan(x))

        kernel_id = 1
        kernel_scale = 1.0_dp
        if (present(kernel)) then
            select case (trim(adjustl(kernel)))
            case ("normal", "dnorm", "gaussian")
                kernel_id = 1
                kernel_scale = 1.0_dp
            case ("epanechnikov", "depan")
                kernel_id = 2
                kernel_scale = sqrt(0.2_dp)
            case ("biweight", "dbiwt")
                kernel_id = 3
                kernel_scale = sqrt(1.0_dp / 7.0_dp)
            case default
                allocate(grid(0), density(0), initial_density(0))
                bandwidth_used = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 2
                return
            end select
        end if

        if (present(bandwidth)) then
            bw = bandwidth
        else
            call bw_nrd0(clean, bw, bandwidth_status)
            if (bandwidth_status /= 0) then
                allocate(grid(0), density(0), initial_density(0))
                bandwidth_used = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 3
                return
            end if
            bw = bw / kernel_scale
        end if

        adjust_value = 1.0_dp
        if (present(adjust)) adjust_value = adjust
        if (bw <= 0.0_dp .or. adjust_value <= 0.0_dp) then
            allocate(grid(0), density(0), initial_density(0))
            bandwidth_used = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 3
            return
        end if
        bw = bw * adjust_value
        bandwidth_used = bw

        n_grid = 500
        if (present(n_points)) n_grid = n_points
        if (n_grid < 2) then
            allocate(grid(0), density(0), initial_density(0))
            status = 4
            return
        end if

        cut_value = 3.0_dp
        if (present(cut)) cut_value = cut
        from_value = minval(clean) - cut_value * bw
        to_value = maxval(clean) + cut_value * bw
        if (present(from)) from_value = from
        if (present(to)) to_value = to
        if (to_value <= from_value) then
            allocate(grid(0), density(0), initial_density(0))
            status = 4
            return
        end if

        allocate(grid(n_grid), density(n_grid), initial_density(n_grid))
        allocate(pilot_at_data(n), factors(n), log_pilot(n), local_density(n_grid))
        do i = 1, n_grid
            grid(i) = from_value + real(i - 1, dp) * (to_value - from_value) / real(n_grid - 1, dp)
        end do

        fac = 1.0_dp / (real(n, dp) * bw)
        do i = 1, n_grid
            initial_density(i) = 0.0_dp
            do j = 1, n
                initial_density(i) = initial_density(i) + kernel_from_id(kernel_id, (clean(j) - grid(i)) / bw)
            end do
            initial_density(i) = fac * initial_density(i)
        end do

        do i = 1, n
            nearest = 1
            best_distance = abs(clean(i) - grid(1))
            do j = 2, n_grid
                distance = abs(clean(i) - grid(j))
                if (distance < best_distance) then
                    best_distance = distance
                    nearest = j
                end if
            end do
            pilot_at_data(i) = initial_density(nearest)
            if (pilot_at_data(i) <= 0.0_dp) then
                density = ieee_value(0.0_dp, ieee_quiet_nan)
                status = 5
                return
            end if
            log_pilot(i) = log(pilot_at_data(i))
        end do

        call mean_valid(log_pilot, pbar, n_valid, mean_status)
        if (mean_status /= 0) then
            density = ieee_value(0.0_dp, ieee_quiet_nan)
            status = 5
            return
        end if
        pbar = exp(pbar)
        factors = (pilot_at_data / pbar)**(-0.5_dp)

        do i = 1, n_grid
            local_density(i) = 0.0_dp
            do j = 1, n
                local_density(i) = local_density(i) + &
                    kernel_from_id(kernel_id, (clean(j) - grid(i)) / (factors(j) * bw)) / factors(j)
            end do
            density(i) = fac * local_density(i)
        end do
        status = 0
    end subroutine adaptive_kernel

end module car_kernels

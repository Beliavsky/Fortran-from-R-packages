module rangen_sampling
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use rangen_kinds, only : dp, i8
    use rangen_pcg32, only : pcg32_state
    implicit none
    private

    type(pcg32_state), save :: sample_rng
    logical, save :: initialized = .false.

    public :: set_sampling_seed, sample_int, sample_real, col_sample, row_sample, nano_time

contains

    subroutine ensure_initialized()
        integer :: count
        if (initialized) return
        call system_clock(count)
        call set_sampling_seed(int(count, i8))
    end subroutine ensure_initialized

    subroutine set_sampling_seed(seed_value)
        integer(i8), intent(in) :: seed_value
        call sample_rng%seed(seed_value, 3935559000370003845_i8)
        initialized = .true.
    end subroutine set_sampling_seed

    function draw_index(n) result(j)
        integer, intent(in) :: n
        integer :: j
        integer(i8) :: u
        call ensure_initialized()
        if (n <= 0) then
            j = 0
            return
        end if
        u = sample_rng%next_uint32()
        j = 1 + int(modulo(u, int(n, i8)))
    end function draw_index

    function sample_int(n, sample_size, replace) result(out)
        integer, intent(in) :: n
        integer, intent(in), optional :: sample_size
        logical, intent(in), optional :: replace
        integer, allocatable :: out(:)
        integer, allocatable :: pool(:)
        integer :: m, i, j, tmp
        logical :: rep

        m = n
        if (present(sample_size)) m = sample_size
        rep = .false.
        if (present(replace)) rep = replace
        if (n < 0 .or. m < 0 .or. (.not. rep .and. m > n)) then
            allocate(out(0))
            return
        end if
        allocate(out(m))
        if (m == 0) return
        if (rep) then
            do i = 1, m
                out(i) = draw_index(n)
            end do
        else
            allocate(pool(n))
            pool = [(i, i=1,n)]
            do i = 1, m
                j = i - 1 + draw_index(n - i + 1)
                tmp = pool(i)
                pool(i) = pool(j)
                pool(j) = tmp
                out(i) = pool(i)
            end do
        end if
    end function sample_int

    function sample_real(x, sample_size, replace) result(out)
        real(dp), intent(in) :: x(:)
        integer, intent(in), optional :: sample_size
        logical, intent(in), optional :: replace
        real(dp), allocatable :: out(:)
        integer, allocatable :: idx(:)
        integer :: m
        logical :: rep

        m = size(x)
        if (present(sample_size)) m = sample_size
        rep = .false.
        if (present(replace)) rep = replace
        idx = sample_int(size(x), m, rep)
        allocate(out(size(idx)))
        if (size(idx) > 0) out = x(idx)
    end function sample_real

    function col_sample(x, sizes, replace, parallel, cores) result(out)
        real(dp), intent(in) :: x(:,:)
        integer, intent(in) :: sizes(:)
        logical, intent(in) :: replace(:)
        logical, intent(in), optional :: parallel
        integer, intent(in), optional :: cores
        real(dp), allocatable :: out(:,:)
        real(dp), allocatable :: tmp(:)
        real(dp) :: nanv
        integer :: j, m

        if (size(sizes) == 0 .or. size(replace) == 0) then
            allocate(out(0, size(x,2)))
            return
        end if
        m = maxval(sizes)
        allocate(out(max(m,0), size(x,2)))
        nanv = ieee_value(0.0_dp, ieee_quiet_nan)
        out = nanv
        if (present(parallel)) then
            if (parallel) continue
        end if
        if (present(cores)) then
            if (cores < 0) continue
        end if
        do j = 1, size(x,2)
            tmp = sample_real(x(:,j), sizes(1 + modulo(j - 1, size(sizes))), &
                replace(1 + modulo(j - 1, size(replace))))
            if (size(tmp) > 0) out(1:size(tmp),j) = tmp
        end do
    end function col_sample

    function row_sample(x, sizes, replace, parallel, cores) result(out)
        real(dp), intent(in) :: x(:,:)
        integer, intent(in) :: sizes(:)
        logical, intent(in) :: replace(:)
        logical, intent(in), optional :: parallel
        integer, intent(in), optional :: cores
        real(dp), allocatable :: out(:,:)
        real(dp), allocatable :: tmp(:)
        real(dp) :: nanv
        integer :: i, m

        if (size(sizes) == 0 .or. size(replace) == 0) then
            allocate(out(size(x,1), 0))
            return
        end if
        m = maxval(sizes)
        allocate(out(size(x,1), max(m,0)))
        nanv = ieee_value(0.0_dp, ieee_quiet_nan)
        out = nanv
        if (present(parallel)) then
            if (parallel) continue
        end if
        if (present(cores)) then
            if (cores < 0) continue
        end if
        do i = 1, size(x,1)
            tmp = sample_real(x(i,:), sizes(1 + modulo(i - 1, size(sizes))), &
                replace(1 + modulo(i - 1, size(replace))))
            if (size(tmp) > 0) out(i,1:size(tmp)) = tmp
        end do
    end function row_sample

    function nano_time() result(ns)
        real(dp) :: ns
        integer(i8) :: count, rate
        call system_clock(count, rate)
        if (rate > 0_i8) then
            ns = real(count, dp) * (1.0e9_dp / real(rate, dp))
        else
            ns = 0.0_dp
        end if
    end function nano_time

end module rangen_sampling

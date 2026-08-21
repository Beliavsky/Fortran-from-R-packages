! SPDX-License-Identifier: GPL-2.0-or-later
module goftest_testing
    use goftest_kinds, only : dp
    use goftest_utils, only : shuffle_int
    use goftest_ad, only : ad_test_uniform
    use goftest_cvm, only : cvm_test_uniform
    implicit none
    private

    type, public :: gof_result
        real(dp) :: statistic = 0.0_dp
        real(dp) :: p_value = 1.0_dp
        integer :: groups = 1
        integer :: status = 0
    end type gof_result

    abstract interface
        function cdf_callback(x) result(p)
            import dp
            real(dp), intent(in) :: x
            real(dp) :: p
        end function cdf_callback
    end interface

    public :: ad_test
    public :: cvm_test
    public :: ad_test_values
    public :: cvm_test_values

contains

    function ad_test(x, cdf, estimated, groups) result(res)
        real(dp), intent(in) :: x(:)
        procedure(cdf_callback) :: cdf
        logical, intent(in), optional :: estimated
        integer, intent(in), optional :: groups(:)
        type(gof_result) :: res
        real(dp), allocatable :: u(:)
        integer :: i

        allocate(u(size(x)))
        do i = 1, size(x)
            u(i) = cdf(x(i))
        end do
        res = ad_test_values(u, estimated, groups)
    end function ad_test

    function cvm_test(x, cdf, estimated, groups) result(res)
        real(dp), intent(in) :: x(:)
        procedure(cdf_callback) :: cdf
        logical, intent(in), optional :: estimated
        integer, intent(in), optional :: groups(:)
        type(gof_result) :: res
        real(dp), allocatable :: u(:)
        integer :: i

        allocate(u(size(x)))
        do i = 1, size(x)
            u(i) = cdf(x(i))
        end do
        res = cvm_test_values(u, estimated, groups)
    end function cvm_test

    function ad_test_values(u, estimated, groups) result(res)
        real(dp), intent(in) :: u(:)
        logical, intent(in), optional :: estimated
        integer, intent(in), optional :: groups(:)
        type(gof_result) :: res
        logical :: est

        if (size(u) == 0 .or. any(u < 0.0_dp) .or. any(u > 1.0_dp)) then
            res%status = 3
            return
        end if
        est = .false.
        if (present(estimated)) est = estimated
        if (.not. est .or. size(u) <= 4) then
            call ad_test_uniform(u, res%statistic, res%p_value)
            return
        end if
        res = braun_test(u, .true., groups)
    end function ad_test_values

    function cvm_test_values(u, estimated, groups) result(res)
        real(dp), intent(in) :: u(:)
        logical, intent(in), optional :: estimated
        integer, intent(in), optional :: groups(:)
        type(gof_result) :: res
        logical :: est

        if (size(u) == 0 .or. any(u < 0.0_dp) .or. any(u > 1.0_dp)) then
            res%status = 3
            return
        end if
        est = .false.
        if (present(estimated)) est = estimated
        if (.not. est .or. size(u) <= 4) then
            call cvm_test_uniform(u, res%statistic, res%p_value)
            return
        end if
        res = braun_test(u, .false., groups)
    end function cvm_test_values

    function braun_test(u, use_ad, supplied_groups) result(res)
        real(dp), intent(in) :: u(:)
        logical, intent(in) :: use_ad
        integer, intent(in), optional :: supplied_groups(:)
        type(gof_result) :: res
        integer, allocatable :: g(:)
        real(dp), allocatable :: sub(:)
        real(dp) :: stat, pval, minp, maxstat
        integer :: n, m, i, j, k, nsub

        n = size(u)
        m = nint(sqrt(real(n, dp)))
        res%groups = m
        if (n < 2 * m) then
            res%status = 1
            return
        end if
        allocate(g(n))
        if (present(supplied_groups)) then
            if (size(supplied_groups) /= n) then
                res%status = 2
                return
            end if
            g = supplied_groups
            m = maxval(g)
            res%groups = m
        else
            do i = 1, n
                g(i) = 1 + mod(i - 1, m)
            end do
            call shuffle_int(g)
        end if

        minp = 1.0_dp
        maxstat = -huge(1.0_dp)
        do j = 1, m
            nsub = count(g == j)
            if (nsub == 0) cycle
            allocate(sub(nsub))
            k = 0
            do i = 1, n
                if (g(i) == j) then
                    k = k + 1
                    sub(k) = u(i)
                end if
            end do
            if (use_ad) then
                call ad_test_uniform(sub, stat, pval)
            else
                call cvm_test_uniform(sub, stat, pval)
            end if
            maxstat = max(maxstat, stat)
            minp = min(minp, pval)
            deallocate(sub)
        end do
        res%statistic = maxstat
        res%p_value = 1.0_dp - (1.0_dp - minp)**real(m, dp)
    end function braun_test

end module goftest_testing

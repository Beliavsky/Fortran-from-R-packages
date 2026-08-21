! Computational translation of the R package bivpois.
! License: GPL-2.0-or-later.
module bivpois_gof
    use bivpois_kinds, only : dp
    use bivpois_math, only : mean_counts, variance_counts, correlation_counts
    use bivpois_distribution, only : rbp
    use bivpois_fit, only : bp_mle2, bp_mle2_result
    implicit none
    private

    type, public :: bp_gof_result
        real(dp) :: statistic = 0.0_dp
        real(dp) :: pvalue = 1.0_dp
        real(dp) :: lambda(3) = 0.0_dp
        integer :: replicates = 0
    end type bp_gof_result

    public :: bp_dispersion_statistic, bp_gof, bp_gof2

contains

    real(dp) function bp_dispersion_statistic(x1, x2) result(stat)
        integer, intent(in) :: x1(:), x2(:)
        real(dp) :: r, m1, m2, v1, v2, den

        if (size(x1) < 2 .or. size(x2) /= size(x1)) error stop "bp_dispersion_statistic: invalid input"
        m1 = mean_counts(x1)
        m2 = mean_counts(x2)
        v1 = variance_counts(x1)
        v2 = variance_counts(x2)
        r = correlation_counts(x1, x2)
        den = 1.0_dp - r * r
        if (m1 <= 0.0_dp .or. m2 <= 0.0_dp .or. den <= tiny(1.0_dp)) then
            stat = huge(1.0_dp)
            return
        end if
        stat = real(size(x1), dp) / den * (v1 / m1 + v2 / m2 - &
               2.0_dp * r * r * sqrt((v1 / m1) * (v2 / m2)))
    end function bp_dispersion_statistic

    function bp_gof(x1, x2, r_replicates) result(res)
        integer, intent(in) :: x1(:), x2(:)
        integer, intent(in), optional :: r_replicates
        type(bp_gof_result) :: res
        type(bp_mle2_result) :: fit
        integer :: n, r, i, exceed
        integer, allocatable :: sim(:, :)
        real(dp) :: stat_i

        if (size(x1) < 2 .or. size(x2) /= size(x1)) error stop "bp_gof: invalid input"
        r = 999
        if (present(r_replicates)) r = r_replicates
        if (r < 1) error stop "bp_gof: replicates must be positive"
        n = size(x1)
        fit = bp_mle2(x1, x2)
        res%lambda = fit%lambda
        res%statistic = bp_dispersion_statistic(x1, x2)
        res%replicates = r
        exceed = 0
        allocate(sim(n, 2))
        do i = 1, r
            call rbp(n, fit%lambda, sim)
            stat_i = bp_dispersion_statistic(sim(:, 1), sim(:, 2))
            if (stat_i > res%statistic) exceed = exceed + 1
        end do
        res%pvalue = real(exceed + 1, dp) / real(r + 1, dp)
    end function bp_gof

    function bp_gof2(x1, x2, r_replicates) result(res)
        integer, intent(in) :: x1(:), x2(:)
        integer, intent(in), optional :: r_replicates
        type(bp_gof_result) :: res

        ! R's bp.gof2 is a vectorized implementation of the same Monte Carlo test.
        ! Fortran's compiled loop already avoids interpreter overhead, so both APIs share
        ! the same memory-efficient numerical implementation.
        if (present(r_replicates)) then
            res = bp_gof(x1, x2, r_replicates)
        else
            res = bp_gof(x1, x2)
        end if
    end function bp_gof2

end module bivpois_gof

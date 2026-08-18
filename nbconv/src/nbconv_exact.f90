! SPDX-License-Identifier: GPL-3.0-or-later
! Furman-series translation from the R package nbconv.
module nbconv_exact
    use nbconv_kinds, only : dp
    use nbconv_math, only : log_add_exp
    implicit none
    private

    real(dp), parameter :: neg_huge = -huge(1.0_dp)

    public :: nb_sum_exact

contains

    function nb_sum_exact(ps, phis, counts, n_terms, tolerance, enforce_tolerance, k_mass) result(pmf)
        real(dp), intent(in) :: ps(:), phis(:)
        integer, intent(in) :: counts(:)
        integer, intent(in), optional :: n_terms
        real(dp), intent(in), optional :: tolerance
        logical, intent(in), optional :: enforce_tolerance
        real(dp), intent(out), optional :: k_mass
        real(dp), allocatable :: pmf(:)

        integer :: nt, j, i, k, x, na, ia
        real(dp) :: tol, pmax, qmax, phisum, log_r, a, term
        real(dp) :: log_ksum, log_mass, lgx
        logical :: enforce
        integer, allocatable :: active(:)
        real(dp), allocatable :: pa(:), qa(:), phia(:), loga(:), xi(:), delta(:)

        call validate_inputs(ps, phis)
        if (size(counts) == 0) then
            allocate(pmf(0))
            if (present(k_mass)) k_mass = 1.0_dp
            return
        end if
        if (any(counts < 0)) error stop "nb_sum_exact: counts must be nonnegative"

        nt = 1000
        if (present(n_terms)) nt = n_terms
        if (nt < 1) error stop "nb_sum_exact: n_terms must be positive"
        tol = 1.0e-3_dp
        if (present(tolerance)) tol = tolerance
        if (tol <= 0.0_dp) error stop "nb_sum_exact: tolerance must be positive"
        enforce = .true.
        if (present(enforce_tolerance)) enforce = enforce_tolerance

        na = count(ps < 1.0_dp)
        allocate(pmf(size(counts)))
        if (na == 0) then
            pmf = 0.0_dp
            do j = 1, size(counts)
                if (counts(j) == 0) pmf(j) = 1.0_dp
            end do
            if (present(k_mass)) k_mass = 1.0_dp
            return
        end if

        allocate(active(na), pa(na), qa(na), phia(na), loga(na))
        ia = 0
        do j = 1, size(ps)
            if (ps(j) < 1.0_dp) then
                ia = ia + 1
                active(ia) = j
                pa(ia) = ps(j)
                qa(ia) = 1.0_dp - ps(j)
                phia(ia) = phis(j)
            end if
        end do

        pmax = maxval(pa)
        qmax = 1.0_dp - pmax
        phisum = sum(phia)
        log_r = 0.0_dp
        do j = 1, na
            log_r = log_r - phia(j) * (log(qa(j) * pmax) - log(pa(j) * qmax))
            a = 1.0_dp - qmax * pa(j) / (qa(j) * pmax)
            if (a > 0.0_dp) then
                loga(j) = log(a)
            else
                loga(j) = neg_huge
            end if
        end do

        allocate(xi(nt), delta(nt))
        do i = 1, nt
            xi(i) = neg_huge
            do j = 1, na
                if (loga(j) > neg_huge / 2.0_dp) then
                    term = log(phia(j)) + real(i, dp) * loga(j) - log(real(i, dp))
                    xi(i) = log_add_exp(xi(i), term)
                end if
            end do
        end do

        delta = neg_huge
        delta(1) = 0.0_dp
        do k = 1, nt - 1
            term = neg_huge
            do i = 1, k
                if (xi(i) > neg_huge / 2.0_dp .and. delta(k + 1 - i) > neg_huge / 2.0_dp) then
                    term = log_add_exp(term, log(real(i, dp)) + xi(i) + delta(k + 1 - i))
                end if
            end do
            if (term > neg_huge / 2.0_dp) delta(k + 1) = term - log(real(k, dp))
        end do

        log_ksum = neg_huge
        do k = 1, nt
            if (delta(k) > neg_huge / 2.0_dp) log_ksum = log_add_exp(log_ksum, log_r + delta(k))
        end do
        if (present(k_mass)) k_mass = exp(log_ksum)
        if (enforce) then
            if (abs(exp(log_ksum) - 1.0_dp) > tol) then
                error stop "nb_sum_exact: K-series mass is outside tolerance; increase n_terms"
            end if
        end if

        do j = 1, size(counts)
            x = counts(j)
            lgx = log_gamma(real(x + 1, dp))
            log_mass = neg_huge
            do k = 0, nt - 1
                if (delta(k + 1) <= neg_huge / 2.0_dp) cycle
                term = delta(k + 1) &
                     + log_gamma(phisum + real(x + k, dp)) &
                     - log_gamma(phisum + real(k, dp)) - lgx &
                     + (phisum + real(k, dp)) * log(pmax)
                if (x > 0) term = term + real(x, dp) * log(qmax)
                log_mass = log_add_exp(log_mass, term)
            end do
            if (log_mass <= neg_huge / 2.0_dp) then
                pmf(j) = 0.0_dp
            else
                pmf(j) = exp(log_r + log_mass)
            end if
        end do
    end function nb_sum_exact

    subroutine validate_inputs(ps, phis)
        real(dp), intent(in) :: ps(:), phis(:)

        if (size(ps) /= size(phis)) error stop "nb_sum_exact: ps and phis must have equal length"
        if (size(ps) < 1) error stop "nb_sum_exact: parameter arrays must be nonempty"
        if (any(ps <= 0.0_dp) .or. any(ps > 1.0_dp)) then
            error stop "nb_sum_exact: ps must satisfy 0 < p <= 1"
        end if
        if (any(phis <= 0.0_dp)) error stop "nb_sum_exact: phis must be positive"
    end subroutine validate_inputs

end module nbconv_exact

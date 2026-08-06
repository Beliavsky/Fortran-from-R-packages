! SPDX-License-Identifier: GPL-3.0-only
module rsdc_viterbi
    use rsdc_kinds, only: dp, rsdc_pi
    use rsdc_types, only: rsdc_model, rsdc_const, rsdc_nox, rsdc_tvtp
    use rsdc_parameters, only: correlations_to_matrix, transition_at
    use rsdc_linalg, only: inverse_spd, logdet_spd, is_positive_definite
    implicit none
    private
    public :: rsdc_viterbi_path

contains

    subroutine rsdc_viterbi_path(model, residuals, path, x, ok)
        type(rsdc_model), intent(in) :: model
        real(dp), intent(in) :: residuals(:, :)
        integer, allocatable, intent(out) :: path(:)
        real(dp), intent(in), optional :: x(:, :)
        logical, intent(out), optional :: ok
        real(dp), allocatable :: logdens(:, :), delta(:, :), pmat(:, :), inv(:, :), r(:, :), pi1(:), vals(:)
        integer, allocatable :: psi(:, :)
        real(dp) :: ld, quad, best
        integer :: t, s, i, best_i, nobs, n, k
        logical :: good
        nobs = size(residuals, 1); k = size(residuals, 2); n = model%n_regimes
        allocate(path(nobs))
        if (model%method == rsdc_const .or. n == 1) then
            path = 1
            if (present(ok)) ok = .true.
            return
        end if
        allocate(logdens(n, nobs), delta(n, nobs), psi(n, nobs), pmat(n, n))
        allocate(inv(k, k), r(k, k), pi1(n), vals(n))
        do s = 1, n
            call correlations_to_matrix(model%correlations(s, :), k, r)
            if (.not. is_positive_definite(r, 1.0e-8_dp)) then
                if (present(ok)) ok = .false.
                return
            end if
            do i = 1, k
                r(i, i) = r(i, i) + 1.0e-8_dp
            end do
            call inverse_spd(r, inv, good)
            if (.not. good) then
                if (present(ok)) ok = .false.
                return
            end if
            ld = logdet_spd(r, good)
            do t = 1, nobs
                quad = dot_product(residuals(t, :), matmul(inv, residuals(t, :)))
                logdens(s, t) = -0.5_dp * (ld + quad + real(k, dp) * log(2.0_dp * rsdc_pi))
            end do
        end do
        call get_p(1, pmat, good)
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if
        pi1 = matmul(transpose(pmat), [(1.0_dp / real(n, dp), i=1,n)])
        delta(:, 1) = log(max(pi1, tiny(1.0_dp))) + logdens(:, 1)
        psi(:, 1) = 0
        do t = 2, nobs
            call get_p(t, pmat, good)
            do s = 1, n
                vals = delta(:, t - 1) + log(max(pmat(:, s), tiny(1.0_dp)))
                best_i = 1; best = vals(1)
                do i = 2, n
                    if (vals(i) > best) then
                        best = vals(i); best_i = i
                    end if
                end do
                psi(s, t) = best_i
                delta(s, t) = best + logdens(s, t)
            end do
        end do
        path(nobs) = maxloc(delta(:, nobs), dim=1)
        do t = nobs - 1, 1, -1
            path(t) = psi(path(t + 1), t + 1)
        end do
        if (present(ok)) ok = .true.

    contains
        subroutine get_p(tt, pp, local_ok)
            integer, intent(in) :: tt
            real(dp), intent(out) :: pp(:, :)
            logical, intent(out) :: local_ok
            local_ok = .true.
            if (model%method == rsdc_nox) then
                pp = model%transition_matrix
            else if (model%method == rsdc_tvtp) then
                if (.not. present(x)) then
                    local_ok = .false.; return
                end if
                call transition_at(model%beta, x(tt, :), pp)
            else
                pp = 1.0_dp
            end if
        end subroutine get_p
    end subroutine rsdc_viterbi_path
end module rsdc_viterbi

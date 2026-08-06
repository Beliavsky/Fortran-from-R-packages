! SPDX-License-Identifier: GPL-3.0-only
module rsdc_inference
    use rsdc_kinds, only: dp
    use rsdc_types, only: rsdc_model, rsdc_diagnostics_result
    use rsdc_types, only: rsdc_const, rsdc_tvtp
    use rsdc_likelihood, only: rsdc_loglik_contributions
    use rsdc_linalg, only: inverse_matrix
    implicit none
    private
    public :: rsdc_scores, rsdc_robust_vcov, rsdc_diagnostics

contains

    subroutine rsdc_scores(model, residuals, scores, x, ok)
        type(rsdc_model), intent(in) :: model
        real(dp), intent(in) :: residuals(:, :)
        real(dp), allocatable, intent(out) :: scores(:, :)
        real(dp), intent(in), optional :: x(:, :)
        logical, intent(out), optional :: ok
        real(dp), allocatable :: pp(:), pm(:), fp(:), fm(:), h(:)
        logical :: good1, good2
        integer :: j, np, t
        np = size(model%parameters); t = size(residuals, 1)
        allocate(scores(t, np), pp(np), pm(np), fp(t), fm(t), h(np))
        h = 1.0e-5_dp * (1.0_dp + abs(model%parameters))
        do j = 1, np
            pp = model%parameters; pm = model%parameters
            pp(j) = pp(j) + h(j); pm(j) = pm(j) - h(j)
            if (model%method == rsdc_tvtp) then
                if (.not. present(x)) then
                    if (present(ok)) ok = .false.
                    return
                end if
                call rsdc_loglik_contributions(pp, residuals, model%method, model%n_regimes, fp, good1, x)
                call rsdc_loglik_contributions(pm, residuals, model%method, model%n_regimes, fm, good2, x)
            else
                call rsdc_loglik_contributions(pp, residuals, model%method, model%n_regimes, fp, good1)
                call rsdc_loglik_contributions(pm, residuals, model%method, model%n_regimes, fm, good2)
            end if
            if (.not. good1 .or. .not. good2) then
                if (present(ok)) ok = .false.
                return
            end if
            scores(:, j) = (fp - fm) / (2.0_dp * h(j))
        end do
        if (present(ok)) ok = .true.
    end subroutine rsdc_scores

    subroutine rsdc_robust_vcov(model, residuals, kind, covariance, x, ok)
        type(rsdc_model), intent(in) :: model
        real(dp), intent(in) :: residuals(:, :)
        character(len=*), intent(in) :: kind
        real(dp), allocatable, intent(out) :: covariance(:, :)
        real(dp), intent(in), optional :: x(:, :)
        logical, intent(out), optional :: ok
        real(dp), allocatable :: scores(:, :), meat(:, :), inv(:, :)
        logical :: good
        call rsdc_scores(model, residuals, scores, x, good)
        if (.not. good) then
            if (present(ok)) ok = .false.
            return
        end if
        allocate(meat(size(scores, 2), size(scores, 2)))
        meat = matmul(transpose(scores), scores)
        if (trim(kind) == 'opg') then
            allocate(inv(size(meat, 1), size(meat, 2)))
            call inverse_matrix(meat, inv, good)
            if (good) then
                allocate(covariance(size(inv, 1), size(inv, 2))); covariance = inv
            end if
        else if (trim(kind) == 'sandwich') then
            good = allocated(model%vcov)
            if (good) then
                allocate(covariance(size(meat, 1), size(meat, 2)))
                covariance = matmul(model%vcov, matmul(meat, model%vcov))
            end if
        else
            good = allocated(model%vcov)
            if (good) then
                allocate(covariance(size(model%vcov, 1), size(model%vcov, 2)))
                covariance = model%vcov
            end if
        end if
        if (present(ok)) ok = good
    end subroutine rsdc_robust_vcov

    subroutine rsdc_diagnostics(model, diagnostics)
        type(rsdc_model), intent(in) :: model
        type(rsdc_diagnostics_result), intent(out) :: diagnostics
        real(dp), allocatable :: pi(:), next(:)
        integer :: i, iter, n
        n = model%n_regimes
        allocate(diagnostics%stay_probability(n), diagnostics%expected_duration(n))
        allocate(diagnostics%ergodic_probability(n), pi(n), next(n))
        do i = 1, n
            diagnostics%stay_probability(i) = model%transition_matrix(i, i)
            if (model%transition_matrix(i, i) < 1.0_dp) then
                diagnostics%expected_duration(i) = 1.0_dp / (1.0_dp - model%transition_matrix(i, i))
            else
                diagnostics%expected_duration(i) = huge(1.0_dp)
            end if
        end do
        pi = 1.0_dp / real(n, dp)
        do iter = 1, 100000
            next = matmul(pi, model%transition_matrix)
            if (maxval(abs(next - pi)) < 1.0e-13_dp) exit
            pi = next
        end do
        diagnostics%ergodic_probability = next / sum(next)
    end subroutine rsdc_diagnostics
end module rsdc_inference

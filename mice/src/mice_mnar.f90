! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Numerical NARFCS kernels derived from mice.impute.mnar.norm.R and mice.impute.mnar.logreg.R.
module mice_mnar
    use r_kinds, only : dp
    use mice_categorical, only : impute_logreg
    use mice_regression, only : mice_norm_draw, norm_draw
    use mice_rng, only : mice_rng_state, rng_normal
    use mice_status, only : mice_ok, mice_invalid_argument, mice_invalid_shape
    implicit none
    private

    public :: impute_mnar_norm
    public :: impute_mnar_logreg

contains

    subroutine impute_mnar_norm(y, observed, x, where, unidentifiable_design, delta, rng, imputed, info, ridge)
        real(dp), intent(in) :: y(:) !! Incomplete numeric response for the identifiable normal regression model.
        logical, intent(in) :: observed(:) !! True for rows used to estimate the identifiable normal regression.
        real(dp), intent(in) :: x(:, :) !! Complete identifiable predictors without an intercept column.
        logical, intent(in) :: where(:) !! True for rows where MNAR imputations are requested.
        real(dp), intent(in) :: unidentifiable_design(:, :) !! Row design for user-specified NARFCS sensitivity terms.
        real(dp), intent(in) :: delta(:) !! Sensitivity coefficients multiplying columns of `unidentifiable_design`.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for posterior and residual Gaussian draws.
        real(dp), allocatable, intent(out) :: imputed(:) !! MNAR normal imputations ordered by true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.
        real(dp), intent(in), optional :: ridge !! Relative identifiable-model ridge; default `1e-5`.

        real(dp), allocatable :: design(:, :), offset(:)
        real(dp) :: ridge_value
        type(mice_norm_draw) :: draw
        integer :: i, j

        if (size(observed) /= size(y) .or. size(where) /= size(y) .or. size(x, 1) /= size(y) .or. &
            size(unidentifiable_design, 1) /= size(y)) then
            info = mice_invalid_shape
            return
        end if
        if (size(unidentifiable_design, 2) /= size(delta)) then
            info = mice_invalid_shape
            return
        end if
        ridge_value = 1.0e-5_dp
        if (present(ridge)) ridge_value = ridge
        if (ridge_value < 0.0_dp) then
            info = mice_invalid_argument
            return
        end if
        allocate(design(size(y), size(x, 2) + 1), offset(size(y)))
        design(:, 1) = 1.0_dp
        if (size(x, 2) > 0) design(:, 2:) = x
        if (size(delta) > 0) then
            offset = matmul(unidentifiable_design, delta)
        else
            offset = 0.0_dp
        end if
        call norm_draw(y, observed, design, rng, ridge_value, draw, info)
        if (info /= mice_ok) return
        allocate(imputed(count(where)))
        j = 0
        do i = 1, size(y)
            if (.not. where(i)) cycle
            j = j + 1
            imputed(j) = dot_product(design(i, :), draw%beta) + offset(i) + draw%sigma * rng_normal(rng)
        end do
        info = mice_ok
    end subroutine impute_mnar_norm

    subroutine impute_mnar_logreg(y, observed, x, where, unidentifiable_design, delta, rng, imputed, info)
        real(dp), intent(in) :: y(:) !! Incomplete binary response coded zero or one.
        logical, intent(in) :: observed(:) !! True for rows used in the identifiable logistic model.
        real(dp), intent(in) :: x(:, :) !! Complete identifiable predictors without an intercept column.
        logical, intent(in) :: where(:) !! True for rows where MNAR binary imputations are requested.
        real(dp), intent(in) :: unidentifiable_design(:, :) !! Row design for user-specified NARFCS sensitivity terms.
        real(dp), intent(in) :: delta(:) !! Sensitivity coefficients multiplying columns of `unidentifiable_design`.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used by posterior logistic and Bernoulli draws.
        real(dp), allocatable, intent(out) :: imputed(:) !! Binary MNAR imputations ordered by true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.

        real(dp), allocatable :: offset(:)

        if (size(observed) /= size(y) .or. size(where) /= size(y) .or. size(x, 1) /= size(y) .or. &
            size(unidentifiable_design, 1) /= size(y)) then
            info = mice_invalid_shape
            return
        end if
        if (size(unidentifiable_design, 2) /= size(delta)) then
            info = mice_invalid_shape
            return
        end if
        allocate(offset(size(y)))
        if (size(delta) > 0) then
            offset = matmul(unidentifiable_design, delta)
        else
            offset = 0.0_dp
        end if
        call impute_logreg(y, observed, x, where, rng, imputed, info, prediction_offset=offset)
    end subroutine impute_mnar_logreg

end module mice_mnar

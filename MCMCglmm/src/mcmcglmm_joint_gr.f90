! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 cs_schur/covu covariance updates; see NOTICE.md and upstream/.
module mcmcglmm_joint_gr
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use mcmcglmm_rng, only : rng_state
    use mcmcglmm_covariance, only : covariance_update_dispatch
    implicit none
    private

    public :: joint_gr_decompose
    public :: joint_gr_compose
    public :: joint_gr_covariance_update

contains

    pure subroutine joint_gr_decompose(joint_covariance, split, leading_covariance, regression, &
                                       conditional_covariance, info)
        real(dp), intent(in) :: joint_covariance(:, :) !! Joint covariance partitioned into leading and trailing blocks.
        integer, intent(in) :: split !! Number of leading coordinates correlated with the trailing random/residual block.
        real(dp), allocatable, intent(out) :: leading_covariance(:, :) !! Allocated leading covariance block S11.
        real(dp), allocatable, intent(out) :: regression(:, :) !! Allocated S21*inverse(S11), matching upstream beta_rr.
        real(dp), allocatable, intent(out) :: conditional_covariance(:, :) !! Allocated Schur complement S22-S21*S11^-1*S12.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid partition or singular leading block.
        real(dp), allocatable :: inverse_leading(:, :)
        real(dp), allocatable :: cross_covariance(:, :)
        integer :: n

        info = 0
        n = size(joint_covariance, 1)
        if (size(joint_covariance, 2) /= n .or. split < 1 .or. split >= n) then
            allocate(leading_covariance(0, 0), regression(0, 0), conditional_covariance(0, 0))
            info = 1
            return
        end if
        leading_covariance = joint_covariance(1:split, 1:split)
        call inverse_matrix(leading_covariance, inverse_leading, info)
        if (info /= 0) then
            allocate(regression(0, 0), conditional_covariance(0, 0))
            return
        end if
        cross_covariance = joint_covariance(split + 1:n, 1:split)
        regression = matmul(cross_covariance, inverse_leading)
        conditional_covariance = joint_covariance(split + 1:n, split + 1:n) - &
            matmul(regression, transpose(cross_covariance))
    end subroutine joint_gr_decompose

    pure subroutine joint_gr_compose(leading_covariance, regression, conditional_covariance, joint_covariance, info)
        real(dp), intent(in) :: leading_covariance(:, :) !! Positive-definite leading covariance block S11.
        real(dp), intent(in) :: regression(:, :) !! Trailing-on-leading regression S21*inverse(S11).
        real(dp), intent(in) :: conditional_covariance(:, :) !! Trailing conditional covariance Schur complement.
        real(dp), allocatable, intent(out) :: joint_covariance(:, :) !! Allocated covariance reconstructed from the blocks.
        integer, intent(out) :: info !! Zero on success; nonzero for incompatible shapes.
        real(dp), allocatable :: lower_cross(:, :)
        integer :: leading
        integer :: trailing

        info = 0
        leading = size(leading_covariance, 1)
        trailing = size(conditional_covariance, 1)
        if (leading < 1 .or. trailing < 1 .or. size(leading_covariance, 2) /= leading .or. &
            size(conditional_covariance, 2) /= trailing .or. size(regression, 1) /= trailing .or. &
            size(regression, 2) /= leading) then
            allocate(joint_covariance(0, 0))
            info = 1
            return
        end if
        lower_cross = matmul(regression, leading_covariance)
        allocate(joint_covariance(leading + trailing, leading + trailing))
        joint_covariance = 0.0_dp
        joint_covariance(1:leading, 1:leading) = leading_covariance
        joint_covariance(leading + 1:, 1:leading) = lower_cross
        joint_covariance(1:leading, leading + 1:) = transpose(lower_cross)
        joint_covariance(leading + 1:, leading + 1:) = conditional_covariance + &
            matmul(lower_cross, transpose(regression))
    end subroutine joint_gr_compose

    pure subroutine joint_gr_covariance_update(state, update_code, posterior_sum, sample_degrees_freedom, &
                                               prior_degrees_freedom, prior_scale, old_joint_covariance, split, &
                                               fixed_block, joint_covariance, conditional_covariance, regression, &
                                               accepted, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the selected covariance update.
        integer, intent(in) :: update_code !! Joint covariance update code 0, 1, 2, 3, 4, or 6.
        real(dp), intent(in) :: posterior_sum(:, :) !! Joint data cross-product plus prior scale.
        real(dp), intent(in) :: sample_degrees_freedom !! Effective number of joint covariance observations/levels.
        real(dp), intent(in) :: prior_degrees_freedom !! Joint covariance prior degrees of freedom.
        real(dp), intent(in) :: prior_scale(:, :) !! Joint covariance prior scale matrix.
        real(dp), intent(in) :: old_joint_covariance(:, :) !! Current joint covariance before the update.
        integer, intent(in) :: split !! Number of leading covu coordinates used by the Schur decomposition.
        real(dp), intent(in) :: fixed_block(:, :) !! Fixed lower-right block for update codes 2/4; zero-size otherwise.
        real(dp), allocatable, intent(out) :: joint_covariance(:, :) !! Allocated updated full joint covariance.
        real(dp), allocatable, intent(out) :: conditional_covariance(:, :) !! Allocated trailing Schur-complement covariance.
        real(dp), allocatable, intent(out) :: regression(:, :) !! Allocated trailing-on-leading regression beta_rr.
        logical, intent(out) :: accepted !! Correlation-MH acceptance flag, true for direct Gibbs/fixed updates.
        integer, intent(out) :: info !! Zero on success; nonzero for covariance update or decomposition failure.
        real(dp), allocatable :: fixed_block_out(:, :)
        real(dp), allocatable :: leading_covariance(:, :)

        call covariance_update_dispatch(state, update_code, posterior_sum, sample_degrees_freedom, &
            prior_degrees_freedom, prior_scale, old_joint_covariance, split, fixed_block, joint_covariance, &
            fixed_block_out, accepted, info)
        if (info /= 0) then
            allocate(conditional_covariance(0, 0), regression(0, 0))
            return
        end if
        call joint_gr_decompose(joint_covariance, split, leading_covariance, regression, conditional_covariance, info)
    end subroutine joint_gr_covariance_update

end module mcmcglmm_joint_gr

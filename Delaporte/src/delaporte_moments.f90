! Copyright (c) 2016, Avraham Adler
! All rights reserved.
! SPDX-License-Identifier: BSD-2-Clause
!
! Method-of-moments algorithm adapted from src/delaporte.f90 upstream.

module delaporte_moments
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, &
        ieee_quiet_nan, ieee_value
    use delaporte_kinds, only : dp
    use delaporte_utils, only : zero, half, one, two, three, three_halves
    implicit none
    private

    public :: momdelap

contains

    pure subroutine momdelap(obs, params, skew_type, status)
        real(dp), intent(in) :: obs(:)
        real(dp), intent(out) :: params(3)
        integer, intent(in), optional :: skew_type
        integer, intent(out), optional :: status
        real(dp) :: nnm1, pcor, mu_d, m2, m3, t1
        real(dp) :: delta, delta_i, nn, var_d, skew_d, vmm_d, ii
        integer :: i, tp, stat

        tp = 2
        if (present(skew_type)) tp = skew_type
        stat = 0

        if (size(obs) < 3 .or. tp < 1 .or. tp > 3) then
            params = ieee_value(0.0_dp, ieee_quiet_nan)
            stat = 1
            if (present(status)) status = stat
            return
        end if

        nn = real(size(obs), dp)
        nnm1 = nn - one
        select case (tp)
        case (1)
            pcor = one
        case (2)
            pcor = sqrt(nn * nnm1) / (nn - two)
        case (3)
            pcor = (nnm1 / nn) ** three_halves
        end select

        mu_d = zero
        m2 = zero
        m3 = zero
        do i = 1, size(obs)
            ii = real(i, dp)
            delta = obs(i) - mu_d
            delta_i = delta / ii
            t1 = delta * delta_i * (ii - one)
            mu_d = mu_d + delta_i
            m3 = m3 + t1 * delta_i * (ii - two) - three * delta_i * m2
            m2 = m2 + t1
        end do

        var_d = m2 / nnm1
        skew_d = pcor * sqrt(nn) * m3 / (m2 ** three_halves)
        vmm_d = var_d - mu_d
        params(2) = half * (skew_d * var_d ** three_halves - mu_d - &
            three * vmm_d) / vmm_d
        params(1) = vmm_d / (params(2) ** 2)
        params(3) = mu_d - params(1) * params(2)

        if (any(params <= zero) .or. .not. all(ieee_is_finite(params))) stat = 2
        if (present(status)) status = stat
    end subroutine momdelap

end module delaporte_moments

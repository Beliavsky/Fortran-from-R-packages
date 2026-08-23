! SPDX-License-Identifier: GPL-3.0-only
module ecpdist
    use ecpdist_kinds, only: dp
    use ecpdist_distribution, only: decp, pecp, qecp, recp, secp, hecp, &
        ecp_cumhaz, ecp_valid_parameters
    use ecpdist_moments, only: ecp_integral_result, ecp_kmoment, &
        ecp_kmoment_cond, ecp_mrl, ecp_shape
    implicit none
    public
end module ecpdist

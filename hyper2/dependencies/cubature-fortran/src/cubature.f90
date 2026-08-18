! SPDX-License-Identifier: GPL-3.0-or-later
module cubature
    use cubature_kinds, only : dp, i8
    use cubature_types, only : cubature_result, cubature_integrand, cubature_integrand_v, &
        ERROR_INDIVIDUAL, ERROR_PAIRED, ERROR_L2, ERROR_L1, ERROR_LINF, &
        CUBATURE_SUCCESS, CUBATURE_MAXEVAL, CUBATURE_BADARG, CUBATURE_FAILURE, &
        cuhre_options, divonne_options, suave_options, vegas_options
    use hcubature_mod, only : hcubature, hcubature_v, adapt_integrate, adaptIntegrate
    use pcubature_mod, only : pcubature, pcubature_v
    use cuba_mod, only : cuhre, cuhre_v, divonne, divonne_v, suave, suave_v, vegas, vegas_v
    use cubintegrate_mod, only : cubintegrate
    implicit none
    public
end module cubature

! SPDX-License-Identifier: GPL-2.0-or-later
module ltsa_types
    use ltsa_kinds, only : dp
    use ltsa_status, only : ltsa_error
    implicit none
    private

    type, public :: dl_ar_result
        real(dp), allocatable :: phi(:)
        real(dp), allocatable :: pacf(:)
        real(dp), allocatable :: prediction_variance(:)
        type(ltsa_error) :: error
    end type dl_ar_result

    type, public :: exact_likelihood_result
        real(dp) :: log_likelihood = 0.0_dp
        real(dp) :: sigma_sq = 0.0_dp
        type(ltsa_error) :: error
    end type exact_likelihood_result

    type, public :: forecast_result
        real(dp), allocatable :: forecasts(:,:)
        real(dp), allocatable :: sd_forecasts(:,:)
        type(ltsa_error) :: error
    end type forecast_result

    type, public :: innovation_variance_result
        real(dp) :: variance = 0.0_dp
        integer :: selected_order = 0
        type(ltsa_error) :: error
    end type innovation_variance_result

end module ltsa_types

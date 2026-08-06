! SPDX-License-Identifier: GPL-2.0-or-later
module ltsa_compat
    use ltsa_kinds, only : dp
    use ltsa_status, only : ltsa_error
    use ltsa_types, only : dl_ar_result, exact_likelihood_result, forecast_result, innovation_variance_result
    use ltsa_linalg, only : is_toeplitz
    use ltsa_durbin_levinson, only : dl_acf_to_ar, dl_residuals, dl_loglikelihood, dl_simulate
    use ltsa_toeplitz, only : trench_inverse, toeplitz_inverse_update, trench_mean
    use ltsa_arma, only : tacvf_arma
    use ltsa_simulation, only : sim_glp, dh_condition, dh_simulate
    use ltsa_likelihood, only : trench_loglikelihood, exact_loglikelihood
    use ltsa_forecast, only : prediction_variance, trench_forecast
    use ltsa_innovation, only : innovation_variance
    implicit none
    private

    public :: dhsimulate, dlacftoar, dlloglikelihood, dlresiduals, dlsimulate
    public :: exactloglikelihood, predictionvariance, innovationvariance, simglp
    public :: toeplitzinverseupdate, trenchforecast, trenchinverse
    public :: trenchloglikelihood, trenchmean, istoeplitz, tacvfarma

contains

    subroutine dhsimulate(n, r, z, error, report_test_only, condition, source_compatible)
        integer, intent(in) :: n
        real(dp), intent(in) :: r(:)
        real(dp), allocatable, intent(out) :: z(:)
        type(ltsa_error), intent(out) :: error
        logical, intent(in), optional :: report_test_only, source_compatible
        logical, intent(out), optional :: condition
        logical :: test_only, valid
        test_only = .false.
        if (present(report_test_only)) test_only = report_test_only
        valid = dh_condition(n,r)
        if (present(condition)) condition = valid
        if (test_only) then
            allocate(z(0))
            error%code = 0
            error%message = ''
        else
            call dh_simulate(n,r,z,error,source_compatible)
        end if
    end subroutine dhsimulate

    function dlacftoar(r) result(value)
        real(dp), intent(in) :: r(:)
        type(dl_ar_result) :: value
        value = dl_acf_to_ar(r)
    end function dlacftoar

    function dlloglikelihood(r,z,error) result(value)
        real(dp), intent(in) :: r(:),z(:)
        type(ltsa_error), intent(out), optional :: error
        real(dp) :: value
        value = dl_loglikelihood(r,z,error)
    end function dlloglikelihood

    subroutine dlresiduals(r,z,residuals,error,standardized)
        real(dp), intent(in) :: r(:),z(:)
        real(dp), allocatable, intent(out) :: residuals(:)
        type(ltsa_error), intent(out) :: error
        logical, intent(in), optional :: standardized
        call dl_residuals(r,z,residuals,error,standardized)
    end subroutine dlresiduals

    subroutine dlsimulate(n,r,z,error,innovations)
        integer, intent(in) :: n
        real(dp), intent(in) :: r(:)
        real(dp), allocatable, intent(out) :: z(:)
        type(ltsa_error), intent(out) :: error
        real(dp), intent(in), optional :: innovations(:)
        call dl_simulate(n,r,z,error,innovations)
    end subroutine dlsimulate

    function exactloglikelihood(r,z,innovation_variance) result(value)
        real(dp), intent(in) :: r(:),z(:)
        logical, intent(in), optional :: innovation_variance
        type(exact_likelihood_result) :: value
        value = exact_loglikelihood(r,z,innovation_variance)
    end function exactloglikelihood

    subroutine predictionvariance(r,max_lead,variances,error,dlq)
        real(dp), intent(in) :: r(:)
        integer, intent(in) :: max_lead
        real(dp), allocatable, intent(out) :: variances(:)
        type(ltsa_error), intent(out) :: error
        logical, intent(in), optional :: dlq
        call prediction_variance(r,max_lead,variances,error,dlq)
    end subroutine predictionvariance

    function innovationvariance(z,method,max_order,smooth_span) result(value)
        real(dp), intent(in) :: z(:)
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_order,smooth_span
        type(innovation_variance_result) :: value
        value = innovation_variance(z,method,max_order,smooth_span)
    end function innovationvariance

    subroutine simglp(psi,innovations,z,error)
        real(dp), intent(in) :: psi(:),innovations(:)
        real(dp), allocatable, intent(out) :: z(:)
        type(ltsa_error), intent(out) :: error
        call sim_glp(psi,innovations,z,error)
    end subroutine simglp

    subroutine toeplitzinverseupdate(gi,r,rnew,updated,error)
        real(dp), intent(in) :: gi(:,:),r(:),rnew
        real(dp), allocatable, intent(out) :: updated(:,:)
        type(ltsa_error), intent(out) :: error
        call toeplitz_inverse_update(gi,r,rnew,updated,error)
    end subroutine toeplitzinverseupdate

    function trenchforecast(z,r,mean_value,origin,max_lead,update_algorithm) result(value)
        real(dp), intent(in) :: z(:),r(:),mean_value
        integer, intent(in) :: origin,max_lead
        logical, intent(in), optional :: update_algorithm
        type(forecast_result) :: value
        value = trench_forecast(z,r,mean_value,origin,max_lead,update_algorithm)
    end function trenchforecast

    subroutine trenchinverse(g,gi,error)
        real(dp), intent(in) :: g(:,:)
        real(dp), allocatable, intent(out) :: gi(:,:)
        type(ltsa_error), intent(out) :: error
        call trench_inverse(g,gi,error)
    end subroutine trenchinverse

    function trenchloglikelihood(r,z,error) result(value)
        real(dp), intent(in) :: r(:),z(:)
        type(ltsa_error), intent(out), optional :: error
        real(dp) :: value
        value = trench_loglikelihood(r,z,error)
    end function trenchloglikelihood

    function trenchmean(r,z,error) result(value)
        real(dp), intent(in) :: r(:),z(:)
        type(ltsa_error), intent(out), optional :: error
        real(dp) :: value
        value = trench_mean(r,z,error)
    end function trenchmean

    logical function istoeplitz(g) result(value)
        real(dp), intent(in) :: g(:,:)
        value = is_toeplitz(g)
    end function istoeplitz

    subroutine tacvfarma(phi,theta,max_lag,sigma2,acvf,error)
        real(dp), intent(in) :: phi(:),theta(:),sigma2
        integer, intent(in) :: max_lag
        real(dp), allocatable, intent(out) :: acvf(:)
        type(ltsa_error), intent(out) :: error
        call tacvf_arma(phi,theta,max_lag,sigma2,acvf,error)
    end subroutine tacvfarma

end module ltsa_compat

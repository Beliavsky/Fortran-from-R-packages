! SPDX-License-Identifier: GPL-2.0-or-later
module ltsa_innovation
    use ltsa_kinds, only : dp, two_pi, pi
    use ltsa_status, only : ltsa_error, ltsa_success, ltsa_invalid_input, set_error
    use ltsa_types, only : innovation_variance_result
    use ltsa_durbin_levinson, only : durbin_levinson_table
    implicit none
    private

    public :: innovation_variance

contains

    function innovation_variance(z, method, max_order, smooth_span) result(result_value)
        real(dp), intent(in) :: z(:)
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: max_order, smooth_span
        type(innovation_variance_result) :: result_value
        character(len=:), allocatable :: selected_method
        integer :: pmax, span
        selected_method = 'ar'
        if (present(method)) selected_method = lower_string(trim(method))
        if (size(z) < 10) then
            call set_error(result_value%error, ltsa_invalid_input, 'series length must be at least 10')
            return
        end if
        pmax = min(size(z)-2, max(1,int(10.0_dp*log10(real(size(z),dp)))))
        if (present(max_order)) pmax = max(0,min(max_order,size(z)-2))
        span = 1
        if (present(smooth_span)) span = max(1,smooth_span)
        select case (selected_method)
        case ('ar')
            call innovation_variance_ar(z,pmax,result_value)
        case ('kolmogoroff','kolmogorov')
            call innovation_variance_kolmogoroff(z,span,result_value)
        case default
            call set_error(result_value%error, ltsa_invalid_input, 'method must be AR or Kolmogoroff')
        end select
    end function innovation_variance

    subroutine innovation_variance_ar(z,pmax,result_value)
        real(dp), intent(in) :: z(:)
        integer, intent(in) :: pmax
        type(innovation_variance_result), intent(inout) :: result_value
        real(dp), allocatable :: centered(:), acvf(:), table(:,:), pacf(:), variances(:)
        real(dp) :: mean_value, criterion, best
        integer :: h, p, n
        n = size(z)
        mean_value = sum(z)/real(n,dp)
        allocate(centered(n))
        centered = z-mean_value
        allocate(acvf(pmax+1))
        do h = 0, pmax
            acvf(h+1) = dot_product(centered(1:n-h),centered(1+h:n))/real(n,dp)
        end do
        if (acvf(1) <= 0.0_dp) then
            call set_error(result_value%error, ltsa_invalid_input, 'series variance is zero')
            return
        end if
        call durbin_levinson_table(acvf,table,pacf,variances,result_value%error)
        if (.not. result_value%error%ok()) return
        best = log(variances(1))
        result_value%variance = variances(1)
        result_value%selected_order = 0
        do p = 1, pmax
            criterion = log(variances(p+1))+2.0_dp*real(p,dp)/real(n,dp)
            if (criterion < best) then
                best = criterion
                result_value%variance = variances(p+1)
                result_value%selected_order = p
            end if
        end do
    end subroutine innovation_variance_ar

    subroutine innovation_variance_kolmogoroff(z,span,result_value)
        real(dp), intent(in) :: z(:)
        integer, intent(in) :: span
        type(innovation_variance_result), intent(inout) :: result_value
        real(dp), allocatable :: x(:), spectrum(:), smoothed(:)
        real(dp) :: mean_value, re, im, angle
        integer :: j, k, lo, hi, m, n
        n = size(z)
        mean_value = sum(z)/real(n,dp)
        allocate(x(n))
        x = z-mean_value
        m = n/2
        allocate(spectrum(m),smoothed(m))
        do k = 1, m
            re = 0.0_dp
            im = 0.0_dp
            do j = 1, n
                angle = two_pi*real(k*(j-1),dp)/real(n,dp)
                re = re+x(j)*cos(angle)
                im = im-x(j)*sin(angle)
            end do
            spectrum(k) = max(tiny(1.0_dp),(re*re+im*im)/(two_pi*real(n,dp)))
        end do
        if (span > 1) then
            do k = 1, m
                lo = max(1,k-span/2)
                hi = min(m,k+span/2)
                smoothed(k) = sum(spectrum(lo:hi))/real(hi-lo+1,dp)
            end do
            spectrum = smoothed
        end if
        result_value%variance = exp(2.0_dp*sum(log(two_pi*spectrum))/real(n,dp))/two_pi
        result_value%selected_order = 0
        result_value%error%code = ltsa_success
        result_value%error%message = ''
    end subroutine innovation_variance_kolmogoroff

    pure function lower_string(s) result(out)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: out
        integer :: i, c
        do i = 1, len(s)
            c = iachar(s(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) then
                out(i:i) = achar(c+32)
            else
                out(i:i) = s(i:i)
            end if
        end do
    end function lower_string

end module ltsa_innovation

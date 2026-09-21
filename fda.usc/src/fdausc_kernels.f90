module fdausc_kernels
    use r_kinds, only : dp
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    public :: ker_norm, ker_cos, ker_epa, ker_tri, ker_quar, ker_unif
    public :: aker_norm, aker_cos, aker_epa, aker_tri, aker_quar, aker_unif
    public :: iker_norm, iker_cos, iker_epa, iker_tri, iker_quar, iker_unif
    public :: kernel_value, asymmetric_kernel_value, integrated_kernel_value

contains

    pure elemental real(dp) function ker_norm(u) result(k)
        real(dp), intent(in) :: u !! Standardized kernel argument on the real line.
        k = exp(-0.5_dp*u*u)/sqrt(2.0_dp*pi)
    end function ker_norm

    pure elemental real(dp) function ker_cos(u) result(k)
        real(dp), intent(in) :: u !! Standardized kernel argument; support is [-1, 1].
        if (abs(u) <= 1.0_dp) then
            k = 0.25_dp*pi*cos(0.5_dp*pi*u)
        else
            k = 0.0_dp
        end if
    end function ker_cos

    pure elemental real(dp) function ker_epa(u) result(k)
        real(dp), intent(in) :: u !! Standardized kernel argument; support is [-1, 1].
        if (abs(u) <= 1.0_dp) then
            k = 0.75_dp*(1.0_dp - u*u)
        else
            k = 0.0_dp
        end if
    end function ker_epa

    pure elemental real(dp) function ker_tri(u) result(k)
        real(dp), intent(in) :: u !! Standardized kernel argument; support is [-1, 1].
        if (abs(u) <= 1.0_dp) then
            k = (35.0_dp/32.0_dp)*(1.0_dp - u*u)**3
        else
            k = 0.0_dp
        end if
    end function ker_tri

    pure elemental real(dp) function ker_quar(u) result(k)
        real(dp), intent(in) :: u !! Standardized kernel argument; support is [-1, 1].
        if (abs(u) <= 1.0_dp) then
            k = (15.0_dp/16.0_dp)*(1.0_dp - u*u)**2
        else
            k = 0.0_dp
        end if
    end function ker_quar

    pure elemental real(dp) function ker_unif(u) result(k)
        real(dp), intent(in) :: u !! Standardized kernel argument; support is [-1, 1].
        if (abs(u) <= 1.0_dp) then
            k = 0.5_dp
        else
            k = 0.0_dp
        end if
    end function ker_unif

    pure elemental real(dp) function aker_norm(u) result(k)
        real(dp), intent(in) :: u !! One-sided standardized kernel argument; support begins at zero.
        if (u >= 0.0_dp) then
            k = 2.0_dp*ker_norm(u)
        else
            k = 0.0_dp
        end if
    end function aker_norm

    pure elemental real(dp) function aker_cos(u) result(k)
        real(dp), intent(in) :: u !! One-sided cosine-kernel argument, matching upstream behavior for u >= 0.
        if (u >= 0.0_dp) then
            k = 0.5_dp*pi*cos(0.5_dp*pi*u)
        else
            k = 0.0_dp
        end if
    end function aker_cos

    pure elemental real(dp) function aker_epa(u) result(k)
        real(dp), intent(in) :: u !! One-sided Epanechnikov argument; support is [0, 1].
        if (u >= 0.0_dp .and. u <= 1.0_dp) then
            k = 1.5_dp*(1.0_dp - u*u)
        else
            k = 0.0_dp
        end if
    end function aker_epa

    pure elemental real(dp) function aker_tri(u) result(k)
        real(dp), intent(in) :: u !! One-sided triweight argument; support is [0, 1].
        if (u >= 0.0_dp .and. u <= 1.0_dp) then
            k = (35.0_dp/16.0_dp)*(1.0_dp - u*u)**3
        else
            k = 0.0_dp
        end if
    end function aker_tri

    pure elemental real(dp) function aker_quar(u) result(k)
        real(dp), intent(in) :: u !! One-sided quartic argument; support is [0, 1].
        if (u >= 0.0_dp .and. u <= 1.0_dp) then
            k = (15.0_dp/8.0_dp)*(1.0_dp - u*u)**2
        else
            k = 0.0_dp
        end if
    end function aker_quar

    pure elemental real(dp) function aker_unif(u) result(k)
        real(dp), intent(in) :: u !! One-sided uniform-kernel argument; support is [0, 1].
        if (u >= 0.0_dp .and. u <= 1.0_dp) then
            k = 1.0_dp
        else
            k = 0.0_dp
        end if
    end function aker_unif

    pure elemental real(dp) function iker_norm(u) result(k)
        real(dp), intent(in) :: u !! Upper integration limit for the normal kernel, with lower limit -1.
        k = 0.5_dp*(1.0_dp + erf(u/sqrt(2.0_dp))) &
            - 0.5_dp*(1.0_dp + erf(-1.0_dp/sqrt(2.0_dp)))
    end function iker_norm

    pure elemental real(dp) function iker_cos(u) result(k)
        real(dp), intent(in) :: u !! Upper integration limit for the cosine kernel, with lower limit -1.
        if (u <= -1.0_dp) then
            k = 0.0_dp
        else if (u >= 1.0_dp) then
            k = 1.0_dp
        else
            k = 0.5_dp*(sin(0.5_dp*pi*u) + 1.0_dp)
        end if
    end function iker_cos

    pure elemental real(dp) function iker_epa(u) result(k)
        real(dp), intent(in) :: u !! Upper integration limit for the Epanechnikov kernel, lower limit -1.
        real(dp) :: z
        z = max(-1.0_dp, min(1.0_dp, u))
        if (u <= -1.0_dp) then
            k = 0.0_dp
        else if (u >= 1.0_dp) then
            k = 1.0_dp
        else
            k = 0.75_dp*(z - z**3/3.0_dp + 2.0_dp/3.0_dp)
        end if
    end function iker_epa

    pure elemental real(dp) function iker_tri(u) result(k)
        real(dp), intent(in) :: u !! Upper integration limit for the triweight kernel, lower limit -1.
        real(dp) :: z
        z = max(-1.0_dp, min(1.0_dp, u))
        if (u <= -1.0_dp) then
            k = 0.0_dp
        else if (u >= 1.0_dp) then
            k = 1.0_dp
        else
            k = (35.0_dp/32.0_dp)*(z - z**3 + 0.6_dp*z**5 - z**7/7.0_dp + 16.0_dp/35.0_dp)
        end if
    end function iker_tri

    pure elemental real(dp) function iker_quar(u) result(k)
        real(dp), intent(in) :: u !! Upper integration limit for the quartic kernel, lower limit -1.
        real(dp) :: z
        z = max(-1.0_dp, min(1.0_dp, u))
        if (u <= -1.0_dp) then
            k = 0.0_dp
        else if (u >= 1.0_dp) then
            k = 1.0_dp
        else
            k = (15.0_dp/16.0_dp)*(z - 2.0_dp*z**3/3.0_dp + z**5/5.0_dp + 8.0_dp/15.0_dp)
        end if
    end function iker_quar

    pure elemental real(dp) function iker_unif(u) result(k)
        real(dp), intent(in) :: u !! Upper integration limit for the uniform kernel, lower limit -1.
        if (u <= -1.0_dp) then
            k = 0.0_dp
        else if (u >= 1.0_dp) then
            k = 1.0_dp
        else
            k = 0.5_dp*(u + 1.0_dp)
        end if
    end function iker_unif

    pure elemental real(dp) function kernel_value(u, kernel_code) result(k)
        real(dp), intent(in) :: u !! Standardized argument at which to evaluate the selected symmetric kernel.
        integer, intent(in) :: kernel_code !! Kernel code: 1 normal, 2 cosine, 3 Epanechnikov, 4 triweight, 5 quartic, 6 uniform.
        select case (kernel_code)
        case (2)
            k = ker_cos(u)
        case (3)
            k = ker_epa(u)
        case (4)
            k = ker_tri(u)
        case (5)
            k = ker_quar(u)
        case (6)
            k = ker_unif(u)
        case default
            k = ker_norm(u)
        end select
    end function kernel_value

    pure elemental real(dp) function asymmetric_kernel_value(u, kernel_code) result(k)
        real(dp), intent(in) :: u !! One-sided standardized argument at which the selected kernel is evaluated.
        integer, intent(in) :: kernel_code !! Kernel code: 1 normal, 2 cosine, 3 Epanechnikov, 4 triweight, 5 quartic, 6 uniform.
        select case (kernel_code)
        case (2)
            k = aker_cos(u)
        case (3)
            k = aker_epa(u)
        case (4)
            k = aker_tri(u)
        case (5)
            k = aker_quar(u)
        case (6)
            k = aker_unif(u)
        case default
            k = aker_norm(u)
        end select
    end function asymmetric_kernel_value

    pure elemental real(dp) function integrated_kernel_value(u, kernel_code) result(k)
        real(dp), intent(in) :: u !! Upper integration limit for the selected symmetric kernel, starting at -1.
        integer, intent(in) :: kernel_code !! Kernel code: 1 normal, 2 cosine, 3 Epanechnikov, 4 triweight, 5 quartic, 6 uniform.
        select case (kernel_code)
        case (2)
            k = iker_cos(u)
        case (3)
            k = iker_epa(u)
        case (4)
            k = iker_tri(u)
        case (5)
            k = iker_quar(u)
        case (6)
            k = iker_unif(u)
        case default
            k = iker_norm(u)
        end select
    end function integrated_kernel_value

end module fdausc_kernels

module skewt
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
        ieee_positive_inf, ieee_negative_inf
    use skewt_special, only : dp, student_t_pdf, student_t_cdf, &
        student_t_quantile, student_t_random
    implicit none
    private

    public :: dp, dskt, pskt, qskt, rskt

    interface dskt
        module procedure dskt_scalar
        module procedure dskt_vector
    end interface dskt

    interface pskt
        module procedure pskt_scalar
        module procedure pskt_vector
    end interface pskt

    interface qskt
        module procedure qskt_scalar
        module procedure qskt_vector
    end interface qskt

contains

    pure real(dp) function invalid_value() result(x)
        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function invalid_value

    real(dp) function dskt_scalar(x, df, gamma) result(f)
        real(dp), intent(in) :: x, df
        real(dp), intent(in), optional :: gamma
        real(dp) :: g, c

        g = 1.0_dp
        if (present(gamma)) g = gamma
        if (df <= 0.0_dp .or. g <= 0.0_dp) then
            f = invalid_value()
            return
        end if

        c = 2.0_dp / (g + 1.0_dp / g)
        if (x < 0.0_dp) then
            f = c * student_t_pdf(g * x, df)
        else
            f = c * student_t_pdf(x / g, df)
        end if
    end function dskt_scalar

    function dskt_vector(x, df, gamma) result(f)
        real(dp), intent(in) :: x(:), df
        real(dp), intent(in), optional :: gamma
        real(dp) :: f(size(x))
        integer :: i

        if (present(gamma)) then
            do i = 1, size(x)
                f(i) = dskt_scalar(x(i), df, gamma)
            end do
        else
            do i = 1, size(x)
                f(i) = dskt_scalar(x(i), df)
            end do
        end if
    end function dskt_vector

    real(dp) function pskt_scalar(x, df, gamma) result(p)
        real(dp), intent(in) :: x, df
        real(dp), intent(in), optional :: gamma
        real(dp) :: g, pzero

        g = 1.0_dp
        if (present(gamma)) g = gamma
        if (df <= 0.0_dp .or. g <= 0.0_dp) then
            p = invalid_value()
            return
        end if

        pzero = 1.0_dp / (g * g + 1.0_dp)
        if (x < 0.0_dp) then
            p = 2.0_dp / (g * g + 1.0_dp) * student_t_cdf(g * x, df)
        else
            p = pzero + 2.0_dp * g * g / (g * g + 1.0_dp) * &
                (student_t_cdf(x / g, df) - 0.5_dp)
        end if
        p = min(1.0_dp, max(0.0_dp, p))
    end function pskt_scalar

    function pskt_vector(x, df, gamma) result(p)
        real(dp), intent(in) :: x(:), df
        real(dp), intent(in), optional :: gamma
        real(dp) :: p(size(x))
        integer :: i

        if (present(gamma)) then
            do i = 1, size(x)
                p(i) = pskt_scalar(x(i), df, gamma)
            end do
        else
            do i = 1, size(x)
                p(i) = pskt_scalar(x(i), df)
            end do
        end if
    end function pskt_vector

    real(dp) function qskt_scalar(p, df, gamma) result(x)
        real(dp), intent(in) :: p, df
        real(dp), intent(in), optional :: gamma
        real(dp) :: g, pzero, ptarg

        g = 1.0_dp
        if (present(gamma)) g = gamma
        if (df <= 0.0_dp .or. g <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
            x = invalid_value()
            return
        end if
        if (p <= 0.0_dp) then
            x = ieee_value(x, ieee_negative_inf)
            return
        else if (p >= 1.0_dp) then
            x = ieee_value(x, ieee_positive_inf)
            return
        end if

        pzero = 1.0_dp / (g * g + 1.0_dp)
        if (p < pzero) then
            ptarg = 0.5_dp * (g * g + 1.0_dp) * p
            x = student_t_quantile(ptarg, df) / g
        else
            ptarg = 0.5_dp + 0.5_dp * (1.0_dp + 1.0_dp / (g * g)) * &
                (p - pzero)
            x = g * student_t_quantile(ptarg, df)
        end if
    end function qskt_scalar

    function qskt_vector(p, df, gamma) result(x)
        real(dp), intent(in) :: p(:), df
        real(dp), intent(in), optional :: gamma
        real(dp) :: x(size(p))
        integer :: i

        if (present(gamma)) then
            do i = 1, size(p)
                x(i) = qskt_scalar(p(i), df, gamma)
            end do
        else
            do i = 1, size(p)
                x(i) = qskt_scalar(p(i), df)
            end do
        end if
    end function qskt_vector

    subroutine rskt(x, df, gamma)
        real(dp), intent(out) :: x(:)
        real(dp), intent(in) :: df
        real(dp), intent(in), optional :: gamma
        real(dp) :: g, u
        integer :: i

        g = 1.0_dp
        if (present(gamma)) g = gamma
        if (df <= 0.0_dp .or. g <= 0.0_dp) then
            x = invalid_value()
            return
        end if

        ! Direct inverse transform exactly mirrors upstream rskt().
        do i = 1, size(x)
            call random_number(u)
            x(i) = qskt_scalar(u, df, g)
        end do
    end subroutine rskt

end module skewt

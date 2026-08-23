module rnanoflann_metrics
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use rnanoflann_kinds, only: dp
    implicit none
    private

    integer, parameter, public :: metric_euclidean = 1
    integer, parameter, public :: metric_hellinger = 2
    integer, parameter, public :: metric_manhattan = 3
    integer, parameter, public :: metric_canberra = 4
    integer, parameter, public :: metric_kullback_leibler = 5
    integer, parameter, public :: metric_jensen_shannon = 6
    integer, parameter, public :: metric_itakura_saito = 7
    integer, parameter, public :: metric_bhattacharyya = 8
    integer, parameter, public :: metric_jeffries_matusita = 9
    integer, parameter, public :: metric_minimum = 10
    integer, parameter, public :: metric_maximum = 11
    integer, parameter, public :: metric_total_variation = 12
    integer, parameter, public :: metric_sorensen = 13
    integer, parameter, public :: metric_cosine = 14
    integer, parameter, public :: metric_gower = 15
    integer, parameter, public :: metric_minkowski = 16
    integer, parameter, public :: metric_soergel = 17
    integer, parameter, public :: metric_kulczynski = 18
    integer, parameter, public :: metric_wave_hedges = 19
    integer, parameter, public :: metric_motyka = 20
    integer, parameter, public :: metric_harmonic_mean = 21

    public :: metric_code, metric_distance, metric_distance_prepared

contains

    pure integer function metric_code(name) result(code)
        character(*), intent(in) :: name
        character(len=:), allocatable :: key

        key = lower_ascii(trim(adjustl(name)))
        select case (key)
        case ("euclidean")
            code = metric_euclidean
        case ("hellinger")
            code = metric_hellinger
        case ("manhattan")
            code = metric_manhattan
        case ("canberra")
            code = metric_canberra
        case ("kullback_leibler")
            code = metric_kullback_leibler
        case ("jensen_shannon")
            code = metric_jensen_shannon
        case ("itakura_saito")
            code = metric_itakura_saito
        case ("bhattacharyya")
            code = metric_bhattacharyya
        case ("jeffries_matusita")
            code = metric_jeffries_matusita
        case ("minimum")
            code = metric_minimum
        case ("maximum")
            code = metric_maximum
        case ("total_variation")
            code = metric_total_variation
        case ("sorensen")
            code = metric_sorensen
        case ("cosine")
            code = metric_cosine
        case ("gower")
            code = metric_gower
        case ("minkowski")
            code = metric_minkowski
        case ("soergel")
            code = metric_soergel
        case ("kulczynski")
            code = metric_kulczynski
        case ("wave_hedges")
            code = metric_wave_hedges
        case ("motyka")
            code = metric_motyka
        case ("harmonic_mean")
            code = metric_harmonic_mean
        case default
            code = 0
        end select
    end function metric_code

    pure real(dp) function metric_distance(x, y, method, square, p) result(d)
        real(dp), intent(in) :: x(:), y(:)
        character(*), intent(in) :: method
        logical, intent(in), optional :: square
        real(dp), intent(in), optional :: p
        logical :: use_square
        real(dp) :: power
        integer :: code
        real(dp), allocatable :: sx(:), sy(:)

        use_square = .false.
        if (present(square)) use_square = square
        power = 0.0_dp
        if (present(p)) power = p
        code = metric_code(method)
        if (code == 0) then
            d = huge(0.0_dp)
            return
        end if
        if (size(x) /= size(y)) then
            d = huge(0.0_dp)
            return
        end if

        if (code == metric_hellinger) then
            if (any(x < 0.0_dp) .or. any(y < 0.0_dp)) then
                d = huge(0.0_dp)
                return
            end if
            allocate(sx(size(x)), sy(size(y)))
            sx = sqrt(x)
            sy = sqrt(y)
            d = metric_distance_prepared(sx, sy, code, use_square, power)
        else
            d = metric_distance_prepared(x, y, code, use_square, power)
        end if
    end function metric_distance

    pure real(dp) function metric_distance_prepared(x, y, code, square, p) result(d)
        real(dp), intent(in) :: x(:), y(:)
        integer, intent(in) :: code
        logical, intent(in) :: square
        real(dp), intent(in) :: p
        real(dp) :: s, denom, bc, v, sx, sy
        integer :: j

        select case (code)
        case (metric_euclidean)
            s = sum((x - y) ** 2)
            if (square) then
                d = s
            else
                d = sqrt(s)
            end if
        case (metric_hellinger)
            s = sum((x - y) ** 2)
            if (square) then
                d = 0.5_dp * s
            else
                d = sqrt(0.5_dp * s)
            end if
        case (metric_manhattan)
            d = sum(abs(x - y))
        case (metric_canberra)
            d = 0.0_dp
            do j = 1, size(x)
                denom = abs(x(j)) + abs(y(j))
                if (denom > 0.0_dp) then
                    d = d + abs(x(j) - y(j)) / denom
                else
                    d = ieee_value(d, ieee_quiet_nan)
                end if
            end do
        case (metric_kullback_leibler)
            d = 0.0_dp
            do j = 1, size(x)
                v = (y(j) - x(j)) * (log(y(j)) - log(x(j)))
                if (ieee_is_finite(v)) d = d + v
            end do
        case (metric_jensen_shannon)
            d = 0.0_dp
            do j = 1, size(x)
                v = xlogx(x(j)) + xlogx(y(j))
                if (x(j) + y(j) > 0.0_dp) then
                    v = v - (x(j) + y(j)) * log(0.5_dp * (x(j) + y(j)))
                end if
                if (ieee_is_finite(v) .and. v > 0.0_dp) d = d + v
            end do
        case (metric_itakura_saito)
            d = 0.0_dp
            do j = 1, size(x)
                v = x(j) / y(j) - (log(x(j)) - log(y(j))) - 1.0_dp
                if (ieee_is_finite(v)) d = d + v
            end do
        case (metric_bhattacharyya)
            bc = sum(sqrt(x * y))
            d = -log(bc)
        case (metric_jeffries_matusita)
            bc = sum(sqrt(x * y))
            d = sqrt(2.0_dp - 2.0_dp * bc)
        case (metric_minimum)
            d = minval(abs(x - y))
        case (metric_maximum)
            d = maxval(abs(x - y))
        case (metric_total_variation)
            d = 0.5_dp * sum(abs(x - y))
        case (metric_sorensen)
            d = sum(abs(x - y) / (x + y))
        case (metric_cosine)
            sx = sqrt(sum(x * x))
            sy = sqrt(sum(y * y))
            d = sum(x * y) / (sx * sy)
        case (metric_gower)
            d = sum(abs(x - y)) / real(size(x), dp)
        case (metric_minkowski)
            if (p <= 0.0_dp) then
                d = huge(0.0_dp)
            else
                d = sum(abs(x - y) ** p) ** (1.0_dp / p)
            end if
        case (metric_soergel)
            d = sum(abs(x - y)) / sum(max(x, y))
        case (metric_kulczynski)
            d = sum(abs(x - y)) / sum(min(x, y))
        case (metric_wave_hedges)
            d = sum(abs(x - y) / max(x, y))
        case (metric_motyka)
            d = 1.0_dp - sum(min(x, y)) / sum(x + y)
        case (metric_harmonic_mean)
            d = 2.0_dp * sum(x * y) / sum(x + y)
        case default
            d = huge(0.0_dp)
        end select
    end function metric_distance_prepared

    pure real(dp) function xlogx(x) result(v)
        real(dp), intent(in) :: x
        if (x > 0.0_dp) then
            v = x * log(x)
        else if (x >= 0.0_dp) then
            v = 0.0_dp
        else
            v = ieee_value(x, ieee_quiet_nan)
        end if
    end function xlogx

    pure function lower_ascii(text) result(out)
        character(*), intent(in) :: text
        character(len=len(text)) :: out
        integer :: i, c
        do i = 1, len(text)
            c = iachar(text(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) then
                out(i:i) = achar(c + 32)
            else
                out(i:i) = text(i:i)
            end if
        end do
    end function lower_ascii

end module rnanoflann_metrics

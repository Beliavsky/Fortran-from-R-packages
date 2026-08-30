module e1071_probplot
    use e1071_kinds, only: dp
    implicit none
    private

    abstract interface
        function quantile_callback(p) result(q)
            import dp
            real(dp), intent(in) :: p !! Probability in (0,1) whose distribution quantile is requested.
            real(dp) :: q
        end function quantile_callback
    end interface

    type, public :: probplot_distribution
        procedure(quantile_callback), pointer, nopass :: quantile => null()
    end type probplot_distribution

    type, public :: probplot_result
        real(dp), allocatable :: ordered(:)
        real(dp), allocatable :: theoretical(:)
        real(dp), allocatable :: probabilities(:)
        real(dp), allocatable :: probability_quantiles(:)
        real(dp) :: intercept = 0.0_dp
        real(dp) :: slope = 0.0_dp
    end type probplot_result

    public :: probplot_normal, probplot_custom

contains

    subroutine probplot_custom(x, distribution, result, probabilities)
        real(dp), intent(in) :: x(:) !! Sample values used to construct the probability plot; order is irrelevant.
        type(probplot_distribution), intent(in) :: distribution !! Quantile callback defining the reference distribution.
        type(probplot_result), intent(out) :: result !! Ordered sample, theoretical quantiles, ticks, and quartile line.
        real(dp), intent(in), optional :: probabilities(:) !! Optional reference probabilities strictly between zero and one.
        real(dp), allocatable :: probs(:)
        real(dp), allocatable :: ppoints(:)
        real(dp) :: qx(2)
        real(dp) :: qy(2)
        real(dp) :: a
        integer :: i
        integer :: n

        if (.not. associated(distribution%quantile)) error stop "probplot_custom: quantile callback is not associated"
        n = size(x)
        if (n < 2) error stop "probplot_custom: at least two observations are required"
        result%ordered = x
        call sort_real(result%ordered)
        call prepare_probabilities(n, probabilities, probs)
        result%probabilities = probs
        allocate(result%probability_quantiles(size(probs)), ppoints(n), result%theoretical(n))
        a = 0.5_dp
        if (n <= 10) a = 0.375_dp
        do i = 1, n
            ppoints(i) = (real(i, dp) - a) / (real(n, dp) + 1.0_dp - 2.0_dp * a)
            result%theoretical(i) = distribution%quantile(ppoints(i))
        end do
        do i = 1, size(probs)
            result%probability_quantiles(i) = distribution%quantile(probs(i))
        end do
        qx(1) = quantile_type7(result%ordered, 0.25_dp)
        qx(2) = quantile_type7(result%ordered, 0.75_dp)
        qy = [distribution%quantile(0.25_dp), distribution%quantile(0.75_dp)]
        if (abs(qx(2) - qx(1)) <= 0.0_dp) error stop "probplot_custom: sample interquartile range is zero"
        result%slope = (qy(2) - qy(1)) / (qx(2) - qx(1))
        result%intercept = qy(1) - result%slope * qx(1)
    end subroutine probplot_custom

    subroutine probplot_normal(x, result, probabilities)
        real(dp), intent(in) :: x(:) !! Sample values used to construct the normal probability plot; order is irrelevant.
        type(probplot_result), intent(out) :: result !! Ordered sample, normal scores, tick probabilities, and quartile
        !! reference line.
        real(dp), intent(in), optional :: probabilities(:) !! Optional probabilities in (0,1) used for labeled reference levels.
        real(dp), allocatable :: probs(:)
        real(dp), allocatable :: ppoints(:)
        real(dp) :: qx(2)
        real(dp) :: qy(2)
        real(dp) :: a
        integer :: i
        integer :: n

        n = size(x)
        if (n < 2) error stop "probplot_normal: at least two observations are required"
        result%ordered = x
        call sort_real(result%ordered)
        call prepare_probabilities(n, probabilities, probs)
        result%probabilities = probs
        allocate(result%probability_quantiles(size(probs)), ppoints(n), result%theoretical(n))
        a = 0.5_dp
        if (n <= 10) a = 0.375_dp
        do i = 1, n
            ppoints(i) = (real(i, dp) - a) / (real(n, dp) + 1.0_dp - 2.0_dp * a)
            result%theoretical(i) = normal_quantile(ppoints(i))
        end do
        do i = 1, size(probs)
            result%probability_quantiles(i) = normal_quantile(probs(i))
        end do
        qx(1) = quantile_type7(result%ordered, 0.25_dp)
        qx(2) = quantile_type7(result%ordered, 0.75_dp)
        qy = [normal_quantile(0.25_dp), normal_quantile(0.75_dp)]
        if (abs(qx(2) - qx(1)) <= 0.0_dp) error stop "probplot_normal: sample interquartile range is zero"
        result%slope = (qy(2) - qy(1)) / (qx(2) - qx(1))
        result%intercept = qy(1) - result%slope * qx(1)
    end subroutine probplot_normal

    subroutine prepare_probabilities(n, probabilities, probs)
        integer, intent(in) :: n !! Sample size controlling whether extreme 0.1% and 99.9% ticks are included.
        real(dp), intent(in), optional :: probabilities(:) !! Optional user-selected tick probabilities in the open unit interval.
        real(dp), allocatable, intent(out) :: probs(:) !! Validated probability ticks used by the probability-plot result.

        if (present(probabilities)) then
            if (any(probabilities <= 0.0_dp) .or. any(probabilities >= 1.0_dp)) then
                error stop "probplot: probabilities must lie strictly between zero and one"
            end if
            allocate(probs(size(probabilities)))
            probs(:) = probabilities
        else if (n >= 1000) then
            allocate(probs(15))
            probs(:) = [0.001_dp, 0.01_dp, 0.05_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, &
                        0.6_dp, 0.7_dp, 0.8_dp, 0.9_dp, 0.95_dp, 0.99_dp, 0.999_dp]
        else
            allocate(probs(13))
            probs(:) = [0.01_dp, 0.05_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.6_dp, &
                        0.7_dp, 0.8_dp, 0.9_dp, 0.95_dp, 0.99_dp]
        end if
    end subroutine prepare_probabilities

    pure function normal_quantile(p) result(x)
        real(dp), intent(in) :: p !! Probability strictly between zero and one whose standard-normal quantile is returned.
        real(dp) :: x
        real(dp), parameter :: a1 = -3.969683028665376e1_dp
        real(dp), parameter :: a2 = 2.209460984245205e2_dp
        real(dp), parameter :: a3 = -2.759285104469687e2_dp
        real(dp), parameter :: a4 = 1.383577518672690e2_dp
        real(dp), parameter :: a5 = -3.066479806614716e1_dp
        real(dp), parameter :: a6 = 2.506628277459239_dp
        real(dp), parameter :: b1 = -5.447609879822406e1_dp
        real(dp), parameter :: b2 = 1.615858368580409e2_dp
        real(dp), parameter :: b3 = -1.556989798598866e2_dp
        real(dp), parameter :: b4 = 6.680131188771972e1_dp
        real(dp), parameter :: b5 = -1.328068155288572e1_dp
        real(dp), parameter :: c1 = -7.784894002430293e-3_dp
        real(dp), parameter :: c2 = -3.223964580411365e-1_dp
        real(dp), parameter :: c3 = -2.400758277161838_dp
        real(dp), parameter :: c4 = -2.549732539343734_dp
        real(dp), parameter :: c5 = 4.374664141464968_dp
        real(dp), parameter :: c6 = 2.938163982698783_dp
        real(dp), parameter :: d1 = 7.784695709041462e-3_dp
        real(dp), parameter :: d2 = 3.224671290700398e-1_dp
        real(dp), parameter :: d3 = 2.445134137142996_dp
        real(dp), parameter :: d4 = 3.754408661907416_dp
        real(dp) :: q
        real(dp) :: r

        if (p <= 0.0_dp .or. p >= 1.0_dp) error stop "normal_quantile: p must lie in (0,1)"
        if (p < 0.02425_dp) then
            q = sqrt(-2.0_dp * log(p))
            x = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / &
                ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0_dp)
        else if (p > 0.97575_dp) then
            q = sqrt(-2.0_dp * log(1.0_dp - p))
            x = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / &
                 ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0_dp)
        else
            q = p - 0.5_dp
            r = q * q
            x = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q / &
                (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0_dp)
        end if
    end function normal_quantile

    pure function quantile_type7(sorted, probability) result(value)
        real(dp), intent(in) :: sorted(:) !! Sample sorted in ascending order before applying R's default type-7 interpolation.
        real(dp), intent(in) :: probability !! Probability in [0,1] at which the sample quantile is evaluated.
        real(dp) :: value
        real(dp) :: h
        real(dp) :: fraction
        integer :: lower

        h = 1.0_dp + (real(size(sorted), dp) - 1.0_dp) * probability
        lower = int(floor(h))
        fraction = h - real(lower, dp)
        if (lower >= size(sorted)) then
            value = sorted(size(sorted))
        else
            value = sorted(lower) + fraction * (sorted(lower + 1) - sorted(lower))
        end if
    end function quantile_type7

    subroutine sort_real(x)
        real(dp), intent(inout) :: x(:) !! Real vector sorted in ascending order in place.
        real(dp) :: key
        integer :: i
        integer :: j

        do i = 2, size(x)
            key = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= key) exit
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = key
        end do
    end subroutine sort_real

end module e1071_probplot

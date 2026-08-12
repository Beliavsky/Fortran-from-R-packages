module ppso_benchmarks
    use ppso_kinds, only : dp
    implicit none
    private
    public :: rastrigin_function, ackley_function, griewank_function
    public :: sample_function, sample_function2

contains

    function rastrigin_function(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        real(dp), parameter :: twopi = 6.2831853071795864769252867665590058_dp
        f = sum(x*x - cos(twopi*x))
    end function rastrigin_function

    function ackley_function(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        real(dp), parameter :: twopi = 6.2831853071795864769252867665590058_dp
        real(dp) :: s1, s2
        if (size(x) == 0) then
            f = 0.0_dp
            return
        end if
        s1 = sum(x*x)
        s2 = sum(cos(twopi*x))
        f = -20.0_dp*exp(-0.2_dp*sqrt(s1/real(size(x),dp))) &
            - exp(s2/real(size(x),dp))
    end function ackley_function

    function griewank_function(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        integer :: i
        f = sum(x*x/4000.0_dp) - product([(cos(x(i)/sqrt(real(i,dp))), i=1,size(x))]) + 1.0_dp
    end function griewank_function

    function sample_function(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = -abs(sum(cos(x) + 1.0_dp))
    end function sample_function

    function sample_function2(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f, rr
        if (size(x) < 2) error stop "sample_function2 requires at least two parameters"
        rr = max(1.0e-320_dp, sqrt((4.0_dp*x(1))**2 + (4.0_dp*x(2))**2))
        f = -5.0_dp*sin(rr)/rr
    end function sample_function2

end module ppso_benchmarks

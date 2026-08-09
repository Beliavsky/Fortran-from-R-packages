! SPDX-License-Identifier: GPL-3.0-or-later
!
! Modern Fortran translation of the computational code in globalOptTests 1.1.
module global_opt_tests
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)
    real(dp), parameter :: pi_obj = 3.14159265359_dp
    real(dp), parameter :: pi_meta = 3.1415926535897932384626433832795_dp
    integer, parameter, public :: n_benchmarks = 50

    character(len=16), parameter, public :: benchmark_names(n_benchmarks) = [character(len=16) :: &
        'Ackleys', 'AluffiPentini', 'BeckerLago', 'Bohachevsky1', 'Bohachevsky2', &
        'Branin', 'Camel3', 'Camel6', 'CosMix2', 'CosMix4', 'DekkersAarts', 'Easom', &
        'EMichalewicz', 'Expo', 'GoldPrice', 'Griewank', 'Gulf', 'Hartman3', 'Hartman6', &
        'Hosaki', 'Kowalik', 'LM1', 'LM2n10', 'LM2n5', 'McCormic', 'MeyerRoth', &
        'MieleCantrell', 'Modlangerman', 'ModRosenbrock', 'MultiGauss', 'Neumaier2', &
        'Neumaier3', 'Paviani', 'Periodic', 'PowellQ', 'PriceTransistor', 'Rastrigin', &
        'Rosenbrock', 'Salomon', 'Schaffer1', 'Schaffer2', 'Schubert', 'Schwefel', &
        'Shekel10', 'Shekel5', 'Shekel7', 'Shekelfox5', 'Wood', 'Zeldasine10', 'Zeldasine20']

    public :: go_test, get_default_bounds, get_problem_dimension, get_global_opt
    public :: ackleys, aluffipentini, beckerlago, bohachevsky1, bohachevsky2
    public :: branin, camel3, camel6, cosmix2, cosmix4, dekkersaarts, easom
    public :: emichalewicz, expo, goldprice, griewank, gulf, hartman3, hartman6
    public :: hosaki, kowalik, lm1, lm2n10, lm2n5, mccormic, meyerroth
    public :: mielecantrell, modlangerman, modrosenbrock, multigauss, neumaier2
    public :: neumaier3, paviani, periodic, powellq, pricetransistor, rastrigin
    public :: rosenbrock, salomon, schaffer1, schaffer2, schubert, schwefel
    public :: shekel10, shekel5, shekel7, shekelfox5, wood, zeldasine10, zeldasine20

contains

    function go_test(x, fn_name, check_dim, status) result(value)
        real(dp), intent(in) :: x(:)
        character(len=*), intent(in) :: fn_name
        logical, intent(in), optional :: check_dim
        integer, intent(out), optional :: status
        real(dp) :: value
        logical :: do_check
        integer :: expected, istat

        do_check = .true.
        if (present(check_dim)) do_check = check_dim
        istat = 0

        expected = get_problem_dimension(fn_name, istat)
        if (istat /= 0) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            if (present(status)) status = 2
            return
        end if
        if (do_check .and. size(x) /= expected) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            if (present(status)) status = 1
            return
        end if

        select case (trim(fn_name))
        case ('Ackleys');         value = ackleys(x)
        case ('AluffiPentini');   value = aluffipentini(x)
        case ('BeckerLago');      value = beckerlago(x)
        case ('Bohachevsky1');    value = bohachevsky1(x)
        case ('Bohachevsky2');    value = bohachevsky2(x)
        case ('Branin');          value = branin(x)
        case ('Camel3');          value = camel3(x)
        case ('Camel6');          value = camel6(x)
        case ('CosMix2');         value = cosmix2(x)
        case ('CosMix4');         value = cosmix4(x)
        case ('DekkersAarts');    value = dekkersaarts(x)
        case ('Easom');           value = easom(x)
        case ('EMichalewicz');    value = emichalewicz(x)
        case ('Expo');            value = expo(x)
        case ('GoldPrice');       value = goldprice(x)
        case ('Griewank');        value = griewank(x)
        case ('Gulf');            value = gulf(x)
        case ('Hartman3');        value = hartman3(x)
        case ('Hartman6');        value = hartman6(x)
        case ('Hosaki');          value = hosaki(x)
        case ('Kowalik');         value = kowalik(x)
        case ('LM1');             value = lm1(x)
        case ('LM2n10');          value = lm2n10(x)
        case ('LM2n5');           value = lm2n5(x)
        case ('McCormic');        value = mccormic(x)
        case ('MeyerRoth');       value = meyerroth(x)
        case ('MieleCantrell');   value = mielecantrell(x)
        case ('Modlangerman');    value = modlangerman(x)
        case ('ModRosenbrock');   value = modrosenbrock(x)
        case ('MultiGauss');      value = multigauss(x)
        case ('Neumaier2');       value = neumaier2(x)
        case ('Neumaier3');       value = neumaier3(x)
        case ('Paviani');         value = paviani(x)
        case ('Periodic');        value = periodic(x)
        case ('PowellQ');         value = powellq(x)
        case ('PriceTransistor'); value = pricetransistor(x)
        case ('Rastrigin');       value = rastrigin(x)
        case ('Rosenbrock');      value = rosenbrock(x)
        case ('Salomon');         value = salomon(x)
        case ('Schaffer1');       value = schaffer1(x)
        case ('Schaffer2');       value = schaffer2(x)
        case ('Schubert');        value = schubert(x)
        case ('Schwefel');        value = schwefel(x)
        case ('Shekel10');        value = shekel10(x)
        case ('Shekel5');         value = shekel5(x)
        case ('Shekel7');         value = shekel7(x)
        case ('Shekelfox5');      value = shekelfox5(x)
        case ('Wood');            value = wood(x)
        case ('Zeldasine10');     value = zeldasine10(x)
        case ('Zeldasine20');     value = zeldasine20(x)
        case default
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            istat = 2
        end select
        if (present(status)) status = istat
    end function go_test

    subroutine get_default_bounds(fn_name, lower, upper, status)
        character(len=*), intent(in) :: fn_name
        real(dp), allocatable, intent(out) :: lower(:), upper(:)
        integer, intent(out), optional :: status
        integer :: istat

        istat = 0
        select case (trim(fn_name))
        case ('Ackleys');         call constant_bounds(10, -35.0_dp, 30.0_dp, lower, upper)
        case ('AluffiPentini');   call constant_bounds(2, -12.0_dp, 10.0_dp, lower, upper)
        case ('BeckerLago');      call constant_bounds(2, -12.0_dp, 10.0_dp, lower, upper)
        case ('Bohachevsky1');    call constant_bounds(2, -55.0_dp, 50.0_dp, lower, upper)
        case ('Bohachevsky2');    call constant_bounds(2, -55.0_dp, 50.0_dp, lower, upper)
        case ('Branin')
            allocate(lower(2), upper(2)); lower = [-5.0_dp, 0.0_dp]; upper = [10.0_dp, 15.0_dp]
        case ('Camel3');          call constant_bounds(2, -8.0_dp, 5.0_dp, lower, upper)
        case ('Camel6');          call constant_bounds(2, -8.0_dp, 5.0_dp, lower, upper)
        case ('CosMix2');         call constant_bounds(2, -2.0_dp, 1.0_dp, lower, upper)
        case ('CosMix4');         call constant_bounds(4, -2.0_dp, 1.0_dp, lower, upper)
        case ('DekkersAarts');    call constant_bounds(2, -25.0_dp, 20.0_dp, lower, upper)
        case ('Easom')
            allocate(lower(2), upper(2)); lower = [-12.0_dp, -12.0_dp]; upper = [10.0_dp, 2.0_dp]
        case ('EMichalewicz');    call constant_bounds(5, 0.0_dp, pi_meta, lower, upper)
        case ('Expo');            call constant_bounds(10, -12.0_dp, 10.0_dp, lower, upper)
        case ('GoldPrice');       call constant_bounds(2, -3.0_dp, 2.0_dp, lower, upper)
        case ('Griewank');        call constant_bounds(10, -550.0_dp, 500.0_dp, lower, upper)
        case ('Gulf')
            allocate(lower(3), upper(3)); lower = [0.1_dp, 0.0_dp, 0.0_dp]; upper = [100.0_dp, 25.6_dp, 5.0_dp]
        case ('Hartman3');        call constant_bounds(3, 0.0_dp, 1.0_dp, lower, upper)
        case ('Hartman6');        call constant_bounds(6, 0.0_dp, 1.0_dp, lower, upper)
        case ('Hosaki')
            allocate(lower(2), upper(2)); lower = [0.0_dp, 0.0_dp]; upper = [5.0_dp, 6.0_dp]
        case ('Kowalik');         call constant_bounds(4, 0.0_dp, 0.42_dp, lower, upper)
        case ('LM1');             call constant_bounds(3, -15.0_dp, 10.0_dp, lower, upper)
        case ('LM2n10');          call constant_bounds(10, -10.0_dp, 5.0_dp, lower, upper)
        case ('LM2n5');           call constant_bounds(5, -10.0_dp, 5.0_dp, lower, upper)
        case ('McCormic')
            allocate(lower(2), upper(2)); lower = [-1.5_dp, -3.0_dp]; upper = [4.0_dp, 3.0_dp]
        case ('MeyerRoth');       call constant_bounds(3, -10.0_dp, 10.0_dp, lower, upper)
        case ('MieleCantrell');   call constant_bounds(4, -1.5_dp, 1.0_dp, lower, upper)
        case ('Modlangerman');    call constant_bounds(10, 0.0_dp, 10.0_dp, lower, upper)
        case ('ModRosenbrock')
            allocate(lower(2), upper(2)); lower = [-7.0_dp, -2.0_dp]; upper = [5.0_dp, 2.0_dp]
        case ('MultiGauss')
            allocate(lower(2), upper(2)); lower = [-3.0_dp, -2.0_dp]; upper = [2.0_dp, 2.0_dp]
        case ('Neumaier2')
            allocate(lower(4), upper(4)); lower = 0.0_dp; upper = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
        case ('Neumaier3');       call constant_bounds(10, -115.0_dp, 100.0_dp, lower, upper)
        case ('Paviani');         call constant_bounds(10, 2.0_dp, 10.0_dp, lower, upper)
        case ('Periodic');        call constant_bounds(2, -15.0_dp, 10.0_dp, lower, upper)
        case ('PowellQ');         call constant_bounds(4, -15.0_dp, 10.0_dp, lower, upper)
        case ('PriceTransistor'); call constant_bounds(9, 0.0_dp, 10.0_dp, lower, upper)
        case ('Rastrigin');       call constant_bounds(10, -525.0_dp, 512.0_dp, lower, upper)
        case ('Rosenbrock');      call constant_bounds(10, -40.0_dp, 30.0_dp, lower, upper)
        case ('Salomon');         call constant_bounds(5, -120.0_dp, 100.0_dp, lower, upper)
        case ('Schaffer1');       call constant_bounds(2, -120.0_dp, 100.0_dp, lower, upper)
        case ('Schaffer2');       call constant_bounds(2, -120.0_dp, 100.0_dp, lower, upper)
        case ('Schubert');        call constant_bounds(2, -15.0_dp, 10.0_dp, lower, upper)
        case ('Schwefel');        call constant_bounds(10, -500.0_dp, 500.0_dp, lower, upper)
        case ('Shekel10');        call constant_bounds(4, 0.0_dp, 10.0_dp, lower, upper)
        case ('Shekel5');         call constant_bounds(4, 0.0_dp, 10.0_dp, lower, upper)
        case ('Shekel7');         call constant_bounds(4, 0.0_dp, 10.0_dp, lower, upper)
        case ('Shekelfox5');      call constant_bounds(5, 0.0_dp, 10.0_dp, lower, upper)
        case ('Wood');            call constant_bounds(4, -14.0_dp, 10.0_dp, lower, upper)
        case ('Zeldasine10');     call constant_bounds(10, 0.0_dp, pi_meta, lower, upper)
        case ('Zeldasine20');     call constant_bounds(20, 0.0_dp, pi_meta, lower, upper)
        case default
            allocate(lower(0), upper(0)); istat = 2
        end select
        if (present(status)) status = istat
    end subroutine get_default_bounds

    integer function get_problem_dimension(fn_name, status) result(n)
        character(len=*), intent(in) :: fn_name
        integer, intent(out), optional :: status
        real(dp), allocatable :: lower(:), upper(:)
        integer :: istat
        call get_default_bounds(fn_name, lower, upper, istat)
        n = size(lower)
        if (present(status)) status = istat
    end function get_problem_dimension

    real(dp) function get_global_opt(fn_name, status) result(value)
        character(len=*), intent(in) :: fn_name
        integer, intent(out), optional :: status
        integer :: istat
        istat = 0
        select case (trim(fn_name))
        case ('Ackleys');         value = 0.0_dp
        case ('AluffiPentini');   value = -0.3523_dp
        case ('BeckerLago');      value = 0.0_dp
        case ('Bohachevsky1');    value = 0.0_dp
        case ('Bohachevsky2');    value = 0.0_dp
        case ('Branin');          value = 0.3979_dp
        case ('Camel3');          value = 0.0_dp
        case ('Camel6');          value = -1.0316_dp
        case ('CosMix2');         value = -0.2_dp
        case ('CosMix4');         value = -0.4_dp
        case ('DekkersAarts');    value = -24776.5183_dp
        case ('Easom');           value = -1.0_dp
        case ('EMichalewicz');    value = -4.6877_dp
        case ('Expo');            value = -1.0_dp
        case ('GoldPrice');       value = 3.0_dp
        case ('Griewank');        value = 0.0_dp
        case ('Gulf');            value = 0.0_dp
        case ('Hartman3');        value = -3.8628_dp
        case ('Hartman6');        value = -3.3224_dp
        case ('Hosaki');          value = -2.3458_dp
        case ('Kowalik');         value = 0.0003_dp
        case ('LM1');             value = 0.0_dp
        case ('LM2n10');          value = 0.0_dp
        case ('LM2n5');           value = 0.0_dp
        case ('McCormic');        value = -1.9133_dp
        case ('MeyerRoth');       value = 4.355628e-5_dp
        case ('MieleCantrell');   value = 0.0_dp
        case ('Modlangerman');    value = -0.9650_dp
        case ('ModRosenbrock');   value = 0.0_dp
        case ('MultiGauss');      value = -1.2970_dp
        case ('Neumaier2');       value = 0.0_dp
        case ('Neumaier3');       value = -210.0_dp
        case ('Paviani');         value = -45.7784_dp
        case ('Periodic');        value = 0.9000_dp
        case ('PowellQ');         value = 0.0_dp
        case ('PriceTransistor'); value = 0.0_dp
        case ('Rastrigin');       value = 0.0_dp
        case ('Rosenbrock');      value = 0.0_dp
        case ('Salomon');         value = 0.0_dp
        case ('Schaffer1');       value = 0.0_dp
        case ('Schaffer2');       value = 0.0_dp
        case ('Schubert');        value = -186.7309_dp
        case ('Schwefel');        value = -4189.8289_dp
        case ('Shekel10');        value = -10.5364_dp
        case ('Shekel5');         value = -10.1532_dp
        case ('Shekel7');         value = -10.4029_dp
        case ('Shekelfox5');      value = -10.4056_dp
        case ('Wood');            value = 0.0_dp
        case ('Zeldasine10');     value = -3.5_dp
        case ('Zeldasine20');     value = -3.5_dp
        case default
            value = ieee_value(0.0_dp, ieee_quiet_nan); istat = 2
        end select
        if (present(status)) status = istat
    end function get_global_opt

    subroutine constant_bounds(n, lo, hi, lower, upper)
        integer, intent(in) :: n
        real(dp), intent(in) :: lo, hi
        real(dp), allocatable, intent(out) :: lower(:), upper(:)
        allocate(lower(n), upper(n))
        lower = lo
        upper = hi
    end subroutine constant_bounds

    pure real(dp) function ackleys(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: s1, s2
        integer :: n
        n = size(x)
        s1 = sum(x*x)
        s2 = sum(cos(2.0_dp*pi_obj*x))
        f = -20.0_dp*exp(-0.2_dp*sqrt(s1/real(n,dp))) - exp(s2/real(n,dp)) + 20.0_dp + exp(1.0_dp)
    end function ackleys

    pure real(dp) function aluffipentini(x) result(f)
        real(dp), intent(in) :: x(:)
        f = 0.25_dp*x(1)**4 - 0.5_dp*x(1)**2 + 0.1_dp*x(1) + 0.5_dp*x(2)**2
    end function aluffipentini

    pure real(dp) function beckerlago(x) result(f)
        real(dp), intent(in) :: x(:)
        f = (abs(x(1))-5.0_dp)**2 + (abs(x(2))-5.0_dp)**2
    end function beckerlago

    pure real(dp) function bohachevsky1(x) result(f)
        real(dp), intent(in) :: x(:)
        f = x(1)**2 + 2.0_dp*x(2)**2 - 0.3_dp*cos(3.0_dp*pi_obj*x(1)) &
            - 0.4_dp*cos(4.0_dp*pi_obj*x(2)) + 0.7_dp
    end function bohachevsky1

    pure real(dp) function bohachevsky2(x) result(f)
        real(dp), intent(in) :: x(:)
        f = x(1)**2 + 2.0_dp*x(2)**2 &
            - 0.3_dp*cos(3.0_dp*pi_obj*x(1))*cos(4.0_dp*pi_obj*x(2)) + 0.3_dp
    end function bohachevsky2

    pure real(dp) function branin(x) result(fval)
        real(dp), intent(in) :: x(:)
        real(dp), parameter :: a=1.0_dp, b=5.1_dp/(4.0_dp*pi_obj*pi_obj)
        real(dp), parameter :: c=5.0_dp/pi_obj, d=6.0_dp, e=10.0_dp
        real(dp), parameter :: ff=1.0_dp/(8.0_dp*pi_obj)
        fval = a*(x(2)-b*x(1)**2+c*x(1)-d)**2 + e*(1.0_dp-ff)*cos(x(1)) + e
    end function branin

    pure real(dp) function camel3(x) result(f)
        real(dp), intent(in) :: x(:)
        f = (2.0_dp-1.05_dp*x(1)**2+x(1)**4/6.0_dp)*x(1)**2 + x(1)*x(2) + x(2)**2
    end function camel3

    pure real(dp) function camel6(x) result(f)
        real(dp), intent(in) :: x(:)
        f = (4.0_dp-2.1_dp*x(1)**2+x(1)**4/3.0_dp)*x(1)**2 + x(1)*x(2) &
            + (-4.0_dp+4.0_dp*x(2)**2)*x(2)**2
    end function camel6

    pure real(dp) function cosmix2(x) result(f)
        real(dp), intent(in) :: x(:)
        f = sum(x*x) - 0.1_dp*sum(cos(5.0_dp*pi_obj*x))
    end function cosmix2

    pure real(dp) function cosmix4(x) result(f)
        real(dp), intent(in) :: x(:)
        f = sum(x*x) - 0.1_dp*sum(cos(5.0_dp*pi_obj*x))
    end function cosmix4

    pure real(dp) function dekkersaarts(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: r2
        r2 = x(1)**2 + x(2)**2
        f = 100000.0_dp*x(1)**2 + x(2)**2 - r2**2 + r2**4/100000.0_dp
    end function dekkersaarts

    pure real(dp) function easom(x) result(f)
        real(dp), intent(in) :: x(:)
        f = -cos(x(1))*cos(x(2))*exp(-(x(1)-pi_obj)**2-(x(2)-pi_obj)**2)
    end function easom

    pure real(dp) function emichalewicz(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: y(size(x)), cost, sint
        integer :: j, n
        n = size(x)
        cost = cos(pi_obj/6.0_dp)
        sint = sin(pi_obj/6.0_dp)
        j = 1
        do while (j < n)
            y(j) = x(j)*cost - x(j+1)*sint
            y(j+1) = x(j)*sint + x(j+1)*cost
            j = j + 2
        end do
        if (j == n) y(j) = x(j)
        f = 0.0_dp
        do j = 1, n
            f = f - sin(y(j))*sin(real(j,dp)*y(j)*y(j)/pi_obj)**20
        end do
    end function emichalewicz

    pure real(dp) function expo(x) result(f)
        real(dp), intent(in) :: x(:)
        f = -exp(-0.5_dp*sum(x*x))
    end function expo

    pure real(dp) function goldprice(x) result(f)
        real(dp), intent(in) :: x(:)
        f = 1.0_dp + (x(1)+x(2)+1.0_dp)**2 * &
            (19.0_dp-14.0_dp*x(1)+3.0_dp*x(1)**2-14.0_dp*x(2)+6.0_dp*x(1)*x(2)+3.0_dp*x(2)**2)
        f = f * (30.0_dp + (2.0_dp*x(1)-3.0_dp*x(2))**2 * &
            (18.0_dp-32.0_dp*x(1)+12.0_dp*x(1)**2+48.0_dp*x(2) &
            -36.0_dp*x(1)*x(2)+27.0_dp*x(2)**2))
    end function goldprice

    pure real(dp) function griewank(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: prod
        integer :: j
        prod = 1.0_dp
        do j = 1, size(x)
            prod = prod*cos(x(j)/sqrt(real(j,dp)))
        end do
        f = sum(x*x)/4000.0_dp - prod + 1.0_dp
    end function griewank

    pure real(dp) function gulf(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: u, r
        integer :: j
        ! The original C starts at j=0 and evaluates log(0); that term tends to zero.
        f = 0.0_dp
        do j = 1, 98
            u = 25.0_dp + (-50.0_dp*log(0.01_dp*real(j,dp)))**0.66666_dp
            r = exp(-(u-x(2))**x(3)/x(1)) - 0.01_dp*real(j,dp)
            f = f + r*r
        end do
    end function gulf

    pure real(dp) function hartman3(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp), parameter :: a(3,4) = reshape([ &
            3.0_dp,10.0_dp,30.0_dp, 0.1_dp,10.0_dp,35.0_dp, &
            3.0_dp,10.0_dp,30.0_dp, 0.1_dp,10.0_dp,35.0_dp], [3,4])
        real(dp), parameter :: c(4) = [1.0_dp,1.2_dp,3.0_dp,3.2_dp]
        real(dp), parameter :: p(3,4) = reshape([ &
            0.3689_dp,0.117_dp,0.2673_dp, 0.4699_dp,0.4387_dp,0.747_dp, &
            0.1091_dp,0.8732_dp,0.5547_dp, 0.03815_dp,0.5743_dp,0.8828_dp], [3,4])
        real(dp) :: s
        integer :: i
        f = 0.0_dp
        do i = 1, 4
            s = -sum(a(:,i)*(x-p(:,i))**2)
            f = f - c(i)*exp(s)
        end do
    end function hartman3

    pure real(dp) function hartman6(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp), parameter :: a(6,4) = reshape([ &
            10.0_dp,3.0_dp,17.0_dp,3.5_dp,1.7_dp,8.0_dp, &
            0.05_dp,10.0_dp,17.0_dp,0.1_dp,8.0_dp,14.0_dp, &
            3.0_dp,3.5_dp,1.7_dp,10.0_dp,17.0_dp,8.0_dp, &
            17.0_dp,8.0_dp,0.05_dp,10.0_dp,0.1_dp,14.0_dp], [6,4])
        real(dp), parameter :: c(4) = [1.0_dp,1.2_dp,3.0_dp,3.2_dp]
        real(dp), parameter :: p(6,4) = reshape([ &
            0.1312_dp,0.1696_dp,0.5569_dp,0.0124_dp,0.8283_dp,0.5886_dp, &
            0.2329_dp,0.4135_dp,0.8307_dp,0.3736_dp,0.1004_dp,0.9991_dp, &
            0.2348_dp,0.1451_dp,0.3522_dp,0.2883_dp,0.3047_dp,0.6650_dp, &
            0.4047_dp,0.8828_dp,0.8732_dp,0.5743_dp,0.1091_dp,0.0381_dp], [6,4])
        real(dp) :: s
        integer :: i
        f = 0.0_dp
        do i = 1, 4
            s = -sum(a(:,i)*(x-p(:,i))**2)
            f = f - c(i)*exp(s)
        end do
    end function hartman6

    pure real(dp) function hosaki(x) result(f)
        real(dp), intent(in) :: x(:)
        f = (1.0_dp-8.0_dp*x(1)+7.0_dp*x(1)**2-(7.0_dp/3.0_dp)*x(1)**3 &
            +0.25_dp*x(1)**4)*x(2)**2*exp(-x(2))
    end function hosaki

    pure real(dp) function kowalik(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp), parameter :: a(11) = [0.1957_dp,0.1947_dp,0.1735_dp,0.16_dp,0.0844_dp, &
            0.0627_dp,0.0456_dp,0.0342_dp,0.0323_dp,0.0235_dp,0.0246_dp]
        real(dp), parameter :: b(11) = [0.25_dp,0.5_dp,1.0_dp,2.0_dp,4.0_dp,6.0_dp, &
            8.0_dp,10.0_dp,12.0_dp,14.0_dp,16.0_dp]
        integer :: i
        f = 0.0_dp
        do i = 1, 11
            f = f + (a(i)-x(1)*(1.0_dp+x(2)*b(i))/(1.0_dp+x(3)*b(i)+x(4)*b(i)**2))**2
        end do
    end function kowalik

    pure real(dp) function lm1(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: z1, z2
        integer :: j, n
        n = size(x)
        z1 = 10.0_dp*sin(pi_obj*(1.0_dp+0.25_dp*(x(1)+1.0_dp)))**2
        z2 = (0.25_dp*(x(n)+1.0_dp))**2
        do j = 1, n-1
            z1 = z1 + (0.25_dp*(x(j)+1.0_dp))**2 * &
                (1.0_dp+sin(pi_obj*(0.25_dp*x(j+1)))**2)
        end do
        f = pi_obj/real(n,dp)*(z1+z2)
    end function lm1

    pure real(dp) function lm2n10(x) result(f)
        real(dp), intent(in) :: x(:)
        f = lm2_common(x)
    end function lm2n10

    pure real(dp) function lm2n5(x) result(f)
        real(dp), intent(in) :: x(:)
        f = lm2_common(x)
    end function lm2n5

    pure real(dp) function lm2_common(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: z1, z2
        integer :: j, n
        n = size(x)
        z1 = sin(3.0_dp*pi_obj*x(1))**2
        z2 = (x(n)-1.0_dp)**2*(1.0_dp+sin(2.0_dp*pi_obj*x(n))**2)
        do j = 1, n-1
            z1 = z1 + (x(j)-1.0_dp)**2*(1.0_dp+sin(3.0_dp*pi_obj*x(j+1))**2)
        end do
        f = 0.1_dp*(z1+z2)
    end function lm2_common

    pure real(dp) function mccormic(x) result(f)
        real(dp), intent(in) :: x(:)
        f = sin(x(1)+x(2)) + (x(1)-x(2))**2 - 1.5_dp*x(1) + 2.5_dp*x(2) + 1.0_dp
    end function mccormic

    pure real(dp) function meyerroth(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp), parameter :: t(5) = [1.0_dp,2.0_dp,1.0_dp,2.0_dp,0.1_dp]
        real(dp), parameter :: v(5) = [1.0_dp,1.0_dp,2.0_dp,2.0_dp,0.0_dp]
        real(dp), parameter :: y(5) = [0.126_dp,0.219_dp,0.076_dp,0.126_dp,0.186_dp]
        real(dp) :: num, den
        integer :: i
        f = 0.0_dp
        do i = 1, 5
            num = x(1)*x(3)*t(i)
            den = 1.0_dp+x(1)*t(i)+x(2)*v(i)
            f = f + (num/den-y(i))**2
        end do
    end function meyerroth

    pure real(dp) function mielecantrell(x) result(f)
        real(dp), intent(in) :: x(:)
        f = (exp(x(1))-x(2))**4 + 100.0_dp*(x(2)-x(3))**6 &
            + tan(x(3)-x(4))**4 + x(1)**8
    end function mielecantrell

    pure real(dp) function modlangerman(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp), parameter :: a(10,5) = reshape([ &
            9.681_dp,0.667_dp,4.783_dp,9.095_dp,3.517_dp,9.325_dp,6.544_dp,0.211_dp,5.122_dp,2.020_dp, &
            9.400_dp,2.041_dp,3.788_dp,7.931_dp,2.882_dp,2.672_dp,3.568_dp,1.284_dp,7.033_dp,7.374_dp, &
            8.025_dp,9.152_dp,5.114_dp,7.621_dp,4.564_dp,4.711_dp,2.996_dp,6.126_dp,0.734_dp,4.982_dp, &
            2.196_dp,0.415_dp,5.649_dp,6.979_dp,9.510_dp,9.166_dp,6.304_dp,6.054_dp,9.377_dp,1.426_dp, &
            8.074_dp,8.777_dp,3.467_dp,1.867_dp,6.708_dp,6.349_dp,4.534_dp,0.276_dp,7.633_dp,1.567_dp], [10,5])
        real(dp), parameter :: c(5) = [0.806_dp,0.517_dp,0.1_dp,0.908_dp,0.965_dp]
        real(dp) :: dist
        integer :: i
        f = 0.0_dp
        do i = 1, 5
            dist = sum((x-a(:,i))**2)
            f = f - c(i)*exp(-dist/pi_obj)*cos(pi_obj*dist)
        end do
    end function modlangerman

    pure real(dp) function modrosenbrock(x) result(f)
        real(dp), intent(in) :: x(:)
        f = 100.0_dp*(x(2)-x(1)**2)**2 + (6.4_dp*(x(2)-0.5_dp)**2-x(1)-0.6_dp)**2
    end function modrosenbrock

    pure real(dp) function multigauss(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp), parameter :: a(5) = [0.5_dp,1.2_dp,1.0_dp,1.0_dp,1.2_dp]
        real(dp), parameter :: b(5) = [0.0_dp,1.0_dp,0.0_dp,-0.5_dp,0.0_dp]
        real(dp), parameter :: c(5) = [0.0_dp,0.0_dp,-0.5_dp,0.0_dp,1.0_dp]
        real(dp), parameter :: d(5) = [0.1_dp,0.5_dp,0.5_dp,0.5_dp,0.5_dp]
        integer :: i
        f = 0.0_dp
        do i = 1, 5
            f = f - a(i)*exp(-((x(1)-b(i))**2+(x(2)-c(i))**2)/d(i)**2)
        end do
    end function multigauss

    pure real(dp) function neumaier2(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp), parameter :: b(4) = [8.0_dp,18.0_dp,44.0_dp,114.0_dp]
        real(dp) :: s
        integer :: i, k
        f = 0.0_dp
        do k = 1, size(x)
            s = 0.0_dp
            do i = 1, size(x)
                s = s + x(i)**k
            end do
            f = f + (b(k)-s)**2
        end do
    end function neumaier2

    pure real(dp) function neumaier3(x) result(f)
        real(dp), intent(in) :: x(:)
        integer :: i
        f = sum((x-1.0_dp)**2)
        do i = 2, size(x)
            f = f - x(i)*x(i-1)
        end do
    end function neumaier3

    pure real(dp) function paviani(x) result(f)
        real(dp), intent(in) :: x(:)
        f = sum(log(x-2.0_dp)**2 + log(10.0_dp-x)**2) - product(x)**0.2_dp
    end function paviani

    pure real(dp) function periodic(x) result(f)
        real(dp), intent(in) :: x(:)
        f = 1.0_dp + sin(x(1))**2 + sin(x(2))**2 - 0.1_dp*exp(-x(1)**2-x(2)**2)
    end function periodic

    pure real(dp) function powellq(x) result(f)
        real(dp), intent(in) :: x(:)
        ! Preserves the upstream expression (x1 + 10*x1)^2 exactly.
        f = (x(1)+10.0_dp*x(1))**2 + 5.0_dp*(x(3)-x(4))**2 &
            + (x(2)-2.0_dp*x(3))**4 + 10.0_dp*(x(1)-x(4))**4
    end function powellq

    pure real(dp) function pricetransistor(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp), parameter :: g(4,5) = reshape([ &
            0.485_dp,0.752_dp,0.869_dp,0.982_dp, &
            0.369_dp,1.254_dp,0.703_dp,1.455_dp, &
            5.2095_dp,10.0677_dp,22.9274_dp,20.2153_dp, &
            23.3037_dp,101.779_dp,111.461_dp,191.267_dp, &
            28.5132_dp,111.8467_dp,134.3884_dp,211.4823_dp], [4,5])
        real(dp) :: alpha, beta, sumsqr
        integer :: k
        sumsqr = 0.0_dp
        do k = 1, 4
            alpha = (1.0_dp-x(1)*x(2))*x(3) * &
                (exp(x(5)*(g(k,1)-0.001_dp*g(k,3)*x(7)-0.001_dp*x(8)*g(k,5)))-1.0_dp) &
                - g(k,5) + g(k,4)*x(2)
            beta = (1.0_dp-x(1)*x(2))*x(4) * &
                (exp(x(6)*(g(k,1)-g(k,2)-0.001_dp*g(k,3)*x(7)+g(k,4)*0.001_dp*x(9)))-1.0_dp) &
                - g(k,5)*x(1) + g(k,4)
            sumsqr = sumsqr + alpha*alpha + beta*beta
        end do
        f = (x(1)*x(3)-x(2)*x(4))**2 + sumsqr
    end function pricetransistor

    pure real(dp) function rastrigin(x) result(f)
        real(dp), intent(in) :: x(:)
        f = sum(x*x - 10.0_dp*cos(2.0_dp*pi_obj*x) + 10.0_dp)
    end function rastrigin

    pure real(dp) function rosenbrock(x) result(f)
        real(dp), intent(in) :: x(:)
        integer :: j
        f = 0.0_dp
        do j = 1, size(x)-1
            f = f + 100.0_dp*(x(j)**2-x(j+1))**2 + (1.0_dp-x(j))**2
        end do
    end function rosenbrock

    pure real(dp) function salomon(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: r
        r = sqrt(sum(x*x))
        f = -cos(2.0_dp*pi_obj*r) + 0.1_dp*r + 1.0_dp
    end function salomon

    pure real(dp) function schaffer1(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: r2
        r2 = x(1)**2+x(2)**2
        f = 0.5_dp + (sin(sqrt(r2))**2-0.5_dp)/(1.0_dp+0.001_dp*r2)**2
    end function schaffer1

    pure real(dp) function schaffer2(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: r2, p1, p2
        r2 = x(1)**2+x(2)**2
        p1 = r2**0.25_dp
        p2 = (50.0_dp*r2)**0.1_dp
        f = p1*(sin(sin(p2))+1.0_dp)
    end function schaffer2

    pure real(dp) function schubert(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: s
        integer :: i, j
        f = 1.0_dp
        do i = 1, size(x)
            s = 0.0_dp
            do j = 1, 5
                s = s + real(j,dp)*cos(real(j+1,dp)*x(i)+real(j,dp))
            end do
            f = f*s
        end do
    end function schubert

    pure real(dp) function schwefel(x) result(f)
        real(dp), intent(in) :: x(:)
        f = -sum(x*sin(sqrt(abs(x))))
    end function schwefel

    pure real(dp) function shekel10(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp), parameter :: a(4,10) = reshape([ &
            4.0_dp,4.0_dp,4.0_dp,4.0_dp, 1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
            8.0_dp,8.0_dp,8.0_dp,8.0_dp, 6.0_dp,6.0_dp,6.0_dp,6.0_dp, &
            3.0_dp,7.0_dp,3.0_dp,7.0_dp, 2.0_dp,9.0_dp,2.0_dp,9.0_dp, &
            5.0_dp,5.0_dp,3.0_dp,3.0_dp, 8.0_dp,1.0_dp,8.0_dp,1.0_dp, &
            6.0_dp,2.0_dp,6.0_dp,2.0_dp, 7.0_dp,3.6_dp,7.0_dp,3.6_dp], [4,10])
        real(dp), parameter :: c(10) = [0.1_dp,0.2_dp,0.2_dp,0.4_dp,0.4_dp, &
            0.6_dp,0.3_dp,0.7_dp,0.5_dp,0.5_dp]
        f = shekel_common(x, a, c)
    end function shekel10

    pure real(dp) function shekel5(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp), parameter :: a(4,5) = reshape([ &
            4.0_dp,4.0_dp,4.0_dp,4.0_dp, 1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
            8.0_dp,8.0_dp,8.0_dp,8.0_dp, 6.0_dp,6.0_dp,6.0_dp,6.0_dp, &
            3.0_dp,7.0_dp,3.0_dp,7.0_dp], [4,5])
        real(dp), parameter :: c(5) = [0.1_dp,0.2_dp,0.2_dp,0.4_dp,0.4_dp]
        f = shekel_common(x, a, c)
    end function shekel5

    pure real(dp) function shekel7(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp), parameter :: a(4,7) = reshape([ &
            4.0_dp,4.0_dp,4.0_dp,4.0_dp, 1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
            8.0_dp,8.0_dp,8.0_dp,8.0_dp, 6.0_dp,6.0_dp,6.0_dp,6.0_dp, &
            3.0_dp,7.0_dp,3.0_dp,7.0_dp, 2.0_dp,9.0_dp,2.0_dp,9.0_dp, &
            5.0_dp,5.0_dp,3.0_dp,3.0_dp], [4,7])
        real(dp), parameter :: c(7) = [0.1_dp,0.2_dp,0.2_dp,0.4_dp,0.4_dp,0.6_dp,0.3_dp]
        f = shekel_common(x, a, c)
    end function shekel7

    pure real(dp) function shekel_common(x, a, c) result(f)
        real(dp), intent(in) :: x(:), a(:,:), c(:)
        integer :: i
        f = 0.0_dp
        do i = 1, size(c)
            f = f - 1.0_dp/(sum((x-a(:,i))**2)+c(i))
        end do
    end function shekel_common

    pure real(dp) function shekelfox5(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp), parameter :: a(10,30) = reshape([ &
            9.681_dp,0.667_dp,4.783_dp,9.095_dp,3.517_dp,9.325_dp,6.544_dp,0.211_dp,5.122_dp,2.020_dp, &
            9.400_dp,2.041_dp,3.788_dp,7.931_dp,2.882_dp,2.672_dp,3.568_dp,1.284_dp,7.033_dp,7.374_dp, &
            8.025_dp,9.152_dp,5.114_dp,7.621_dp,4.564_dp,4.711_dp,2.996_dp,6.126_dp,0.734_dp,4.982_dp, &
            2.196_dp,0.415_dp,5.649_dp,6.979_dp,9.510_dp,9.166_dp,6.304_dp,6.054_dp,9.377_dp,1.426_dp, &
            8.074_dp,8.777_dp,3.467_dp,1.863_dp,6.708_dp,6.349_dp,4.534_dp,0.276_dp,7.633_dp,1.567_dp, &
            7.650_dp,5.658_dp,0.720_dp,2.764_dp,3.278_dp,5.283_dp,7.474_dp,6.274_dp,1.409_dp,8.208_dp, &
            1.256_dp,3.605_dp,8.623_dp,6.905_dp,4.584_dp,8.133_dp,6.071_dp,6.888_dp,4.187_dp,5.448_dp, &
            8.314_dp,2.261_dp,4.224_dp,1.781_dp,4.124_dp,0.932_dp,8.129_dp,8.658_dp,1.208_dp,5.762_dp, &
            0.226_dp,8.858_dp,1.420_dp,0.945_dp,1.622_dp,4.698_dp,6.228_dp,9.096_dp,0.972_dp,7.637_dp, &
            7.305_dp,2.228_dp,1.242_dp,5.928_dp,9.133_dp,1.826_dp,4.060_dp,5.204_dp,8.713_dp,8.247_dp, &
            0.652_dp,7.027_dp,0.508_dp,4.876_dp,8.807_dp,4.632_dp,5.808_dp,6.937_dp,3.291_dp,7.016_dp, &
            2.699_dp,3.516_dp,5.874_dp,4.119_dp,4.461_dp,7.496_dp,8.817_dp,0.690_dp,6.593_dp,9.789_dp, &
            8.327_dp,3.897_dp,2.017_dp,9.570_dp,9.825_dp,1.150_dp,1.395_dp,3.885_dp,6.354_dp,0.109_dp, &
            2.132_dp,7.006_dp,7.136_dp,2.641_dp,1.882_dp,5.943_dp,7.273_dp,7.691_dp,2.880_dp,0.564_dp, &
            4.707_dp,5.579_dp,4.080_dp,0.581_dp,9.698_dp,8.542_dp,8.077_dp,8.515_dp,9.231_dp,4.670_dp, &
            8.304_dp,7.559_dp,8.567_dp,0.322_dp,7.128_dp,8.392_dp,1.472_dp,8.524_dp,2.277_dp,7.826_dp, &
            8.632_dp,4.409_dp,4.832_dp,5.768_dp,7.050_dp,6.715_dp,1.711_dp,4.323_dp,4.405_dp,4.591_dp, &
            4.887_dp,9.112_dp,0.170_dp,8.967_dp,9.693_dp,9.867_dp,7.508_dp,7.770_dp,8.382_dp,6.740_dp, &
            2.440_dp,6.686_dp,4.299_dp,1.007_dp,7.008_dp,1.427_dp,9.398_dp,8.480_dp,9.950_dp,1.675_dp, &
            6.306_dp,8.583_dp,6.084_dp,1.138_dp,4.350_dp,3.134_dp,7.853_dp,6.061_dp,7.457_dp,2.258_dp, &
            0.652_dp,2.343_dp,1.370_dp,0.821_dp,1.310_dp,1.063_dp,0.689_dp,8.819_dp,8.833_dp,9.070_dp, &
            5.558_dp,1.272_dp,5.756_dp,9.857_dp,2.279_dp,2.764_dp,1.284_dp,1.677_dp,1.244_dp,1.234_dp, &
            3.352_dp,7.549_dp,9.817_dp,9.437_dp,8.687_dp,4.167_dp,2.570_dp,6.540_dp,0.228_dp,0.027_dp, &
            8.798_dp,0.880_dp,2.370_dp,0.168_dp,1.701_dp,3.680_dp,1.231_dp,2.390_dp,2.499_dp,0.064_dp, &
            1.460_dp,8.057_dp,1.336_dp,7.217_dp,7.914_dp,3.615_dp,9.981_dp,9.198_dp,5.292_dp,1.224_dp, &
            0.432_dp,8.645_dp,8.774_dp,0.249_dp,8.081_dp,7.461_dp,4.416_dp,0.652_dp,4.002_dp,4.644_dp, &
            0.679_dp,2.800_dp,5.523_dp,3.049_dp,2.968_dp,7.225_dp,6.730_dp,4.199_dp,9.614_dp,9.229_dp, &
            4.263_dp,1.074_dp,7.286_dp,5.599_dp,8.291_dp,5.200_dp,9.214_dp,8.272_dp,4.398_dp,4.506_dp, &
            9.496_dp,4.830_dp,3.150_dp,8.270_dp,5.079_dp,1.231_dp,5.731_dp,9.494_dp,1.883_dp,9.732_dp, &
            4.138_dp,2.562_dp,2.532_dp,9.661_dp,5.611_dp,5.500_dp,6.886_dp,2.341_dp,9.699_dp,6.500_dp], [10,30])
        real(dp), parameter :: c(30) = [ &
            0.806_dp,0.517_dp,0.1_dp,0.908_dp,0.965_dp,0.669_dp,0.524_dp,0.902_dp,0.531_dp,0.876_dp, &
            0.462_dp,0.491_dp,0.463_dp,0.714_dp,0.352_dp,0.869_dp,0.813_dp,0.811_dp,0.828_dp,0.964_dp, &
            0.789_dp,0.360_dp,0.369_dp,0.992_dp,0.332_dp,0.817_dp,0.632_dp,0.883_dp,0.608_dp,0.326_dp]
        integer :: j
        f = 0.0_dp
        do j = 1, 30
            f = f - 1.0_dp/(sum((x-a(1:size(x),j))**2)+c(j))
        end do
    end function shekelfox5

    pure real(dp) function wood(x) result(f)
        real(dp), intent(in) :: x(:)
        f = 100.0_dp*(x(2)-x(1)**2)**2 + (1.0_dp-x(1))**2 &
            + 90.0_dp*(x(4)-x(3)**2)**2 + (1.0_dp-x(3))**2 &
            + 10.1_dp*((x(2)-1.0_dp)**2+(x(4)-1.0_dp)**2) &
            + 19.8_dp*(x(2)-1.0_dp)*(x(4)-1.0_dp)
    end function wood

    pure real(dp) function zeldasine10(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: p1, p2, z
        z = pi_obj/6.0_dp
        p1 = product(sin(x-z))
        p2 = product(sin(5.0_dp*(x-z)))
        f = -(2.5_dp*p1+p2)
    end function zeldasine10

    pure real(dp) function zeldasine20(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: p1, p2, z
        z = pi_obj/6.0_dp
        p1 = product(sin(x-z))
        p2 = product(sin(5.0_dp*(x-z)))
        f = -(2.5_dp*p1+p2)
    end function zeldasine20

end module global_opt_tests

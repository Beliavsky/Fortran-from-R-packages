module zigg
    use iso_fortran_env, only : int32, int64, real64
    implicit none
    private

    integer, parameter, public :: dp = real64
    integer(int64), parameter :: u32_mod = 4294967296_int64
    integer(int64), parameter :: u32_mask = 4294967295_int64
    integer(int64), parameter :: u31 = 2147483648_int64
    real(dp), parameter :: uni_scale = 0.2328306e-9_dp

    type, public :: ziggurat_rng
        private
        real(dp) :: fn(0:127) = 0.0_dp
        real(dp) :: fe(0:255) = 0.0_dp
        integer(int64) :: kn(0:127) = 0_int64
        integer(int64) :: ke(0:255) = 0_int64
        real(dp) :: wn(0:127) = 0.0_dp
        real(dp) :: we(0:255) = 0.0_dp
        integer(int64) :: jcong = 380116160_int64
        integer(int64) :: jsr = 123456789_int64
        integer(int64) :: w = 521288629_int64
        integer(int64) :: z = 362436069_int64
        integer(int64) :: hz = 0_int64
        integer(int64) :: iz = 0_int64
        integer(int64) :: jz = 0_int64
        logical :: initialized = .false.
    contains
        procedure, private :: ensure_initialized
        procedure, private :: initialize_tables
        procedure, private :: znew
        procedure, private :: wnew
        procedure, private :: mwc
        procedure, private :: shr3
        procedure, private :: cong
        procedure, public :: kiss
        procedure, private :: set_seed_i32
        procedure, private :: set_seed_i64
        generic, public :: set_seed => set_seed_i32, set_seed_i64
        procedure, public :: get_seed
        procedure, public :: get_state
        procedure, public :: set_state
        procedure, public :: runi
        procedure, public :: rnorm
        procedure, public :: rexp
        procedure, private :: nfix
        procedure, private :: efix
        procedure, public :: fill_uniform
        procedure, public :: fill_normal
        procedure, public :: fill_exponential
    end type ziggurat_rng

    type(ziggurat_rng), save :: default_rng

    public :: zrnorm, zrexp, zrunif, zsetseed, zgetstate, zsetstate

    interface zsetseed
        module procedure zsetseed_i32
        module procedure zsetseed_i64
    end interface zsetseed

contains

    pure integer(int64) function u32(x) result(y)
        integer(int64), intent(in) :: x
        y = iand(x, u32_mask)
    end function u32

    pure integer(int64) function signed32(x) result(y)
        integer(int64), intent(in) :: x
        integer(int64) :: ux
        ux = u32(x)
        if (ux >= u31) then
            y = ux - u32_mod
        else
            y = ux
        end if
    end function signed32

    subroutine ensure_initialized(self)
        class(ziggurat_rng), intent(inout) :: self
        integer(int64) :: count
        if (self%initialized) return
        call self%initialize_tables()
        call system_clock(count=count)
        call self%set_seed_i64(u32(count))
    end subroutine ensure_initialized

    subroutine initialize_tables(self)
        class(ziggurat_rng), intent(inout) :: self
        real(dp) :: dn, de, tn, te, q
        real(dp), parameter :: m1 = 2147483648.0_dp
        real(dp), parameter :: m2 = 4294967296.0_dp
        real(dp), parameter :: vn = 9.91256303526217e-3_dp
        real(dp), parameter :: ve = 3.949659822581572e-3_dp
        integer :: i

        dn = 3.442619855899_dp
        de = 7.697117470131487_dp
        tn = dn
        te = de

        q = vn / exp(-0.5_dp * dn * dn)
        self%kn(0) = int((dn / q) * m1, int64)
        self%kn(1) = 0_int64
        self%wn(0) = q / m1
        self%wn(127) = dn / m1
        self%fn(0) = 1.0_dp
        self%fn(127) = exp(-0.5_dp * dn * dn)
        do i = 126, 1, -1
            dn = sqrt(-2.0_dp * log(vn / dn + exp(-0.5_dp * dn * dn)))
            self%kn(i + 1) = int((dn / tn) * m1, int64)
            tn = dn
            self%fn(i) = exp(-0.5_dp * dn * dn)
            self%wn(i) = dn / m1
        end do

        q = ve / exp(-de)
        ! Preserve the upstream C++ cast placement exactly:
        ! ke[0] = (uint32_t) (de/q) * m2;  Since de/q < 1, this is zero.
        self%ke(0) = int(real(int(de / q, int64), dp) * m2, int64)
        self%ke(1) = 0_int64
        self%we(0) = q / m2
        self%we(255) = de / m2
        self%fe(0) = 1.0_dp
        self%fe(255) = exp(-de)
        do i = 254, 1, -1
            de = -log(ve / de + exp(-de))
            self%ke(i + 1) = int((de / te) * m2, int64)
            te = de
            self%fe(i) = exp(-de)
            self%we(i) = de / m2
        end do

        self%initialized = .true.
    end subroutine initialize_tables

    integer(int64) function znew(self) result(v)
        class(ziggurat_rng), intent(inout) :: self
        self%z = u32(36969_int64 * iand(self%z, 65535_int64) + shiftr(self%z, 16))
        v = self%z
    end function znew

    integer(int64) function wnew(self) result(v)
        class(ziggurat_rng), intent(inout) :: self
        self%w = u32(18000_int64 * iand(self%w, 65535_int64) + shiftr(self%w, 16))
        v = self%w
    end function wnew

    integer(int64) function mwc(self) result(v)
        class(ziggurat_rng), intent(inout) :: self
        integer(int64) :: zv, wv
        zv = self%znew()
        wv = self%wnew()
        v = u32(u32(shiftl(zv, 16)) + wv)
    end function mwc

    integer(int64) function shr3(self) result(v)
        class(ziggurat_rng), intent(inout) :: self
        integer(int64) :: old
        old = self%jsr
        self%jsr = u32(ieor(self%jsr, u32(shiftl(self%jsr, 13))))
        self%jsr = u32(ieor(self%jsr, shiftr(self%jsr, 17)))
        self%jsr = u32(ieor(self%jsr, u32(shiftl(self%jsr, 5))))
        v = u32(old + self%jsr)
    end function shr3

    integer(int64) function cong(self) result(v)
        class(ziggurat_rng), intent(inout) :: self
        self%jcong = u32(69069_int64 * self%jcong + 1234567_int64)
        v = self%jcong
    end function cong

    integer(int64) function kiss(self) result(v)
        class(ziggurat_rng), intent(inout) :: self
        integer(int64) :: a, b, c
        call self%ensure_initialized()
        a = self%mwc()
        b = self%cong()
        c = self%shr3()
        v = u32(ieor(a, b) + c)
    end function kiss

    subroutine set_seed_i32(self, seed)
        class(ziggurat_rng), intent(inout) :: self
        integer(int32), intent(in) :: seed
        call self%set_seed_i64(int(seed, int64))
    end subroutine set_seed_i32

    subroutine set_seed_i64(self, seed)
        class(ziggurat_rng), intent(inout) :: self
        integer(int64), intent(in) :: seed
        integer(int64) :: s
        if (.not. self%initialized) call self%initialize_tables()
        s = u32(seed)
        self%jsr = 123456789_int64
        if (self%jsr /= s) self%jsr = u32(ieor(self%jsr, s))
        self%z = 362436069_int64
        self%w = 521288629_int64
        self%jcong = 380116160_int64
        self%hz = 0_int64
        self%iz = 0_int64
        self%jz = 0_int64
    end subroutine set_seed_i64

    integer(int64) function get_seed(self) result(seed)
        class(ziggurat_rng), intent(inout) :: self
        call self%ensure_initialized()
        seed = self%jsr
    end function get_seed

    subroutine get_state(self, state)
        class(ziggurat_rng), intent(inout) :: self
        integer(int64), intent(out) :: state(4)
        call self%ensure_initialized()
        state = [self%jsr, self%z, self%w, self%jcong]
    end subroutine get_state

    subroutine set_state(self, state)
        class(ziggurat_rng), intent(inout) :: self
        integer(int64), intent(in) :: state(4)
        if (.not. self%initialized) call self%initialize_tables()
        self%jsr = u32(state(1))
        self%z = u32(state(2))
        self%w = u32(state(3))
        self%jcong = u32(state(4))
    end subroutine set_state

    real(dp) function runi(self) result(x)
        class(ziggurat_rng), intent(inout) :: self
        x = 0.5_dp + real(signed32(self%kiss()), dp) * uni_scale
    end function runi

    real(dp) function rnorm(self) result(x)
        class(ziggurat_rng), intent(inout) :: self
        integer(int64) :: k
        call self%ensure_initialized()
        k = self%kiss()
        self%hz = signed32(k)
        self%iz = iand(k, 127_int64)
        if (abs(self%hz) < self%kn(self%iz)) then
            x = real(self%hz, dp) * self%wn(self%iz)
        else
            x = self%nfix()
        end if
    end function rnorm

    real(dp) function rexp(self) result(x)
        class(ziggurat_rng), intent(inout) :: self
        call self%ensure_initialized()
        self%jz = self%kiss()
        self%iz = iand(self%jz, 255_int64)
        if (self%jz < self%ke(self%iz)) then
            x = real(self%jz, dp) * self%we(self%iz)
        else
            x = self%efix()
        end if
    end function rexp

    real(dp) function nfix(self) result(value)
        class(ziggurat_rng), intent(inout) :: self
        real(dp), parameter :: r = 3.442620_dp
        real(dp) :: x, y
        integer(int64) :: k

        do
            x = real(self%hz, dp) * self%wn(self%iz)
            if (self%iz == 0_int64) then
                do
                    x = -log(self%runi()) * 0.2904764_dp
                    y = -log(self%runi())
                    if (y + y >= x * x) exit
                end do
                if (self%hz > 0_int64) then
                    value = r + x
                else
                    value = -r - x
                end if
                return
            end if

            if (self%fn(self%iz) + self%runi() * (self%fn(self%iz - 1) - self%fn(self%iz)) < &
                    exp(-0.5_dp * x * x)) then
                value = x
                return
            end if

            k = self%shr3()
            self%hz = signed32(k)
            self%iz = iand(k, 127_int64)
            if (abs(self%hz) < self%kn(self%iz)) then
                value = real(self%hz, dp) * self%wn(self%iz)
                return
            end if
        end do
    end function nfix

    real(dp) function efix(self) result(value)
        class(ziggurat_rng), intent(inout) :: self
        real(dp) :: x
        do
            if (self%iz == 0_int64) then
                value = 7.69711_dp - log(self%runi())
                return
            end if
            x = real(self%jz, dp) * self%we(self%iz)
            if (self%fe(self%iz) + self%runi() * (self%fe(self%iz - 1) - self%fe(self%iz)) < exp(-x)) then
                value = x
                return
            end if
            self%jz = self%shr3()
            self%iz = iand(self%jz, 255_int64)
            if (self%jz < self%ke(self%iz)) then
                value = real(self%jz, dp) * self%we(self%iz)
                return
            end if
        end do
    end function efix

    subroutine fill_uniform(self, x)
        class(ziggurat_rng), intent(inout) :: self
        real(dp), intent(out) :: x(:)
        integer :: i
        do i = 1, size(x)
            x(i) = self%runi()
        end do
    end subroutine fill_uniform

    subroutine fill_normal(self, x)
        class(ziggurat_rng), intent(inout) :: self
        real(dp), intent(out) :: x(:)
        integer :: i
        do i = 1, size(x)
            x(i) = self%rnorm()
        end do
    end subroutine fill_normal

    subroutine fill_exponential(self, x)
        class(ziggurat_rng), intent(inout) :: self
        real(dp), intent(out) :: x(:)
        integer :: i
        do i = 1, size(x)
            x(i) = self%rexp()
        end do
    end subroutine fill_exponential

    function zrnorm(n) result(x)
        integer, intent(in) :: n
        real(dp), allocatable :: x(:)
        allocate(x(max(n, 0)))
        call default_rng%fill_normal(x)
    end function zrnorm

    function zrexp(n) result(x)
        integer, intent(in) :: n
        real(dp), allocatable :: x(:)
        allocate(x(max(n, 0)))
        call default_rng%fill_exponential(x)
    end function zrexp

    function zrunif(n) result(x)
        integer, intent(in) :: n
        real(dp), allocatable :: x(:)
        allocate(x(max(n, 0)))
        call default_rng%fill_uniform(x)
    end function zrunif

    subroutine zsetseed_i32(seed)
        integer(int32), intent(in) :: seed
        call default_rng%set_seed(seed)
    end subroutine zsetseed_i32

    subroutine zsetseed_i64(seed)
        integer(int64), intent(in) :: seed
        call default_rng%set_seed(seed)
    end subroutine zsetseed_i64

    subroutine zgetstate(state)
        integer(int64), intent(out) :: state(4)
        call default_rng%get_state(state)
    end subroutine zgetstate

    subroutine zsetstate(state)
        integer(int64), intent(in) :: state(4)
        call default_rng%set_state(state)
    end subroutine zsetstate

end module zigg

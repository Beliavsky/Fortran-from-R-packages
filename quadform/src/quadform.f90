module quadform
    use quadform_kinds, only : dp
    use quadform_linalg, only : solve_linear
    implicit none
    private

    public :: dp
    public :: ht, cprod, tcprod
    public :: quad_form_chol, quad_form, quad_form_inv
    public :: quad3_form_ab, quad3_form_bc, quad3_form, quad3_form_inv
    public :: quad3_tform_ab, quad3_tform_bc, quad3_tform
    public :: quad_tform, quad_tform_inv
    public :: quad_diag, quad_tdiag, quad3_diag, quad3_tdiag
    public :: quad_trace, quad_ttrace
    public :: cp, tcp, qf, qfi, q3, q3i, q3t, qt, qti, qd, qtd, q3d, q3td, qtr, qttr

    interface ht
        module procedure ht_real
        module procedure ht_complex
    end interface
    interface cprod
        module procedure cprod_real
        module procedure cprod_complex
    end interface
    interface tcprod
        module procedure tcprod_real
        module procedure tcprod_complex
    end interface
    interface quad_form_chol
        module procedure quad_form_chol_real
        module procedure quad_form_chol_complex
    end interface
    interface quad_form
        module procedure quad_form_real
        module procedure quad_form_complex
    end interface
    interface quad_form_inv
        module procedure quad_form_inv_real
        module procedure quad_form_inv_complex
    end interface
    interface quad3_form_ab
        module procedure quad3_form_ab_real
        module procedure quad3_form_ab_complex
    end interface
    interface quad3_form_bc
        module procedure quad3_form_bc_real
        module procedure quad3_form_bc_complex
    end interface
    interface quad3_form
        module procedure quad3_form_real
        module procedure quad3_form_complex
    end interface
    interface quad3_form_inv
        module procedure quad3_form_inv_real
        module procedure quad3_form_inv_complex
    end interface
    interface quad3_tform_ab
        module procedure quad3_tform_ab_real
        module procedure quad3_tform_ab_complex
    end interface
    interface quad3_tform_bc
        module procedure quad3_tform_bc_real
        module procedure quad3_tform_bc_complex
    end interface
    interface quad3_tform
        module procedure quad3_tform_real
        module procedure quad3_tform_complex
    end interface
    interface quad_tform
        module procedure quad_tform_real
        module procedure quad_tform_complex
    end interface
    interface quad_tform_inv
        module procedure quad_tform_inv_real
        module procedure quad_tform_inv_complex
    end interface
    interface quad_diag
        module procedure quad_diag_real
        module procedure quad_diag_complex
    end interface
    interface quad_tdiag
        module procedure quad_tdiag_real
        module procedure quad_tdiag_complex
    end interface
    interface quad3_diag
        module procedure quad3_diag_real
        module procedure quad3_diag_complex
    end interface
    interface quad3_tdiag
        module procedure quad3_tdiag_real
        module procedure quad3_tdiag_complex
    end interface
    interface quad_trace
        module procedure quad_trace_real
        module procedure quad_trace_complex
    end interface
    interface quad_ttrace
        module procedure quad_ttrace_real
        module procedure quad_ttrace_complex
    end interface

    interface cp
        module procedure cprod_real
        module procedure cprod_complex
    end interface
    interface tcp
        module procedure tcprod_real
        module procedure tcprod_complex
    end interface
    interface qf
        module procedure quad_form_real
        module procedure quad_form_complex
    end interface
    interface qfi
        module procedure quad_form_inv_real
        module procedure quad_form_inv_complex
    end interface
    interface q3
        module procedure quad3_form_real
        module procedure quad3_form_complex
    end interface
    interface q3i
        module procedure quad3_form_inv_real
        module procedure quad3_form_inv_complex
    end interface
    interface q3t
        module procedure quad3_tform_real
        module procedure quad3_tform_complex
    end interface
    interface qt
        module procedure quad_tform_real
        module procedure quad_tform_complex
    end interface
    interface qti
        module procedure quad_tform_inv_real
        module procedure quad_tform_inv_complex
    end interface
    interface qd
        module procedure quad_diag_real
        module procedure quad_diag_complex
    end interface
    interface qtd
        module procedure quad_tdiag_real
        module procedure quad_tdiag_complex
    end interface
    interface q3d
        module procedure quad3_diag_real
        module procedure quad3_diag_complex
    end interface
    interface q3td
        module procedure quad3_tdiag_real
        module procedure quad3_tdiag_complex
    end interface
    interface qtr
        module procedure quad_trace_real
        module procedure quad_trace_complex
    end interface
    interface qttr
        module procedure quad_ttrace_real
        module procedure quad_ttrace_complex
    end interface

contains

    pure function ht_real(x) result(y)
        real(dp), intent(in) :: x(:, :)
        real(dp) :: y(size(x, 2), size(x, 1))
        y = transpose(x)
    end function ht_real

    pure function ht_complex(x) result(y)
        complex(dp), intent(in) :: x(:, :)
        complex(dp) :: y(size(x, 2), size(x, 1))
        y = transpose(conjg(x))
    end function ht_complex

    pure function cprod_real(x, y) result(z)
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(in), optional :: y(:, :)
        real(dp), allocatable :: z(:, :)
        if (present(y)) then
            z = matmul(transpose(x), y)
        else
            z = matmul(transpose(x), x)
        end if
    end function cprod_real

    pure function cprod_complex(x, y) result(z)
        complex(dp), intent(in) :: x(:, :)
        complex(dp), intent(in), optional :: y(:, :)
        complex(dp), allocatable :: z(:, :)
        if (present(y)) then
            z = matmul(transpose(conjg(x)), y)
        else
            z = matmul(transpose(conjg(x)), x)
        end if
    end function cprod_complex

    pure function tcprod_real(x, y) result(z)
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(in), optional :: y(:, :)
        real(dp), allocatable :: z(:, :)
        if (present(y)) then
            z = matmul(x, transpose(y))
        else
            z = matmul(x, transpose(x))
        end if
    end function tcprod_real

    pure function tcprod_complex(x, y) result(z)
        complex(dp), intent(in) :: x(:, :)
        complex(dp), intent(in), optional :: y(:, :)
        complex(dp), allocatable :: z(:, :)
        if (present(y)) then
            z = matmul(x, transpose(conjg(y)))
        else
            z = matmul(x, transpose(conjg(x)))
        end if
    end function tcprod_complex

    pure function quad_form_chol_real(chol, x) result(z)
        real(dp), intent(in) :: chol(:, :), x(:, :)
        real(dp), allocatable :: z(:, :), jj(:, :)
        jj = matmul(transpose(chol), x)
        z = matmul(transpose(jj), jj)
    end function quad_form_chol_real

    pure function quad_form_chol_complex(chol, x) result(z)
        complex(dp), intent(in) :: chol(:, :), x(:, :)
        complex(dp), allocatable :: z(:, :), jj(:, :)
        jj = matmul(transpose(conjg(chol)), x)
        z = matmul(transpose(conjg(jj)), jj)
    end function quad_form_chol_complex

    pure function quad_form_real(m, x) result(z)
        real(dp), intent(in) :: m(:, :), x(:, :)
        real(dp), allocatable :: z(:, :)
        z = matmul(transpose(x), matmul(m, x))
    end function quad_form_real

    pure function quad_form_complex(m, x) result(z)
        complex(dp), intent(in) :: m(:, :), x(:, :)
        complex(dp), allocatable :: z(:, :)
        z = matmul(transpose(conjg(x)), matmul(m, x))
    end function quad_form_complex

    function quad_form_inv_real(m, x, info) result(z)
        real(dp), intent(in) :: m(:, :), x(:, :)
        integer, intent(out), optional :: info
        real(dp), allocatable :: z(:, :), sol(:, :)
        integer :: istat
        sol = solve_linear(m, x, istat)
        if (istat /= 0) then
            allocate(z(0, 0))
        else
            z = matmul(transpose(x), sol)
        end if
        if (present(info)) info = istat
    end function quad_form_inv_real

    function quad_form_inv_complex(m, x, info) result(z)
        complex(dp), intent(in) :: m(:, :), x(:, :)
        integer, intent(out), optional :: info
        complex(dp), allocatable :: z(:, :), sol(:, :)
        integer :: istat
        sol = solve_linear(m, x, istat)
        if (istat /= 0) then
            allocate(z(0, 0))
        else
            z = matmul(transpose(conjg(x)), sol)
        end if
        if (present(info)) info = istat
    end function quad_form_inv_complex

    pure function quad3_form_ab_real(m, left, right) result(z)
        real(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        real(dp), allocatable :: z(:, :)
        z = matmul(transpose(matmul(transpose(m), left)), right)
    end function quad3_form_ab_real

    pure function quad3_form_ab_complex(m, left, right) result(z)
        complex(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        complex(dp), allocatable :: z(:, :)
        z = matmul(transpose(matmul(transpose(m), conjg(left))), right)
    end function quad3_form_ab_complex

    pure function quad3_form_bc_real(m, left, right) result(z)
        real(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        real(dp), allocatable :: z(:, :)
        z = matmul(transpose(left), matmul(m, right))
    end function quad3_form_bc_real

    pure function quad3_form_bc_complex(m, left, right) result(z)
        complex(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        complex(dp), allocatable :: z(:, :)
        z = matmul(transpose(conjg(left)), matmul(m, right))
    end function quad3_form_bc_complex

    pure function quad3_form_real(m, left, right) result(z)
        real(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        real(dp), allocatable :: z(:, :)
        if (size(left, 2) < size(right, 2)) then
            z = quad3_form_ab_real(m, left, right)
        else
            z = quad3_form_bc_real(m, left, right)
        end if
    end function quad3_form_real

    pure function quad3_form_complex(m, left, right) result(z)
        complex(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        complex(dp), allocatable :: z(:, :)
        if (size(left, 2) < size(right, 2)) then
            z = quad3_form_ab_complex(m, left, right)
        else
            z = quad3_form_bc_complex(m, left, right)
        end if
    end function quad3_form_complex

    function quad3_form_inv_real(m, left, right, info) result(z)
        real(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        integer, intent(out), optional :: info
        real(dp), allocatable :: z(:, :), sol(:, :)
        integer :: istat
        sol = solve_linear(m, right, istat)
        if (istat /= 0) then
            allocate(z(0, 0))
        else
            z = matmul(transpose(left), sol)
        end if
        if (present(info)) info = istat
    end function quad3_form_inv_real

    function quad3_form_inv_complex(m, left, right, info) result(z)
        complex(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        integer, intent(out), optional :: info
        complex(dp), allocatable :: z(:, :), sol(:, :)
        integer :: istat
        sol = solve_linear(m, right, istat)
        if (istat /= 0) then
            allocate(z(0, 0))
        else
            z = matmul(transpose(conjg(left)), sol)
        end if
        if (present(info)) info = istat
    end function quad3_form_inv_complex

    pure function quad3_tform_ab_real(m, left, right) result(z)
        real(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        real(dp), allocatable :: z(:, :)
        z = matmul(matmul(left, m), transpose(right))
    end function quad3_tform_ab_real

    pure function quad3_tform_ab_complex(m, left, right) result(z)
        complex(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        complex(dp), allocatable :: z(:, :)
        z = matmul(matmul(left, m), transpose(conjg(right)))
    end function quad3_tform_ab_complex

    pure function quad3_tform_bc_real(m, left, right) result(z)
        real(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        real(dp), allocatable :: z(:, :)
        z = matmul(left, matmul(m, transpose(right)))
    end function quad3_tform_bc_real

    pure function quad3_tform_bc_complex(m, left, right) result(z)
        complex(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        complex(dp), allocatable :: z(:, :)
        z = matmul(left, matmul(m, transpose(conjg(right))))
    end function quad3_tform_bc_complex

    pure function quad3_tform_real(m, left, right) result(z)
        real(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        real(dp), allocatable :: z(:, :)
        if (size(left, 1) < size(right, 1)) then
            z = quad3_tform_ab_real(m, left, right)
        else
            z = quad3_tform_bc_real(m, left, right)
        end if
    end function quad3_tform_real

    pure function quad3_tform_complex(m, left, right) result(z)
        complex(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        complex(dp), allocatable :: z(:, :)
        if (size(left, 1) < size(right, 1)) then
            z = quad3_tform_ab_complex(m, left, right)
        else
            z = quad3_tform_bc_complex(m, left, right)
        end if
    end function quad3_tform_complex

    pure function quad_tform_real(m, x) result(z)
        real(dp), intent(in) :: m(:, :), x(:, :)
        real(dp), allocatable :: z(:, :)
        z = matmul(x, matmul(m, transpose(x)))
    end function quad_tform_real

    pure function quad_tform_complex(m, x) result(z)
        complex(dp), intent(in) :: m(:, :), x(:, :)
        complex(dp), allocatable :: z(:, :)
        z = matmul(x, matmul(m, transpose(conjg(x))))
    end function quad_tform_complex

    function quad_tform_inv_real(m, x, info) result(z)
        real(dp), intent(in) :: m(:, :), x(:, :)
        integer, intent(out), optional :: info
        real(dp), allocatable :: z(:, :), hx(:, :)
        integer :: istat
        hx = transpose(x)
        z = quad_form_inv_real(m, hx, istat)
        if (present(info)) info = istat
    end function quad_tform_inv_real

    function quad_tform_inv_complex(m, x, info) result(z)
        complex(dp), intent(in) :: m(:, :), x(:, :)
        integer, intent(out), optional :: info
        complex(dp), allocatable :: z(:, :), hx(:, :)
        integer :: istat
        hx = transpose(conjg(x))
        z = quad_form_inv_complex(m, hx, istat)
        if (present(info)) info = istat
    end function quad_tform_inv_complex

    pure function quad_diag_real(m, x) result(d)
        real(dp), intent(in) :: m(:, :), x(:, :)
        real(dp) :: d(size(x, 2))
        real(dp) :: mx(size(m, 1), size(x, 2))
        integer :: j
        mx = matmul(m, x)
        do j = 1, size(x, 2)
            d(j) = sum(x(:, j) * mx(:, j))
        end do
    end function quad_diag_real

    pure function quad_diag_complex(m, x) result(d)
        complex(dp), intent(in) :: m(:, :), x(:, :)
        complex(dp) :: d(size(x, 2))
        complex(dp) :: mx(size(m, 1), size(x, 2))
        integer :: j
        mx = matmul(m, x)
        do j = 1, size(x, 2)
            d(j) = sum(conjg(x(:, j)) * mx(:, j))
        end do
    end function quad_diag_complex

    pure function quad_tdiag_real(m, x) result(d)
        real(dp), intent(in) :: m(:, :), x(:, :)
        real(dp) :: d(size(x, 1))
        real(dp) :: xm(size(x, 1), size(m, 2))
        integer :: i
        xm = matmul(x, m)
        do i = 1, size(x, 1)
            d(i) = sum(xm(i, :) * x(i, :))
        end do
    end function quad_tdiag_real

    pure function quad_tdiag_complex(m, x) result(d)
        complex(dp), intent(in) :: m(:, :), x(:, :)
        complex(dp) :: d(size(x, 1))
        complex(dp) :: xm(size(x, 1), size(m, 2))
        integer :: i
        xm = matmul(x, m)
        do i = 1, size(x, 1)
            d(i) = sum(xm(i, :) * conjg(x(i, :)))
        end do
    end function quad_tdiag_complex

    pure function quad3_diag_real(m, left, right) result(d)
        real(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        real(dp) :: d(min(size(left, 2), size(right, 2)))
        real(dp) :: mr(size(m, 1), size(right, 2))
        integer :: j
        mr = matmul(m, right)
        do j = 1, size(d)
            d(j) = sum(left(:, j) * mr(:, j))
        end do
    end function quad3_diag_real

    pure function quad3_diag_complex(m, left, right) result(d)
        complex(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        complex(dp) :: d(min(size(left, 2), size(right, 2)))
        complex(dp) :: mr(size(m, 1), size(right, 2))
        integer :: j
        mr = matmul(m, right)
        do j = 1, size(d)
            d(j) = sum(conjg(left(:, j)) * mr(:, j))
        end do
    end function quad3_diag_complex

    pure function quad3_tdiag_real(m, left, right) result(d)
        real(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        real(dp) :: d(min(size(left, 1), size(right, 1)))
        real(dp) :: lm(size(left, 1), size(m, 2))
        integer :: i
        lm = matmul(left, m)
        do i = 1, size(d)
            d(i) = sum(lm(i, :) * right(i, :))
        end do
    end function quad3_tdiag_real

    pure function quad3_tdiag_complex(m, left, right) result(d)
        complex(dp), intent(in) :: m(:, :), left(:, :), right(:, :)
        complex(dp) :: d(min(size(left, 1), size(right, 1)))
        complex(dp) :: lm(size(left, 1), size(m, 2))
        integer :: i
        lm = matmul(left, m)
        do i = 1, size(d)
            d(i) = sum(lm(i, :) * conjg(right(i, :)))
        end do
    end function quad3_tdiag_complex

    pure function quad_trace_real(m, x) result(v)
        real(dp), intent(in) :: m(:, :), x(:, :)
        real(dp) :: v
        v = sum(quad_diag_real(m, x))
    end function quad_trace_real

    pure function quad_trace_complex(m, x) result(v)
        complex(dp), intent(in) :: m(:, :), x(:, :)
        complex(dp) :: v
        v = sum(quad_diag_complex(m, x))
    end function quad_trace_complex

    pure function quad_ttrace_real(m, x) result(v)
        real(dp), intent(in) :: m(:, :), x(:, :)
        real(dp) :: v
        v = sum(quad_tdiag_real(m, x))
    end function quad_ttrace_real

    pure function quad_ttrace_complex(m, x) result(v)
        complex(dp), intent(in) :: m(:, :), x(:, :)
        complex(dp) :: v
        v = sum(quad_tdiag_complex(m, x))
    end function quad_ttrace_complex

end module quadform

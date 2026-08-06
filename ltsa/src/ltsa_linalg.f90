! SPDX-License-Identifier: GPL-2.0-or-later
module ltsa_linalg
    use ltsa_kinds, only : dp
    use ltsa_status, only : ltsa_error, ltsa_success, ltsa_invalid_input, ltsa_not_positive_definite, &
                            ltsa_singular, set_error
    implicit none
    private

    public :: toeplitz_matrix, is_toeplitz, cholesky_factor, solve_spd, inverse_spd, solve_linear

contains

    function toeplitz_matrix(r) result(a)
        real(dp), intent(in) :: r(:)
        real(dp), allocatable :: a(:,:)
        integer :: i, j, n
        n = size(r)
        allocate(a(n,n))
        do j = 1, n
            do i = 1, n
                a(i,j) = r(abs(i-j)+1)
            end do
        end do
    end function toeplitz_matrix

    logical function is_toeplitz(a, tolerance) result(answer)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(in), optional :: tolerance
        real(dp) :: tol, scale
        integer :: i, j, n
        answer = .false.
        if (size(a,1) /= size(a,2)) return
        n = size(a,1)
        if (n == 0) return
        tol = epsilon(1.0_dp)
        if (present(tolerance)) tol = max(0.0_dp, tolerance)
        scale = max(1.0_dp, maxval(abs(a)))
        do j = 1, n
            do i = 1, n
                if (abs(a(i,j)-a(1,abs(i-j)+1)) > tol*scale) return
            end do
        end do
        answer = .true.
    end function is_toeplitz

    subroutine cholesky_factor(a, l, error)
        real(dp), intent(in) :: a(:,:)
        real(dp), allocatable, intent(out) :: l(:,:)
        type(ltsa_error), intent(out) :: error
        integer :: i, j, k, n
        real(dp) :: s, tol
        error%code = ltsa_success
        error%message = ''
        if (size(a,1) /= size(a,2) .or. size(a,1) < 1) then
            allocate(l(0,0))
            call set_error(error, ltsa_invalid_input, 'matrix must be nonempty and square')
            return
        end if
        n = size(a,1)
        allocate(l(n,n), source=0.0_dp)
        tol = epsilon(1.0_dp)*max(1.0_dp, maxval(abs(a)))*real(n,dp)
        do i = 1, n
            do j = 1, i
                s = a(i,j)
                do k = 1, j-1
                    s = s-l(i,k)*l(j,k)
                end do
                if (i == j) then
                    if (s <= tol) then
                        call set_error(error, ltsa_not_positive_definite, 'matrix is not positive definite')
                        return
                    end if
                    l(i,j) = sqrt(s)
                else
                    l(i,j) = s/l(j,j)
                end if
            end do
        end do
    end subroutine cholesky_factor

    subroutine solve_spd(a, b, x, error)
        real(dp), intent(in) :: a(:,:), b(:)
        real(dp), allocatable, intent(out) :: x(:)
        type(ltsa_error), intent(out) :: error
        real(dp), allocatable :: l(:,:), y(:)
        integer :: i, j, n
        if (size(a,1) /= size(a,2) .or. size(b) /= size(a,1)) then
            allocate(x(0))
            call set_error(error, ltsa_invalid_input, 'incompatible solve dimensions')
            return
        end if
        call cholesky_factor(a, l, error)
        if (.not. error%ok()) then
            allocate(x(0))
            return
        end if
        n = size(b)
        allocate(y(n), x(n))
        do i = 1, n
            y(i) = b(i)
            do j = 1, i-1
                y(i) = y(i)-l(i,j)*y(j)
            end do
            y(i) = y(i)/l(i,i)
        end do
        do i = n, 1, -1
            x(i) = y(i)
            do j = i+1, n
                x(i) = x(i)-l(j,i)*x(j)
            end do
            x(i) = x(i)/l(i,i)
        end do
    end subroutine solve_spd

    subroutine inverse_spd(a, ainv, logdet, error)
        real(dp), intent(in) :: a(:,:)
        real(dp), allocatable, intent(out) :: ainv(:,:)
        real(dp), intent(out) :: logdet
        type(ltsa_error), intent(out) :: error
        real(dp), allocatable :: l(:,:), y(:), x(:)
        integer :: i, j, k, n
        call cholesky_factor(a, l, error)
        if (.not. error%ok()) then
            allocate(ainv(0,0))
            logdet = 0.0_dp
            return
        end if
        n = size(a,1)
        allocate(ainv(n,n), y(n), x(n))
        logdet = 2.0_dp*sum(log([(l(i,i), i=1,n)]))
        do k = 1, n
            y = 0.0_dp
            do i = 1, n
                if (i == k) y(i) = 1.0_dp
                do j = 1, i-1
                    y(i) = y(i)-l(i,j)*y(j)
                end do
                y(i) = y(i)/l(i,i)
            end do
            do i = n, 1, -1
                x(i) = y(i)
                do j = i+1, n
                    x(i) = x(i)-l(j,i)*x(j)
                end do
                x(i) = x(i)/l(i,i)
            end do
            ainv(:,k) = x
        end do
        ainv = 0.5_dp*(ainv+transpose(ainv))
    end subroutine inverse_spd

    subroutine solve_linear(a, b, x, error)
        real(dp), intent(in) :: a(:,:), b(:)
        real(dp), allocatable, intent(out) :: x(:)
        type(ltsa_error), intent(out) :: error
        real(dp), allocatable :: aug(:,:), rowtmp(:)
        real(dp) :: pivot, factor, tol
        integer :: i, k, p, n
        error%code = ltsa_success
        error%message = ''
        if (size(a,1) /= size(a,2) .or. size(b) /= size(a,1) .or. size(b) < 1) then
            allocate(x(0))
            call set_error(error, ltsa_invalid_input, 'invalid linear-system dimensions')
            return
        end if
        n = size(b)
        allocate(aug(n,n+1), rowtmp(n+1), x(n))
        aug(:,1:n) = a
        aug(:,n+1) = b
        tol = epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))*real(n,dp)
        do k = 1, n
            p = k-1+maxloc(abs(aug(k:n,k)), dim=1)
            if (abs(aug(p,k)) <= tol) then
                call set_error(error, ltsa_singular, 'linear system is singular')
                x = 0.0_dp
                return
            end if
            if (p /= k) then
                rowtmp = aug(k,:)
                aug(k,:) = aug(p,:)
                aug(p,:) = rowtmp
            end if
            pivot = aug(k,k)
            aug(k,k:n+1) = aug(k,k:n+1)/pivot
            do i = 1, n
                if (i == k) cycle
                factor = aug(i,k)
                if (factor /= 0.0_dp) aug(i,k:n+1) = aug(i,k:n+1)-factor*aug(k,k:n+1)
            end do
        end do
        do i = 1, n
            x(i) = aug(i,n+1)
        end do
    end subroutine solve_linear

end module ltsa_linalg

! SPDX-License-Identifier: GPL-2.0-or-later
module msm_linalg
    use msm_kinds, only : dp
    implicit none
    private
    public :: eye, solve_linear, inverse_matrix, expm, expm_frechet

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*)
            integer, intent(out) :: info
            real(dp), intent(inout) :: a(lda,*), b(ldb,*)
        end subroutine dgesv
    end interface
contains
    function eye(n) result(a)
        integer, intent(in) :: n
        real(dp) :: a(n,n)
        integer :: i
        a = 0.0_dp
        do i = 1, n
            a(i,i) = 1.0_dp
        end do
    end function eye

    function solve_linear(a, b) result(x)
        real(dp), intent(in) :: a(:,:), b(:,:)
        real(dp), allocatable :: x(:,:), ac(:,:)
        integer, allocatable :: ipiv(:)
        integer :: n, nrhs, info
        n = size(a,1)
        if (size(a,2) /= n .or. size(b,1) /= n) error stop "solve_linear: nonconformable arrays"
        nrhs = size(b,2)
        ac = a
        x = b
        allocate(ipiv(n))
        call dgesv(n, nrhs, ac, n, ipiv, x, n, info)
        if (info /= 0) error stop "solve_linear: dgesv failed"
    end function solve_linear

    function inverse_matrix(a) result(ainv)
        real(dp), intent(in) :: a(:,:)
        real(dp), allocatable :: ainv(:,:)
        ainv = solve_linear(a, eye(size(a,1)))
    end function inverse_matrix

    pure function norm1(a) result(v)
        real(dp), intent(in) :: a(:,:)
        real(dp) :: v
        integer :: j
        v = 0.0_dp
        do j = 1, size(a,2)
            v = max(v, sum(abs(a(:,j))))
        end do
    end function norm1

    function expm(a) result(x)
        ! Higham-style scaling and squaring with diagonal Pade approximants.
        real(dp), intent(in) :: a(:,:)
        real(dp), allocatable :: x(:,:), aa(:,:), a2(:,:), a4(:,:), a6(:,:), u(:,:), v(:,:), id(:,:)
        real(dp) :: na
        real(dp), parameter :: c(14) = [ &
            64764752532480000.0_dp, 32382376266240000.0_dp, 7771770303897600.0_dp, &
            1187353796428800.0_dp, 129060195264000.0_dp, 10559470521600.0_dp, &
            670442572800.0_dp, 33522128640.0_dp, 1323241920.0_dp, 40840800.0_dp, &
            960960.0_dp, 16380.0_dp, 182.0_dp, 1.0_dp ]
        real(dp), parameter :: theta13 = 5.4_dp
        integer :: n, s, k
        n = size(a,1)
        if (size(a,2) /= n) error stop "expm: matrix must be square"
        if (n == 0) then
            allocate(x(0,0))
            return
        end if
        if (n == 1) then
            allocate(x(1,1))
            x(1,1) = exp(a(1,1))
            return
        end if
        na = norm1(a)
        s = 0
        if (na > theta13) s = max(0, ceiling(log(na/theta13)/log(2.0_dp)))
        aa = a / (2.0_dp**s)
        id = eye(n)
        a2 = matmul(aa,aa)
        a4 = matmul(a2,a2)
        a6 = matmul(a2,a4)
        u = matmul(aa, matmul(a6, c(14)*a6 + c(12)*a4 + c(10)*a2) + &
            c(8)*a6 + c(6)*a4 + c(4)*a2 + c(2)*id)
        v = matmul(a6, c(13)*a6 + c(11)*a4 + c(9)*a2) + &
            c(7)*a6 + c(5)*a4 + c(3)*a2 + c(1)*id
        x = solve_linear(v-u, v+u)
        do k = 1, s
            x = matmul(x,x)
        end do
    end function expm

    subroutine expm_frechet(a, e, f, l)
        ! Block identity: exp([A,E;0,A]) = [exp(A), L_exp(A,E); 0, exp(A)].
        real(dp), intent(in) :: a(:,:), e(:,:)
        real(dp), allocatable, intent(out) :: f(:,:), l(:,:)
        real(dp), allocatable :: z(:,:), ez(:,:)
        integer :: n
        n = size(a,1)
        if (size(a,2) /= n .or. any(shape(e) /= [n,n])) error stop "expm_frechet: shape mismatch"
        allocate(z(2*n,2*n))
        z = 0.0_dp
        z(1:n,1:n) = a
        z(1:n,n+1:2*n) = e
        z(n+1:2*n,n+1:2*n) = a
        ez = expm(z)
        f = ez(1:n,1:n)
        l = ez(1:n,n+1:2*n)
    end subroutine expm_frechet
end module msm_linalg

! SPDX-License-Identifier: GPL-2.0-or-later
module ltsa_toeplitz
    use ltsa_kinds, only : dp
    use ltsa_status, only : ltsa_error, ltsa_success, ltsa_invalid_input, ltsa_not_positive_definite, set_error
    use ltsa_linalg, only : is_toeplitz, toeplitz_matrix, inverse_spd
    use ltsa_durbin_levinson, only : durbin_levinson_table
    implicit none
    private

    public :: trench_inverse, toeplitz_inverse_update, trench_quadratic_logdet, trench_mean

contains

    subroutine trench_inverse(g, gi, error)
        real(dp), intent(in) :: g(:,:)
        real(dp), allocatable, intent(out) :: gi(:,:)
        type(ltsa_error), intent(out) :: error
        real(dp), allocatable :: r(:), rn(:), table(:,:), pacf(:), variances(:), u(:), b(:,:)
        real(dp) :: scale, sigsq, unused_logdet
        integer :: i, j, n, n1
        error%code = ltsa_success
        error%message = ''
        if (size(g,1) /= size(g,2) .or. size(g,1) < 1) then
            allocate(gi(0,0))
            call set_error(error, ltsa_invalid_input, 'G must be a nonempty square matrix')
            return
        end if
        if (.not. is_toeplitz(g, 64.0_dp*epsilon(1.0_dp))) then
            allocate(gi(0,0))
            call set_error(error, ltsa_invalid_input, 'G must be symmetric Toeplitz')
            return
        end if
        n = size(g,1)
        scale = g(1,1)
        if (scale <= 0.0_dp) then
            allocate(gi(0,0))
            call set_error(error, ltsa_invalid_input, 'G(1,1) must be positive')
            return
        end if
        if (n == 1) then
            allocate(gi(1,1))
            gi(1,1) = 1.0_dp/scale
            return
        end if
        allocate(r(n), rn(n))
        r = g(1,:)
        rn = r/scale
        call durbin_levinson_table(rn, table, pacf, variances, error)
        if (.not. error%ok()) then
            allocate(gi(0,0))
            return
        end if
        sigsq = variances(n)
        allocate(u(n-1), b(n,n), source=0.0_dp)
        do i = 1, n-1
            u(i) = -table(n-1,n-i)/sigsq
        end do
        b(1,1) = 1.0_dp/sigsq
        do j = 2, n
            b(1,j) = u(n-j+1)
        end do
        n1 = (n-1)/2
        do i = 2, n1+1
            do j = i, n-i+1
                b(i,j) = b(i-1,j-1)+(u(n-j+1)*u(n-i+1)-u(i-1)*u(j-1))*sigsq
            end do
        end do
        do j = 1, n
            do i = 1, j
                if (i > (n+1-(j-i))/2) b(i,j) = b(n-j+1,n-i+1)
            end do
        end do
        do i = 1, n
            do j = 1, i-1
                b(i,j) = b(j,i)
            end do
        end do
        if (maxval(abs(matmul(g,b/scale)-identity_matrix(n))) > 1.0e-7_dp*real(n,dp)) then
            call inverse_spd(g, gi, unused_logdet, error)
            return
        end if
        allocate(gi(n,n))
        gi = b/scale
    end subroutine trench_inverse

    function identity_matrix(n) result(a)
        integer, intent(in) :: n
        real(dp) :: a(n,n)
        integer :: i
        a = 0.0_dp
        do i = 1, n
            a(i,i) = 1.0_dp
        end do
    end function identity_matrix

    subroutine toeplitz_inverse_update(gi, r, rnew, updated, error)
        real(dp), intent(in) :: gi(:,:), r(:), rnew
        real(dp), allocatable, intent(out) :: updated(:,:)
        type(ltsa_error), intent(out) :: error
        real(dp), allocatable :: g(:), gig(:)
        real(dp) :: denominator, e
        integer :: n
        error%code = ltsa_success
        error%message = ''
        n = size(r)
        if (size(gi,1) /= n .or. size(gi,2) /= n .or. n < 1) then
            allocate(updated(0,0))
            call set_error(error, ltsa_invalid_input, 'GI and r have incompatible dimensions')
            return
        end if
        allocate(g(n), gig(n))
        g(1) = rnew
        if (n > 1) g(2:n) = r(n:2:-1)
        gig = matmul(gi,g)
        denominator = r(1)-dot_product(g,gig)
        if (denominator <= epsilon(1.0_dp)*max(1.0_dp,abs(r(1)))) then
            allocate(updated(0,0))
            call set_error(error, ltsa_not_positive_definite, 'updated Toeplitz matrix is not positive definite')
            return
        end if
        e = 1.0_dp/denominator
        allocate(updated(n+1,n+1))
        updated(1:n,1:n) = gi+e*spread(gig,2,n)*spread(gig,1,n)
        updated(1:n,n+1) = -e*gig
        updated(n+1,1:n) = -e*gig
        updated(n+1,n+1) = e
    end subroutine toeplitz_inverse_update

    subroutine trench_quadratic_logdet(r, z, quadratic, logdet, error)
        real(dp), intent(in) :: r(:), z(:)
        real(dp), intent(out) :: quadratic, logdet
        type(ltsa_error), intent(out) :: error
        real(dp), allocatable :: rn(:), g(:,:), gi(:,:), table(:,:), pacf(:), variances(:)
        if (size(r) /= size(z) .or. size(r) < 1 .or. r(1) <= 0.0_dp) then
            quadratic = 0.0_dp
            logdet = 0.0_dp
            call set_error(error, ltsa_invalid_input, 'r and z must be equal-length and r(1) positive')
            return
        end if
        allocate(rn(size(r)))
        rn = r/r(1)
        g = toeplitz_matrix(rn)
        call trench_inverse(g, gi, error)
        if (.not. error%ok()) then
            quadratic = 0.0_dp
            logdet = 0.0_dp
            return
        end if
        call durbin_levinson_table(rn, table, pacf, variances, error)
        if (.not. error%ok()) then
            quadratic = 0.0_dp
            logdet = 0.0_dp
            return
        end if
        quadratic = dot_product(z,matmul(gi,z))
        logdet = sum(log(variances))
    end subroutine trench_quadratic_logdet

    function trench_mean(r, z, error) result(mean_value)
        real(dp), intent(in) :: r(:), z(:)
        type(ltsa_error), intent(out), optional :: error
        real(dp) :: mean_value
        type(ltsa_error) :: local_error
        real(dp), allocatable :: g(:,:), gi(:,:), one(:), weights(:)
        if (size(r) /= size(z) .or. size(r) < 1) then
            mean_value = 0.0_dp
            call set_error(local_error, ltsa_invalid_input, 'r and z must be nonempty and equal-length')
            if (present(error)) error = local_error
            return
        end if
        g = toeplitz_matrix(r)
        call trench_inverse(g, gi, local_error)
        if (.not. local_error%ok()) then
            mean_value = 0.0_dp
            if (present(error)) error = local_error
            return
        end if
        allocate(one(size(z)), source=1.0_dp)
        weights = matmul(gi,one)
        mean_value = dot_product(weights,z)/sum(weights)
        if (present(error)) error = local_error
    end function trench_mean

end module ltsa_toeplitz

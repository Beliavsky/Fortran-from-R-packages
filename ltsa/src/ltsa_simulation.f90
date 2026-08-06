! SPDX-License-Identifier: GPL-2.0-or-later
module ltsa_simulation
    use ltsa_kinds, only : dp
    use ltsa_status, only : ltsa_error, ltsa_success, ltsa_invalid_input, ltsa_dh_condition_failed, set_error
    use ltsa_random, only : ltsa_normal
    use ltsa_fft, only : fft_inplace, next_power_of_two
    implicit none
    private

    public :: sim_glp, dh_condition, dh_simulate

contains

    subroutine sim_glp(psi, innovations, z, error)
        real(dp), intent(in) :: psi(:), innovations(:)
        real(dp), allocatable, intent(out) :: z(:)
        type(ltsa_error), intent(out) :: error
        integer :: i, k, n, q
        error%code = ltsa_success
        error%message = ''
        q = size(psi)-1
        n = size(innovations)-q
        if (size(psi) < 1 .or. n < 1) then
            allocate(z(0))
            call set_error(error, ltsa_invalid_input, 'innovations must contain at least length(psi) observations')
            return
        end if
        allocate(z(n), source=0.0_dp)
        do i = 1, n
            do k = 0, q
                z(i) = z(i)+psi(k+1)*innovations(q+i-k)
            end do
        end do
    end subroutine sim_glp

    logical function dh_condition(n, r, eigenvalues) result(valid)
        integer, intent(in) :: n
        real(dp), intent(in) :: r(:)
        real(dp), allocatable, intent(out), optional :: eigenvalues(:)
        complex(dp), allocatable :: d(:)
        real(dp), allocatable :: acvf(:), lambda(:)
        integer :: emb, i, l, m
        valid = .false.
        if (n < 1 .or. size(r) < 1) then
            if (present(eigenvalues)) allocate(eigenvalues(0))
            return
        end if
        if (n == 1) then
            valid = r(1) >= 0.0_dp
            if (present(eigenvalues)) then
                allocate(eigenvalues(1))
                eigenvalues(1) = r(1)
            end if
            return
        end if
        emb = next_power_of_two(n-1)
        l = 2*emb
        m = min(size(r),emb+1)
        allocate(acvf(emb+1), source=0.0_dp)
        acvf(1:m) = r(1:m)
        allocate(d(l), source=cmplx(0.0_dp,0.0_dp,dp))
        do i = 1, emb+1
            d(i) = cmplx(acvf(i),0.0_dp,dp)
        end do
        do i = 2, emb
            d(l-i+2) = cmplx(acvf(i),0.0_dp,dp)
        end do
        call fft_inplace(d)
        allocate(lambda(l))
        lambda = real(d,dp)
        valid = minval(lambda) >= -128.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(lambda)))
        where (lambda < 0.0_dp .and. lambda > -128.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(lambda)))) lambda = 0.0_dp
        if (present(eigenvalues)) then
            allocate(eigenvalues(l))
            eigenvalues = lambda
        end if
    end function dh_condition

    subroutine dh_simulate(n, r, z, error, source_compatible)
        integer, intent(in) :: n
        real(dp), intent(in) :: r(:)
        real(dp), allocatable, intent(out) :: z(:)
        type(ltsa_error), intent(out) :: error
        logical, intent(in), optional :: source_compatible
        logical :: source_mode, valid
        real(dp), allocatable :: lambda(:)
        complex(dp), allocatable :: w(:)
        integer :: emb, l, k
        source_mode = .false.
        if (present(source_compatible)) source_mode = source_compatible
        error%code = ltsa_success
        error%message = ''
        if (n < 1 .or. size(r) < 1) then
            allocate(z(0))
            call set_error(error, ltsa_invalid_input, 'n and r must be nonempty')
            return
        end if
        if (n == 1) then
            if (r(1) < 0.0_dp) then
                allocate(z(0))
                call set_error(error, ltsa_dh_condition_failed, 'negative variance')
                return
            end if
            allocate(z(1))
            z(1) = sqrt(r(1))*ltsa_normal()
            return
        end if
        valid = dh_condition(n, r, lambda)
        if (.not. valid) then
            allocate(z(0))
            call set_error(error, ltsa_dh_condition_failed, 'Davies-Harte nonnegativity condition failed')
            return
        end if
        l = size(lambda)
        emb = l/2
        allocate(w(l), source=cmplx(0.0_dp,0.0_dp,dp))
        if (source_mode) then
            w(1) = sqrt(lambda(1))*cmplx(2.0_dp+sqrt(2.0_dp)*ltsa_normal(),0.0_dp,dp)
            w(emb+1) = sqrt(lambda(emb+1))*cmplx(2.0_dp+sqrt(2.0_dp)*ltsa_normal(),0.0_dp,dp)
            do k = 2, emb
                w(k) = sqrt(lambda(k))*cmplx(ltsa_normal(),ltsa_normal(),dp)
                w(l-k+2) = conjg(w(k))
            end do
            call fft_inplace(w, inverse=.true.)
            allocate(z(n))
            z = real(w(1:n),dp)*sqrt(real(l,dp)/2.0_dp)
        else
            w(1) = sqrt(lambda(1))*cmplx(ltsa_normal(),0.0_dp,dp)
            w(emb+1) = sqrt(lambda(emb+1))*cmplx(ltsa_normal(),0.0_dp,dp)
            do k = 2, emb
                w(k) = sqrt(0.5_dp*lambda(k))*cmplx(ltsa_normal(),ltsa_normal(),dp)
                w(l-k+2) = conjg(w(k))
            end do
            call fft_inplace(w, inverse=.true.)
            allocate(z(n))
            z = real(w(1:n),dp)*sqrt(real(l,dp))
        end if
    end subroutine dh_simulate

end module ltsa_simulation

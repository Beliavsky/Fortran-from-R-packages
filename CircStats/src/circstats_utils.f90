module circstats_utils
    use circstats_kinds, only: dp, pi, twopi
    implicit none
    private
    public :: wrap_2pi, wrap_pi, sort_real, sample_variance, sample_covariance
    public :: quantile_type7, solve_linear, inverse_matrix, ols_fit
    public :: randu, randn, randexp, randcauchy

contains

    pure elemental real(dp) function wrap_2pi(x) result(y)
        real(dp), intent(in) :: x
        y = modulo(x, twopi)
    end function wrap_2pi

    pure elemental real(dp) function wrap_pi(x) result(y)
        real(dp), intent(in) :: x
        y = modulo(x + pi, twopi) - pi
    end function wrap_pi

    subroutine sort_real(x)
        real(dp), intent(inout) :: x(:)
        integer :: i, j
        real(dp) :: key
        do i = 2, size(x)
            key = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= key) exit
                x(j+1) = x(j)
                j = j - 1
            end do
            x(j+1) = key
        end do
    end subroutine sort_real

    pure real(dp) function sample_variance(x) result(v)
        real(dp), intent(in) :: x(:)
        real(dp) :: m
        integer :: n
        n = size(x)
        if (n < 2) then
            v = 0.0_dp
            return
        end if
        m = sum(x)/real(n,dp)
        v = sum((x-m)**2)/real(n-1,dp)
    end function sample_variance

    pure real(dp) function sample_covariance(x,y) result(v)
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: mx, my
        integer :: n
        n = min(size(x),size(y))
        if (n < 2) then
            v = 0.0_dp
            return
        end if
        mx = sum(x(1:n))/real(n,dp)
        my = sum(y(1:n))/real(n,dp)
        v = sum((x(1:n)-mx)*(y(1:n)-my))/real(n-1,dp)
    end function sample_covariance

    function quantile_type7(x, p) result(q)
        real(dp), intent(in) :: x(:), p
        real(dp) :: q, h, frac
        real(dp), allocatable :: w(:)
        integer :: n, j
        n = size(x)
        if (n == 0) then
            q = 0.0_dp
            return
        end if
        allocate(w(n))
        w = x
        call sort_real(w)
        if (p <= 0.0_dp) then
            q = w(1)
        else if (p >= 1.0_dp) then
            q = w(n)
        else
            h = 1.0_dp + real(n-1,dp)*p
            j = floor(h)
            frac = h-real(j,dp)
            if (j >= n) then
                q = w(n)
            else
                q = (1.0_dp-frac)*w(j) + frac*w(j+1)
            end if
        end if
    end function quantile_type7

    subroutine solve_linear(a, b, x, info)
        real(dp), intent(in) :: a(:,:), b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: info
        real(dp), allocatable :: aug(:,:)
        real(dp) :: pivot, factor, tmp
        integer :: n, i, j, k, imax
        n = size(b)
        info = 0
        allocate(aug(n,n+1))
        aug(:,1:n) = a
        aug(:,n+1) = b
        do k = 1, n
            imax = k
            do i = k+1, n
                if (abs(aug(i,k)) > abs(aug(imax,k))) imax = i
            end do
            if (abs(aug(imax,k)) <= 100.0_dp*epsilon(1.0_dp)) then
                info = k
                x = 0.0_dp
                return
            end if
            if (imax /= k) then
                do j = k, n+1
                    tmp = aug(k,j)
                    aug(k,j) = aug(imax,j)
                    aug(imax,j) = tmp
                end do
            end if
            pivot = aug(k,k)
            aug(k,k:n+1) = aug(k,k:n+1)/pivot
            do i = 1, n
                if (i == k) cycle
                factor = aug(i,k)
                aug(i,k:n+1) = aug(i,k:n+1) - factor*aug(k,k:n+1)
            end do
        end do
        x = aug(:,n+1)
    end subroutine solve_linear

    subroutine inverse_matrix(a, ainv, info)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: ainv(:,:)
        integer, intent(out) :: info
        real(dp), allocatable :: e(:), x(:)
        integer :: n, j, stat
        n = size(a,1)
        allocate(e(n),x(n))
        ainv = 0.0_dp
        info = 0
        do j = 1, n
            e = 0.0_dp
            e(j) = 1.0_dp
            call solve_linear(a,e,x,stat)
            if (stat /= 0) then
                info = stat
                return
            end if
            ainv(:,j) = x
        end do
    end subroutine inverse_matrix

    subroutine ols_fit(x, y, beta, fitted, info)
        real(dp), intent(in) :: x(:,:), y(:)
        real(dp), intent(out) :: beta(:), fitted(:)
        integer, intent(out) :: info
        real(dp), allocatable :: xtx(:,:), xty(:)
        allocate(xtx(size(x,2),size(x,2)),xty(size(x,2)))
        xtx = matmul(transpose(x),x)
        xty = matmul(transpose(x),y)
        call solve_linear(xtx,xty,beta,info)
        if (info == 0) then
            fitted = matmul(x,beta)
        else
            fitted = 0.0_dp
        end if
    end subroutine ols_fit

    real(dp) function randu() result(u)
        call random_number(u)
        do while (u <= 0.0_dp .or. u >= 1.0_dp)
            call random_number(u)
        end do
    end function randu

    real(dp) function randn() result(z)
        real(dp) :: u1, u2
        u1 = randu()
        u2 = randu()
        z = sqrt(-2.0_dp*log(u1))*cos(twopi*u2)
    end function randn

    real(dp) function randexp() result(x)
        x = -log(randu())
    end function randexp

    real(dp) function randcauchy(location, scale) result(x)
        real(dp), intent(in) :: location, scale
        x = location + scale*tan(pi*(randu()-0.5_dp))
    end function randcauchy
end module circstats_utils

module fdausc_simulation
    use r_kinds, only : dp
    use r_linalg, only : symmetric_eigen
    implicit none
    private

    public :: simulate_ou_process, process_covariance, simulate_gaussian_process
    public :: wild_residuals, random_basis_combinations, grid_basis_combinations

contains

    subroutine simulate_ou_process(n, t, mu, alpha, sigma, x0, curves, info)
        integer, intent(in) :: n !! Number of independent Ornstein-Uhlenbeck trajectories to simulate.
        real(dp), intent(in) :: t(:) !! Nonnegative observation times for every trajectory.
        real(dp), intent(in) :: mu !! Long-run process mean.
        real(dp), intent(in) :: alpha !! Positive mean-reversion rate.
        real(dp), intent(in) :: sigma !! Positive diffusion scale.
        real(dp), intent(in) :: x0(:) !! Initial process value for each simulated trajectory; size must be n.
        real(dp), intent(out) :: curves(:, :) !! Simulated trajectories with shape [n,size(t)].
        integer, intent(out), optional :: info !! Zero on successful covariance eigendecomposition; nonzero otherwise.
        real(dp), allocatable :: covariance(:, :)
        real(dp), allocatable :: noise(:, :)
        integer :: i
        integer :: j
        integer :: stat

        allocate(covariance(size(t), size(t)), noise(n, size(t)))
        do i = 1, size(t)
            do j = 1, size(t)
                covariance(i, j) = sigma*sigma/(2.0_dp*alpha)*( &
                    exp(alpha*(2.0_dp*min(t(i), t(j)) - t(i) - t(j))) - exp(-alpha*(t(i) + t(j))))
            end do
        end do
        call simulate_gaussian_process(n, covariance, noise, stat)
        if (present(info)) info = stat
        if (stat /= 0) then
            curves = 0.0_dp
            return
        end if
        do i = 1, n
            curves(i, :) = noise(i, :) + mu + (x0(i) - mu)*exp(-alpha*t)
        end do
    end subroutine simulate_ou_process

    pure subroutine process_covariance(t, covariance_code, scale, theta, hurst, covariance)
        real(dp), intent(in) :: t(:) !! Observation times at which the requested process covariance is evaluated.
        integer, intent(in) :: covariance_code !! Covariance code: 1 Brownian, 2 OU, 3 exponential, 4 fractional Brownian.
        real(dp), intent(in) :: scale !! Positive process scale parameter.
        real(dp), intent(in) :: theta !! Positive OU rate or exponential correlation range parameter.
        real(dp), intent(in) :: hurst !! Fractional-Brownian Hurst exponent for covariance_code=4.
        real(dp), intent(out) :: covariance(:, :) !! Symmetric process covariance on t.
        integer :: i
        integer :: j

        do i = 1, size(t)
            do j = i, size(t)
                select case (covariance_code)
                case (2)
                    covariance(i, j) = scale/(2.0_dp*theta)*exp(-theta*(t(i) + t(j)))*( &
                        exp(2.0_dp*theta*min(t(i), t(j))) - 1.0_dp)
                case (3)
                    covariance(i, j) = scale*exp(-abs(t(i) - t(j))/theta)
                case (4)
                    covariance(i, j) = 0.5_dp*abs(scale)**hurst*(abs(t(i))**(2.0_dp*hurst) &
                        + abs(t(j))**(2.0_dp*hurst) - abs(t(i) - t(j))**(2.0_dp*hurst))
                case default
                    covariance(i, j) = scale*min(t(i), t(j))
                end select
                covariance(j, i) = covariance(i, j)
            end do
        end do
    end subroutine process_covariance

    subroutine simulate_gaussian_process(n, covariance, draws, info, mean_curve)
        integer, intent(in) :: n !! Number of independent Gaussian-process draws requested.
        real(dp), intent(in) :: covariance(:, :) !! Symmetric positive-semidefinite covariance matrix on the observation grid.
        real(dp), intent(out) :: draws(:, :) !! Simulated rows with shape [n,size(covariance,1)].
        integer, intent(out) :: info !! Zero on successful symmetric eigendecomposition; nonzero otherwise.
        real(dp), intent(in), optional :: mean_curve(:) !! Optional pointwise mean added to every generated draw.
        real(dp), allocatable :: values(:)
        real(dp), allocatable :: vectors(:, :)
        real(dp), allocatable :: root(:, :)
        real(dp), allocatable :: z(:)
        integer :: i
        integer :: j
        integer :: p

        call symmetric_eigen(covariance, values, vectors, info, descending=.true.)
        if (info /= 0) then
            draws = 0.0_dp
            return
        end if
        p = size(covariance, 1)
        allocate(root(p, p), z(p))
        root = 0.0_dp
        do j = 1, p
            root(:, j) = vectors(:, j)*sqrt(max(0.0_dp, values(j)))
        end do
        do i = 1, n
            call random_normal_vector(z)
            draws(i, :) = matmul(root, z)
            if (present(mean_curve)) draws(i, :) = draws(i, :) + mean_curve
        end do
    end subroutine simulate_gaussian_process

    subroutine wild_residuals(residuals, type_code, out)
        real(dp), intent(in) :: residuals(:) !! Residual vector multiplied by independent wild-bootstrap weights.
        integer, intent(in) :: type_code !! Weight selector: 1 golden two-point, 2 Rademacher, 3 standard normal.
        real(dp), intent(out) :: out(:) !! Wild-bootstrap residuals aligned with residuals.
        real(dp) :: u
        real(dp) :: z
        real(dp) :: a
        real(dp) :: b
        real(dp) :: prob_a
        integer :: i

        a = (1.0_dp - sqrt(5.0_dp))/2.0_dp
        b = (1.0_dp + sqrt(5.0_dp))/2.0_dp
        prob_a = (5.0_dp + sqrt(5.0_dp))/10.0_dp
        do i = 1, size(residuals)
            select case (type_code)
            case (2)
                call random_number(u)
                if (u < 0.5_dp) then
                    z = -1.0_dp
                else
                    z = 1.0_dp
                end if
            case (3)
                call random_normal_scalar(z)
            case default
                call random_number(u)
                if (u < prob_a) then
                    z = a
                else
                    z = b
                end if
            end select
            out(i) = residuals(i)*z
        end do
    end subroutine wild_residuals

    subroutine random_basis_combinations(n, basis_curves, sdarg, norm_target, mean_curve, curves)
        integer, intent(in) :: n !! Number of random functional combinations to generate.
        real(dp), intent(in) :: basis_curves(:, :) !! Basis or component curves in rows.
        real(dp), intent(in) :: sdarg(:) !! Standard-deviation multiplier for each basis coefficient.
        real(dp), intent(in) :: norm_target !! Euclidean norm imposed on each coefficient vector before sdarg scaling.
        real(dp), intent(in) :: mean_curve(:) !! Mean function added to every generated combination.
        real(dp), intent(out) :: curves(:, :) !! Generated functional observations with shape [n,size(basis_curves,2)].
        real(dp), allocatable :: coef(:)
        real(dp) :: cnorm
        integer :: i

        allocate(coef(size(basis_curves, 1)))
        do i = 1, n
            call random_normal_vector(coef)
            cnorm = sqrt(dot_product(coef, coef))
            if (cnorm > sqrt(tiny(1.0_dp))) coef = coef*norm_target/cnorm
            coef = coef*sdarg
            curves(i, :) = matmul(coef, basis_curves) + mean_curve
        end do
    end subroutine random_basis_combinations

    pure subroutine grid_basis_combinations(coefficients, basis_curves, mean_curve, curves)
        real(dp), intent(in) :: coefficients(:, :) !! User-supplied coefficient rows, one coefficient per basis curve.
        real(dp), intent(in) :: basis_curves(:, :) !! Basis or component curves in rows.
        real(dp), intent(in) :: mean_curve(:) !! Mean function added to every generated combination.
        real(dp), intent(out) :: curves(:, :) !! Functional combinations coefficients*basis_curves plus mean_curve.
        integer :: i

        curves = matmul(coefficients, basis_curves)
        do i = 1, size(curves, 1)
            curves(i, :) = curves(i, :) + mean_curve
        end do
    end subroutine grid_basis_combinations

    subroutine random_normal_vector(z)
        real(dp), intent(out) :: z(:) !! Independent standard normal variates generated with the Box-Muller transform.
        integer :: i
        real(dp) :: z1
        real(dp) :: z2

        i = 1
        do while (i <= size(z))
            call random_normal_pair(z1, z2)
            z(i) = z1
            if (i + 1 <= size(z)) z(i + 1) = z2
            i = i + 2
        end do
    end subroutine random_normal_vector

    subroutine random_normal_scalar(z)
        real(dp), intent(out) :: z !! One standard normal variate generated with the Box-Muller transform.
        real(dp) :: z2
        call random_normal_pair(z, z2)
    end subroutine random_normal_scalar

    subroutine random_normal_pair(z1, z2)
        real(dp), intent(out) :: z1 !! First independent standard normal variate from a Box-Muller pair.
        real(dp), intent(out) :: z2 !! Second independent standard normal variate from a Box-Muller pair.
        real(dp) :: u1
        real(dp) :: u2
        real(dp) :: radius
        real(dp), parameter :: pi = acos(-1.0_dp)

        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, tiny(1.0_dp))
        radius = sqrt(-2.0_dp*log(u1))
        z1 = radius*cos(2.0_dp*pi*u2)
        z2 = radius*sin(2.0_dp*pi*u2)
    end subroutine random_normal_pair

end module fdausc_simulation

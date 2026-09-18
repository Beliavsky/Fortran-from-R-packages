module tclust_simulation
    use tclust_kinds, only : dp
    use tclust_types, only : simulation_result
    use tclust_rng, only : set_tclust_seed, rand_uniform, rand_normal, random_orthonormal
    use tclust_linalg, only : mvn_sample
    implicit none
    private

    public :: simulate_tclust
    public :: simulate_rlg

contains

    subroutine simulate_tclust(n, result, p, k, type_id, balanced, seed)
        integer, intent(in) :: n !! Total number of observations, including approximately ten percent contamination.
        type(simulation_result), intent(out) :: result !! Simulated data matrix and true labels, with zero marking contamination.
        integer, intent(in), optional :: p !! Ambient dimension, at least two; default 4.
        integer, intent(in), optional :: k !! Number of clusters, either three or six; default 3.
        integer, intent(in), optional :: type_id !! Scenario one for spherical or two for elliptical clusters; default 2.
        integer, intent(in), optional :: balanced !! One for equal cluster sizes or two for 25/30/35-percent proportions.
        integer, intent(in), optional :: seed !! Optional deterministic seed for the intrinsic Fortran RNG.
        integer :: pp
        integer :: kk
        integer :: typ
        integer :: bal
        integer :: nn
        integer :: n1
        integer :: n2
        integer :: n3
        integer :: nh
        integer :: row
        integer :: i
        integer :: j
        real(dp) :: angle1
        real(dp) :: angle2
        real(dp) :: rot1(2, 2)
        real(dp) :: rot2(2, 2)

        if (n < 10) error stop 'simulate_tclust: n must be at least 10'
        pp = 4
        if (present(p)) pp = p
        kk = 3
        if (present(k)) kk = k
        typ = 2
        if (present(type_id)) typ = type_id
        bal = 1
        if (present(balanced)) bal = balanced
        if (pp < 2) error stop 'simulate_tclust: p must be at least two'
        if (kk /= 3 .and. kk /= 6) error stop 'simulate_tclust: k must be three or six'
        if (typ /= 1 .and. typ /= 2) error stop 'simulate_tclust: type_id must be one or two'
        if (bal /= 1 .and. bal /= 2) error stop 'simulate_tclust: balanced must be one or two'
        if (present(seed)) call set_tclust_seed(seed)

        if (kk == 6) then
            nn = n / 2
        else
            nn = n
        end if
        if (bal == 1) then
            n1 = ceiling(0.30_dp * real(nn, dp))
            n2 = n1
            n3 = n1
        else
            n1 = ceiling(0.25_dp * real(nn, dp))
            n2 = ceiling(0.30_dp * real(nn, dp))
            n3 = ceiling(0.35_dp * real(nn, dp))
        end if
        nh = n1 + n2 + n3
        if ((kk == 3 .and. nh > n) .or. (kk == 6 .and. 2 * nh > n)) then
            error stop 'simulate_tclust: cluster sizes exceed n for this small input'
        end if
        allocate(result%x(n, pp), result%true_cluster(n))
        result%x = 0.0_dp
        result%true_cluster = 0
        angle1 = rand_uniform(0.0_dp, 2.0_dp * acos(-1.0_dp))
        angle2 = rand_uniform(0.0_dp, 2.0_dp * acos(-1.0_dp))
        call rotation2(angle1, rot1)
        call rotation2(angle2, rot2)

        row = 1
        if (typ == 1) then
            call fill_cluster2(result%x, row, n1, [4.0_dp, 4.0_dp], identity2(), 1)
            call fill_cluster2(result%x, row, n2, [-4.0_dp, 4.0_dp], identity2(), 2)
            call fill_cluster2(result%x, row, n3, [0.0_dp, 0.0_dp], identity2(), 3)
            if (kk == 6) then
                call fill_cluster2(result%x, row, n1, [24.0_dp, 0.0_dp], identity2(), 4)
                call fill_cluster2(result%x, row, n2, [16.0_dp, 0.0_dp], identity2(), 5)
                call fill_cluster2(result%x, row, n3, [20.0_dp, 4.0_dp], identity2(), 6)
            end if
        else
            call fill_cluster2(result%x, row, n1, [20.0_dp, 20.0_dp], &
                               matmul(transpose(rot1), matmul(diag2(1.0_dp, 81.0_dp), rot1)), 1)
            call fill_cluster2(result%x, row, n2, [-20.0_dp, -20.0_dp], &
                               matmul(transpose(rot1), matmul(diag2(81.0_dp, 1.0_dp), rot1)), 2)
            call fill_cluster2(result%x, row, n3, [0.0_dp, 0.0_dp], diag2(9.0_dp, 9.0_dp), 3)
            if (kk == 6) then
                call fill_cluster2(result%x, row, n1, [40.0_dp, 20.0_dp], &
                                   matmul(transpose(rot2), matmul(diag2(1.0_dp, 81.0_dp), rot2)), 4)
                call fill_cluster2(result%x, row, n2, [0.0_dp, -20.0_dp], &
                                   matmul(transpose(rot2), matmul(diag2(81.0_dp, 1.0_dp), rot2)), 5)
                call fill_cluster2(result%x, row, n3, [20.0_dp, 0.0_dp], diag2(9.0_dp, 9.0_dp), 6)
            end if
        end if
        do i = row, n
            do j = 1, pp
                if (typ == 1) then
                    if (kk == 3) then
                        result%x(i, j) = rand_uniform(-20.0_dp, 20.0_dp)
                    else
                        result%x(i, j) = rand_uniform(-20.0_dp, 40.0_dp)
                    end if
                else
                    if (kk == 3) then
                        result%x(i, j) = rand_uniform(-50.0_dp, 50.0_dp)
                    else
                        result%x(i, j) = rand_uniform(-60.0_dp, 80.0_dp)
                    end if
                end if
            end do
        end do
        if (pp > 2) then
            do i = 1, n
                do j = 3, pp
                    result%x(i, j) = rand_normal()
                end do
            end do
        end if
        if (kk == 3) then
            result%true_cluster(1:n1) = 1
            result%true_cluster(n1 + 1:n1 + n2) = 2
            result%true_cluster(n1 + n2 + 1:nh) = 3
        else
            result%true_cluster(1:n1) = 1
            result%true_cluster(n1 + 1:n1 + n2) = 2
            result%true_cluster(n1 + n2 + 1:nh) = 3
            result%true_cluster(nh + 1:nh + n1) = 4
            result%true_cluster(nh + n1 + 1:nh + n1 + n2) = 5
            result%true_cluster(nh + n1 + n2 + 1:2 * nh) = 6
        end if

    contains

        subroutine fill_cluster2(y, pos, m, mean2, cov2, label)
            real(dp), intent(inout) :: y(:, :) !! Full simulation matrix receiving generated cluster rows.
            integer, intent(inout) :: pos !! One-based insertion row, advanced by m on return.
            integer, intent(in) :: m !! Number of observations generated for this cluster.
            real(dp), intent(in) :: mean2(2) !! Mean vector for the first two coordinates.
            real(dp), intent(in) :: cov2(2, 2) !! Covariance matrix for the first two coordinates.
            integer, intent(in) :: label !! Positive true-cluster label assigned to generated rows.
            real(dp) :: draw(2)
            integer :: ii

            do ii = 1, m
                call mvn_sample(mean2, cov2, draw)
                y(pos, 1:2) = draw
                result%true_cluster(pos) = label
                pos = pos + 1
            end do
        end subroutine fill_cluster2

    end subroutine simulate_tclust

    subroutine simulate_rlg(result, q, p, n, variance, sep_means, alpha, seed)
        type(simulation_result), intent(out) :: result !! Simulated robust-linear-grouping data and true labels.
        integer, intent(in), optional :: q !! Common intrinsic cluster dimension; default 2 and strictly less than p.
        integer, intent(in), optional :: p !! Ambient dimension; default 10.
        integer, intent(in), optional :: n !! Total number of observations; default 200.
        real(dp), intent(in), optional :: variance !! Isotropic noise variance around each affine subspace; default 0.01.
        real(dp), intent(in), optional :: sep_means !! Standard deviation of random cluster-location shifts. Default zero.
        real(dp), intent(in), optional :: alpha !! Contamination fraction; default 0.05.
        integer, intent(in), optional :: seed !! Optional deterministic seed for the intrinsic Fortran RNG.
        real(dp), allocatable :: r1(:, :)
        real(dp), allocatable :: r2(:, :)
        real(dp), allocatable :: base(:, :)
        real(dp), allocatable :: u1(:, :)
        real(dp), allocatable :: u2(:, :)
        real(dp), allocatable :: u3(:, :)
        real(dp), allocatable :: mu(:)
        real(dp), allocatable :: z(:)
        real(dp) :: var
        real(dp) :: sep
        real(dp) :: a
        real(dp) :: sd
        real(dp) :: rr
        integer :: qq
        integer :: pp
        integer :: nn
        integer :: partn
        integer :: ncont
        integer :: n3
        integer :: row
        integer :: i
        integer :: j

        qq = 2
        if (present(q)) qq = q
        pp = 10
        if (present(p)) pp = p
        nn = 200
        if (present(n)) nn = n
        var = 0.01_dp
        if (present(variance)) var = variance
        sep = 0.0_dp
        if (present(sep_means)) sep = sep_means
        a = 0.05_dp
        if (present(alpha)) a = alpha
        if (qq < 0 .or. qq >= pp) error stop 'simulate_rlg: q must satisfy 0 <= q < p'
        if (nn < 6 .or. var < 0.0_dp .or. a < 0.0_dp .or. a >= 1.0_dp) error stop 'simulate_rlg: invalid parameters'
        if (present(seed)) call set_tclust_seed(seed)
        partn = max(2, floor(real(nn, dp) * (1.0_dp - a) / 3.0_dp))
        ncont = floor(real(nn, dp) * a)
        n3 = nn - 2 * partn - ncont
        if (n3 < 2) error stop 'simulate_rlg: n is too small for requested contamination'
        sd = sqrt(var)
        allocate(r1(pp, pp), r2(pp, pp), base(pp, pp), u1(pp, qq), u2(pp, qq), u3(pp, qq), mu(pp), z(qq))
        call random_orthonormal(r1)
        call random_orthonormal(r2)
        call random_orthonormal(base)
        if (qq > 0) then
            u1 = base(:, 1:qq)
            u2 = matmul(r1, u1)
            u3 = matmul(r2, u2)
        end if
        allocate(result%x(nn, pp), result%true_cluster(nn))
        result%x = 0.0_dp
        result%true_cluster = 0
        row = 1
        call fill_affine_cluster(partn, u1, 1)
        call fill_affine_cluster(partn, u2, 2)
        call fill_affine_cluster(n3, u3, 3)
        rr = max(1.0_dp, maxval(abs(result%x(1:row - 1, :))))
        do i = row, nn
            do j = 1, pp
                result%x(i, j) = rand_uniform(-rr, rr)
            end do
        end do

    contains

        subroutine fill_affine_cluster(m, u, label)
            integer, intent(in) :: m !! Number of observations generated for this affine cluster.
            real(dp), intent(in) :: u(:, :) !! p-by-q orthonormal basis of the affine component.
            integer, intent(in) :: label !! Positive true-cluster label written for generated observations.
            integer :: ii
            integer :: jj

            do jj = 1, pp
                mu(jj) = sep * rand_normal()
            end do
            do ii = 1, m
                do jj = 1, qq
                    z(jj) = rand_uniform(0.0_dp, 1.0_dp)
                end do
                result%x(row, :) = mu
                if (qq > 0) result%x(row, :) = result%x(row, :) + matmul(u, z)
                do jj = 1, pp
                    result%x(row, jj) = result%x(row, jj) + sd * rand_normal()
                end do
                result%true_cluster(row) = label
                row = row + 1
            end do
        end subroutine fill_affine_cluster

    end subroutine simulate_rlg

    pure subroutine rotation2(angle, r)
        real(dp), intent(in) :: angle !! Rotation angle in radians.
        real(dp), intent(out) :: r(2, 2) !! Two-dimensional orthogonal rotation matrix.

        r(1, 1) = cos(angle)
        r(1, 2) = -sin(angle)
        r(2, 1) = sin(angle)
        r(2, 2) = cos(angle)
    end subroutine rotation2

    pure function identity2() result(a)
        real(dp) :: a(2, 2)

        a = 0.0_dp
        a(1, 1) = 1.0_dp
        a(2, 2) = 1.0_dp
    end function identity2

    pure function diag2(a1, a2) result(a)
        real(dp), intent(in) :: a1 !! First diagonal element.
        real(dp), intent(in) :: a2 !! Second diagonal element.
        real(dp) :: a(2, 2)

        a = 0.0_dp
        a(1, 1) = a1
        a(2, 2) = a2
    end function diag2

end module tclust_simulation

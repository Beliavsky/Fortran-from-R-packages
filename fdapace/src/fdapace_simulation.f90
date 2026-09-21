module fdapace_simulation
    use fdapace_basis, only : create_basis
    use fdapace_kinds, only : dp
    use fdapace_math, only : mean_value, normal_random, pi_dp, sample_sd
    use fdapace_types, only : fpca_inputs, gp_functional_data, sparse_gp_result, sparse_sample
    implicit none
    private

    public :: make_fpca_inputs_dense
    public :: make_gp_functional_data
    public :: make_sparse_gp
    public :: sparsify
    public :: wiener

contains

    function make_fpca_inputs_dense(t_vec, y_vec) result(inputs)
        real(dp), intent(in) :: t_vec(:) !! Common observation grid corresponding to columns of y_vec.
        real(dp), intent(in) :: y_vec(:,:) !! Dense observations with subjects in rows and common times in columns.
        type(fpca_inputs) :: inputs
        integer :: i

        if (size(y_vec, 2) /= size(t_vec)) error stop "make_fpca_inputs_dense: columns of y_vec do not match t_vec"
        allocate(inputs%lid(size(y_vec, 1)), inputs%ly(size(y_vec, 1)), inputs%lt(size(y_vec, 1)))
        do i = 1, size(y_vec, 1)
            inputs%lid(i) = i
            allocate(inputs%ly(i)%v(size(t_vec)), inputs%lt(i)%v(size(t_vec)))
            inputs%ly(i)%v = y_vec(i, :)
            inputs%lt(i)%v = t_vec
        end do
    end function make_fpca_inputs_dense

    function make_gp_functional_data(n, m, mu, k, lambda, sigma, basis_type) result(data)
        integer, intent(in) :: n !! Number of simulated curves; must be at least two.
        integer, intent(in) :: m !! Number of equidistant readings per curve; the R interface expects at least twenty.
        real(dp), intent(in) :: mu(:) !! Mean curve evaluated at the m equidistant grid points.
        integer, intent(in) :: k !! Number of basis functions and score components.
        real(dp), intent(in) :: lambda(:) !! Nonnegative component variances, one for each of the k basis functions.
        real(dp), intent(in), optional :: sigma !! Gaussian measurement-noise standard deviation; defaults to zero.
        character(len=*), intent(in), optional :: basis_type !! Basis family passed to create_basis; defaults to cos.
        type(gp_functional_data) :: data
        character(len=16) :: btype
        real(dp) :: noise_sd
        real(dp) :: sdv
        integer :: i
        integer :: j

        if (n < 2) error stop "make_gp_functional_data: n must be at least two"
        if (m < 20) error stop "make_gp_functional_data: m must be at least twenty"
        if (size(mu) /= m) error stop "make_gp_functional_data: mu length does not match m"
        if (size(lambda) /= k) error stop "make_gp_functional_data: lambda length does not match k"
        if (any(lambda < 0.0_dp)) error stop "make_gp_functional_data: lambda must be nonnegative"
        noise_sd = 0.0_dp
        if (present(sigma)) noise_sd = sigma
        if (noise_sd < 0.0_dp) error stop "make_gp_functional_data: sigma must be nonnegative"
        btype = "cos"
        if (present(basis_type)) btype = trim(basis_type)

        allocate(data%pts(m), data%xi(n, k), data%y(n, m))
        if (m == 1) then
            data%pts = 0.0_dp
        else
            do i = 1, m
                data%pts(i) = real(i - 1, dp) / real(m - 1, dp)
            end do
        end if
        data%phi = create_basis(k, data%pts, btype)
        do j = 1, k
            do i = 1, n
                data%xi(i, j) = normal_random()
            end do
            data%xi(:, j) = data%xi(:, j) - mean_value(data%xi(:, j))
            sdv = sample_sd(data%xi(:, j))
            if (sdv <= 0.0_dp) error stop "make_gp_functional_data: degenerate random score column"
            data%xi(:, j) = data%xi(:, j) / sdv * sqrt(lambda(j))
        end do
        data%y = matmul(data%xi, transpose(data%phi))
        do i = 1, n
            data%y(i, :) = data%y(i, :) + mu
        end do
        if (noise_sd > 0.0_dp) then
            allocate(data%yn(n, m))
            do j = 1, m
                do i = 1, n
                    data%yn(i, j) = data%y(i, j) + noise_sd * normal_random()
                end do
            end do
            data%has_noisy = .true.
        end if
    end function make_gp_functional_data


    function make_sparse_gp(n, sparsity, k, lambda, sigma, basis_type) result(data)
        integer, intent(in) :: n !! Number of sparse Gaussian-process subjects; must be at least two.
        integer, intent(in) :: sparsity(:) !! Candidate positive observation counts sampled with equal probability for each subject.
        integer, intent(in) :: k !! Number of basis-score components in the finite-rank covariance model.
        real(dp), intent(in) :: lambda(:) !! Nonnegative component variances, one for each basis function.
        real(dp), intent(in), optional :: sigma !! Optional nonnegative Gaussian measurement-noise standard deviation; defaults to zero.
        character(len=*), intent(in), optional :: basis_type !! Basis family passed to create_basis; defaults to cos.
        type(sparse_gp_result) :: data
        character(len=16) :: btype
        real(dp), allocatable :: phi(:,:)
        real(dp), allocatable :: truth(:)
        real(dp) :: noise_sd
        real(dp) :: sdv
        real(dp) :: u
        integer :: i
        integer :: j
        integer :: ni
        integer :: pick

        if (n < 2) error stop "make_sparse_gp: n must be at least two"
        if (size(sparsity) < 1 .or. any(sparsity < 1)) error stop "make_sparse_gp: sparsity values must be positive"
        if (k < 1 .or. size(lambda) /= k) error stop "make_sparse_gp: lambda length must equal k"
        if (any(lambda < 0.0_dp)) error stop "make_sparse_gp: lambda must be nonnegative"
        noise_sd = 0.0_dp
        if (present(sigma)) noise_sd = sigma
        if (noise_sd < 0.0_dp) error stop "make_sparse_gp: sigma must be nonnegative"
        btype = "cos"
        if (present(basis_type)) btype = trim(basis_type)
        allocate(data%xi(n, k), data%ni(n), data%sample%ly(n), data%sample%lt(n))
        if (noise_sd > 0.0_dp) then
            allocate(data%true_sample%ly(n), data%true_sample%lt(n))
            data%has_true = .true.
        end if
        do j = 1, k
            do i = 1, n
                data%xi(i, j) = normal_random()
            end do
            data%xi(:, j) = data%xi(:, j) - mean_value(data%xi(:, j))
            sdv = sample_sd(data%xi(:, j))
            if (sdv <= 0.0_dp) error stop "make_sparse_gp: degenerate random score column"
            data%xi(:, j) = data%xi(:, j) / sdv * sqrt(lambda(j))
        end do
        do i = 1, n
            call random_number(u)
            pick = min(size(sparsity), int(u * real(size(sparsity), dp)) + 1)
            ni = sparsity(pick)
            data%ni(i) = ni
            allocate(data%sample%lt(i)%v(ni), data%sample%ly(i)%v(ni))
            do j = 1, ni
                call random_number(data%sample%lt(i)%v(j))
            end do
            call sort_pairs_by_time(data%sample%lt(i)%v)
            phi = create_basis(k, data%sample%lt(i)%v, btype)
            truth = matmul(phi, data%xi(i, :))
            data%sample%ly(i)%v = truth
            if (noise_sd > 0.0_dp) then
                allocate(data%true_sample%lt(i)%v(ni), data%true_sample%ly(i)%v(ni))
                data%true_sample%lt(i)%v = data%sample%lt(i)%v
                data%true_sample%ly(i)%v = truth
                do j = 1, ni
                    data%sample%ly(i)%v(j) = truth(j) + noise_sd * normal_random()
                end do
            end if
        end do
    end function make_sparse_gp

    pure subroutine sort_pairs_by_time(t)
        real(dp), intent(inout) :: t(:) !! Observation times sorted into ascending order in place.
        real(dp) :: key
        integer :: i
        integer :: j

        do i = 2, size(t)
            key = t(i)
            j = i - 1
            do while (j >= 1)
                if (t(j) <= key) exit
                t(j + 1) = t(j)
                j = j - 1
            end do
            t(j + 1) = key
        end do
    end subroutine sort_pairs_by_time

    function wiener(n, pts, k) result(sample)
        integer, intent(in) :: n !! Number of independent standard Wiener-process sample paths.
        real(dp), intent(in) :: pts(:) !! Support points, normally within [0,1], at which paths are evaluated.
        integer, intent(in), optional :: k !! Number of Karhunen-Loeve components; defaults to fifty.
        real(dp), allocatable :: sample(:,:)
        real(dp), allocatable :: basis(:,:)
        real(dp), allocatable :: coeff(:)
        real(dp), allocatable :: normals(:,:)
        integer :: kk
        integer :: i
        integer :: j

        if (n < 1) error stop "wiener: n must be positive"
        kk = 50
        if (present(k)) kk = k
        if (kk < 1) error stop "wiener: k must be positive"
        allocate(basis(size(pts), kk), coeff(kk), normals(kk, n), sample(n, size(pts)))
        do j = 1, kk
            coeff(j) = 1.0_dp / ((real(j, dp) - 0.5_dp) * pi_dp)
            basis(:, j) = sqrt(2.0_dp) * sin(pts * (real(j, dp) - 0.5_dp) * pi_dp)
            do i = 1, n
                normals(j, i) = normal_random()
            end do
        end do
        do j = 1, kk
            basis(:, j) = basis(:, j) * coeff(j)
        end do
        sample = transpose(matmul(basis, normals))
    end function wiener

    function sparsify(sample, pts, sparsity, aggressive, fragment) result(out)
        real(dp), intent(in) :: sample(:,:) !! Dense functional observations with one curve per row.
        real(dp), intent(in) :: pts(:) !! Observation grid corresponding to sample columns.
        integer, intent(in) :: sparsity(:) !! Candidate positive numbers of retained observations, selected with equal probability per curve.
        logical, intent(in), optional :: aggressive !! If true, repeatedly sample until retained indices satisfy the upstream remote-spacing criterion.
        real(dp), intent(in), optional :: fragment !! Positive approximate fragment fraction; mutually exclusive with aggressive mode.
        type(sparse_sample) :: out
        logical :: remote
        real(dp) :: frag
        integer, allocatable :: candidate(:)
        integer, allocatable :: indices(:)
        integer :: i
        integer :: n_keep
        integer :: ns
        integer :: nuse
        integer :: tries
        real(dp) :: mid
        real(dp) :: ran
        real(dp) :: span

        if (size(sample, 2) /= size(pts)) error stop "sparsify: sample columns do not match pts"
        if (size(sparsity) < 1 .or. any(sparsity < 1)) error stop "sparsify: sparsity values must be positive"
        remote = .false.
        if (present(aggressive)) remote = aggressive
        frag = 0.0_dp
        if (present(fragment)) frag = fragment
        if (remote .and. frag > 0.0_dp) error stop "sparsify: aggressive and fragment modes are mutually exclusive"
        if (frag < 0.0_dp) error stop "sparsify: fragment must be nonnegative"

        allocate(out%lt(size(sample, 1)), out%ly(size(sample, 1)))
        span = 0.0_dp
        if (size(pts) > 0) span = maxval(pts) - minval(pts)
        do i = 1, size(sample, 1)
            call random_number(ran)
            ns = min(size(sparsity), int(ran * real(size(sparsity), dp)) + 1)
            n_keep = sparsity(ns)
            if (frag > 0.0_dp) then
                call random_number(ran)
                mid = minval(pts) + ran * span
                candidate = pack([(ns, ns=1, size(pts))], &
                    pts >= mid - 0.5_dp * span * frag .and. pts <= mid + 0.5_dp * span * frag)
                if (n_keep >= size(candidate)) then
                    indices = candidate
                else
                    call sample_without_replacement(candidate, n_keep, indices)
                    call sort_int(indices)
                end if
            else if (remote) then
                if (n_keep > size(pts)) error stop "sparsify: requested sparsity exceeds grid size"
                tries = 0
                do
                    tries = tries + 1
                    allocate(candidate(size(pts)))
                    candidate = [(ns, ns=1, size(pts))]
                    call sample_without_replacement(candidate, n_keep, indices)
                    deallocate(candidate)
                    call sort_int(indices)
                    if (n_keep <= 1) exit
                    if (real(minval(indices(2:) - indices(:n_keep - 1)), dp) >= &
                        (1.0_dp / real(n_keep, dp))**1.5_dp * real(size(pts), dp)) exit
                    if (tries >= 100000) error stop "sparsify: aggressive spacing criterion could not be met"
                end do
            else
                if (n_keep > size(pts)) error stop "sparsify: requested sparsity exceeds grid size"
                allocate(candidate(size(pts)))
                candidate = [(ns, ns=1, size(pts))]
                call sample_without_replacement(candidate, n_keep, indices)
                deallocate(candidate)
                call sort_int(indices)
            end if
            nuse = size(indices)
            allocate(out%lt(i)%v(nuse), out%ly(i)%v(nuse))
            if (nuse > 0) then
                out%lt(i)%v = pts(indices)
                out%ly(i)%v = sample(i, indices)
            end if
            if (allocated(indices)) deallocate(indices)
            if (allocated(candidate)) deallocate(candidate)
        end do
    end function sparsify

    subroutine sample_without_replacement(population, n_keep, sample_indices)
        integer, intent(in) :: population(:) !! Distinct integer population values from which to sample.
        integer, intent(in) :: n_keep !! Number of values to select without replacement.
        integer, allocatable, intent(out) :: sample_indices(:) !! Selected population values in random order.
        integer, allocatable :: work(:)
        integer :: i
        integer :: j
        integer :: tmp
        real(dp) :: u

        if (n_keep < 0 .or. n_keep > size(population)) error stop "sample_without_replacement: invalid sample size"
        allocate(work(size(population)))
        work = population
        do i = 1, n_keep
            call random_number(u)
            j = i + min(size(population) - i, int(u * real(size(population) - i + 1, dp)))
            tmp = work(i)
            work(i) = work(j)
            work(j) = tmp
        end do
        allocate(sample_indices(n_keep))
        if (n_keep > 0) sample_indices = work(1:n_keep)
    end subroutine sample_without_replacement

    pure subroutine sort_int(x)
        integer, intent(inout) :: x(:) !! Integer vector sorted into ascending order in place.
        integer :: i
        integer :: j
        integer :: key

        do i = 2, size(x)
            key = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= key) exit
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = key
        end do
    end subroutine sort_int

end module fdapace_simulation

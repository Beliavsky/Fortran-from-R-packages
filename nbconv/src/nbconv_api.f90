! SPDX-License-Identifier: GPL-3.0-or-later
module nbconv_api
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use nbconv_kinds, only : dp
    use nbconv_exact, only : nb_sum_exact
    use nbconv_approximations, only : nb_sum_moments, nb_sum_saddlepoint
    use nbconv_rng, only : negbin_rng_mu, nbconv_seed
    implicit none
    private

    type, public :: nbconv_summary
        real(dp) :: mean = 0.0_dp
        real(dp) :: variance = 0.0_dp
        real(dp) :: skewness = 0.0_dp
        real(dp) :: excess_kurtosis = 0.0_dp
        real(dp) :: k_mean = 0.0_dp
    end type nbconv_summary

    public :: dnbconv, dnbconv_mu, dnbconv_p
    public :: pnbconv, pnbconv_mu, pnbconv_p
    public :: qnbconv, qnbconv_mu, qnbconv_p
    public :: rnbconv, rnbconv_mu, rnbconv_p
    public :: nbconv_params, nbconv_params_mu, nbconv_params_p
    public :: nbconv_seed

contains

    function dnbconv(counts, values, phis, parameterization, method, n_terms, tolerance, normalize) result(pmf)
        integer, intent(in) :: counts(:)
        real(dp), intent(in) :: values(:), phis(:)
        character(len=*), intent(in), optional :: parameterization, method
        integer, intent(in), optional :: n_terms
        real(dp), intent(in), optional :: tolerance
        logical, intent(in), optional :: normalize
        real(dp), allocatable :: pmf(:)
        character(len=8) :: par

        par = get_parameterization(parameterization)
        if (trim(par) == "mu") then
            pmf = dnbconv_mu(counts, values, phis, method, n_terms, tolerance, normalize)
        else
            pmf = dnbconv_p(counts, values, phis, method, n_terms, tolerance, normalize)
        end if
    end function dnbconv

    function pnbconv(quants, values, phis, parameterization, method, n_terms, tolerance, normalize) result(cdfq)
        integer, intent(in) :: quants(:)
        real(dp), intent(in) :: values(:), phis(:)
        character(len=*), intent(in), optional :: parameterization, method
        integer, intent(in), optional :: n_terms
        real(dp), intent(in), optional :: tolerance
        logical, intent(in), optional :: normalize
        real(dp), allocatable :: cdfq(:)
        character(len=8) :: par

        par = get_parameterization(parameterization)
        if (trim(par) == "mu") then
            cdfq = pnbconv_mu(quants, values, phis, method, n_terms, tolerance, normalize)
        else
            cdfq = pnbconv_p(quants, values, phis, method, n_terms, tolerance, normalize)
        end if
    end function pnbconv

    function qnbconv(probs, max_count, values, phis, parameterization, method, n_terms, tolerance, normalize) result(quants)
        real(dp), intent(in) :: probs(:)
        integer, intent(in) :: max_count
        real(dp), intent(in) :: values(:), phis(:)
        character(len=*), intent(in), optional :: parameterization, method
        integer, intent(in), optional :: n_terms
        real(dp), intent(in), optional :: tolerance
        logical, intent(in), optional :: normalize
        integer, allocatable :: quants(:)
        character(len=8) :: par

        par = get_parameterization(parameterization)
        if (trim(par) == "mu") then
            quants = qnbconv_mu(probs, max_count, values, phis, method, n_terms, tolerance, normalize)
        else
            quants = qnbconv_p(probs, max_count, values, phis, method, n_terms, tolerance, normalize)
        end if
    end function qnbconv

    function rnbconv(n_samp, values, phis, parameterization) result(samples)
        integer, intent(in) :: n_samp
        real(dp), intent(in) :: values(:), phis(:)
        character(len=*), intent(in), optional :: parameterization
        integer, allocatable :: samples(:)
        character(len=8) :: par

        par = get_parameterization(parameterization)
        if (trim(par) == "mu") then
            samples = rnbconv_mu(n_samp, values, phis)
        else
            samples = rnbconv_p(n_samp, values, phis)
        end if
    end function rnbconv

    function nbconv_params(values, phis, parameterization) result(params)
        real(dp), intent(in) :: values(:), phis(:)
        character(len=*), intent(in), optional :: parameterization
        type(nbconv_summary) :: params
        character(len=8) :: par

        par = get_parameterization(parameterization)
        if (trim(par) == "mu") then
            params = nbconv_params_mu(values, phis)
        else
            params = nbconv_params_p(values, phis)
        end if
    end function nbconv_params

    function dnbconv_mu(counts, mus, phis, method, n_terms, tolerance, normalize) result(pmf)
        integer, intent(in) :: counts(:)
        real(dp), intent(in) :: mus(:), phis(:)
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: n_terms
        real(dp), intent(in), optional :: tolerance
        logical, intent(in), optional :: normalize
        real(dp), allocatable :: pmf(:)
        real(dp), allocatable :: ps(:)
        character(len=16) :: meth
        integer :: nt
        real(dp) :: tol
        logical :: norm

        call validate_mu(mus, phis)
        meth = get_method(method)
        nt = 1000
        if (present(n_terms)) nt = n_terms
        tol = 1.0e-3_dp
        if (present(tolerance)) tol = tolerance
        norm = .true.
        if (present(normalize)) norm = normalize

        select case (trim(meth))
        case ("exact")
            allocate(ps(size(mus)))
            ps = phis / (phis + mus)
            pmf = nb_sum_exact(ps, phis, counts, n_terms=nt, tolerance=tol)
        case ("moments")
            pmf = nb_sum_moments(mus, phis, counts)
        case ("saddlepoint")
            pmf = nb_sum_saddlepoint(mus, phis, counts, normalize=norm)
        case default
            error stop "dnbconv_mu: unknown method"
        end select
    end function dnbconv_mu

    function dnbconv_p(counts, ps, phis, method, n_terms, tolerance, normalize) result(pmf)
        integer, intent(in) :: counts(:)
        real(dp), intent(in) :: ps(:), phis(:)
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: n_terms
        real(dp), intent(in), optional :: tolerance
        logical, intent(in), optional :: normalize
        real(dp), allocatable :: pmf(:)
        real(dp), allocatable :: mus(:)
        character(len=16) :: meth
        integer :: nt
        real(dp) :: tol
        logical :: norm

        call validate_p(ps, phis)
        meth = get_method(method)
        nt = 1000
        if (present(n_terms)) nt = n_terms
        tol = 1.0e-3_dp
        if (present(tolerance)) tol = tolerance
        norm = .true.
        if (present(normalize)) norm = normalize

        select case (trim(meth))
        case ("exact")
            pmf = nb_sum_exact(ps, phis, counts, n_terms=nt, tolerance=tol)
        case ("moments", "saddlepoint")
            allocate(mus(size(ps)))
            mus = phis * (1.0_dp - ps) / ps
            if (trim(meth) == "moments") then
                pmf = nb_sum_moments(mus, phis, counts)
            else
                pmf = nb_sum_saddlepoint(mus, phis, counts, normalize=norm)
            end if
        case default
            error stop "dnbconv_p: unknown method"
        end select
    end function dnbconv_p

    function pnbconv_mu(quants, mus, phis, method, n_terms, tolerance, normalize) result(cdfq)
        integer, intent(in) :: quants(:)
        real(dp), intent(in) :: mus(:), phis(:)
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: n_terms
        real(dp), intent(in), optional :: tolerance
        logical, intent(in), optional :: normalize
        real(dp), allocatable :: cdfq(:)
        real(dp), allocatable :: pmf(:), cdf(:)
        integer, allocatable :: counts(:)
        integer :: i, maxq

        if (size(quants) == 0) then
            allocate(cdfq(0))
            return
        end if
        if (any(quants < 0)) error stop "pnbconv_mu: quants must be nonnegative"
        maxq = maxval(quants)
        allocate(counts(0:maxq), cdf(0:maxq), cdfq(size(quants)))
        counts = [(i, i=0,maxq)]
        pmf = dnbconv_mu(counts, mus, phis, method, n_terms, tolerance, normalize)
        cdf(0) = pmf(1)
        do i = 1, maxq
            cdf(i) = cdf(i - 1) + pmf(i + 1)
        end do
        do i = 1, size(quants)
            cdfq(i) = cdf(quants(i))
        end do
    end function pnbconv_mu

    function pnbconv_p(quants, ps, phis, method, n_terms, tolerance, normalize) result(cdfq)
        integer, intent(in) :: quants(:)
        real(dp), intent(in) :: ps(:), phis(:)
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: n_terms
        real(dp), intent(in), optional :: tolerance
        logical, intent(in), optional :: normalize
        real(dp), allocatable :: cdfq(:)
        real(dp), allocatable :: pmf(:), cdf(:)
        integer, allocatable :: counts(:)
        integer :: i, maxq

        if (size(quants) == 0) then
            allocate(cdfq(0))
            return
        end if
        if (any(quants < 0)) error stop "pnbconv_p: quants must be nonnegative"
        maxq = maxval(quants)
        allocate(counts(0:maxq), cdf(0:maxq), cdfq(size(quants)))
        counts = [(i, i=0,maxq)]
        pmf = dnbconv_p(counts, ps, phis, method, n_terms, tolerance, normalize)
        cdf(0) = pmf(1)
        do i = 1, maxq
            cdf(i) = cdf(i - 1) + pmf(i + 1)
        end do
        do i = 1, size(quants)
            cdfq(i) = cdf(quants(i))
        end do
    end function pnbconv_p

    function qnbconv_mu(probs, max_count, mus, phis, method, n_terms, tolerance, normalize) result(quants)
        real(dp), intent(in) :: probs(:)
        integer, intent(in) :: max_count
        real(dp), intent(in) :: mus(:), phis(:)
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: n_terms
        real(dp), intent(in), optional :: tolerance
        logical, intent(in), optional :: normalize
        integer, allocatable :: quants(:)
        real(dp), allocatable :: pmf(:), cdf(:)
        integer, allocatable :: counts(:)
        integer :: i, j

        if (max_count < 0) error stop "qnbconv_mu: max_count must be nonnegative"
        if (any(probs < 0.0_dp) .or. any(probs > 1.0_dp)) then
            error stop "qnbconv_mu: probabilities must lie in [0,1]"
        end if
        allocate(counts(0:max_count), cdf(0:max_count), quants(size(probs)))
        counts = [(i, i=0,max_count)]
        pmf = dnbconv_mu(counts, mus, phis, method, n_terms, tolerance, normalize)
        cdf(0) = pmf(1)
        do i = 1, max_count
            cdf(i) = cdf(i - 1) + pmf(i + 1)
        end do
        if (size(probs) > 0) then
            if (cdf(max_count) + 64.0_dp * epsilon(1.0_dp) < maxval(probs)) then
                error stop "qnbconv_mu: requested probability exceeds evaluated count range"
            end if
        end if
        do j = 1, size(probs)
            quants(j) = max_count
            do i = 0, max_count
                if (cdf(i) >= probs(j)) then
                    quants(j) = i
                    exit
                end if
            end do
        end do
    end function qnbconv_mu

    function qnbconv_p(probs, max_count, ps, phis, method, n_terms, tolerance, normalize) result(quants)
        real(dp), intent(in) :: probs(:)
        integer, intent(in) :: max_count
        real(dp), intent(in) :: ps(:), phis(:)
        character(len=*), intent(in), optional :: method
        integer, intent(in), optional :: n_terms
        real(dp), intent(in), optional :: tolerance
        logical, intent(in), optional :: normalize
        integer, allocatable :: quants(:)
        real(dp), allocatable :: pmf(:), cdf(:)
        integer, allocatable :: counts(:)
        integer :: i, j

        if (max_count < 0) error stop "qnbconv_p: max_count must be nonnegative"
        if (any(probs < 0.0_dp) .or. any(probs > 1.0_dp)) then
            error stop "qnbconv_p: probabilities must lie in [0,1]"
        end if
        allocate(counts(0:max_count), cdf(0:max_count), quants(size(probs)))
        counts = [(i, i=0,max_count)]
        pmf = dnbconv_p(counts, ps, phis, method, n_terms, tolerance, normalize)
        cdf(0) = pmf(1)
        do i = 1, max_count
            cdf(i) = cdf(i - 1) + pmf(i + 1)
        end do
        if (size(probs) > 0) then
            if (cdf(max_count) + 64.0_dp * epsilon(1.0_dp) < maxval(probs)) then
                error stop "qnbconv_p: requested probability exceeds evaluated count range"
            end if
        end if
        do j = 1, size(probs)
            quants(j) = max_count
            do i = 0, max_count
                if (cdf(i) >= probs(j)) then
                    quants(j) = i
                    exit
                end if
            end do
        end do
    end function qnbconv_p

    function rnbconv_mu(n_samp, mus, phis) result(samples)
        integer, intent(in) :: n_samp
        real(dp), intent(in) :: mus(:), phis(:)
        integer, allocatable :: samples(:)
        integer :: i, j

        call validate_mu(mus, phis)
        if (n_samp < 0) error stop "rnbconv_mu: n_samp must be nonnegative"
        allocate(samples(n_samp))
        samples = 0
        do i = 1, n_samp
            do j = 1, size(mus)
                samples(i) = samples(i) + negbin_rng_mu(mus(j), phis(j))
            end do
        end do
    end function rnbconv_mu

    function rnbconv_p(n_samp, ps, phis) result(samples)
        integer, intent(in) :: n_samp
        real(dp), intent(in) :: ps(:), phis(:)
        integer, allocatable :: samples(:)
        real(dp), allocatable :: mus(:)

        call validate_p(ps, phis)
        allocate(mus(size(ps)))
        mus = phis * (1.0_dp - ps) / ps
        samples = rnbconv_mu(n_samp, mus, phis)
    end function rnbconv_p

    function nbconv_params_mu(mus, phis) result(params)
        real(dp), intent(in) :: mus(:), phis(:)
        type(nbconv_summary) :: params
        real(dp), allocatable :: ps(:)

        call validate_mu(mus, phis)
        allocate(ps(size(mus)))
        ps = phis / (phis + mus)
        params = calculate_params(mus, phis, ps)
    end function nbconv_params_mu

    function nbconv_params_p(ps, phis) result(params)
        real(dp), intent(in) :: ps(:), phis(:)
        type(nbconv_summary) :: params
        real(dp), allocatable :: mus(:)

        call validate_p(ps, phis)
        allocate(mus(size(ps)))
        mus = phis * (1.0_dp - ps) / ps
        params = calculate_params(mus, phis, ps)
    end function nbconv_params_p

    function calculate_params(mus, phis, ps) result(params)
        real(dp), intent(in) :: mus(:), phis(:), ps(:)
        type(nbconv_summary) :: params
        real(dp) :: k1, k2, k3, k4, pmax, qmax
        integer :: i

        k1 = sum(mus)
        k2 = sum(mus + mus * mus / phis)
        k3 = sum((2.0_dp * mus + phis) * (mus + phis) * mus / (phis * phis))
        k4 = sum((6.0_dp * mus * mus + 6.0_dp * mus * phis + phis * phis) &
             * (mus + phis) * mus / (phis * phis * phis))
        params%mean = k1
        params%variance = k2
        if (k2 > 0.0_dp) then
            params%skewness = k3 / k2**1.5_dp
            params%excess_kurtosis = k4 / (k2 * k2)
        else
            params%skewness = ieee_value(0.0_dp, ieee_quiet_nan)
            params%excess_kurtosis = ieee_value(0.0_dp, ieee_quiet_nan)
        end if

        pmax = 0.0_dp
        do i = 1, size(ps)
            if (mus(i) > 0.0_dp) pmax = max(pmax, ps(i))
        end do
        if (pmax > 0.0_dp) then
            qmax = 1.0_dp - pmax
            params%k_mean = params%mean * pmax / qmax - sum(phis, mask=mus > 0.0_dp)
        else
            params%k_mean = 0.0_dp
        end if
    end function calculate_params

    function get_parameterization(parameterization) result(par)
        character(len=*), intent(in), optional :: parameterization
        character(len=8) :: par

        par = "mu"
        if (present(parameterization)) par = adjustl(parameterization)
        select case (trim(par))
        case ("mu", "mus", "mean", "means")
            par = "mu"
        case ("p", "ps", "prob", "probability")
            par = "p"
        case default
            error stop "nbconv: parameterization must be mu or p"
        end select
    end function get_parameterization

    function get_method(method) result(meth)
        character(len=*), intent(in), optional :: method
        character(len=16) :: meth

        meth = "exact"
        if (present(method)) meth = adjustl(method)
        select case (trim(meth))
        case ("exact", "moments", "saddlepoint")
        case default
            error stop "nbconv: method must be exact, moments, or saddlepoint"
        end select
    end function get_method

    subroutine validate_mu(mus, phis)
        real(dp), intent(in) :: mus(:), phis(:)

        if (size(mus) /= size(phis)) error stop "nbconv: mus and phis must have equal length"
        if (size(mus) < 1) error stop "nbconv: parameter arrays must be nonempty"
        if (any(mus < 0.0_dp)) error stop "nbconv: mus must be nonnegative"
        if (any(phis <= 0.0_dp)) error stop "nbconv: phis must be positive"
    end subroutine validate_mu

    subroutine validate_p(ps, phis)
        real(dp), intent(in) :: ps(:), phis(:)

        if (size(ps) /= size(phis)) error stop "nbconv: ps and phis must have equal length"
        if (size(ps) < 1) error stop "nbconv: parameter arrays must be nonempty"
        if (any(ps <= 0.0_dp) .or. any(ps > 1.0_dp)) error stop "nbconv: ps must satisfy 0 < p <= 1"
        if (any(phis <= 0.0_dp)) error stop "nbconv: phis must be positive"
    end subroutine validate_p

end module nbconv_api

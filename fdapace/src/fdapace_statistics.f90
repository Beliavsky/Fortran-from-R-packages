module fdapace_statistics
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use fdapace_kinds, only : dp
    use fdapace_math, only : fill_missing_linear, linear_interp, mean_value, sample_sd, sort_real
    use fdapace_smoothing, only : lwls1d, trapz_rcpp
    use fdapace_types, only : bwnn_result, dyn_test_result, fccor_result, real_vector
    implicit none
    private

    public :: bw_nn
    public :: dyn_corr
    public :: dyn_test
    public :: fc_cor

contains

    function dyn_corr(x, y, t) result(z)
        real(dp), intent(in) :: x(:,:) !! First dense functional sample with subjects in rows; NaNs are linearly imputed.
        real(dp), intent(in) :: y(:,:) !! Second paired dense functional sample with the same shape as x; NaNs are linearly imputed.
        real(dp), intent(in) :: t(:) !! Common increasing observation grid corresponding to sample columns.
        real(dp), allocatable :: z(:)
        real(dp), allocatable :: xfill(:,:)
        real(dp), allocatable :: yfill(:,:)
        real(dp), allocatable :: temp1_x(:,:)
        real(dp), allocatable :: temp1_y(:,:)
        real(dp), allocatable :: mx(:)
        real(dp), allocatable :: my(:)
        real(dp), allocatable :: rx(:)
        real(dp), allocatable :: ry(:)
        real(dp) :: aver_x
        real(dp) :: aver_y
        real(dp) :: denom_x
        real(dp) :: denom_y
        real(dp) :: span
        integer :: i
        integer :: info

        if (any(shape(x) /= shape(y))) error stop "dyn_corr: x and y shapes differ"
        if (size(x, 2) /= size(t)) error stop "dyn_corr: sample columns do not match t"
        if (size(t) < 2) error stop "dyn_corr: at least two time points are required"
        span = t(size(t)) - t(1)
        if (span <= 0.0_dp) error stop "dyn_corr: t must span a positive interval"

        allocate(xfill(size(x, 1), size(x, 2)), yfill(size(y, 1), size(y, 2)))
        do i = 1, size(x, 1)
            call fill_missing_linear(t, x(i, :), xfill(i, :), info)
            if (info /= 0) error stop "dyn_corr: a row of x cannot be imputed"
            call fill_missing_linear(t, y(i, :), yfill(i, :), info)
            if (info /= 0) error stop "dyn_corr: a row of y cannot be imputed"
        end do

        allocate(temp1_x(size(x, 1), size(x, 2)), temp1_y(size(y, 1), size(y, 2)))
        do i = 1, size(x, 1)
            aver_x = trapz_rcpp(t, xfill(i, :)) / span
            aver_y = trapz_rcpp(t, yfill(i, :)) / span
            temp1_x(i, :) = xfill(i, :) - aver_x
            temp1_y(i, :) = yfill(i, :) - aver_y
        end do
        allocate(mx(size(t)), my(size(t)))
        mx = sum(temp1_x, dim=1) / real(size(x, 1), dp)
        my = sum(temp1_y, dim=1) / real(size(y, 1), dp)
        allocate(rx(size(t)), ry(size(t)), z(size(x, 1)))
        do i = 1, size(x, 1)
            rx = temp1_x(i, :) - mx
            ry = temp1_y(i, :) - my
            denom_x = sqrt(trapz_rcpp(t, rx * rx) / span)
            denom_y = sqrt(trapz_rcpp(t, ry * ry) / span)
            if (denom_x <= 0.0_dp .or. denom_y <= 0.0_dp) then
                z(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            else
                z(i) = trapz_rcpp(t, (rx / denom_x) * (ry / denom_y)) / span
            end if
        end do
    end function dyn_corr

    subroutine dyn_test(x1, y1, t1, result, b, x2, y2, t2)
        real(dp), intent(in) :: x1(:,:) !! First process of sample one, with subjects in rows and optional NaNs.
        real(dp), intent(in) :: y1(:,:) !! Second paired process of sample one, with the same shape as x1.
        real(dp), intent(in) :: t1(:) !! Common observation grid for sample one.
        type(dyn_test_result), intent(out) :: result !! Studentized observed statistic and two-sided bootstrap p-value.
        integer, intent(in), optional :: b !! Number of paired bootstrap replicates; defaults to 1000 and must be positive.
        real(dp), intent(in), optional :: x2(:,:) !! Optional first process of sample two for the paired two-sample test.
        real(dp), intent(in), optional :: y2(:,:) !! Optional second process of sample two, required when x2 is present.
        real(dp), intent(in), optional :: t2(:) !! Optional observation grid for sample two, required when x2 is present.
        real(dp), allocatable :: boot1(:)
        real(dp), allocatable :: boot2(:)
        real(dp), allocatable :: boot_diff(:)
        real(dp), allocatable :: bx1(:,:)
        real(dp), allocatable :: by1(:,:)
        real(dp), allocatable :: bx2(:,:)
        real(dp), allocatable :: by2(:,:)
        real(dp), allocatable :: d1(:)
        real(dp), allocatable :: d2(:)
        real(dp), allocatable :: db1(:)
        real(dp), allocatable :: db2(:)
        integer, allocatable :: idx(:)
        integer :: bb
        integer :: ib
        integer :: i
        integer :: n
        real(dp) :: u
        real(dp) :: sdv

        bb = 1000
        if (present(b)) bb = b
        if (bb < 1) error stop "dyn_test: b must be positive"
        n = size(x1, 1)
        if (n < 2) error stop "dyn_test: at least two subjects are required"
        if (any(shape(x1) /= shape(y1))) error stop "dyn_test: x1 and y1 shapes differ"
        if (size(x1, 2) /= size(t1)) error stop "dyn_test: x1 columns do not match t1"

        d1 = dyn_corr(x1, y1, t1)
        allocate(boot1(bb), idx(n), bx1(n, size(t1)), by1(n, size(t1)))
        sdv = sample_sd(d1)
        if (sdv <= 0.0_dp .or. ieee_is_nan(sdv)) then
            result%stats = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            result%stats = mean_value(d1) * sqrt(real(n, dp)) / sdv
        end if

        if (.not. present(x2)) then
            do ib = 1, bb
                do i = 1, n
                    call random_number(u)
                    idx(i) = min(n, int(u * real(n, dp)) + 1)
                    bx1(i, :) = x1(idx(i), :)
                    by1(i, :) = y1(idx(i), :)
                end do
                db1 = dyn_corr(bx1, by1, t1)
                sdv = sample_sd(db1)
                if (sdv <= 0.0_dp .or. ieee_is_nan(sdv)) then
                    boot1(ib) = ieee_value(0.0_dp, ieee_quiet_nan)
                else
                    boot1(ib) = mean_value(db1 - d1) * sqrt(real(n, dp)) / sdv
                end if
            end do
            result%pval = bootstrap_pvalue(boot1, result%stats)
            return
        end if

        if (.not. present(y2) .or. .not. present(t2)) error stop "dyn_test: x2, y2, and t2 must be supplied together"
        if (size(x2, 1) /= n .or. any(shape(x2) /= shape(y2))) error stop "dyn_test: paired sample-two dimensions are invalid"
        if (size(x2, 2) /= size(t2)) error stop "dyn_test: x2 columns do not match t2"
        d2 = dyn_corr(x2, y2, t2)
        sdv = sample_sd(d2 - d1)
        if (sdv <= 0.0_dp .or. ieee_is_nan(sdv)) then
            result%stats = ieee_value(0.0_dp, ieee_quiet_nan)
        else
            result%stats = mean_value(d2 - d1) * sqrt(real(n, dp)) / sdv
        end if
        allocate(boot2(bb), boot_diff(bb), bx2(n, size(t2)), by2(n, size(t2)))
        do ib = 1, bb
            do i = 1, n
                call random_number(u)
                idx(i) = min(n, int(u * real(n, dp)) + 1)
                bx1(i, :) = x1(idx(i), :)
                by1(i, :) = y1(idx(i), :)
                bx2(i, :) = x2(idx(i), :)
                by2(i, :) = y2(idx(i), :)
            end do
            db1 = dyn_corr(bx1, by1, t1)
            db2 = dyn_corr(bx2, by2, t2)
            sdv = sample_sd(db1)
            if (sdv > 0.0_dp .and. .not. ieee_is_nan(sdv)) then
                boot1(ib) = mean_value(db1 - d1) * sqrt(real(n, dp)) / sdv
            else
                boot1(ib) = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
            sdv = sample_sd(db2)
            if (sdv > 0.0_dp .and. .not. ieee_is_nan(sdv)) then
                boot2(ib) = mean_value(db2 - d2) * sqrt(real(n, dp)) / sdv
            else
                boot2(ib) = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
            sdv = sample_sd(db2 - db1)
            if (sdv > 0.0_dp .and. .not. ieee_is_nan(sdv)) then
                boot_diff(ib) = mean_value(db2 - db1) * sqrt(real(n, dp)) / sdv
            else
                boot_diff(ib) = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
        end do
        result%pval = bootstrap_pvalue(boot_diff, result%stats)
    end subroutine dyn_test

    pure real(dp) function bootstrap_pvalue(boot, observed) result(value)
        real(dp), intent(in) :: boot(:) !! Bootstrap test statistics; NaN replicates are ignored in the numerator but retain the upstream denominator.
        real(dp), intent(in) :: observed !! Observed studentized statistic whose absolute magnitude defines the two-sided tail.
        integer :: i
        integer :: count_extreme

        if (size(boot) == 0 .or. ieee_is_nan(observed)) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        count_extreme = 0
        do i = 1, size(boot)
            if (ieee_is_nan(boot(i))) cycle
            if (boot(i) > abs(observed) .or. boot(i) < -abs(observed)) count_extreme = count_extreme + 1
        end do
        value = real(count_extreme, dp) / real(size(boot), dp)
    end function bootstrap_pvalue

    function fc_cor(x, y, lt, bw, kernel_type, tout) result(res)
        type(real_vector), intent(in) :: x(:) !! Per-subject observations of the first process.
        type(real_vector), intent(in) :: y(:) !! Per-subject observations of the second process, aligned with x and lt.
        type(real_vector), intent(in) :: lt(:) !! Per-subject observation times shared by x and y.
        real(dp), intent(in) :: bw(:) !! One bandwidth recycled five times, or five bandwidths for both means, both variances, and covariance.
        character(len=*), intent(in), optional :: kernel_type !! Smoothing kernel; defaults to epan.
        real(dp), intent(in), optional :: tout(:) !! Output grid; defaults to sorted unique observed time points.
        type(fccor_result) :: res
        real(dp), allocatable :: xvec(:)
        real(dp), allocatable :: yvec(:)
        real(dp), allocatable :: tvec(:)
        real(dp), allocatable :: tall(:)
        real(dp), allocatable :: mux(:)
        real(dp), allocatable :: muy(:)
        real(dp), allocatable :: xcent(:)
        real(dp), allocatable :: ycent(:)
        real(dp), allocatable :: varx(:)
        real(dp), allocatable :: vary(:)
        real(dp), allocatable :: covxy(:)
        real(dp), allocatable :: ones(:)
        character(len=16) :: kern
        integer :: i
        integer :: nall
        integer :: pos

        if (size(x) /= size(y) .or. size(x) /= size(lt)) error stop "fc_cor: list lengths differ"
        nall = 0
        do i = 1, size(lt)
            if (size(x(i)%v) /= size(lt(i)%v) .or. size(y(i)%v) /= size(lt(i)%v)) then
                error stop "fc_cor: within-subject lengths differ"
            end if
            nall = nall + size(lt(i)%v)
        end do
        allocate(xvec(nall), yvec(nall), tvec(nall))
        pos = 0
        do i = 1, size(lt)
            xvec(pos + 1:pos + size(lt(i)%v)) = x(i)%v
            yvec(pos + 1:pos + size(lt(i)%v)) = y(i)%v
            tvec(pos + 1:pos + size(lt(i)%v)) = lt(i)%v
            pos = pos + size(lt(i)%v)
        end do
        call sort_triplets(tvec, xvec, yvec)
        tall = unique_sorted(tvec)
        if (size(bw) == 1) then
            res%bw = bw(1)
        else if (size(bw) == 5) then
            res%bw = bw
        else
            error stop "fc_cor: bw must have length one or five"
        end if
        kern = "epan"
        if (present(kernel_type)) kern = trim(kernel_type)
        if (present(tout)) then
            res%tout = tout
        else
            res%tout = tall
        end if
        allocate(ones(nall))
        ones = 1.0_dp
        mux = lwls1d(res%bw(1), kern, tvec, xvec, tall, ones, 1, 0)
        muy = lwls1d(res%bw(2), kern, tvec, yvec, tall, ones, 1, 0)
        allocate(xcent(nall), ycent(nall))
        do i = 1, nall
            xcent(i) = xvec(i) - linear_interp(tall, mux, tvec(i))
            ycent(i) = yvec(i) - linear_interp(tall, muy, tvec(i))
        end do
        varx = lwls1d(res%bw(3), kern, tvec, xcent**2, res%tout, ones, 1, 0)
        vary = lwls1d(res%bw(4), kern, tvec, ycent**2, res%tout, ones, 1, 0)
        covxy = lwls1d(res%bw(5), kern, tvec, xcent * ycent, res%tout, ones, 1, 0)
        allocate(res%corr(size(res%tout)))
        do i = 1, size(res%tout)
            if (varx(i) <= 0.0_dp .or. vary(i) <= 0.0_dp .or. ieee_is_nan(varx(i)) .or. ieee_is_nan(vary(i))) then
                res%corr(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            else
                res%corr(i) = covxy(i) / sqrt(varx(i) * vary(i))
            end if
        end do
    end function fc_cor

    function bw_nn(lt, k, only_mean, only_cov) result(res)
        type(real_vector), intent(in) :: lt(:) !! Per-subject observation-time vectors used to construct unique one- and two-dimensional supports.
        integer, intent(in), optional :: k !! Number of unique neighbors required; defaults to three and must be positive.
        logical, intent(in), optional :: only_mean !! If true, compute only the mean bandwidth.
        logical, intent(in), optional :: only_cov !! If true, compute only the covariance bandwidth.
        type(bwnn_result) :: res
        real(dp), allocatable :: pair1(:)
        real(dp), allocatable :: pair2(:)
        real(dp), allocatable :: up1(:)
        real(dp), allocatable :: up2(:)
        real(dp), allocatable :: distances(:)
        real(dp), allocatable :: times(:)
        real(dp) :: kth
        integer :: kk
        integer :: i
        integer :: j
        integer :: q
        integer :: npairs
        integer :: nuniq
        integer :: pos
        logical :: mean_only
        logical :: cov_only

        kk = 3
        if (present(k)) kk = k
        if (kk < 1) error stop "bw_nn: k must be at least one"
        mean_only = .false.
        if (present(only_mean)) mean_only = only_mean
        cov_only = .false.
        if (present(only_cov)) cov_only = only_cov
        if (mean_only .and. cov_only) error stop "bw_nn: both only_mean and only_cov are true"

        npairs = 0
        do i = 1, size(lt)
            npairs = npairs + size(lt(i)%v) * size(lt(i)%v)
        end do
        allocate(pair1(npairs), pair2(npairs))
        pos = 0
        do i = 1, size(lt)
            do j = 1, size(lt(i)%v)
                do q = 1, size(lt(i)%v)
                    pos = pos + 1
                    pair1(pos) = lt(i)%v(j)
                    pair2(pos) = lt(i)%v(q)
                end do
            end do
        end do
        call unique_pairs(pair1, pair2, up1, up2)
        nuniq = size(up1)

        if (.not. mean_only) then
            if (nuniq < kk + 1) error stop "bw_nn: too few unique covariance pairs for requested k"
            allocate(distances(nuniq))
            res%cov_bw = 0.0_dp
            do i = 1, nuniq
                do j = 1, nuniq
                    distances(j) = max(abs(up1(j) - up1(i)), abs(up2(j) - up2(i)))
                end do
                call sort_real(distances)
                kth = distances(kk + 1)
                res%cov_bw = max(res%cov_bw, kth)
            end do
            res%has_cov = .true.
        end if

        if (.not. cov_only) then
            times = unique_sorted(up1)
            if (size(times) < kk + 1) error stop "bw_nn: too few unique mean-grid points for requested k"
            res%mu_bw = 0.0_dp
            do i = 1, size(times) - kk
                res%mu_bw = max(res%mu_bw, times(i + kk) - times(i))
            end do
            res%has_mean = .true.
        end if
    end function bw_nn

    subroutine sort_triplets(key, a, b)
        real(dp), intent(inout) :: key(:) !! Sort key reordered into ascending order.
        real(dp), intent(inout) :: a(:) !! First companion vector permuted identically with key.
        real(dp), intent(inout) :: b(:) !! Second companion vector permuted identically with key.
        integer :: i
        integer :: j
        real(dp) :: ka
        real(dp) :: aa
        real(dp) :: bb

        do i = 2, size(key)
            ka = key(i)
            aa = a(i)
            bb = b(i)
            j = i - 1
            do while (j >= 1)
                if (key(j) <= ka) exit
                key(j + 1) = key(j)
                a(j + 1) = a(j)
                b(j + 1) = b(j)
                j = j - 1
            end do
            key(j + 1) = ka
            a(j + 1) = aa
            b(j + 1) = bb
        end do
    end subroutine sort_triplets

    function unique_sorted(x) result(values)
        real(dp), intent(in) :: x(:) !! Numeric values from which sorted exact duplicates are removed.
        real(dp), allocatable :: values(:)
        real(dp), allocatable :: work(:)
        integer :: i
        integer :: nuniq

        if (size(x) == 0) then
            allocate(values(0))
            return
        end if
        allocate(work(size(x)))
        work = x
        call sort_real(work)
        nuniq = 1
        do i = 2, size(work)
            if (work(i) /= work(nuniq)) then
                nuniq = nuniq + 1
                work(nuniq) = work(i)
            end if
        end do
        allocate(values(nuniq))
        values = work(1:nuniq)
    end function unique_sorted

    subroutine unique_pairs(a, b, ua, ub)
        real(dp), intent(in) :: a(:) !! First coordinates of candidate pairs.
        real(dp), intent(in) :: b(:) !! Second coordinates of candidate pairs.
        real(dp), allocatable, intent(out) :: ua(:) !! First coordinates of exact unique pairs in first-occurrence order.
        real(dp), allocatable, intent(out) :: ub(:) !! Second coordinates of exact unique pairs in first-occurrence order.
        real(dp), allocatable :: ta(:)
        real(dp), allocatable :: tb(:)
        integer :: i
        integer :: j
        integer :: nuniq
        logical :: seen

        if (size(a) /= size(b)) error stop "unique_pairs: coordinate lengths differ"
        allocate(ta(size(a)), tb(size(b)))
        nuniq = 0
        do i = 1, size(a)
            seen = .false.
            do j = 1, nuniq
                if (a(i) == ta(j) .and. b(i) == tb(j)) then
                    seen = .true.
                    exit
                end if
            end do
            if (.not. seen) then
                nuniq = nuniq + 1
                ta(nuniq) = a(i)
                tb(nuniq) = b(i)
            end if
        end do
        allocate(ua(nuniq), ub(nuniq))
        ua = ta(1:nuniq)
        ub = tb(1:nuniq)
    end subroutine unique_pairs

end module fdapace_statistics

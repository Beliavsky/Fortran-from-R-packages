module tclust_metrics
    use tclust_kinds, only : dp
    use tclust_types, only : rand_index_result, fm_index_result
    implicit none
    private

    public :: rand_index_labels
    public :: rand_index_table
    public :: fowlkes_mallows_labels
    public :: fowlkes_mallows_table

contains

    subroutine rand_index_labels(c1, c2, result, noisecluster)
        integer, intent(in) :: c1(:) !! Labels for the first partition; length n.
        integer, intent(in) :: c2(:) !! Labels for the second partition; must have the same length as c1.
        type(rand_index_result), intent(out) :: result !! Adjusted Rand, Rand, Mirkin, and Hubert indices.
        integer, intent(in), optional :: noisecluster !! Optional label excluded whenever it occurs in either partition.
        integer, allocatable :: tab(:, :)

        call contingency_from_labels(c1, c2, tab, noisecluster)
        call rand_index_table(tab, result)
    end subroutine rand_index_labels

    subroutine rand_index_table(table, result)
        integer, intent(in) :: table(:, :) !! Nonnegative contingency table comparing two partitions.
        type(rand_index_result), intent(out) :: result !! Adjusted Rand, Rand, Mirkin, and Hubert indices.
        real(dp), allocatable :: rows(:)
        real(dp), allocatable :: cols(:)
        real(dp) :: n
        real(dp) :: nis
        real(dp) :: njs
        real(dp) :: totcomp
        real(dp) :: t2
        real(dp) :: t3
        real(dp) :: nc
        real(dp) :: agree
        real(dp) :: disagree

        if (any(table < 0)) error stop 'rand_index_table: contingency counts must be nonnegative'
        n = real(sum(table), dp)
        if (n < 2.0_dp) error stop 'rand_index_table: at least two classified observations are required'
        allocate(rows(size(table, 1)), cols(size(table, 2)))
        rows = real(sum(table, dim=2), dp)
        cols = real(sum(table, dim=1), dp)
        nis = sum(rows**2)
        njs = sum(cols**2)
        totcomp = n * (n - 1.0_dp) / 2.0_dp
        t2 = sum(real(table, dp)**2)
        t3 = 0.5_dp * (nis + njs)
        nc = (n * (n**2 + 1.0_dp) - (n + 1.0_dp) * nis - (n + 1.0_dp) * njs + &
              2.0_dp * nis * njs / n) / (2.0_dp * (n - 1.0_dp))
        agree = totcomp + t2 - t3
        disagree = -t2 + t3
        if (abs(totcomp - nc) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(totcomp))) then
            result%adjusted_rand = 0.0_dp
        else
            result%adjusted_rand = (agree - nc) / (totcomp - nc)
        end if
        result%rand = agree / totcomp
        result%mirkin = disagree / totcomp
        result%hubert = (agree - disagree) / totcomp
    end subroutine rand_index_table

    subroutine fowlkes_mallows_labels(c1, c2, result, noisecluster)
        integer, intent(in) :: c1(:) !! Labels for the first partition; length n.
        integer, intent(in) :: c2(:) !! Labels for the second partition; must have the same length as c1.
        type(fm_index_result), intent(out) :: result !! Adjusted/raw Fowlkes-Mallows index, expectation, and null variance.
        integer, intent(in), optional :: noisecluster !! Optional label excluded whenever it occurs in either partition.
        integer, allocatable :: tab(:, :)

        call contingency_from_labels(c1, c2, tab, noisecluster)
        call fowlkes_mallows_table(tab, result)
    end subroutine fowlkes_mallows_labels

    subroutine fowlkes_mallows_table(table, result)
        integer, intent(in) :: table(:, :) !! Nonnegative contingency table comparing two partitions.
        type(fm_index_result), intent(out) :: result !! Adjusted/raw Fowlkes-Mallows index, expectation, and null variance.
        real(dp), allocatable :: rows(:)
        real(dp), allocatable :: cols(:)
        real(dp) :: n
        real(dp) :: tk
        real(dp) :: pk
        real(dp) :: qk
        real(dp) :: pk2
        real(dp) :: qk2
        real(dp) :: denom

        if (any(table < 0)) error stop 'fowlkes_mallows_table: contingency counts must be nonnegative'
        n = real(sum(table), dp)
        if (n < 4.0_dp) error stop 'fowlkes_mallows_table: at least four observations are required for the variance formula'
        allocate(rows(size(table, 1)), cols(size(table, 2)))
        rows = real(sum(table, dim=2), dp)
        cols = real(sum(table, dim=1), dp)
        tk = sum(real(table, dp)**2) - n
        pk = sum(rows**2) - n
        qk = sum(cols**2) - n
        if (pk <= 0.0_dp .or. qk <= 0.0_dp) error stop 'fowlkes_mallows_table: degenerate partition'
        denom = sqrt(pk * qk)
        result%raw = tk / denom
        result%expected = denom / (n * (n - 1.0_dp))
        pk2 = sum(rows * (rows - 1.0_dp) * (rows - 2.0_dp))
        qk2 = sum(cols * (cols - 1.0_dp) * (cols - 2.0_dp))
        result%variance = 2.0_dp / (n * (n - 1.0_dp)) + &
            4.0_dp * pk2 * qk2 / ((n * (n - 1.0_dp) * (n - 2.0_dp)) * pk * qk) + &
            (pk - 2.0_dp - 4.0_dp * pk2 / pk) * (qk - 2.0_dp - 4.0_dp * qk2 / qk) / &
            (n * (n - 1.0_dp) * (n - 2.0_dp) * (n - 3.0_dp)) - &
            pk * qk / (n**2 * (n - 1.0_dp)**2)
        if (abs(1.0_dp - result%expected) <= 100.0_dp * epsilon(1.0_dp)) then
            result%adjusted = 0.0_dp
        else
            result%adjusted = (result%raw - result%expected) / (1.0_dp - result%expected)
        end if
    end subroutine fowlkes_mallows_table

    subroutine contingency_from_labels(c1, c2, table, noisecluster)
        integer, intent(in) :: c1(:) !! Labels for the first partition.
        integer, intent(in) :: c2(:) !! Labels for the second partition.
        integer, allocatable, intent(out) :: table(:, :) !! Dense contingency table over sorted unique retained labels.
        integer, intent(in), optional :: noisecluster !! Optional label excluded from either partition.
        integer, allocatable :: u1(:)
        integer, allocatable :: u2(:)
        integer, allocatable :: w1(:)
        integer, allocatable :: w2(:)
        logical, allocatable :: keep(:)
        integer :: n
        integer :: i
        integer :: m
        integer :: r
        integer :: c

        if (size(c1) /= size(c2)) error stop 'contingency_from_labels: label lengths differ'
        n = size(c1)
        allocate(keep(n))
        keep = .true.
        if (present(noisecluster)) keep = (c1 /= noisecluster .and. c2 /= noisecluster)
        m = count(keep)
        if (m < 1) error stop 'contingency_from_labels: no retained observations'
        allocate(w1(m), w2(m))
        m = 0
        do i = 1, n
            if (keep(i)) then
                m = m + 1
                w1(m) = c1(i)
                w2(m) = c2(i)
            end if
        end do
        call unique_sorted(w1, u1)
        call unique_sorted(w2, u2)
        allocate(table(size(u1), size(u2)))
        table = 0
        do i = 1, size(w1)
            r = find_label(u1, w1(i))
            c = find_label(u2, w2(i))
            table(r, c) = table(r, c) + 1
        end do
    end subroutine contingency_from_labels

    subroutine unique_sorted(x, u)
        integer, intent(in) :: x(:) !! Integer labels whose sorted unique values are requested.
        integer, allocatable, intent(out) :: u(:) !! Sorted vector of unique integer labels.
        integer, allocatable :: tmp(:)
        integer :: i
        integer :: j
        integer :: t
        integer :: m

        allocate(tmp(size(x)))
        tmp = x
        do i = 1, size(tmp) - 1
            do j = i + 1, size(tmp)
                if (tmp(j) < tmp(i)) then
                    t = tmp(i)
                    tmp(i) = tmp(j)
                    tmp(j) = t
                end if
            end do
        end do
        m = 1
        do i = 2, size(tmp)
            if (tmp(i) /= tmp(m)) then
                m = m + 1
                tmp(m) = tmp(i)
            end if
        end do
        allocate(u(m))
        u = tmp(:m)
    end subroutine unique_sorted

    pure integer function find_label(labels, value) result(pos)
        integer, intent(in) :: labels(:) !! Sorted unique label vector to search.
        integer, intent(in) :: value !! Label value whose one-based position is requested.
        integer :: i

        pos = 0
        do i = 1, size(labels)
            if (labels(i) == value) then
                pos = i
                return
            end if
        end do
    end function find_label

end module tclust_metrics

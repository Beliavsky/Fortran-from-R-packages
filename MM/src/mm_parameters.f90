! Computational translation of the R package MM 1.7-0.
! Upstream license: GPL-2. This translation is GPL-2.
module mm_parameters
    use mm_kinds, only : dp
    use mm_types, only : paras_type, mb_type
    implicit none
    private

    public :: paras, paras_from_values, paras_from_p_theta
    public :: paras_dimension, p, theta, set_p, set_theta

    interface set_theta
        module procedure set_theta_matrix
        module procedure set_theta_scalar
    end interface set_theta
    public :: valid_paras, make_mb, valid_mb

contains

    function paras(k) result(out)
        integer, intent(in) :: k
        type(paras_type) :: out
        integer :: nv

        if (k < 2) error stop "paras: k must be at least 2"
        nv = k - 1 + k * (k - 1) / 2
        allocate(out%vals(nv), out%pnames(k))
        out%vals = 1.0_dp
        out%vals(1:k - 1) = 1.0_dp / real(k, dp)
        out%pnames = ''
    end function paras

    function paras_from_values(vals, names) result(out)
        real(dp), intent(in) :: vals(:)
        character(len=*), intent(in), optional :: names(:)
        type(paras_type) :: out
        integer :: k

        k = dimension_from_nvals(size(vals))
        if (k < 2) error stop "paras_from_values: invalid value-vector length"
        allocate(out%vals(size(vals)), out%pnames(k))
        out%vals = vals
        out%pnames = ''
        if (present(names)) then
            if (size(names) /= k) error stop "paras_from_values: names have wrong length"
            out%pnames = names
        end if
    end function paras_from_values

    function paras_from_p_theta(prob, th, names) result(out)
        real(dp), intent(in) :: prob(:), th(:,:)
        character(len=*), intent(in), optional :: names(:)
        type(paras_type) :: out
        integer :: k, i, j, q

        k = size(th, 1)
        if (k < 2 .or. size(th, 2) /= k) then
            error stop "paras_from_p_theta: theta must be square"
        end if
        if (size(prob) /= k .and. size(prob) /= k - 1) then
            error stop "paras_from_p_theta: probability vector has wrong length"
        end if
        out = paras(k)
        out%vals(1:k - 1) = prob(1:k - 1)
        q = k
        do j = 2, k
            do i = 1, j - 1
                out%vals(q) = th(i, j)
                q = q + 1
            end do
        end do
        if (present(names)) then
            if (size(names) /= k) error stop "paras_from_p_theta: names have wrong length"
            out%pnames = names
        end if
    end function paras_from_p_theta

    integer function paras_dimension(x) result(k)
        type(paras_type), intent(in) :: x
        if (.not. allocated(x%vals)) then
            k = 0
        else
            k = dimension_from_nvals(size(x%vals))
        end if
    end function paras_dimension

    function p(x) result(prob)
        type(paras_type), intent(in) :: x
        real(dp), allocatable :: prob(:)
        integer :: k

        k = paras_dimension(x)
        if (k < 2) error stop "p: invalid paras object"
        allocate(prob(k))
        prob(1:k - 1) = x%vals(1:k - 1)
        prob(k) = 1.0_dp - sum(prob(1:k - 1))
    end function p

    function theta(x) result(th)
        type(paras_type), intent(in) :: x
        real(dp), allocatable :: th(:,:)
        integer :: k, i, j, q

        k = paras_dimension(x)
        if (k < 2) error stop "theta: invalid paras object"
        allocate(th(k, k))
        th = 1.0_dp
        q = k
        do j = 2, k
            do i = 1, j - 1
                th(i, j) = x%vals(q)
                q = q + 1
            end do
        end do
    end function theta

    subroutine set_p(x, prob)
        type(paras_type), intent(inout) :: x
        real(dp), intent(in) :: prob(:)
        integer :: k

        k = paras_dimension(x)
        if (size(prob) /= k .and. size(prob) /= k - 1) then
            error stop "set_p: wrong probability-vector length"
        end if
        x%vals(1:k - 1) = prob(1:k - 1)
    end subroutine set_p

    subroutine set_theta_matrix(x, th)
        type(paras_type), intent(inout) :: x
        real(dp), intent(in) :: th(:,:)
        integer :: k, i, j, q

        k = paras_dimension(x)
        if (size(th, 1) /= k .or. size(th, 2) /= k) error stop "set_theta: wrong matrix size"
        q = k
        do j = 2, k
            do i = 1, j - 1
                x%vals(q) = th(i, j)
                q = q + 1
            end do
        end do
    end subroutine set_theta_matrix

    subroutine set_theta_scalar(x, value)
        type(paras_type), intent(inout) :: x
        real(dp), intent(in) :: value
        integer :: k

        k = paras_dimension(x)
        x%vals(k:) = value
    end subroutine set_theta_scalar

    logical function valid_paras(x) result(ok)
        type(paras_type), intent(in) :: x
        real(dp), allocatable :: prob(:), th(:,:)
        integer :: k

        k = paras_dimension(x)
        if (k < 2) then
            ok = .false.
            return
        end if
        prob = p(x)
        th = theta(x)
        ok = all(prob >= 0.0_dp) .and. all(th > 0.0_dp)
    end function valid_paras

    function make_mb(dep, m, names) result(out)
        integer, intent(in) :: dep(:,:), m(:)
        character(len=*), intent(in), optional :: names(:)
        type(mb_type) :: out
        integer :: k

        k = size(dep, 2)
        if (size(m) /= k) error stop "make_mb: m has wrong length"
        allocate(out%counts(size(dep, 1), k), out%m(k), out%pnames(k))
        out%counts = dep
        out%m = m
        out%pnames = ''
        if (present(names)) then
            if (size(names) /= k) error stop "make_mb: names have wrong length"
            out%pnames = names
        end if
        if (.not. valid_mb(out)) error stop "make_mb: invalid counts or m"
    end function make_mb

    logical function valid_mb(x) result(ok)
        type(mb_type), intent(in) :: x
        integer :: j

        ok = allocated(x%counts) .and. allocated(x%m)
        if (.not. ok) return
        if (size(x%counts, 2) /= size(x%m)) then
            ok = .false.
            return
        end if
        if (any(x%m < 0) .or. any(x%counts < 0)) then
            ok = .false.
            return
        end if
        do j = 1, size(x%m)
            if (any(x%counts(:, j) > x%m(j))) then
                ok = .false.
                return
            end if
        end do
    end function valid_mb

    integer function dimension_from_nvals(nv) result(k)
        integer, intent(in) :: nv
        real(dp) :: kr

        if (nv < 2) then
            k = -1
            return
        end if
        kr = 0.5_dp * (-1.0_dp + sqrt(9.0_dp + 8.0_dp * real(nv, dp)))
        k = nint(kr)
        if (k - 1 + k * (k - 1) / 2 /= nv) k = -1
    end function dimension_from_nvals

end module mm_parameters

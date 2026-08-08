! SPDX-License-Identifier: GPL-2.0-or-later
module rceim_utils
    use rceim_kinds, only : dp
    implicit none
    private
    public :: enforce_domain, sort_population_by_score, sample_mean, sample_sd
contains
    pure real(dp) function sample_mean(x) result(m)
        real(dp), intent(in) :: x(:)
        if (size(x) == 0) then
            m = 0.0_dp
        else
            m = sum(x)/real(size(x), dp)
        end if
    end function sample_mean

    pure real(dp) function sample_sd(x) result(s)
        real(dp), intent(in) :: x(:)
        real(dp) :: m
        integer :: n
        n = size(x)
        if (n <= 1) then
            s = 0.0_dp
        else
            m = sample_mean(x)
            s = sqrt(max(0.0_dp, sum((x-m)**2)/real(n-1, dp)))
        end if
    end function sample_sd

    subroutine enforce_domain(params, lower, upper)
        real(dp), intent(inout) :: params(:,:)
        real(dp), intent(in) :: lower(:), upper(:)
        integer :: j
        if (size(params,2) /= size(lower) .or. size(lower) /= size(upper)) then
            error stop "enforce_domain: incompatible dimensions"
        end if
        do j = 1, size(params,2)
            params(:,j) = max(lower(j), min(upper(j), params(:,j)))
        end do
    end subroutine enforce_domain

    subroutine sort_population_by_score(params, scores)
        real(dp), intent(inout) :: params(:,:)
        real(dp), intent(inout) :: scores(:)
        real(dp) :: key_score
        real(dp), allocatable :: key_row(:)
        integer :: i, j, n

        n = size(scores)
        if (size(params,1) /= n) error stop "sort_population_by_score: incompatible dimensions"
        allocate(key_row(size(params,2)))
        do i = 2, n
            key_score = scores(i)
            key_row = params(i,:)
            j = i - 1
            do while (j >= 1)
                if (scores(j) <= key_score) exit
                scores(j+1) = scores(j)
                params(j+1,:) = params(j,:)
                j = j - 1
            end do
            scores(j+1) = key_score
            params(j+1,:) = key_row
        end do
    end subroutine sort_population_by_score
end module rceim_utils

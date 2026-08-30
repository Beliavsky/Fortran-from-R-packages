! SPDX-License-Identifier: GPL-2.0-only
module multcomp_contrasts
  use multcomp_kinds, only : dp
  use multcomp_types, only : contrast_matrix_type
  implicit none
  private

  public :: contr_mat
  public :: tukey_pairs

contains

  subroutine contr_mat(n, contrast_type, result, base)
    real(dp), intent(in) :: n(:) !! Group sample sizes or weights used by weighted contrast families.
    character(len=*), intent(in) :: contrast_type !! Contrast family name matching multcomp::contrMat.
    type(contrast_matrix_type), intent(out) :: result !! Constructed contrast matrix with one group per column.
    integer, intent(in), optional :: base !! One-based control group for Dunnett contrasts; default 1.

    character(len=:), allocatable :: key
    integer :: b
    integer :: i
    integer :: j
    integer :: k
    integer :: row
    real(dp) :: denom

    k = size(n)
    result%ok = .false.
    result%message = ''
    if (k < 2) then
      result%message = 'at least two groups are required'
      return
    end if
    if (any(n < 0.0_dp)) then
      result%message = 'group sizes must be nonnegative'
      return
    end if

    b = 1
    if (present(base)) b = base
    if (b < 1 .or. b > k) then
      result%message = 'base group is outside the valid range'
      return
    end if

    key = lowercase(trim(adjustl(contrast_type)))
    result%contrast_type = trim(contrast_type)

    select case (key)
    case ('dunnett')
      allocate(result%value(k - 1, k))
      result%value = 0.0_dp
      row = 0
      do i = 1, k
        if (i == b) cycle
        row = row + 1
        result%value(row, i) = 1.0_dp
        result%value(row, b) = -1.0_dp
      end do

    case ('tukey')
      allocate(result%value(k * (k - 1) / 2, k))
      result%value = 0.0_dp
      row = 0
      do i = 1, k - 1
        do j = i + 1, k
          row = row + 1
          result%value(row, j) = 1.0_dp
          result%value(row, i) = -1.0_dp
        end do
      end do

    case ('sequen')
      allocate(result%value(k - 1, k))
      result%value = 0.0_dp
      do i = 2, k
        result%value(i - 1, i) = 1.0_dp
        result%value(i - 1, i - 1) = -1.0_dp
      end do

    case ('ave')
      if (k < 3) then
        result%message = 'AVE contrasts require at least three groups'
        return
      end if
      allocate(result%value(k, k))
      result%value = 0.0_dp
      denom = sum(n(2:k))
      if (denom <= 0.0_dp) then
        result%message = 'AVE contrast denominator is zero'
        return
      end if
      result%value(1, 1) = 1.0_dp
      result%value(1, 2:k) = -n(2:k) / denom
      do i = 2, k - 1
        denom = sum(n(1:i - 1)) + sum(n(i + 1:k))
        if (denom <= 0.0_dp) then
          result%message = 'AVE contrast denominator is zero'
          return
        end if
        result%value(i, 1:i - 1) = -n(1:i - 1) / denom
        result%value(i, i) = 1.0_dp
        result%value(i, i + 1:k) = -n(i + 1:k) / denom
      end do
      denom = sum(n(1:k - 1))
      if (denom <= 0.0_dp) then
        result%message = 'AVE contrast denominator is zero'
        return
      end if
      result%value(k, 1:k - 1) = -n(1:k - 1) / denom
      result%value(k, k) = 1.0_dp

    case ('changepoint')
      allocate(result%value(k - 1, k))
      result%value = 0.0_dp
      do i = 1, k - 1
        if (sum(n(1:i)) <= 0.0_dp .or. sum(n(i + 1:k)) <= 0.0_dp) then
          result%message = 'Changepoint contrast denominator is zero'
          return
        end if
        result%value(i, 1:i) = -n(1:i) / sum(n(1:i))
        result%value(i, i + 1:k) = n(i + 1:k) / sum(n(i + 1:k))
      end do

    case ('williams')
      if (k < 3) then
        result%message = 'Williams contrasts require at least three groups'
        return
      end if
      allocate(result%value(k - 1, k))
      result%value = 0.0_dp
      do i = 1, k - 2
        result%value(i, 1) = -1.0_dp
        denom = sum(n(k - i + 1:k))
        if (denom <= 0.0_dp) then
          result%message = 'Williams contrast denominator is zero'
          return
        end if
        result%value(i, k - i + 1:k) = n(k - i + 1:k) / denom
      end do
      result%value(k - 1, 1) = -1.0_dp
      denom = sum(n(2:k))
      if (denom <= 0.0_dp) then
        result%message = 'Williams contrast denominator is zero'
        return
      end if
      result%value(k - 1, 2:k) = n(2:k) / denom

    case ('marcus')
      allocate(result%value(k * (k - 1) / 2, k))
      result%value = 0.0_dp
      row = 0
      do i = 1, k - 1
        do j = 1, i
          row = row + 1
          if (sum(n(i + 1:k)) <= 0.0_dp .or. sum(n(1:j)) <= 0.0_dp) then
            result%message = 'Marcus contrast denominator is zero'
            return
          end if
          result%value(row, i + 1:k) = n(i + 1:k) / sum(n(i + 1:k))
          result%value(row, 1:j) = result%value(row, 1:j) - &
            n(1:j) / sum(n(1:j))
        end do
      end do

    case ('mcdermott')
      if (k < 3) then
        result%message = 'McDermott contrasts require at least three groups'
        return
      end if
      allocate(result%value(k - 1, k))
      result%value = 0.0_dp
      do i = 1, k - 2
        denom = sum(n(1:i))
        if (denom <= 0.0_dp) then
          result%message = 'McDermott contrast denominator is zero'
          return
        end if
        result%value(i, 1:i) = -n(1:i) / denom
        result%value(i, i + 1) = 1.0_dp
      end do
      denom = sum(n(1:k - 1))
      if (denom <= 0.0_dp) then
        result%message = 'McDermott contrast denominator is zero'
        return
      end if
      result%value(k - 1, 1:k - 1) = -n(1:k - 1) / denom
      result%value(k - 1, k) = 1.0_dp

    case ('umbrellawilliams')
      allocate(result%value(k * (k - 1) / 2, k))
      result%value = 0.0_dp
      row = 0
      do j = 1, k - 1
        do i = 1, k - j
          row = row + 1
          result%value(row, 1) = -1.0_dp
          denom = sum(n(k - i - j + 2:k - j + 1))
          if (denom <= 0.0_dp) then
            result%message = 'UmbrellaWilliams contrast denominator is zero'
            return
          end if
          result%value(row, k - i - j + 2:k - j + 1) = &
            n(k - i - j + 2:k - j + 1) / denom
        end do
      end do

    case ('grandmean')
      if (sum(n) <= 0.0_dp) then
        result%message = 'GrandMean contrast denominator is zero'
        return
      end if
      allocate(result%value(k, k))
      do i = 1, k
        result%value(i, :) = -n / sum(n)
        result%value(i, i) = result%value(i, i) + 1.0_dp
      end do

    case default
      result%message = 'unknown contrast type'
      return
    end select

    result%ok = .true.
  end subroutine contr_mat

  subroutine tukey_pairs(k, lower_group, upper_group)
    integer, intent(in) :: k !! Number of factor levels represented by the Tukey contrast family.
    integer, allocatable, intent(out) :: lower_group(:) !! Group carrying coefficient -1 in each Tukey contrast.
    integer, allocatable, intent(out) :: upper_group(:) !! Group carrying coefficient +1 in each Tukey contrast.

    integer :: i
    integer :: j
    integer :: row

    allocate(lower_group(k * (k - 1) / 2), upper_group(k * (k - 1) / 2))
    row = 0
    do i = 1, k - 1
      do j = i + 1, k
        row = row + 1
        lower_group(row) = i
        upper_group(row) = j
      end do
    end do
  end subroutine tukey_pairs

  pure function lowercase(text) result(out)
    character(len=*), intent(in) :: text !! Contrast family text to normalize for matching.
    character(len=len(text)) :: out

    integer :: c
    integer :: i

    out = text
    do i = 1, len(text)
      c = iachar(out(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) out(i:i) = achar(c + 32)
    end do
  end function lowercase

end module multcomp_contrasts

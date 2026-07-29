! SPDX-License-Identifier: GPL-2.0-or-later
module actuar_types
  use actuar_kinds, only : dp
  implicit none
  private

  type, public :: aggregate_distribution
    real(dp), allocatable :: support(:)
    real(dp), allocatable :: probability(:)
    real(dp), allocatable :: cdf(:)
    logical :: ok = .false.
    character(len=160) :: message = ''
  contains
    procedure :: mean => aggregate_mean
    procedure :: variance => aggregate_variance
    procedure :: quantile => aggregate_quantile
  end type aggregate_distribution

  type, public :: credibility_result
    real(dp), allocatable :: means(:)
    real(dp), allocatable :: weights(:)
    real(dp), allocatable :: estimates(:)
    real(dp) :: collective_mean = 0.0_dp
    real(dp) :: process_variance = 0.0_dp
    real(dp) :: structural_variance = 0.0_dp
    real(dp) :: k = 0.0_dp
    logical :: ok = .false.
    character(len=160) :: message = ''
  end type credibility_result

  type, public :: phase_type_result
    real(dp) :: value = 0.0_dp
    logical :: ok = .false.
    character(len=160) :: message = ''
  end type phase_type_result

contains

  pure function aggregate_mean(this) result(x)
    class(aggregate_distribution), intent(in) :: this
    real(dp) :: x
    if (.not. allocated(this%support) .or. .not. allocated(this%probability)) then
      x = 0.0_dp
    else
      x = sum(this%support*this%probability)
    end if
  end function aggregate_mean

  pure function aggregate_variance(this) result(x)
    class(aggregate_distribution), intent(in) :: this
    real(dp) :: x, mu
    if (.not. allocated(this%support) .or. .not. allocated(this%probability)) then
      x = 0.0_dp
    else
      mu = sum(this%support*this%probability)
      x = sum((this%support-mu)**2*this%probability)
    end if
  end function aggregate_variance

  pure function aggregate_quantile(this, p) result(x)
    class(aggregate_distribution), intent(in) :: this
    real(dp), intent(in) :: p
    real(dp) :: x
    integer :: i
    x = 0.0_dp
    if (.not. allocated(this%support) .or. .not. allocated(this%cdf)) return
    do i = 1, size(this%cdf)
      if (this%cdf(i) >= p) then
        x = this%support(i)
        return
      end if
    end do
    x = this%support(size(this%support))
  end function aggregate_quantile

end module actuar_types

! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
module highfrequency_leadlag
  use highfrequency_kinds, only: dp
  use highfrequency_types, only: lead_lag_result
  use highfrequency_realized, only: hayashi_yoshida_covariance
  implicit none
  private
  public :: lead_lag

contains

  function lead_lag(times1, prices1, times2, prices2, lags, normalize) result(result)
    integer, intent(in) :: times1(:), times2(:), lags(:)
    real(dp), intent(in) :: prices1(:), prices2(:)
    logical, intent(in), optional :: normalize
    type(lead_lag_result) :: result
    integer, allocatable :: shifted(:)
    real(dp) :: c, scale
    logical :: norm
    integer :: i, loc
    norm=.true.
    if(present(normalize)) norm=normalize
    allocate(result%lags(size(lags)),result%contrast(size(lags)),shifted(size(times2)))
    result%lags=lags
    scale=sqrt(max(tiny(1.0_dp), &
      hayashi_yoshida_covariance(times1,prices1,times1,prices1)* &
      hayashi_yoshida_covariance(times2,prices2,times2,prices2)))
    do i=1,size(lags)
      shifted=times2+lags(i)
      c=hayashi_yoshida_covariance(times1,prices1,shifted,prices2)
      if(norm) c=c/scale
      result%contrast(i)=c
    end do
    if(size(lags)>0)then
      loc=maxloc(abs(result%contrast),dim=1)
      result%optimal_lag=result%lags(loc)
      result%maximum_contrast=result%contrast(loc)
    end if
  end function lead_lag

end module highfrequency_leadlag

! SPDX-License-Identifier: GPL-2.0-or-later
module compositions_zero
  use compositions_kinds, only: dp
  implicit none
  private
  public :: zero_replace, zero_replace_rows, detection_limits_from_negative
  public :: missing_summary_counts, simulate_missings_mcar

  type, public :: missing_summary_counts
    integer :: observed=0
    integer :: below_detection=0
    integer :: structural_zero=0
    integer :: missing=0
  end type
contains
  function zero_replace(x,d,a,mask) result(y)
    real(dp), intent(in) :: x(:),d(:)
    real(dp), intent(in), optional :: a
    logical, intent(in), optional :: mask(:)
    real(dp) :: y(size(x)),fac
    logical :: m(size(x))
    if(size(d)/=size(x)) error stop 'zero_replace: detection limit mismatch'
    fac=2.0_dp/3.0_dp; if(present(a)) fac=a
    if(present(mask)) then; m=mask; else; m=(x<=0.0_dp); end if
    y=x; where(m) y=fac*d
  end function zero_replace

  function zero_replace_rows(x,d,a,mask) result(y)
    real(dp), intent(in) :: x(:,:),d(:,:)
    real(dp), intent(in), optional :: a
    logical, intent(in), optional :: mask(:,:)
    real(dp) :: y(size(x,1),size(x,2)),fac
    logical :: m(size(x,1),size(x,2))
    if(any(shape(d)/=shape(x))) error stop 'zero_replace_rows: shape mismatch'
    fac=2.0_dp/3.0_dp; if(present(a)) fac=a
    if(present(mask)) then; m=mask; else; m=(x<=0.0_dp); end if
    y=x; where(m) y=fac*d
  end function zero_replace_rows

  function detection_limits_from_negative(x) result(d)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: d(size(x,1),size(x,2))
    d=0.0_dp; where(x<0.0_dp) d=-x
  end function detection_limits_from_negative

  function simulate_missings_mcar(x,prob,seed) result(y)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(in) :: prob
    integer, intent(in), optional :: seed
    real(dp) :: y(size(x,1),size(x,2)),u
    integer :: i,j,nseed,k
    integer, allocatable :: put(:)
    if(prob<0.0_dp.or.prob>1.0_dp) error stop 'simulate_missings_mcar: prob outside [0,1]'
    if(present(seed)) then
      call random_seed(size=nseed); allocate(put(nseed)); do k=1,nseed; put(k)=mod(abs(seed)+104729*k,huge(1)-1)+1; end do
      call random_seed(put=put)
    end if
    y=x
    do j=1,size(x,2); do i=1,size(x,1)
      call random_number(u); if(u<prob) y(i,j)=0.0_dp
    end do; end do
  end function simulate_missings_mcar
end module compositions_zero

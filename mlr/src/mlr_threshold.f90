module mlr_threshold
  use mlr_kinds, only : dp
  use mlr_utils, only : argsort_real
  implicit none
  private
  public :: threshold_measure, tune_binary_threshold
  abstract interface
    real(dp) function threshold_measure(truth,response)
      import dp
      integer,intent(in)::truth(:),response(:)
    end function
  end interface
contains
  subroutine tune_binary_threshold(prob,truth,negative,positive,measure,minimize,threshold,performance)
    real(dp),intent(in)::prob(:);integer,intent(in)::truth(:),negative,positive
    procedure(threshold_measure)::measure;logical,intent(in)::minimize
    real(dp),intent(out)::threshold,performance
    integer,allocatable::idx(:),response(:);integer::i,n;real(dp)::th,v,best
    if(size(prob)/=size(truth).or.size(prob)==0)error stop 'tune_binary_threshold: invalid data'
    n=size(prob);call argsort_real(prob,idx);allocate(response(n))
    best=merge(huge(1.0_dp),-huge(1.0_dp),minimize);threshold=0.5_dp
    do i=0,n
      if(i==0)then
        th=prob(idx(1))-epsilon(1.0_dp)*max(1.0_dp,abs(prob(idx(1))))
      else if(i==n)then
        th=prob(idx(n))+epsilon(1.0_dp)*max(1.0_dp,abs(prob(idx(n))))
      else
        th=0.5_dp*(prob(idx(i))+prob(idx(i+1)))
      end if
      where(prob>th);response=positive;elsewhere;response=negative;end where
      v=measure(truth,response)
      if(i==0.or.(minimize.and.v<best).or.((.not.minimize).and.v>best))then
        best=v;threshold=th
      end if
    end do
    performance=best
  end subroutine
end module mlr_threshold

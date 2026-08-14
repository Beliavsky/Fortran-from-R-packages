module mlr_tuning
  use mlr_kinds, only : dp
  use mlr_rng, only : rng_state, rng_uniform
  use mlr_types, only : tune_result
  implicit none
  private
  public :: tuning_objective, grid_search, random_search
  abstract interface
    real(dp) function tuning_objective(par)
      import dp
      real(dp),intent(in)::par(:)
    end function
  end interface
contains
  subroutine grid_search(grid,objective,minimize,result)
    real(dp),intent(in)::grid(:,:);procedure(tuning_objective)::objective;logical,intent(in)::minimize
    type(tune_result),intent(out)::result
    integer::i,best;real(dp)::v,bv
    if(size(grid,1)<1)error stop 'grid_search: empty grid'
    best=1;bv=objective(grid(1,:))
    do i=2,size(grid,1)
      v=objective(grid(i,:))
      if((minimize.and.v<bv).or.((.not.minimize).and.v>bv))then;best=i;bv=v;end if
    end do
    result%par=grid(best,:);result%objective=bv;result%evaluations=size(grid,1)
  end subroutine

  subroutine random_search(lower,upper,nevals,objective,minimize,rng,result)
    real(dp),intent(in)::lower(:),upper(:);integer,intent(in)::nevals;procedure(tuning_objective)::objective
    logical,intent(in)::minimize;type(rng_state),intent(inout)::rng;type(tune_result),intent(out)::result
    real(dp),allocatable::p(:);real(dp)::v,bv;integer::i,j
    if(size(lower)/=size(upper).or.nevals<1)error stop 'random_search: invalid arguments'
    allocate(p(size(lower)));bv=merge(huge(1.0_dp),-huge(1.0_dp),minimize)
    do i=1,nevals
      do j=1,size(p);p(j)=lower(j)+(upper(j)-lower(j))*rng_uniform(rng);end do
      v=objective(p)
      if(i==1.or.(minimize.and.v<bv).or.((.not.minimize).and.v>bv))then;result%par=p;bv=v;end if
    end do
    result%objective=bv;result%evaluations=nevals
  end subroutine
end module mlr_tuning

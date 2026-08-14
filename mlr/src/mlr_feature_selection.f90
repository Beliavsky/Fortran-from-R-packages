module mlr_feature_selection
  use mlr_kinds, only : dp
  use mlr_rng, only : rng_state, rng_uniform
  use mlr_types, only : feature_select_result
  implicit none
  private
  public :: subset_objective, feature_select_exhaustive, feature_select_forward, feature_select_random
  abstract interface
    real(dp) function subset_objective(selected)
      import dp
      logical,intent(in)::selected(:)
    end function
  end interface
contains
  subroutine feature_select_exhaustive(p,objective,minimize,result,min_features,max_features)
    integer,intent(in)::p;procedure(subset_objective)::objective;logical,intent(in)::minimize
    type(feature_select_result),intent(out)::result;integer,intent(in),optional::min_features,max_features
    integer::mask,i,lo,hi,nsel,limit;logical,allocatable::sel(:);real(dp)::v,bv
    if(p<1.or.p>24)error stop 'feature_select_exhaustive: p must be 1..24'
    lo=1;if(present(min_features))lo=min_features;hi=p;if(present(max_features))hi=max_features
    allocate(sel(p));limit=2**p-1;result%evaluations=0;bv=merge(huge(1.0_dp),-huge(1.0_dp),minimize)
    do mask=1,limit
      do i=1,p;sel(i)=btest(mask,i-1);end do;nsel=count(sel);if(nsel<lo.or.nsel>hi)cycle
      v=objective(sel);result%evaluations=result%evaluations+1
      if(result%evaluations==1.or.(minimize.and.v<bv).or.((.not.minimize).and.v>bv))then
        result%selected=sel;bv=v
      end if
    end do
    result%objective=bv
  end subroutine

  subroutine feature_select_forward(p,objective,minimize,result,max_features)
    integer,intent(in)::p;procedure(subset_objective)::objective;logical,intent(in)::minimize
    type(feature_select_result),intent(out)::result;integer,intent(in),optional::max_features
    logical,allocatable::sel(:),trial(:);integer::step,j,bestj,mx;real(dp)::v,beststep,current
    allocate(sel(p),trial(p));sel=.false.;mx=p;if(present(max_features))mx=min(p,max_features)
    current=merge(huge(1.0_dp),-huge(1.0_dp),minimize);result%evaluations=0
    do step=1,mx
      bestj=0;beststep=merge(huge(1.0_dp),-huge(1.0_dp),minimize)
      do j=1,p
        if(sel(j))cycle;trial=sel;trial(j)=.true.;v=objective(trial);result%evaluations=result%evaluations+1
        if(bestj==0.or.(minimize.and.v<beststep).or.((.not.minimize).and.v>beststep))then;bestj=j;beststep=v;end if
      end do
      if(step>1)then
        if(minimize.and.beststep>=current)exit
        if((.not.minimize).and.beststep<=current)exit
      end if
      sel(bestj)=.true.;current=beststep;result%selected=sel;result%objective=current
    end do
  end subroutine

  subroutine feature_select_random(p,nevals,objective,minimize,rng,result,prob)
    integer,intent(in)::p,nevals;procedure(subset_objective)::objective;logical,intent(in)::minimize
    type(rng_state),intent(inout)::rng;type(feature_select_result),intent(out)::result;real(dp),intent(in),optional::prob
    logical,allocatable::sel(:);real(dp)::pr,v,bv;integer::i,j
    pr=0.5_dp;if(present(prob))pr=prob;allocate(sel(p));bv=merge(huge(1.0_dp),-huge(1.0_dp),minimize)
    do i=1,nevals
      do j=1,p;sel(j)=rng_uniform(rng)<pr;end do;if(.not.any(sel))sel(1)=.true.;v=objective(sel)
      if(i==1.or.(minimize.and.v<bv).or.((.not.minimize).and.v>bv))then;result%selected=sel;bv=v;end if
    end do
    result%objective=bv;result%evaluations=nevals
  end subroutine
end module mlr_feature_selection

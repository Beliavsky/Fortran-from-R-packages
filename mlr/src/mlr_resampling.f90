module mlr_resampling
  use mlr_kinds, only : dp
  use mlr_rng, only : rng_state, rng_shuffle, rng_integer
  use mlr_types, only : resample_plan
  implicit none
  private
  public :: make_holdout, make_kfold, make_repeated_kfold, make_subsample, make_bootstrap_oob
contains
  subroutine make_holdout(n,ratio,rng,plan)
    integer,intent(in)::n;real(dp),intent(in)::ratio;type(rng_state),intent(inout)::rng
    type(resample_plan),intent(out)::plan
    integer,allocatable::idx(:);integer::nt
    if(n<2.or.ratio<=0.0_dp.or.ratio>=1.0_dp)error stop "make_holdout: invalid arguments"
    allocate(idx(n));idx=[(nt,nt=1,n)];call rng_shuffle(rng,idx);nt=max(1,min(n-1,nint(ratio*real(n,dp))))
    allocate(plan%train(1),plan%test(1),plan%group(1));plan%group=1
    plan%train(1)%idx=idx(1:nt);plan%test(1)%idx=idx(nt+1:n)
  end subroutine

  subroutine make_kfold(n,k,rng,plan)
    integer,intent(in)::n,k;type(rng_state),intent(inout)::rng;type(resample_plan),intent(out)::plan
    integer,allocatable::idx(:);integer::i,j,lo,hi,nt
    if(k<2.or.k>n)error stop "make_kfold: invalid k"
    allocate(idx(n));idx=[(i,i=1,n)];call rng_shuffle(rng,idx)
    allocate(plan%train(k),plan%test(k),plan%group(k));plan%group=1
    do i=1,k
      lo=1+(i-1)*n/k;hi=i*n/k;plan%test(i)%idx=idx(lo:hi);nt=n-(hi-lo+1)
      allocate(plan%train(i)%idx(nt));j=0
      if(lo>1)then;plan%train(i)%idx(1:lo-1)=idx(1:lo-1);j=lo-1;end if
      if(hi<n)plan%train(i)%idx(j+1:nt)=idx(hi+1:n)
    end do
  end subroutine

  subroutine make_repeated_kfold(n,k,repeats,rng,plan)
    integer,intent(in)::n,k,repeats;type(rng_state),intent(inout)::rng;type(resample_plan),intent(out)::plan
    type(resample_plan)::tmp;integer::r,i,pos
    allocate(plan%train(k*repeats),plan%test(k*repeats),plan%group(k*repeats));pos=0
    do r=1,repeats
      call make_kfold(n,k,rng,tmp)
      do i=1,k;pos=pos+1;plan%train(pos)%idx=tmp%train(i)%idx;plan%test(pos)%idx=tmp%test(i)%idx;plan%group(pos)=r;end do
    end do
  end subroutine

  subroutine make_subsample(n,iters,ratio,rng,plan)
    integer,intent(in)::n,iters;real(dp),intent(in)::ratio;type(rng_state),intent(inout)::rng
    type(resample_plan),intent(out)::plan
    integer::r,i,nt;integer,allocatable::idx(:)
    nt=max(1,min(n-1,nint(ratio*real(n,dp))));allocate(plan%train(iters),plan%test(iters),plan%group(iters));plan%group=1
    do r=1,iters
      allocate(idx(n));idx=[(i,i=1,n)];call rng_shuffle(rng,idx)
      plan%train(r)%idx=idx(1:nt);plan%test(r)%idx=idx(nt+1:n);deallocate(idx)
    end do
  end subroutine

  subroutine make_bootstrap_oob(n,iters,rng,plan)
    integer,intent(in)::n,iters;type(rng_state),intent(inout)::rng;type(resample_plan),intent(out)::plan
    integer::r,i,j,noob;logical,allocatable::used(:)
    allocate(plan%train(iters),plan%test(iters),plan%group(iters));plan%group=1
    do r=1,iters
      allocate(plan%train(r)%idx(n),used(n));used=.false.
      do i=1,n;j=rng_integer(rng,1,n);plan%train(r)%idx(i)=j;used(j)=.true.;end do
      noob=count(.not.used);allocate(plan%test(r)%idx(noob));j=0
      do i=1,n;if(.not.used(i))then;j=j+1;plan%test(r)%idx(j)=i;end if;end do
      deallocate(used)
    end do
  end subroutine
end module mlr_resampling

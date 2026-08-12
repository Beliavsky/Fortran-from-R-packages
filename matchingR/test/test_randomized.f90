program test_randomized
   use matchingr
   implicit none
   integer :: n, trial, i, j, r
   real(dp), allocatable :: ua(:,:), ub(:,:), ur(:,:), us(:,:)
   integer, allocatable :: pref(:,:), tmatch(:), prop(:,:), eng(:,:)
   type(marriage_result_t) :: mr
   type(roommate_result_t) :: rr
   integer, allocatable :: seed(:)

   call random_seed(size=n)
   allocate(seed(n)); seed=[(104729+37*i,i=1,n)]; call random_seed(put=seed)

   do trial=1,50
      n=4+mod(trial,5)
      allocate(ua(n,n),ub(n,n))
      call random_number(ua); call random_number(ub)
      mr=marriage_market(ua,ub)
      allocate(prop(n,1),eng(n,1)); prop(:,1)=mr%proposals; eng(:,1)=mr%engagements
      if(.not.gale_shapley_stable(ua,ub,prop,eng)) error stop "random Gale-Shapley instability"
      deallocate(ua,ub,prop,eng)
   end do

   do trial=1,40
      n=2*(2+mod(trial,4))
      allocate(ur(n,n)); call random_number(ur)
      do i=1,n
         ur(i,i)=-1.0_dp
      end do
      allocate(us(n-1,n))
      do j=1,n
         r=0
         do i=1,n
            if(i==j) cycle
            r=r+1; us(r,j)=ur(i,j)
         end do
      end do
      pref=sort_index_one_sided(us)
      rr=stable_roommates_preferences(pref)
      if(rr%stable_exists) then
         do i=1,n
            if(rr%matching(i)<1 .or. rr%matching(i)>n) error stop "roommate index"
            if(rr%matching(rr%matching(i))/=i) error stop "roommate not reciprocal"
         end do
         if(.not.roommate_stable(pref,rr%matching,.false.)) error stop "random roommate instability"
      end if
      deallocate(ur,us)
   end do

   do trial=1,40
      n=3+mod(trial,8)
      allocate(ua(n,n)); call random_number(ua)
      pref=sort_index(ua)
      tmatch=top_trading_cycles_preferences(pref)
      if(.not.top_trading_stable(pref,tmatch)) error stop "random TTC pair instability"
      if(.not.is_permutation(tmatch)) error stop "TTC not permutation"
      deallocate(ua)
   end do

   print *, "test_randomized: PASS"
contains
   logical function is_permutation(x) result(ok)
      integer,intent(in)::x(:)
      logical :: seen(size(x))
      integer :: k
      seen=.false.; ok=.false.
      do k=1,size(x)
         if(x(k)<1 .or. x(k)>size(x)) return
         if(seen(x(k))) return
         seen(x(k))=.true.
      end do
      ok=.true.
   end function is_permutation
end program test_randomized

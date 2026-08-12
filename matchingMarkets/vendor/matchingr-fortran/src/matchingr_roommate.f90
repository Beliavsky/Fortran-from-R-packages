module matchingr_roommate
   use matchingr_kinds, only : dp
   use matchingr_types, only : roommate_result_t
   use matchingr_utils, only : sort_index_one_sided, check_roommate_preferences
   implicit none
   private
   public :: stable_roommates, stable_roommates_preferences, roommate_stable

   type :: pref_list_t
      integer, allocatable :: a(:)
      integer :: n = 0
   end type pref_list_t

contains

   function stable_roommates(utils) result(res)
      real(dp), intent(in) :: utils(:,:)
      type(roommate_result_t) :: res
      real(dp), allocatable :: u(:,:)
      integer, allocatable :: pref(:,:)
      integer :: n, i, j, r
      n = size(utils,2)
      allocate(u(0,0),pref(0,0))
      if (size(utils,1) == n) then
         allocate(u(n-1,n))
         do j=1,n
            r=0
            do i=1,n
               if(i==j) cycle
               r=r+1; u(r,j)=utils(i,j)
            end do
         end do
      else
         if(size(utils,1)/=n-1) error stop "stable_roommates: utilities must be n x n or (n-1) x n"
         u=utils
      end if
      pref=sort_index_one_sided(u)
      res=stable_roommates_preferences(pref)
   end function stable_roommates

   function stable_roommates_preferences(pref_in) result(res)
      integer, intent(in) :: pref_in(:,:)
      type(roommate_result_t) :: res
      integer, allocatable :: pref(:,:), aug(:,:), match(:)
      integer :: n, j
      logical :: zb
      allocate(pref(0,0),match(0))
      if(.not.check_roommate_preferences(pref_in,zb)) error stop "stable_roommates_preferences: incomplete preferences"
      pref=pref_in; if(zb) pref=pref+1
      n=size(pref,2)
      if(mod(n,2)==1) then
         allocate(aug(n,n+1))
         aug(1:n-1,1:n)=pref
         aug(n,1:n)=n+1
         aug(:,n+1)=[(j,j=1,n)]
         match=irving_core(aug)
         allocate(res%matching(n)); res%matching=0
         if(any(match/=0)) then
            do j=1,n
               if(match(j)<=n) res%matching(j)=match(j)
            end do
            res%stable_exists=.true.
         end if
      else
         match=irving_core(pref)
         allocate(res%matching(n)); res%matching=match
         res%stable_exists=any(match/=0)
      end if
   end function stable_roommates_preferences

   function irving_core(pref) result(matchings)
      integer, intent(in) :: pref(:,:)
      integer, allocatable :: matchings(:)
      integer :: n, person, proposee, op, op_curr, i, stable_count
      integer :: new_index, new_x, rot_tail, ii, target, old_size
      integer, allocatable :: proposal_to(:), proposal_from(:), proposed_to(:)
      integer, allocatable :: x(:), index(:)
      type(pref_list_t), allocatable :: table(:)
      logical :: stable, erased

      n=size(pref,2)
      allocate(matchings(n)); matchings=0
      if(size(pref,1)/=n-1) return
      allocate(proposal_to(n),proposal_from(n),proposed_to(n))
      proposal_to=0; proposal_from=0; proposed_to=1

      stable=.false.
      do while(.not.stable)
         stable=.true.
         do person=1,n
            if(proposal_to(person)/=0) cycle
            if(proposed_to(person)>n-1) then
               matchings=0; return
            end if
            proposee=pref(proposed_to(person),person)
            op=position_in_pref(pref(:,proposee),person)
            if(op==0) then; matchings=0; return; end if
            if(proposal_from(proposee)==0) then
               op_curr=n
            else
               op_curr=position_in_pref(pref(:,proposee),proposal_from(proposee))
               if(op_curr==0) op_curr=n
            end if
            if(op<op_curr) then
               proposal_to(person)=proposee
               if(proposal_from(proposee)/=0) then
                  proposal_to(proposal_from(proposee))=0
                  stable=.false.
               end if
               proposal_from(proposee)=person
            else
               stable=.false.
            end if
            proposed_to(person)=proposed_to(person)+1
         end do
      end do

      allocate(table(n))
      do person=1,n
         allocate(table(person)%a(n-1)); table(person)%a=pref(:,person); table(person)%n=n-1
      end do

      ! Phase 1 reduction: in person's list keep entries through the proposer
      ! currently held; delete each discarded pair reciprocally.
      do person=1,n
         do
            if(table(person)%n<=0) then; matchings=0; return; end if
            if(table(person)%a(table(person)%n)==proposal_from(person)) exit
            target=table(person)%a(table(person)%n)
            erased=list_remove(table(target),person)
            if(.not.erased) then; matchings=0; return; end if
            call list_pop_back(table(person))
         end do
      end do

      do
         stable=.true.
         do person=1,n
            if(table(person)%n<=1) cycle
            stable=.false.
            allocate(x(n+1),index(n+1)); x=0; index=0
            new_index=person
            stable_count=0
            do
               if(table(new_index)%n<2) then
                  matchings=0; deallocate(x,index); return
               end if
               new_x=table(new_index)%a(2)
               if(table(new_x)%n<1) then
                  matchings=0; deallocate(x,index); return
               end if
               new_index=table(new_x)%a(table(new_x)%n)
               rot_tail=0
               do i=1,stable_count
                  if(index(i)==new_index) then; rot_tail=i; exit; end if
               end do
               stable_count=stable_count+1
               x(stable_count)=new_x; index(stable_count)=new_index
               if(rot_tail/=0) exit
            end do

            ! C++ uses elements strictly after the first repeated index.
            do ii=rot_tail+1,stable_count
               do
                  if(table(x(ii))%n<1) then; matchings=0; deallocate(x,index); return; end if
                  if(table(x(ii))%a(table(x(ii))%n)==index(ii-1)) exit
                  target=table(x(ii))%a(table(x(ii))%n)
                  old_size=table(target)%n
                  erased=list_remove(table(target),x(ii))
                  if(.not.erased .or. table(target)%n==old_size) then
                     matchings=0; deallocate(x,index); return
                  end if
                  if(table(x(ii))%n==1) then
                     matchings=0; deallocate(x,index); return
                  end if
                  call list_pop_back(table(x(ii)))
               end do
            end do
            deallocate(x,index)
         end do
         if(stable) exit
      end do

      do person=1,n
         if(table(person)%n<1) then; matchings=0; return; end if
         matchings(person)=table(person)%a(1)
      end do
   end function irving_core

   integer function position_in_pref(col, value) result(pos)
      integer, intent(in) :: col(:), value
      integer :: i
      pos=0
      do i=1,size(col)
         if(col(i)==value) then; pos=i; return; end if
      end do
   end function position_in_pref

   logical function list_remove(list,value) result(erased)
      type(pref_list_t), intent(inout) :: list
      integer, intent(in) :: value
      integer :: i
      erased=.false.
      do i=1,list%n
         if(list%a(i)==value) then
            if(i<list%n) list%a(i:list%n-1)=list%a(i+1:list%n)
            list%n=list%n-1; erased=.true.; return
         end if
      end do
   end function list_remove

   subroutine list_pop_back(list)
      type(pref_list_t), intent(inout) :: list
      if(list%n>0) list%n=list%n-1
   end subroutine list_pop_back

   logical function roommate_stable(pref_in,matching,legacy_cpp_bug) result(stable)
      integer, intent(in) :: pref_in(:,:), matching(:)
      logical, intent(in), optional :: legacy_cpp_bug
      integer, allocatable :: pmat(:,:)
      logical :: zb, legacy, ip, jp
      integer :: n,i,j,k,target_j
      if(.not.check_roommate_preferences(pref_in,zb)) then; stable=.false.; return; end if
      pmat=pref_in; if(zb) pmat=pmat+1
      n=size(pmat,2)
      if(size(matching)/=n) then; stable=.false.; return; end if
      legacy=.false.; if(present(legacy_cpp_bug)) legacy=legacy_cpp_bug
      stable=.true.
      do i=1,n
         do j=i+1,n
            ip=.false.; jp=.false.
            if(legacy) then
               ! Preserve matchingR C++ source order: the candidate test occurs
               ! before the current-match break, and j's scan compares with j.
               do k=1,n-1
                  if(pmat(k,i)==j) ip=.true.
                  if(pmat(k,i)==matching(i)) exit
               end do
               target_j=j
               do k=1,n-1
                  if(pmat(k,j)==target_j) jp=.true.
                  if(pmat(k,j)==matching(j)) exit
               end do
            else
               do k=1,n-1
                  if(pmat(k,i)==matching(i)) exit
                  if(pmat(k,i)==j) ip=.true.
               end do
               target_j=i
               do k=1,n-1
                  if(pmat(k,j)==matching(j)) exit
                  if(pmat(k,j)==target_j) jp=.true.
               end do
            end if
            if(ip .and. jp) then; stable=.false.; return; end if
         end do
      end do
   end function roommate_stable

end module matchingr_roommate

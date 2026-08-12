module matchingr_ttc
   use matchingr_kinds, only : dp
   use matchingr_utils, only : sort_index, check_preferences
   implicit none
   private
   public :: top_trading_cycles, top_trading_cycles_preferences, top_trading_stable

contains

   function top_trading_cycles(utils) result(matching)
      real(dp), intent(in) :: utils(:,:)
      integer, allocatable :: matching(:)
      integer, allocatable :: pref(:,:)
      if(size(utils,1)/=size(utils,2)) error stop "top_trading_cycles: utilities must be square"
      pref=sort_index(utils)
      matching=top_trading_cycles_preferences(pref)
   end function top_trading_cycles

   function top_trading_cycles_preferences(pref_in) result(matching)
      integer, intent(in) :: pref_in(:,:)
      integer, allocatable :: matching(:)
      integer, allocatable :: pref(:,:), points(:), state(:), pos(:), path(:)
      logical, allocatable :: active(:), in_cycle(:)
      logical :: zb
      integer :: n, i, j, k, p, cur, plen, cycle_start, matched_count
      if(.not.check_preferences(pref_in,zb)) error stop "top_trading_cycles_preferences: incomplete preferences"
      pref=pref_in; if(zb) pref=pref+1
      n=size(pref,2)
      if(size(pref,1)/=n) error stop "top_trading_cycles_preferences: preference matrix must be square"
      allocate(matching(n),points(n),state(n),pos(n),path(n),active(n),in_cycle(n))
      matching=0; active=.true.; matched_count=0
      do while(matched_count<n)
         do i=1,n
            if(.not.active(i)) cycle
            points(i)=0
            do k=1,n
               p=pref(k,i)
               if(active(p)) then; points(i)=p; exit; end if
            end do
            if(points(i)==0) error stop "top_trading_cycles_preferences: no active preference"
         end do
         state=0; pos=0; in_cycle=.false.
         do i=1,n
            if(.not.active(i) .or. state(i)/=0) cycle
            plen=0; cur=i
            do
               if(state(cur)==0) then
                  plen=plen+1; path(plen)=cur; state(cur)=i; pos(cur)=plen
                  cur=points(cur)
               else if(state(cur)==i) then
                  cycle_start=pos(cur)
                  do j=cycle_start,plen
                     in_cycle(path(j))=.true.
                  end do
                  exit
               else
                  exit
               end if
            end do
         end do
         if(.not.any(in_cycle .and. active)) error stop "top_trading_cycles_preferences: cycle detection failure"
         do i=1,n
            if(active(i) .and. in_cycle(i)) then
               matching(i)=points(i)
            end if
         end do
         do i=1,n
            if(active(i) .and. in_cycle(i)) then
               active(i)=.false.; matched_count=matched_count+1
            end if
         end do
      end do
   end function top_trading_cycles_preferences

   logical function top_trading_stable(pref_in,matching) result(stable)
      integer, intent(in) :: pref_in(:,:), matching(:)
      integer, allocatable :: pref(:,:)
      logical :: zb, ip, jp
      integer :: n,i,j,k
      if(.not.check_preferences(pref_in,zb)) then; stable=.false.; return; end if
      pref=pref_in; if(zb) pref=pref+1
      n=size(pref,2)
      if(size(pref,1)/=n .or. size(matching)/=n) then; stable=.false.; return; end if
      stable=.true.
      do i=1,n
         do j=i,n
            ip=.false.; jp=.false.
            do k=1,n
               if(pref(k,i)==matching(i)) exit
               if(pref(k,i)==j) ip=.true.
            end do
            do k=1,n
               if(pref(k,j)==matching(j)) exit
               if(pref(k,j)==i) jp=.true.
            end do
            if(ip .and. jp) then; stable=.false.; return; end if
         end do
      end do
   end function top_trading_stable

end module matchingr_ttc

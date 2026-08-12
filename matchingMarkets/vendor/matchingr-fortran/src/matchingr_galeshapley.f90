module matchingr_galeshapley
   use matchingr_kinds, only : dp
   use matchingr_types, only : marriage_result_t, college_result_t
   use matchingr_utils, only : sort_index, rank_index, check_preferences, repeat_rows_real, repeat_cols_real
   implicit none
   private
   public :: marriage_market, marriage_market_preferences
   public :: college_admissions, college_admissions_preferences
   public :: gale_shapley_stable

contains

   function marriage_market(proposer_utils, reviewer_utils) result(res)
      real(dp), intent(in) :: proposer_utils(:,:), reviewer_utils(:,:)
      type(marriage_result_t) :: res
      integer, allocatable :: pref(:,:)
      if (size(proposer_utils,1) /= size(reviewer_utils,2) .or. &
          size(proposer_utils,2) /= size(reviewer_utils,1)) then
         error stop "marriage_market: incompatible utility dimensions"
      end if
      pref = sort_index(proposer_utils)
      res = gale_shapley_core(pref, reviewer_utils)
   end function marriage_market

   function marriage_market_preferences(proposer_pref, reviewer_pref) result(res)
      integer, intent(in) :: proposer_pref(:,:), reviewer_pref(:,:)
      type(marriage_result_t) :: res
      integer, allocatable :: pp(:,:), rp(:,:), rr(:,:)
      real(dp), allocatable :: ru(:,:)
      logical :: zb
      if (.not. check_preferences(proposer_pref, zb)) error stop "incomplete proposer preferences"
      pp = proposer_pref
      if (zb) pp = pp + 1
      if (.not. check_preferences(reviewer_pref, zb)) error stop "incomplete reviewer preferences"
      rp = reviewer_pref
      if (zb) rp = rp + 1
      if (size(pp,1) /= size(rp,2) .or. size(pp,2) /= size(rp,1)) then
         error stop "marriage_market_preferences: incompatible dimensions"
      end if
      rr = rank_index(rp)
      allocate(ru(size(rr,1),size(rr,2)))
      ru = -real(rr,dp)
      res = gale_shapley_core(pp, ru)
   end function marriage_market_preferences

   function gale_shapley_core(proposer_pref, reviewer_utils) result(res)
      integer, intent(in) :: proposer_pref(:,:)
      real(dp), intent(in) :: reviewer_utils(:,:)
      type(marriage_result_t) :: res
      integer :: m, n, cap, head, tail, proposer, j, w, old, i, ns1, ns2
      integer, allocatable :: queue(:)
      m = size(proposer_pref,2)
      n = size(proposer_pref,1)
      if (size(reviewer_utils,1) /= m .or. size(reviewer_utils,2) /= n) then
         error stop "gale_shapley_core: incompatible dimensions"
      end if
      allocate(res%proposals(m),res%engagements(n))
      res%proposals = 0; res%engagements = 0
      cap = max(1,m*(n+2)+m)
      allocate(queue(cap)); head = 1; tail = 0
      do i = m, 1, -1
         tail = tail + 1; queue(tail) = i
      end do
      do while (head <= tail)
         proposer = queue(head); head = head + 1
         do j = 1, n
            w = proposer_pref(j,proposer)
            if (w < 1 .or. w > n) error stop "gale_shapley_core: preference out of range"
            if (res%engagements(w) == 0) then
               res%engagements(w) = proposer
               res%proposals(proposer) = w
               exit
            end if
            old = res%engagements(w)
            if (reviewer_utils(proposer,w) > reviewer_utils(old,w)) then
               res%proposals(old) = 0
               if (tail >= size(queue)) call grow_queue(queue)
               tail = tail + 1; queue(tail) = old
               res%engagements(w) = proposer
               res%proposals(proposer) = w
               exit
            end if
         end do
      end do
      ns1 = count(res%proposals == 0); ns2 = count(res%engagements == 0)
      allocate(res%single_proposers(ns1),res%single_reviewers(ns2))
      ns1 = 0
      do i = 1,m
         if (res%proposals(i)==0) then; ns1=ns1+1; res%single_proposers(ns1)=i; end if
      end do
      ns2 = 0
      do i = 1,n
         if (res%engagements(i)==0) then; ns2=ns2+1; res%single_reviewers(ns2)=i; end if
      end do
   end function gale_shapley_core

   subroutine grow_queue(q)
      integer, allocatable, intent(inout) :: q(:)
      integer, allocatable :: tmp(:)
      allocate(tmp(max(2,2*size(q))))
      tmp(1:size(q)) = q
      call move_alloc(tmp,q)
   end subroutine grow_queue

   logical function gale_shapley_stable(proposer_utils, reviewer_utils, proposals, engagements) result(stable)
      real(dp), intent(in) :: proposer_utils(:,:), reviewer_utils(:,:)
      integer, intent(in) :: proposals(:,:), engagements(:,:)
      integer :: m, n, sp, sr, w, f, ip, ir, pcur, rcur
      real(dp) :: upcur, urcur
      m = size(proposer_utils,2); n = size(proposer_utils,1)
      if (size(reviewer_utils,1)/=m .or. size(reviewer_utils,2)/=n) then
         stable=.false.; return
      end if
      if (size(proposals,1)/=m .or. size(engagements,1)/=n) then
         stable=.false.; return
      end if
      sp=size(proposals,2); sr=size(engagements,2)
      stable=.true.
      do w=1,m
         do f=1,n
            do ip=1,sp
               pcur=proposals(w,ip)
               if (pcur==0) then; upcur=-huge(1.0_dp); else; upcur=proposer_utils(pcur,w); end if
               do ir=1,sr
                  rcur=engagements(f,ir)
                  if (rcur==0) then; urcur=-huge(1.0_dp); else; urcur=reviewer_utils(rcur,f); end if
                  if (reviewer_utils(w,f)>urcur .and. proposer_utils(f,w)>upcur) then
                     stable=.false.; return
                  end if
               end do
            end do
         end do
      end do
   end function gale_shapley_stable

   function college_admissions(student_utils, college_utils, slots, student_optimal) result(res)
      real(dp), intent(in) :: student_utils(:,:), college_utils(:,:)
      integer, intent(in) :: slots(:)
      logical, intent(in), optional :: student_optimal
      type(college_result_t) :: res
      logical :: so
      integer :: ns, nc, c, s, pos, k, nv
      real(dp), allocatable :: pu(:,:), ru(:,:)
      integer, allocatable :: pref(:,:)
      type(marriage_result_t) :: mr
      ns=size(student_utils,2); nc=size(student_utils,1)
      if (size(college_utils,1)/=ns .or. size(college_utils,2)/=nc) error stop "college_admissions: dimension mismatch"
      if (size(slots)/=nc .or. any(slots<0)) error stop "college_admissions: invalid slots"
      so=.true.; if (present(student_optimal)) so=student_optimal
      allocate(res%slots(nc)); res%slots=slots
      allocate(res%matched_students(ns)); res%matched_students=0
      allocate(res%matched_colleges(nc,max(1,maxval(slots)))); res%matched_colleges=0
      if (so) then
         pu = repeat_rows_real(student_utils,slots)
         ru = repeat_cols_real(college_utils,slots)
         pref = sort_index(pu)
         mr = gale_shapley_core(pref,ru)
         pos=0; nv=0
         do c=1,nc
            do s=1,slots(c)
               pos=pos+1
               if (mr%engagements(pos)>0) then
                  res%matched_colleges(c,s)=mr%engagements(pos)
                  res%matched_students(mr%engagements(pos))=c
               else
                  nv=nv+1
               end if
            end do
         end do
      else
         pu = repeat_cols_real(college_utils,slots)
         ru = repeat_rows_real(student_utils,slots)
         pref = sort_index(pu)
         mr = gale_shapley_core(pref,ru)
         pos=0; nv=0
         do c=1,nc
            do s=1,slots(c)
               pos=pos+1
               if (mr%proposals(pos)>0) then
                  res%matched_colleges(c,s)=mr%proposals(pos)
                  res%matched_students(mr%proposals(pos))=c
               else
                  nv=nv+1
               end if
            end do
         end do
      end if
      allocate(res%unmatched_students(count(res%matched_students==0)))
      k=0
      do s=1,ns
         if(res%matched_students(s)==0) then;k=k+1;res%unmatched_students(k)=s;end if
      end do
      allocate(res%unmatched_colleges(nv)); k=0
      do c=1,nc
         do s=1,slots(c)
            if(res%matched_colleges(c,s)==0) then;k=k+1;res%unmatched_colleges(k)=c;end if
         end do
      end do
   end function college_admissions

   function college_admissions_preferences(student_pref, college_pref, slots, student_optimal) result(res)
      integer, intent(in) :: student_pref(:,:), college_pref(:,:), slots(:)
      logical, intent(in), optional :: student_optimal
      type(college_result_t) :: res
      integer, allocatable :: sp(:,:), cp(:,:), sr(:,:), cr(:,:)
      real(dp), allocatable :: su(:,:), cu(:,:)
      logical :: zb, so
      if(.not.check_preferences(student_pref,zb)) error stop "invalid student preferences"
      sp=student_pref; if(zb) sp=sp+1
      if(.not.check_preferences(college_pref,zb)) error stop "invalid college preferences"
      cp=college_pref; if(zb) cp=cp+1
      sr=rank_index(sp); cr=rank_index(cp)
      allocate(su(size(sr,1),size(sr,2)),cu(size(cr,1),size(cr,2)))
      su=-real(sr,dp); cu=-real(cr,dp)
      so=.true.; if(present(student_optimal)) so=student_optimal
      res=college_admissions(su,cu,slots,so)
   end function college_admissions_preferences

end module matchingr_galeshapley

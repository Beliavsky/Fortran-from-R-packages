module matchingmarkets_mechanisms
   use matchingmarkets_kinds, only : dp, i8
   use matchingmarkets_types
   use matchingmarkets_rng, only : rng_t
   use matchingmarkets_utils
   implicit none
   private
   public :: iaa, hri3_eadam, rsd, ttc_tenants, ttc_school, ttcc_kidney
   public :: stability_check, hri_all, sri_all, hri2_couples_exact

contains

   function iaa(s_pref,c_pref,slots,acceptance) result(res)
      integer,intent(in)::s_pref(:,:),c_pref(:,:),slots(:)
      character(len=*),intent(in),optional::acceptance
      type(assignment_result_t)::res
      character(len=16)::mode
      integer :: ns,nc,s,c,i,worst,worstr,iter,nactive,nrej
      integer,allocatable::next(:),active(:),rejected(:),offer(:),assign(:)
      logical,allocatable::accepted(:,:)
      logical :: deferred
      ns=size(s_pref,2); nc=size(s_pref,1)
      if(size(c_pref,1)/=ns .or. size(c_pref,2)/=nc) error stop 'iaa: dimension mismatch'
      if(size(slots)/=nc) error stop 'iaa: slots mismatch'
      mode='immediate'; if(present(acceptance)) mode=adjustl(acceptance)
      deferred=(trim(mode)=='deferred')
      allocate(next(ns),active(ns),rejected(ns),offer(ns),assign(ns),accepted(ns,nc))
      next=0;assign=0;accepted=.false.;active=[(s,s=1,ns)];nactive=ns;iter=0
      do while(nactive>0)
         iter=iter+1;offer=0
         do i=1,nactive
            s=active(i)
            c=next_valid(s_pref(:,s),next(s)+1)
            if(c>0) then
               next(s)=position_from(s_pref(:,s),next(s)+1,c);offer(s)=c
            end if
         end do
         nrej=0
         do c=1,nc
            if(.not.any(offer(active(:nactive))==c))cycle
            if(deferred) then
               do i=1,nactive
                  s=active(i);if(offer(s)/=c)cycle
                  if(rank_of(c_pref(:,c),s)<huge(1)) then
                     accepted(s,c)=.true.
                  else
                     nrej=nrej+1;rejected(nrej)=s
                  end if
               end do
               do while(count(accepted(:,c))>slots(c))
                  worst=0;worstr=-1
                  do s=1,ns
                     if(.not.accepted(s,c))cycle
                     if(rank_of(c_pref(:,c),s)>worstr) then
                        worst=s;worstr=rank_of(c_pref(:,c),s)
                     end if
                  end do
                  accepted(worst,c)=.false.;nrej=nrej+1;rejected(nrej)=worst
               end do
            else
               do s=1,ns
                  if(offer(s)/=c)cycle
                  if(rank_of(c_pref(:,c),s)>=huge(1)) then
                     nrej=nrej+1;rejected(nrej)=s;offer(s)=0
                  end if
               end do
               do i=1,min(slots(c)-count(accepted(:,c)),ns)
                  worst=0;worstr=huge(1)
                  do s=1,ns
                     if(offer(s)/=c .or. accepted(s,c))cycle
                     if(rank_of(c_pref(:,c),s)<worstr) then
                        worst=s;worstr=rank_of(c_pref(:,c),s)
                     end if
                  end do
                  if(worst==0 .or. worstr==huge(1))exit
                  accepted(worst,c)=.true.;offer(worst)=0
               end do
               do s=1,ns
                  if(offer(s)==c) then;nrej=nrej+1;rejected(nrej)=s;end if
               end do
            end if
         end do
         nactive=0
         do i=1,nrej
            s=rejected(i)
            if(next_valid(s_pref(:,s),next(s)+1)>0) then
               nactive=nactive+1;active(nactive)=s
            end if
         end do
      end do
      do c=1,nc;do s=1,ns;if(accepted(s,c))assign(s)=c;end do;end do
      allocate(res%assignment(ns));res%assignment=assign;res%iterations=iter
      res%pairs=pairs_from_assignment(assign);res%singles=singles_from_assignment(assign)
      allocate(res%free_capacity(nc));res%free_capacity=slots
      do s=1,ns;if(assign(s)>0)res%free_capacity(assign(s))=res%free_capacity(assign(s))-1;end do
   end function iaa

   function hri3_eadam(s_pref_in,c_pref,slots,acceptance,consent,max_ea) result(out)
      integer,intent(in)::s_pref_in(:,:),c_pref(:,:),slots(:)
      character(len=*),intent(in),optional::acceptance
      logical,intent(in),optional::consent(:)
      integer,intent(in),optional::max_ea
      type(eadam_result_t)::out
      integer :: ns,nc,ea,lim,iter,s,c,i,k,worst,worstr,nactive,nrej,nfinal
      integer,allocatable :: sp(:,:),next(:),active(:),rejected(:),finals(:),offer(:),assign(:)
      integer,allocatable :: ip_college(:),current_ip(:)
      logical,allocatable :: cons(:),was_reject(:),accepted(:,:),new_assign(:,:)
      logical :: deferred,interrupting_flag
      character(len=16)::mode
      ns=size(s_pref_in,2);nc=size(s_pref_in,1)
      if(size(c_pref,1)/=ns .or. size(c_pref,2)/=nc) error stop 'hri3_eadam: dimensions'
      if(size(slots)/=nc) error stop 'hri3_eadam: slots mismatch'
      sp=s_pref_in
      allocate(cons(ns));cons=.false.;if(present(consent)) cons=consent
      lim=1000000;if(present(max_ea))lim=max_ea
      mode='deferred';if(present(acceptance))mode=adjustl(acceptance)
      deferred=(trim(mode)=='deferred')
      allocate(ip_college(ns));ip_college=0
      ea=-1
      do
         ea=ea+1
         do s=1,ns
            if(ip_college(s)>0) call remove_pref(sp(:,s),ip_college(s))
         end do
         ip_college=0
         allocate(next(ns),active(ns),rejected(ns),finals(ns),offer(ns),assign(ns), &
            current_ip(ns),was_reject(nc),accepted(ns,nc),new_assign(ns,nc))
         next=0;assign=0;current_ip=0;was_reject=.false.;accepted=.false.
         nactive=0;nfinal=0
         do s=1,ns
            if(first_valid_from(sp(:,s),1)>0) then
               nactive=nactive+1;active(nactive)=s
            else
               nfinal=nfinal+1;finals(nfinal)=s
            end if
         end do
         iter=0
         do while(nactive>0)
            iter=iter+1;offer=0;new_assign=.false.
            do i=1,nactive
               s=active(i)
               c=next_valid(sp(:,s),next(s)+1)
               if(c>0) then
                  next(s)=position_from(sp(:,s),next(s)+1,c)
                  offer(s)=c
               end if
            end do
            nrej=0;current_ip=0
            do c=1,nc
               if(.not.any(offer(active(:nactive))==c)) cycle
               if(deferred) then
                  do i=1,nactive
                     s=active(i);if(offer(s)/=c)cycle
                     if(rank_of(c_pref(:,c),s)<huge(1)) then
                        accepted(s,c)=.true.;new_assign(s,c)=.true.
                     else
                        nrej=nrej+1;rejected(nrej)=s
                     end if
                  end do
                  interrupting_flag=was_reject(c)
                  do while(count(accepted(:,c))>slots(c))
                     worst=0;worstr=-1
                     do s=1,ns
                        if(.not.accepted(s,c))cycle
                        if(rank_of(c_pref(:,c),s)>worstr) then
                           worst=s;worstr=rank_of(c_pref(:,c),s)
                        end if
                     end do
                     if(worst==0)exit
                     accepted(worst,c)=.false.;nrej=nrej+1;rejected(nrej)=worst
                     was_reject(c)=.true.
                     if(interrupting_flag .and. .not.new_assign(worst,c) .and. cons(worst)) then
                        current_ip(worst)=c
                     end if
                  end do
               else
                  do i=1,nactive
                     s=active(i);if(offer(s)/=c)cycle
                     if(rank_of(c_pref(:,c),s)>=huge(1)) then
                        nrej=nrej+1;rejected(nrej)=s
                     else if(count(accepted(:,c))<slots(c)) then
                        accepted(s,c)=.true.
                     else
                        nrej=nrej+1;rejected(nrej)=s
                     end if
                  end do
               end if
            end do
            if(any(current_ip>0)) ip_college=current_ip
            nactive=0
            do i=1,nrej
               s=rejected(i)
               if(next_valid(sp(:,s),next(s)+1)>0) then
                  nactive=nactive+1;active(nactive)=s
               else
                  nfinal=nfinal+1;finals(nfinal)=s
               end if
            end do
         end do
         assign=0
         do c=1,nc
            do s=1,ns
               if(accepted(s,c)) assign(s)=c
            end do
         end do
         if(.not.any(ip_college>0) .or. ea>=lim) exit
         deallocate(next,active,rejected,finals,offer,assign,current_ip,was_reject,accepted,new_assign)
      end do
      allocate(out%assignment(ns));out%assignment=assign
      out%pairs=pairs_from_assignment(assign);out%singles=singles_from_assignment(assign)
      out%student_prefs=sp;out%iterations=iter;out%ea_iterations=ea
      allocate(out%interrupting_pairs(2,count(ip_college>0)));k=0
      do s=1,ns
         if(ip_college(s)>0) then;k=k+1;out%interrupting_pairs(:,k)=[s,ip_college(s)];end if
      end do
   end function hri3_eadam

   function rsd(prefs,slots,priority,seed) result(res)
      integer,intent(in)::prefs(:,:),slots(:)
      integer,intent(in),optional::priority(:)
      integer(i8),intent(in),optional::seed
      type(assignment_result_t)::res
      integer :: ni,no,i,j,s,o
      integer,allocatable::prio(:),cap(:)
      type(rng_t)::rng
      ni=size(prefs,2);no=size(prefs,1)
      allocate(prio(ni));prio=[(i,i=1,ni)]
      if(present(priority)) then;prio=priority;else;if(present(seed))call rng%seed(seed);call rng%shuffle(prio);end if
      cap=slots;allocate(res%assignment(ni));res%assignment=0
      do i=1,ni
         s=prio(i)
         do j=1,no
            o=prefs(j,s)
            if(o>=1 .and. o<=size(cap)) then
               if(cap(o)>0) then;res%assignment(s)=o;cap(o)=cap(o)-1;exit;end if
            end if
         end do
      end do
      res%pairs=pairs_from_assignment(res%assignment);res%singles=singles_from_assignment(res%assignment)
      res%free_capacity=cap
   end function rsd

   function ttc_tenants(prefs,houses,priority) result(res)
      integer,intent(in)::prefs(:,:),houses(:)
      integer,intent(in),optional::priority(:)
      type(assignment_result_t)::res
      integer :: n,nh,i,s,h,owner,start,cur,k,pos
      integer,allocatable::prio(:),path(:),phouse(:)
      logical,allocatable::done_s(:),done_h(:)
      n=size(prefs,2);nh=size(houses)
      allocate(prio(n));prio=[(i,i=1,n)];if(present(priority))prio=priority
      allocate(done_s(n),done_h(nh),res%assignment(n),path(n+1),phouse(n+1))
      done_s=.false.;done_h=.false.;res%assignment=0
      do while(.not.all(done_s) .and. .not.all(done_h))
         start=0
         do i=1,n
            if(.not.done_s(prio(i))) then;start=prio(i);exit;end if
         end do
         if(start==0) exit
         cur=start;k=0
         do
            h=best_available(prefs(:,cur),done_h)
            if(h==0) then;done_s(cur)=.true.;exit;end if
            k=k+1;path(k)=cur;phouse(k)=h
            owner=0;if(h<=size(houses))owner=houses(h)
            if(owner>=1 .and. owner<=n .and. .not.done_s(owner)) then
               cur=owner
            else
               cur=start
            end if
            pos=find_in(path(:k),cur)
            if(pos>0) then
               do i=pos,k
                  s=path(i);h=phouse(i);res%assignment(s)=h;done_s(s)=.true.;done_h(h)=.true.
               end do
               exit
            end if
         end do
      end do
      res%pairs=pairs_from_assignment(res%assignment);res%singles=singles_from_assignment(res%assignment)
   end function ttc_tenants

   function ttc_school(s_pref,c_pref,slots,priority) result(res)
      integer,intent(in)::s_pref(:,:),c_pref(:,:),slots(:)
      integer,intent(in),optional::priority(:)
      type(assignment_result_t)::res
      integer :: ns,nc,i,j,s,c,k,pos,nexts
      integer,allocatable::prio(:),cap(:),path_s(:),path_c(:)
      logical,allocatable::done(:)
      ns=size(s_pref,2);nc=size(s_pref,1)
      allocate(prio(ns));prio=[(i,i=1,ns)];if(present(priority))prio=priority
      cap=slots;allocate(done(ns),res%assignment(ns),path_s(ns+1),path_c(ns+1));done=.false.;res%assignment=0
      do while(.not.all(done))
         s=0;do i=1,ns;if(.not.done(prio(i))) then;s=prio(i);exit;end if;end do
         if(s==0) exit;k=0
         do
            c=best_capacity(s_pref(:,s),cap)
            if(c==0) then;done(s)=.true.;exit;end if
            k=k+1;path_s(k)=s;path_c(k)=c
            nexts=first_unmatched_ranked(c_pref(:,c),done)
            if(nexts==0) then
               pos=k
            else
               s=nexts;pos=find_in(path_s(:k),s)
            end if
            if(pos>0) then
               do j=pos,k
                  if(cap(path_c(j))>0) then
                     res%assignment(path_s(j))=path_c(j);cap(path_c(j))=cap(path_c(j))-1;done(path_s(j))=.true.
                  end if
               end do
               exit
            end if
         end do
      end do
      res%pairs=pairs_from_assignment(res%assignment);res%singles=singles_from_assignment(res%assignment);res%free_capacity=cap
   end function ttc_school

   function ttcc_kidney(prefs,priority) result(res)
      integer,intent(in)::prefs(:,:)
      integer,intent(in),optional::priority(:)
      type(assignment_result_t)::res
      integer :: n,i,k,s,w,pos,start,bestlen,j
      integer,allocatable::prio(:),path(:),obj(:),bestp(:),besto(:)
      logical,allocatable::matched(:),wait(:)
      n=size(prefs,2);allocate(prio(n));prio=[(i,i=1,n)];if(present(priority))prio=priority
      allocate(matched(n+1),wait(n+1),res%assignment(n),path(n+1),obj(n+1),bestp(n+1),besto(n+1))
      matched=.false.;wait=.false.;matched(n+1)=.false.;wait(n+1)=.false.;res%assignment=0
      do while(any(.not.(matched(1:n).or.wait(1:n))))
         bestlen=0
         do i=1,n
            start=prio(i);if(matched(start).or.wait(start))cycle
            s=start;k=0
            do
               w=best_kidney(prefs(:,s),matched,wait)
               if(w==0) exit
               k=k+1;path(k)=s;obj(k)=w
               if(w==n+1) then
                  if(k>=bestlen) then;bestlen=k;bestp(:k)=path(:k);besto(:k)=obj(:k);end if
                  exit
               end if
               s=w;pos=find_in(path(:k),s)
               if(pos>0) then
                  do j=pos,k
                     res%assignment(path(j))=obj(j);matched(path(j))=.true.
                  end do
                  bestlen=-1;exit
               end if
            end do
            if(bestlen==-1) exit
         end do
         if(bestlen==-1) cycle
         if(bestlen<=0) exit
         wait(bestp(bestlen))=.true.
         if(bestlen>1) then
            do j=1,bestlen-1
               res%assignment(bestp(j))=besto(j);matched(bestp(j))=.true.
            end do
         end if
      end do
      res%pairs=pairs_from_assignment(res%assignment);res%singles=singles_from_assignment(res%assignment)
   end function ttcc_kidney

   function stability_check(assign,c_pref,s_pref,slots) result(bp)
      integer,intent(in)::assign(:),c_pref(:,:),s_pref(:,:),slots(:)
      integer,allocatable::bp(:,:)
      integer :: ns,nc,c,s,k,cur,rs,rc,worst,used
      integer,allocatable::tmp(:,:)
      ns=size(assign);nc=size(c_pref,2);allocate(tmp(2,ns*max(1,nc)));k=0
      do c=1,nc
         used=count(assign==c);worst=0
         if(used>=slots(c).and.slots(c)>0) then
            worst=1
            do s=1,ns
               if(assign(s)==c) worst=max(worst,rank_of(c_pref(:,c),s))
            end do
         end if
         do s=1,ns
            if(assign(s)==c) cycle
            rc=rank_of(c_pref(:,c),s);if(rc==huge(1))cycle
            if(used>=slots(c).and.rc>=worst)cycle
            cur=assign(s)
            rs=rank_of(s_pref(:,s),c)
            if(rs==huge(1))cycle
            if(cur==0 .or. rs<rank_of(s_pref(:,s),cur)) then
               k=k+1;tmp(:,k)=[c,s]
            end if
         end do
      end do
      allocate(bp(2,k));if(k>0)bp=tmp(:,:k)
   end function stability_check

   function hri_all(s_pref,c_pref,slots,max_solutions) result(set)
      integer,intent(in)::s_pref(:,:),c_pref(:,:),slots(:)
      integer,intent(in),optional::max_solutions
      type(stable_set_t)::set
      integer :: ns,nc,lim,countsol
      integer,allocatable::a(:),cap(:),store(:,:)
      ns=size(s_pref,2);nc=size(s_pref,1);lim=10000;if(present(max_solutions))lim=max_solutions
      allocate(a(ns),cap(nc),store(ns,lim));a=0;cap=0;countsol=0
      call rec(1)
      set%count=countsol;allocate(set%assignments(ns,countsol));if(countsol>0)set%assignments=store(:,:countsol)
   contains
      recursive subroutine rec(s)
         integer,intent(in)::s
         integer :: c
         integer,allocatable::bp(:,:)
         if(countsol>=lim)return
         if(s>ns) then
            bp=stability_check(a,c_pref,s_pref,slots)
            if(size(bp,2)==0) then;countsol=countsol+1;store(:,countsol)=a;end if
            return
         end if
         a(s)=0;call rec(s+1)
         do c=1,nc
            if(cap(c)>=slots(c))cycle
            if(rank_of(s_pref(:,s),c)==huge(1).or.rank_of(c_pref(:,c),s)==huge(1))cycle
            a(s)=c;cap(c)=cap(c)+1;call rec(s+1);cap(c)=cap(c)-1
         end do
         a(s)=0
      end subroutine rec
   end function hri_all

   function sri_all(pref,max_solutions) result(set)
      integer,intent(in)::pref(:,:)
      integer,intent(in),optional::max_solutions
      type(stable_set_t)::set
      integer :: n,lim,countsol
      integer,allocatable::mate(:),store(:,:)
      n=size(pref,2);lim=10000;if(present(max_solutions))lim=max_solutions
      allocate(mate(n),store(n,lim));mate=0;countsol=0
      call pair_next()
      set%count=countsol;allocate(set%assignments(n,countsol));if(countsol>0)set%assignments=store(:,:countsol)
   contains
      recursive subroutine pair_next()
         integer::i,j
         if(countsol>=lim)return
         i=0;do j=1,n;if(mate(j)==0)then;i=j;exit;end if;end do
         if(i==0) then
            if(roommate_is_stable(pref,mate)) then;countsol=countsol+1;store(:,countsol)=mate;end if
            return
         end if
         do j=i+1,n
            if(mate(j)/=0)cycle
            if(rank_of(pref(:,i),j)==huge(1).or.rank_of(pref(:,j),i)==huge(1))cycle
            mate(i)=j;mate(j)=i;call pair_next();mate(i)=0;mate(j)=0
         end do
         if(mod(n,2)==1) then;mate(i)=-1;call pair_next();mate(i)=0;end if
      end subroutine pair_next
   end function sri_all

   function hri2_couples_exact(s_pref,c_pref,slots,couple_members,couple_choice,max_nodes) result(res)
      integer,intent(in)::s_pref(:,:),c_pref(:,:),slots(:)
      integer,intent(in)::couple_members(:,:),couple_choice(:,:,:)
      integer,intent(in),optional::max_nodes
      type(assignment_result_t)::res
      ! Exact exhaustive fallback: couples are assigned jointly according to couple_choice(:,1:2,k).
      ! Stability is checked for individual resident/college blocking pairs; joint-couple blocking is
      ! enforced against listed joint alternatives in the recursive complete-state check.
      integer :: ns,nc,ncou,limit,nodes,bestscore,k
      integer,allocatable::a(:),cap(:),best(:)
      logical,allocatable::in_couple(:)
      ns=size(c_pref,1);nc=size(c_pref,2);ncou=size(couple_members,2)
      allocate(res%pairs(2,0),res%singles(0),res%free_capacity(0))
      limit=1000000;if(present(max_nodes))limit=max_nodes
      allocate(a(ns),cap(nc),best(ns),in_couple(ns));a=0;cap=0;best=0;in_couple=.false.;nodes=0;bestscore=huge(1)
      do k=1,ncou
         if(all(couple_members(:,k)>=1).and.all(couple_members(:,k)<=ns)) in_couple(couple_members(:,k))=.true.
      end do
      call rec_single(1)
      allocate(res%assignment(ns));res%assignment=best
      res%pairs=pairs_from_assignment(best);res%singles=singles_from_assignment(best)
   contains
      recursive subroutine rec_single(s)
         integer,intent(in)::s
         integer::c
         if(nodes>=limit)return
         if(s>ns) then;call rec_couple(1);return;end if
         if(in_couple(s)) then;call rec_single(s+1);return;end if
         nodes=nodes+1;a(s)=0;call rec_single(s+1)
         do c=1,nc
            if(cap(c)<slots(c)) then;a(s)=c;cap(c)=cap(c)+1;call rec_single(s+1);cap(c)=cap(c)-1;end if
         end do
         a(s)=0
      end subroutine rec_single
      recursive subroutine rec_couple(kc)
         integer,intent(in)::kc
         integer::r,c1,c2,u,v,score
         integer,allocatable::bp(:,:)
         if(nodes>=limit)return
         if(kc>ncou) then
            nodes=nodes+1;bp=stability_check(a,c_pref,s_pref,slots);if(size(bp,2)>0)return
            if(.not.couples_stable())return
            score=0
            do r=1,min(ns,size(s_pref,2))
               if(a(r)>0) score=score+rank_of(s_pref(:,r),a(r))
            end do
            if(score<bestscore)then;bestscore=score;best=a;end if
            return
         end if
         u=couple_members(1,kc);v=couple_members(2,kc);a(u)=0;a(v)=0;call rec_couple(kc+1)
         do r=1,size(couple_choice,1)
            c1=couple_choice(r,1,kc);c2=couple_choice(r,2,kc)
            if(c1<1.or.c1>nc.or.c2<1.or.c2>nc)cycle
            if(c1==c2)then;if(cap(c1)+2>slots(c1))cycle;else;if(cap(c1)>=slots(c1).or.cap(c2)>=slots(c2))cycle;end if
            a(u)=c1;a(v)=c2;cap(c1)=cap(c1)+1;cap(c2)=cap(c2)+1
            call rec_couple(kc+1);cap(c1)=cap(c1)-1;cap(c2)=cap(c2)-1
         end do
         a(u)=0;a(v)=0
      end subroutine rec_couple
      logical function couples_stable() result(ok)
         integer::kc,r,u,v,c1,c2,cur_rank
         ok=.true.
         do kc=1,ncou
            u=couple_members(1,kc);v=couple_members(2,kc);cur_rank=huge(1)
            do r=1,size(couple_choice,1)
               if(couple_choice(r,1,kc)==a(u).and.couple_choice(r,2,kc)==a(v))then
                  cur_rank=r;exit
               end if
            end do
            do r=1,min(cur_rank-1,size(couple_choice,1))
               c1=couple_choice(r,1,kc);c2=couple_choice(r,2,kc)
               if(c1<1.or.c1>nc.or.c2<1.or.c2>nc)cycle
               if(joint_college_accepts(u,v,c1,c2))then
                  ok=.false.;return
               end if
            end do
         end do
      end function couples_stable

      logical function joint_college_accepts(u,v,c1,c2) result(ok)
         integer,intent(in)::u,v,c1,c2
         integer::c,s,ncand,nnew,j,t,ru,rv
         integer,allocatable::ranks(:),ids(:)
         logical::uok,vok
         uok=.true.;vok=.true.
         do c=1,nc
            if(c/=c1 .and. c/=c2)cycle
            ncand=count(a==c)
            if(a(u)==c)ncand=ncand-1
            if(a(v)==c)ncand=ncand-1
            nnew=merge(2,1,c1==c2)
            if(c/=c1 .or. c/=c2)nnew=1
            allocate(ranks(max(1,ncand+nnew)),ids(max(1,ncand+nnew)))
            j=0
            do s=1,ns
               if(s==u.or.s==v)cycle
               if(a(s)==c)then;j=j+1;ids(j)=s;ranks(j)=rank_of(c_pref(:,c),s);end if
            end do
            if(c==c1)then;j=j+1;ids(j)=u;ranks(j)=rank_of(c_pref(:,c),u);end if
            if(c==c2)then;j=j+1;ids(j)=v;ranks(j)=rank_of(c_pref(:,c),v);end if
            if(any(ranks(:j)>=huge(1)))then;ok=.false.;return;end if
            do s=1,j-1
               do t=s+1,j
                  if(ranks(t)<ranks(s))then
                     ru=ranks(s);ranks(s)=ranks(t);ranks(t)=ru
                     rv=ids(s);ids(s)=ids(t);ids(t)=rv
                  end if
               end do
            end do
            if(c==c1)uok=any(ids(:min(j,slots(c)))==u)
            if(c==c2)vok=any(ids(:min(j,slots(c)))==v)
            deallocate(ranks,ids)
         end do
         ok=uok.and.vok
      end function joint_college_accepts
   end function hri2_couples_exact

   subroutine remove_pref(college_list,c)
      integer,intent(inout)::college_list(:);integer,intent(in)::c
      integer::i,k
      k=0
      do i=1,size(college_list)
         if(college_list(i)==c)cycle
         if(college_list(i)>0)then;k=k+1;college_list(k)=college_list(i);end if
      end do
      if(k<size(college_list))college_list(k+1:)=0
   end subroutine remove_pref
   integer function first_valid_from(x,start) result(v)
      integer,intent(in)::x(:),start;integer::i
      v=0;do i=max(1,start),size(x);if(x(i)>0)then;v=x(i);return;end if;end do
   end function first_valid_from
   integer function next_valid(x,start) result(v)
      integer,intent(in)::x(:),start;v=first_valid_from(x,start)
   end function next_valid
   integer function position_from(x,start,val) result(p)
      integer,intent(in)::x(:),start,val;integer::i
      p=size(x);do i=max(1,start),size(x);if(x(i)==val)then;p=i;return;end if;end do
   end function position_from
   integer function best_available(pref,done) result(v)
      integer,intent(in)::pref(:);logical,intent(in)::done(:);integer::i,x
      v=0;do i=1,size(pref);x=pref(i);if(x>=1.and.x<=size(done))then;if(.not.done(x))then;v=x;return;end if;end if;end do
   end function best_available
   integer function best_capacity(pref,cap) result(v)
      integer,intent(in)::pref(:),cap(:);integer::i,x
      v=0;do i=1,size(pref);x=pref(i);if(x>=1.and.x<=size(cap))then;if(cap(x)>0)then;v=x;return;end if;end if;end do
   end function best_capacity
   integer function first_unmatched_ranked(pref,done) result(v)
      integer,intent(in)::pref(:);logical,intent(in)::done(:);integer::i,x
      v=0;do i=1,size(pref);x=pref(i);if(x>=1.and.x<=size(done))then;if(.not.done(x))then;v=x;return;end if;end if;end do
   end function first_unmatched_ranked
   integer function best_kidney(pref,matched,wait) result(v)
      integer,intent(in)::pref(:);logical,intent(in)::matched(:),wait(:);integer::i,x
      v=0
      do i=1,size(pref)
         x=pref(i)
         if(x>=1 .and. x<=size(matched))then
            if(.not.matched(x) .and. .not.wait(x))then
               v=x;return
            end if
         end if
      end do
   end function best_kidney
   integer function find_in(x,v) result(p)
      integer,intent(in)::x(:),v;integer::i
      p=0;do i=1,size(x);if(x(i)==v)then;p=i;return;end if;end do
   end function find_in
   logical function roommate_is_stable(pref,mate) result(ok)
      integer,intent(in)::pref(:,:),mate(:);integer::i,j
      ok=.true.
      do i=1,size(mate)
         do j=i+1,size(mate)
            if(mate(i)==j)cycle
            if(rank_of(pref(:,i),j)<merge(rank_of(pref(:,i),mate(i)),huge(1),mate(i)>0) .and. &
               rank_of(pref(:,j),i)<merge(rank_of(pref(:,j),mate(j)),huge(1),mate(j)>0))then
               ok=.false.;return
            end if
         end do
      end do
   end function roommate_is_stable
end module matchingmarkets_mechanisms

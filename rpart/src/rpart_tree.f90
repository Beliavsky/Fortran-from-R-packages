module rpart_tree
   use rpart_kinds, only : dp
   use rpart_types
   use rpart_utils, only : finite_dp, argsort_real
   use rpart_methods, only : evaluate_node, choose_split
   implicit none
   private
   public :: grow_tree, fix_complexity, route_direction, route_observation
   public :: count_nodes, count_splits, compute_final_counts

contains

   recursive subroutine grow_tree(model, x, y1, y2, wt, obs, node, depth, sumrisk, nsplit)
      type(rpart_model), intent(in) :: model
      real(dp), intent(in) :: x(:,:), y1(:), y2(:), wt(:)
      integer, intent(in) :: obs(:), depth
      type(rpart_node), intent(inout) :: node
      real(dp), intent(out) :: sumrisk
      integer, intent(out) :: nsplit
      type(rpart_split), allocatable :: candidates(:)
      integer, allocatable :: leftobs(:), rightobs(:), stayobs(:)
      integer :: j, ncand, nleft, nright, nstay, ls, rs
      real(dp) :: alpha, iscale, tempcp, tempcp2, lrisk, rrisk

      alpha = model%control%cp * model%root_risk
      node%depth = depth
      node%nobs = size(obs)
      node%sum_wt = sum(wt(obs))
      call evaluate_node(model,y1,y2,obs,wt,node%response,node%risk)
      if (depth == 0) node%complexity = node%risk
      tempcp = min(node%risk,node%complexity)

      if (node%nobs < model%control%minsplit .or. tempcp <= alpha .or. depth >= model%control%maxdepth) then
         node%complexity = alpha
         sumrisk = node%risk; nsplit=0; return
      end if

      allocate(candidates(model%nvar)); ncand=0; iscale=0.0_dp
      do j=1,model%nvar
         call choose_split(model,x(:,j),y1,y2,obs,wt,model%ncat(j),node%risk,candidates(j))
         candidates(j)%var=j
         if(candidates(j)%improve>iscale) iscale=candidates(j)%improve
      end do
      if(iscale<=0.0_dp) then
         node%complexity=alpha;sumrisk=node%risk;nsplit=0;return
      end if
      do j=1,model%nvar
         if(candidates(j)%improve>iscale*1.0e-10_dp)then
            candidates(j)%improve=candidates(j)%improve/model%vcost(j)
            ncand=ncand+1
         else
            candidates(j)%improve=0.0_dp
         end if
      end do
      if(ncand==0)then
         node%complexity=alpha;sumrisk=node%risk;nsplit=0;return
      end if
      call retain_best(candidates, min(model%control%maxcompete+1,ncand), node%primary)
      if(.not.allocated(node%primary))then
         node%complexity=alpha;sumrisk=node%risk;nsplit=0;return
      end if
      if(size(node%primary)==0)then
         node%complexity=alpha;sumrisk=node%risk;nsplit=0;return
      end if

      if(model%control%maxsurrogate>0) call find_surrogates(model,x,obs,wt,node%primary(1),node%surrogate,node%lastsurrogate)
      call partition_obs(model,x,obs,node%primary(1),node%surrogate,node%lastsurrogate,leftobs,rightobs,stayobs)
      nleft=size(leftobs);nright=size(rightobs);nstay=size(stayobs);node%nfinal=nstay
      if(nleft==0.or.nright==0)then
         if(allocated(node%primary))deallocate(node%primary)
         if(allocated(node%surrogate))deallocate(node%surrogate)
         node%complexity=alpha;sumrisk=node%risk;nsplit=0;return
      end if

      allocate(node%left,node%right)
      node%left%id=2*node%id
      node%right%id=2*node%id+1
      node%left%complexity=tempcp-alpha
      call grow_tree(model,x,y1,y2,wt,leftobs,node%left,depth+1,lrisk,ls)

      tempcp=(node%risk-lrisk)/real(ls+1,dp)
      tempcp2=node%risk-node%left%risk
      if(tempcp<tempcp2)tempcp=tempcp2
      if(tempcp>node%complexity)tempcp=node%complexity
      node%right%complexity=tempcp-alpha
      call grow_tree(model,x,y1,y2,wt,rightobs,node%right,depth+1,rrisk,rs)

      tempcp=(node%risk-(lrisk+rrisk))/real(ls+rs+1,dp)
      if(node%right%complexity>node%left%complexity)then
         if(tempcp>node%left%complexity)then
            lrisk=node%left%risk;ls=0
            tempcp=(node%risk-(lrisk+rrisk))/real(ls+rs+1,dp)
            if(tempcp>node%right%complexity)then
               rrisk=node%right%risk;rs=0
            end if
         end if
      else if(tempcp>node%right%complexity)then
         rs=0;rrisk=node%right%risk
         tempcp=(node%risk-(lrisk+rrisk))/real(ls+rs+1,dp)
         if(tempcp>node%left%complexity)then
            lrisk=node%left%risk;ls=0
         end if
      end if
      node%complexity=(node%risk-(lrisk+rrisk))/real(ls+rs+1,dp)
      if(node%complexity<=alpha)then
         call make_leaf(node,alpha)
         sumrisk=node%risk;nsplit=0
      else
         sumrisk=lrisk+rrisk;nsplit=ls+rs+1
      end if
   end subroutine grow_tree

   subroutine make_leaf(node,alpha)
      type(rpart_node),intent(inout)::node
      real(dp),intent(in)::alpha
      if(allocated(node%left))deallocate(node%left)
      if(allocated(node%right))deallocate(node%right)
      if(allocated(node%primary))deallocate(node%primary)
      if(allocated(node%surrogate))deallocate(node%surrogate)
      node%complexity=alpha;node%nfinal=node%nobs
   end subroutine make_leaf

   subroutine retain_best(candidates,nkeep,best)
      type(rpart_split),intent(in)::candidates(:)
      integer,intent(in)::nkeep
      type(rpart_split),allocatable,intent(out)::best(:)
      integer,allocatable::idx(:)
      real(dp),allocatable::score(:)
      integer::i,k,npos
      npos=count([(candidates(i)%improve>0.0_dp,i=1,size(candidates))])
      k=min(nkeep,npos)
      if(k<=0)then;allocate(best(0));return;end if
      allocate(score(size(candidates)))
      do i=1,size(candidates);score(i)=-candidates(i)%improve;end do
      call argsort_real(score,idx)
      allocate(best(k));do i=1,k;best(i)=candidates(idx(i));end do
   end subroutine retain_best

   subroutine find_surrogates(model,x,obs,wt,primary,surrogates,lastdir)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::x(:,:),wt(:)
      integer,intent(in)::obs(:)
      type(rpart_split),intent(in)::primary
      type(rpart_split),allocatable,intent(out)::surrogates(:)
      integer,intent(out)::lastdir
      integer,allocatable::ydir(:)
      type(rpart_split),allocatable::cand(:)
      integer::k,i,j,npos,nkeep
      real(dp)::lw,rw
      allocate(ydir(size(obs)));ydir=0;lw=0.0_dp;rw=0.0_dp
      do k=1,size(obs)
         i=obs(k);ydir(k)=split_dir_value(primary,x(i,primary%var))
         if(ydir(k)==RPART_LEFT)lw=lw+wt(i)
         if(ydir(k)==RPART_RIGHT)rw=rw+wt(i)
      end do
      if(lw>rw)then;lastdir=RPART_LEFT;else if(rw>lw)then;lastdir=RPART_RIGHT;else;lastdir=0;end if
      allocate(cand(model%nvar));npos=0
      do j=1,model%nvar
         if(j==primary%var)cycle
         cand(j)%var=j;cand(j)%ncat=model%ncat(j)
         call choose_surrogate(model,x(:,j),obs,wt,ydir,lw,rw,cand(j))
         if(cand(j)%adj>1.0e-10_dp)npos=npos+1
      end do
      nkeep=min(model%control%maxsurrogate,npos)
      if(nkeep<=0)then;allocate(surrogates(0));return;end if
      call retain_best_surrogate(cand,nkeep,surrogates)
   end subroutine find_surrogates

   subroutine retain_best_surrogate(cand,nkeep,best)
      type(rpart_split),intent(in)::cand(:)
      integer,intent(in)::nkeep
      type(rpart_split),allocatable,intent(out)::best(:)
      integer,allocatable::idx(:);real(dp),allocatable::score(:)
      integer::i,k
      allocate(score(size(cand)));do i=1,size(cand);score(i)=-cand(i)%improve;end do
      call argsort_real(score,idx);allocate(best(nkeep));k=0
      do i=1,size(idx)
         if(cand(idx(i))%adj>1.0e-10_dp)then;k=k+1;best(k)=cand(idx(i));if(k==nkeep)exit;end if
      end do
   end subroutine retain_best_surrogate

   subroutine choose_surrogate(model,xcol,obs,wt,ydir,tleft,tright,s)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::xcol(:),wt(:),tleft,tright
      integer,intent(in)::obs(:),ydir(:)
      type(rpart_split),intent(inout)::s
      integer,allocatable::pos(:),ord(:),lc(:),rc(:)
      real(dp),allocatable::xx(:),lwt(:),rwt(:)
      integer::i,j,k,n,ll,lr,rl,rr,defdir,lcount,rcount,where,direction
      real(dp)::llw,lrw,rlw,rrw,agree,majority,total,common
      n=0;allocate(pos(size(obs)))
      do k=1,size(obs);i=obs(k);if(finite_dp(xcol(i)))then;n=n+1;pos(n)=k;end if;end do
      if(n==0)return
      if(s%ncat==0)then
         allocate(xx(n));do k=1,n;xx(k)=xcol(obs(pos(k)));end do;call argsort_real(xx,ord)
         ll=0;rl=0;llw=0.0_dp;rlw=0.0_dp
         do k = 1, n
            i = pos(ord(k)); j = obs(i)
            if (ydir(i) == RPART_LEFT) then
               if (wt(j) > 0.0_dp) ll = ll + 1
               llw = llw + wt(j)
            else if (ydir(i) == RPART_RIGHT) then
               if (wt(j) > 0.0_dp) rl = rl + 1
               rlw = rlw + wt(j)
            end if
         end do
         majority=max(llw,rlw);agree=majority;total=llw+rlw;lr=0;rr=0;lrw=0.0_dp;rrw=0.0_dp;where=0;direction=RPART_LEFT
         do k=1,n-1
            i=pos(ord(k));j=obs(i)
            if(ydir(i)==RPART_LEFT)then;if(wt(j)>0.0_dp)then;ll=ll-1;lr=lr+1;end if;llw=llw-wt(j);lrw=lrw+wt(j)
            else if(ydir(i)==RPART_RIGHT)then;if(wt(j)>0.0_dp)then;rl=rl-1;rr=rr+1;end if;rlw=rlw-wt(j);rrw=rrw+wt(j);end if
            if(xx(ord(k+1))>xx(ord(k)).and.(ll+rl)>=2.and.(lr+rr)>=2)then
               if(llw+rrw>agree)then;agree=llw+rrw;where=k;direction=RPART_RIGHT
               else if(lrw+rlw>agree)then;agree=lrw+rlw;where=k;direction=RPART_LEFT;end if
            end if
         end do
         if(where==0)return
         s%spoint=(xx(ord(where))+xx(ord(where+1)))/2.0_dp;s%direction=direction
         common=total
      else
         allocate(lc(s%ncat),rc(s%ncat),lwt(s%ncat),rwt(s%ncat),s%csplit(s%ncat));lc=0;rc=0;lwt=0.0_dp;rwt=0.0_dp;s%csplit=0
         do k=1,n;i=pos(k);j=obs(i);if(nint(xcol(j))<1.or.nint(xcol(j))>s%ncat)cycle
            if (ydir(i) == RPART_LEFT) then
               if (wt(j) > 0.0_dp) lc(nint(xcol(j))) = lc(nint(xcol(j))) + 1
               lwt(nint(xcol(j))) = lwt(nint(xcol(j))) + wt(j)
            else if (ydir(i) == RPART_RIGHT) then
               if (wt(j) > 0.0_dp) rc(nint(xcol(j))) = rc(nint(xcol(j))) + 1
               rwt(nint(xcol(j))) = rwt(nint(xcol(j))) + wt(j)
            end if
         end do
         llw=sum(lwt);rrw=sum(rwt);if(llw>rrw)then;defdir=RPART_LEFT;majority=llw;else;defdir=RPART_RIGHT;majority=rrw;end if
         total=llw+rrw;agree=0.0_dp;lcount=0;rcount=0
         do i=1,s%ncat
            if(lc(i)==0.and.rc(i)==0)then;s%csplit(i)=0
            else if(lwt(i)<rwt(i).or.(abs(lwt(i)-rwt(i))<=epsilon(1.0_dp)*max(1.0_dp,lwt(i),rwt(i)).and.defdir==RPART_RIGHT))then
               agree=agree+rwt(i);s%csplit(i)=RPART_RIGHT;lcount=lcount+lc(i);rcount=rcount+rc(i)
            else
               agree=agree+lwt(i);s%csplit(i)=RPART_LEFT;lcount=lcount+rc(i);rcount=rcount+lc(i)
            end if
         end do
         if(lcount<=1.or.rcount<=1)return
         common=total
      end if
      if (model%control%surrogatestyle == 0) then
         total = tleft + tright
         majority = max(tleft,tright)
      else
         total = common
         majority = min(majority,total)
      end if
      if(total<=0.0_dp.or.total<=majority)return
      s%improve=agree/total;s%adj=(s%improve-majority/total)/(1.0_dp-majority/total)
      s%count=0
   end subroutine choose_surrogate

   subroutine partition_obs(model,x,obs,primary,surrogate,lastdir,leftobs,rightobs,stayobs)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::obs(:)
      type(rpart_split),intent(in)::primary
      type(rpart_split),intent(in),optional::surrogate(:)
      integer,intent(in)::lastdir
      integer,allocatable,intent(out)::leftobs(:),rightobs(:),stayobs(:)
      integer,allocatable::l(:),r(:),z(:)
      integer::i,k,j,d,nl,nr,nz
      allocate(l(size(obs)),r(size(obs)),z(size(obs)));nl=0;nr=0;nz=0
      do k=1,size(obs)
         i=obs(k);d=split_dir_value(primary,x(i,primary%var))
         if(d==0.and.model%control%usesurrogate>0.and.present(surrogate))then
            do j=1,size(surrogate)
               d=split_dir_value(surrogate(j),x(i,surrogate(j)%var));if(d/=0)exit
            end do
         end if
         if(d==0.and.model%control%usesurrogate==2)d=lastdir
         select case(d);case(RPART_LEFT);nl=nl+1;l(nl)=i;case(RPART_RIGHT);nr=nr+1;r(nr)=i;case default;nz=nz+1;z(nz)=i;end select
      end do
      allocate(leftobs(nl),rightobs(nr),stayobs(nz));if(nl>0)leftobs=l(1:nl);if(nr>0)rightobs=r(1:nr);if(nz>0)stayobs=z(1:nz)
   end subroutine partition_obs

   integer function split_dir_value(split,x) result(d)
      type(rpart_split),intent(in)::split
      real(dp),intent(in)::x
      integer::cat
      d=0;if(.not.finite_dp(x))return
      if(split%ncat==0)then
         if(x<split%spoint)then;d=split%direction;else;d=-split%direction;end if
      else
         cat=nint(x);if(cat>=1.and.allocated(split%csplit))then;if(cat<=size(split%csplit))d=split%csplit(cat);end if
      end if
   end function split_dir_value

   integer function route_direction(model,node,xrow) result(d)
      type(rpart_model),intent(in)::model
      type(rpart_node),intent(in)::node
      real(dp),intent(in)::xrow(:)
      integer::j
      d=0
      if(.not.allocated(node%primary))return
      if(size(node%primary)==0)return
      d=split_dir_value(node%primary(1),xrow(node%primary(1)%var))
      if(d==0.and.model%control%usesurrogate>0.and.allocated(node%surrogate))then
         do j=1,size(node%surrogate);d=split_dir_value(node%surrogate(j),xrow(node%surrogate(j)%var));if(d/=0)exit;end do
      end if
      if(d==0.and.model%control%usesurrogate==2)d=node%lastsurrogate
   end function route_direction

   recursive subroutine route_observation(model,node,xrow,cp,node_out)
      type(rpart_model),intent(in)::model
      type(rpart_node),intent(in),target::node
      real(dp),intent(in)::xrow(:),cp
      type(rpart_node),pointer,intent(out)::node_out
      integer::d
      if(.not.allocated(node%left))then;node_out=>node;return;end if
      if(model%root_risk<=0.0_dp)then;node_out=>node;return;end if
      if(node%complexity/model%root_risk<=cp)then;node_out=>node;return;end if
      d=route_direction(model,node,xrow)
      if(d==RPART_LEFT)then;call route_observation(model,node%left,xrow,cp,node_out)
      else if(d==RPART_RIGHT)then;call route_observation(model,node%right,xrow,cp,node_out)
      else;node_out=>node;end if
   end subroutine route_observation

   recursive subroutine fix_complexity(node,parent_cp)
      type(rpart_node),intent(inout)::node
      real(dp),intent(in)::parent_cp
      if(node%complexity>parent_cp)node%complexity=parent_cp
      if (allocated(node%left)) then
         call fix_complexity(node%left,node%complexity)
         call fix_complexity(node%right,node%complexity)
      end if
   end subroutine fix_complexity

   recursive integer function count_nodes(node) result(n)
      type(rpart_node),intent(in)::node
      n=1;if(allocated(node%left))n=n+count_nodes(node%left)+count_nodes(node%right)
   end function count_nodes

   recursive integer function count_splits(node,cp,rootrisk) result(n)
      type(rpart_node),intent(in)::node
      real(dp),intent(in)::cp,rootrisk
      n = 0
      if(.not.allocated(node%left))return
      if(rootrisk<=0.0_dp)return
      if(node%complexity/rootrisk>cp)then
         n = 1 + count_splits(node%left,cp,rootrisk) + count_splits(node%right,cp,rootrisk)
      end if
   end function count_splits

   subroutine compute_final_counts(model,x)
      type(rpart_model),intent(inout)::model
      real(dp),intent(in)::x(:,:)
      integer::i
      call zero_final(model%root)
      do i=1,size(x,1);call increment_final(model,model%root,x(i,:));end do
   end subroutine compute_final_counts

   recursive subroutine zero_final(node)
      type(rpart_node),intent(inout)::node
      node%nfinal=0;if(allocated(node%left))then;call zero_final(node%left);call zero_final(node%right);end if
   end subroutine zero_final

   recursive subroutine increment_final(model,node,xrow)
      type(rpart_model),intent(in)::model
      type(rpart_node),intent(inout)::node
      real(dp),intent(in)::xrow(:)
      integer::d
      if(.not.allocated(node%left))then;node%nfinal=node%nfinal+1;return;end if
      d=route_direction(model,node,xrow)
      if(d==RPART_LEFT)then;call increment_final(model,node%left,xrow)
      else if(d==RPART_RIGHT)then;call increment_final(model,node%right,xrow)
      else;node%nfinal=node%nfinal+1;end if
   end subroutine increment_final

end module rpart_tree

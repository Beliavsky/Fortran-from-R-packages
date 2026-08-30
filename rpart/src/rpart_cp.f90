module rpart_cp
   use rpart_kinds, only : dp
   use rpart_types
   implicit none
   private
   public :: build_cp_table, compute_variable_importance, prune_model

contains

   subroutine build_cp_table(model)
      type(rpart_model),intent(inout)::model
      real(dp),allocatable::vals(:),uniq(:)
      integer::n,i,k
      real(dp)::risk
      n=count_nodes_local(model%root)
      allocate(vals(n+1));k=0
      call collect_cp(model%root,model%root_risk,vals,k)
      k=k+1;vals(k)=model%control%cp
      call sort_desc_unique(vals(1:k),uniq)
      if(allocated(model%cptable))deallocate(model%cptable)
      allocate(model%cptable(size(uniq)))
      do i=1,size(uniq)
         model%cptable(i)%cp=uniq(i)
         call pruned_risk_count(model%root,uniq(i),model%root_risk,risk,n)
         model%cptable(i)%nsplit=n
         if (model%root_risk > 0.0_dp) then
            model%cptable(i)%rel_error = risk/model%root_risk
         else
            model%cptable(i)%rel_error = 0.0_dp
         end if
         model%cptable(i)%xerror=0.0_dp;model%cptable(i)%xstd=0.0_dp
      end do
   end subroutine build_cp_table

   recursive subroutine collect_cp(node,rootrisk,vals,k)
      type(rpart_node),intent(in)::node
      real(dp),intent(in)::rootrisk
      real(dp),intent(inout)::vals(:)
      integer,intent(inout)::k
      if(allocated(node%left))then
         k=k+1
         if(rootrisk>0.0_dp)then;vals(k)=node%complexity/rootrisk;else;vals(k)=0.0_dp;end if
         call collect_cp(node%left,rootrisk,vals,k);call collect_cp(node%right,rootrisk,vals,k)
      end if
   end subroutine collect_cp

   subroutine sort_desc_unique(x,out)
      real(dp),intent(in)::x(:)
      real(dp),allocatable,intent(out)::out(:)
      real(dp),allocatable::z(:),tmp(:)
      integer::i,j,k
      real(dp)::v
      allocate(z(size(x)));z=x
      do i=2,size(z)
         v=z(i);j=i-1
         do while(j>=1)
            if(z(j)>=v)exit
            z(j+1)=z(j);j=j-1
         end do
         z(j+1)=v
      end do
      allocate(tmp(size(z)));k=0
      do i=1,size(z)
         if(k==0)then;k=1;tmp(k)=z(i)
         else if(abs(z(i)-tmp(k))>epsilon(1.0_dp)*max(1.0_dp,abs(z(i)),abs(tmp(k))))then;k=k+1;tmp(k)=z(i);end if
      end do
      allocate(out(k));if(k>0)out=tmp(1:k)
   end subroutine sort_desc_unique

   recursive subroutine pruned_risk_count(node,cp,rootrisk,risk,nsplit)
      type(rpart_node),intent(in)::node
      real(dp),intent(in)::cp,rootrisk
      real(dp),intent(out)::risk
      integer,intent(out)::nsplit
      real(dp)::lr,rr;integer::ls,rs
      if(.not.allocated(node%left))then
         risk=node%risk;nsplit=0
      else if(rootrisk<=0.0_dp)then
         risk=node%risk;nsplit=0
      else if(node%complexity/rootrisk<=cp)then
         risk=node%risk;nsplit=0
      else
         call pruned_risk_count(node%left,cp,rootrisk,lr,ls);call pruned_risk_count(node%right,cp,rootrisk,rr,rs)
         risk=lr+rr;nsplit=1+ls+rs
      end if
   end subroutine pruned_risk_count

   subroutine compute_variable_importance(model)
      type(rpart_model),intent(inout)::model
      if(allocated(model%variable_importance))deallocate(model%variable_importance)
      allocate(model%variable_importance(model%nvar));model%variable_importance=0.0_dp
      call add_importance(model,model%root,model%variable_importance)
   end subroutine compute_variable_importance

   recursive subroutine add_importance(model,node,imp)
      type(rpart_model),intent(in)::model
      type(rpart_node),intent(in)::node
      real(dp),intent(inout)::imp(:)
      real(dp)::base
      integer::j
      if(.not.allocated(node%left).or..not.allocated(node%primary))return
      base=node%primary(1)%improve
      if(model%method==RPART_ANOVA)base=base*node%risk
      imp(node%primary(1)%var)=imp(node%primary(1)%var)+base
      if(allocated(node%surrogate))then
         do j=1,size(node%surrogate)
            imp(node%surrogate(j)%var)=imp(node%surrogate(j)%var)+base*node%surrogate(j)%adj
         end do
      end if
      call add_importance(model,node%left,imp);call add_importance(model,node%right,imp)
   end subroutine add_importance

   subroutine prune_model(model,cp)
      type(rpart_model),intent(inout)::model
      real(dp),intent(in)::cp
      call prune_node(model%root,cp,model%root_risk)
      call build_cp_table(model)
      call compute_variable_importance(model)
      if(allocated(model%where))deallocate(model%where)
      if(allocated(model%fitted))deallocate(model%fitted)
   end subroutine prune_model

   recursive subroutine prune_node(node,cp,rootrisk)
      type(rpart_node),intent(inout)::node
      real(dp),intent(in)::cp,rootrisk
      if(.not.allocated(node%left))return
      if(node%complexity/rootrisk<=cp)then
         deallocate(node%left,node%right)
         if(allocated(node%primary))deallocate(node%primary)
         if(allocated(node%surrogate))deallocate(node%surrogate)
         node%nfinal=node%nobs
      else
         call prune_node(node%left,cp,rootrisk);call prune_node(node%right,cp,rootrisk)
      end if
   end subroutine prune_node

   recursive integer function count_nodes_local(node) result(n)
      type(rpart_node),intent(in)::node
      n=1;if(allocated(node%left))n=n+count_nodes_local(node%left)+count_nodes_local(node%right)
   end function count_nodes_local

end module rpart_cp

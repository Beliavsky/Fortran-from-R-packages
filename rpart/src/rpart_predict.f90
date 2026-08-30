module rpart_predict
   use rpart_kinds, only : dp
   use rpart_types
   use rpart_tree, only : route_direction
   implicit none
   private
   public :: rpart_predict_values, rpart_predict_class, rpart_predict_proba
   public :: rpart_predict_where, rpart_predict_full, rpart_node_path, rpart_predict_one

contains

   subroutine rpart_predict_values(model,x,pred,cp)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::x(:,:)
      real(dp),intent(out)::pred(:)
      real(dp),intent(in),optional::cp
      real(dp)::c
      real(dp),allocatable::resp(:)
      integer::i,id
      c=0.0_dp;if(present(cp))c=cp
      if(size(pred)/=size(x,1))error stop 'rpart_predict_values: wrong output size'
      do i=1,size(x,1)
         call terminal_response(model,model%root,x(i,:),c,resp,id)
         pred(i)=resp(1)
      end do
   end subroutine rpart_predict_values

   subroutine rpart_predict_class(model,x,pred,cp)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::x(:,:)
      integer,intent(out)::pred(:)
      real(dp),intent(in),optional::cp
      real(dp),allocatable::p(:)
      if(model%method/=RPART_CLASS)error stop 'rpart_predict_class: model is not classification'
      allocate(p(size(x,1)));call rpart_predict_values(model,x,p,cp);pred=nint(p)
   end subroutine rpart_predict_class

   subroutine rpart_predict_proba(model,x,prob,cp)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::x(:,:)
      real(dp),intent(out)::prob(:,:)
      real(dp),intent(in),optional::cp
      real(dp)::c,s
      real(dp),allocatable::resp(:)
      integer::i,j,id
      if(model%method/=RPART_CLASS)error stop 'rpart_predict_proba: model is not classification'
      if(size(prob,1)/=size(x,1).or.size(prob,2)/=model%nclass)error stop 'rpart_predict_proba: wrong output shape'
      c=0.0_dp;if(present(cp))c=cp
      do i=1,size(x,1)
         call terminal_response(model,model%root,x(i,:),c,resp,id)
         s=0.0_dp
         do j=1,model%nclass
            if(model%class_freq(j)>0.0_dp)then
               prob(i,j)=resp(1+j)*model%prior(j)/model%class_freq(j)
            else
               prob(i,j)=0.0_dp
            end if
            s=s+prob(i,j)
         end do
         if(s>0.0_dp)prob(i,:)=prob(i,:)/s
      end do
   end subroutine rpart_predict_proba

   subroutine rpart_predict_where(model,x,where,cp)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::x(:,:)
      integer,intent(out)::where(:)
      real(dp),intent(in),optional::cp
      real(dp)::c
      real(dp),allocatable::resp(:)
      integer::i
      if(size(where)/=size(x,1))error stop 'rpart_predict_where: wrong output size'
      c=0.0_dp;if(present(cp))c=cp
      do i=1,size(x,1);call terminal_response(model,model%root,x(i,:),c,resp,where(i));end do
   end subroutine rpart_predict_where

   subroutine rpart_predict_full(model,x,response,where,cp)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::x(:,:)
      real(dp),intent(out)::response(:,:)
      integer,intent(out),optional::where(:)
      real(dp),intent(in),optional::cp
      real(dp)::c
      real(dp),allocatable::resp(:)
      integer::i,id
      if(size(response,1)/=size(x,1).or.size(response,2)/=model%nresp)error stop 'rpart_predict_full: wrong output shape'
      c=0.0_dp;if(present(cp))c=cp
      do i=1,size(x,1)
         call terminal_response(model,model%root,x(i,:),c,resp,id);response(i,:)=resp
         if(present(where))where(i)=id
      end do
   end subroutine rpart_predict_full

   subroutine rpart_predict_one(model,xrow,cp,response,id)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::xrow(:),cp
      real(dp),allocatable,intent(out)::response(:)
      integer,intent(out),optional::id
      integer::nid
      call terminal_response(model,model%root,xrow,cp,response,nid)
      if(present(id))id=nid
   end subroutine rpart_predict_one

   recursive subroutine terminal_response(model,node,xrow,cp,response,id)
      type(rpart_model),intent(in)::model
      type(rpart_node),intent(in)::node
      real(dp),intent(in)::xrow(:),cp
      real(dp),allocatable,intent(out)::response(:)
      integer,intent(out)::id
      integer::d
      if(.not.allocated(node%left))then
         allocate(response(size(node%response)));response=node%response;id=node%id;return
      end if
      if(model%root_risk<=0.0_dp.or.node%complexity/model%root_risk<=cp)then
         allocate(response(size(node%response)));response=node%response;id=node%id;return
      end if
      d=route_direction(model,node,xrow)
      if(d==RPART_LEFT)then
         call terminal_response(model,node%left,xrow,cp,response,id)
      else if(d==RPART_RIGHT)then
         call terminal_response(model,node%right,xrow,cp,response,id)
      else
         allocate(response(size(node%response)));response=node%response;id=node%id
      end if
   end subroutine terminal_response

   subroutine rpart_node_path(model,xrow,path,cp)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::xrow(:)
      integer,allocatable,intent(out)::path(:)
      real(dp),intent(in),optional::cp
      integer,allocatable::tmp(:)
      integer::n
      real(dp)::c
      c=0.0_dp;if(present(cp))c=cp
      allocate(tmp(model%control%maxdepth+1));n=0
      call collect_path(model,model%root,xrow,c,tmp,n)
      allocate(path(n));path=tmp(1:n)
   end subroutine rpart_node_path

   recursive subroutine collect_path(model,node,xrow,cp,tmp,n)
      type(rpart_model),intent(in)::model
      type(rpart_node),intent(in)::node
      real(dp),intent(in)::xrow(:),cp
      integer,intent(inout)::tmp(:),n
      integer::d
      n=n+1;tmp(n)=node%id
      if(.not.allocated(node%left))return
      if(model%root_risk<=0.0_dp)return
      if(node%complexity/model%root_risk<=cp)return
      d=route_direction(model,node,xrow)
      if(d==RPART_LEFT)call collect_path(model,node%left,xrow,cp,tmp,n)
      if(d==RPART_RIGHT)call collect_path(model,node%right,xrow,cp,tmp,n)
   end subroutine collect_path

end module rpart_predict

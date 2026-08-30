module rpart_xpred
   use rpart_kinds, only : dp
   use rpart_types
   use rpart_fit, only : rpart_fit_subset_from_template
   use rpart_predict, only : rpart_predict_one
   use rpart_survival, only : rpart_exp_transform_right, rpart_exp_transform_startstop
   implicit none
   private
   public :: rpart_default_xpred_cp
   public :: rpart_xpred_regression, rpart_xpred_classification, rpart_xpred_poisson
   public :: rpart_xpred_survival, rpart_xpred_survival_startstop
   public :: rpart_xpred_full

contains

   subroutine rpart_default_xpred_cp(model,cp)
      type(rpart_model),intent(in)::model
      real(dp),allocatable,intent(out)::cp(:)
      integer::j,n
      n=size(model%cptable)
      allocate(cp(n))
      if(n==0)return
      cp(1)=(1.0_dp+model%cptable(1)%cp)/2.0_dp
      do j=2,n
         cp(j)=sqrt(max(0.0_dp,model%cptable(j-1)%cp*model%cptable(j)%cp))
      end do
   end subroutine rpart_default_xpred_cp

   subroutine rpart_xpred_regression(model,x,y,groups,pred,cp,weights,stat)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::x(:,:),y(:)
      integer,intent(in)::groups(:)
      real(dp),intent(out)::pred(:,:)
      real(dp),intent(in),optional::cp(:),weights(:)
      integer,intent(out),optional::stat
      real(dp),allocatable::y2(:),wt(:),cpuse(:),full(:,:,:)
      integer::s
      if(model%method/=RPART_ANOVA.or.size(y)/=size(x,1))then
         if(present(stat))stat=1;return
      end if
      allocate(y2(size(y)));y2=0.0_dp
      call prepare_xpred(model,x,groups,pred,cp,weights,wt,cpuse,s)
      if(s/=0)then;if(present(stat))stat=s;return;end if
      allocate(full(size(x,1),size(cpuse),model%nresp))
      call xpred_core(model,x,y,y2,wt,groups,cpuse,full,s)
      if(s==0)pred=full(:,:,1)
      if(present(stat))stat=s
   end subroutine rpart_xpred_regression

   subroutine rpart_xpred_classification(model,x,y,groups,pred,cp,weights,stat)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::y(:),groups(:)
      real(dp),intent(out)::pred(:,:)
      real(dp),intent(in),optional::cp(:),weights(:)
      integer,intent(out),optional::stat
      real(dp),allocatable::y1(:),y2(:),wt(:),cpuse(:),full(:,:,:)
      integer::s
      if(model%method/=RPART_CLASS.or.size(y)/=size(x,1).or.any(y<1).or.any(y>model%nclass))then
         if(present(stat))stat=1;return
      end if
      allocate(y1(size(y)),y2(size(y)));y1=real(y,dp);y2=0.0_dp
      call prepare_xpred(model,x,groups,pred,cp,weights,wt,cpuse,s)
      if(s/=0)then;if(present(stat))stat=s;return;end if
      allocate(full(size(x,1),size(cpuse),model%nresp))
      call xpred_core(model,x,y1,y2,wt,groups,cpuse,full,s)
      if(s==0)pred=full(:,:,1)
      if(present(stat))stat=s
   end subroutine rpart_xpred_classification

   subroutine rpart_xpred_poisson(model,x,time,event,groups,pred,cp,weights,stat)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::x(:,:),time(:),event(:)
      integer,intent(in)::groups(:)
      real(dp),intent(out)::pred(:,:)
      real(dp),intent(in),optional::cp(:),weights(:)
      integer,intent(out),optional::stat
      real(dp),allocatable::wt(:),cpuse(:),full(:,:,:)
      integer::s
      if(model%method/=RPART_POISSON.or.size(time)/=size(x,1).or.size(event)/=size(x,1))then
         if(present(stat))stat=1;return
      end if
      call prepare_xpred(model,x,groups,pred,cp,weights,wt,cpuse,s)
      if(s/=0)then;if(present(stat))stat=s;return;end if
      allocate(full(size(x,1),size(cpuse),model%nresp))
      call xpred_core(model,x,time,event,wt,groups,cpuse,full,s)
      if(s==0)pred=full(:,:,1)
      if(present(stat))stat=s
   end subroutine rpart_xpred_poisson

   subroutine rpart_xpred_survival(model,x,time,status,groups,pred,cp,weights,offset,stat)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::x(:,:),time(:),status(:)
      integer,intent(in)::groups(:)
      real(dp),intent(out)::pred(:,:)
      real(dp),intent(in),optional::cp(:),weights(:),offset(:)
      integer,intent(out),optional::stat
      real(dp),allocatable::newtime(:),wt(:),cpuse(:),full(:,:,:)
      integer::s
      if(model%method/=RPART_EXP.or.size(time)/=size(x,1).or.size(status)/=size(x,1))then
         if(present(stat))stat=1;return
      end if
      allocate(newtime(size(time)))
      call rpart_exp_transform_right(time,status,newtime,offset,s)
      if(s/=0)then;if(present(stat))stat=10+s;return;end if
      call prepare_xpred(model,x,groups,pred,cp,weights,wt,cpuse,s)
      if(s/=0)then;if(present(stat))stat=s;return;end if
      allocate(full(size(x,1),size(cpuse),model%nresp))
      call xpred_core(model,x,newtime,status,wt,groups,cpuse,full,s)
      if(s==0)pred=full(:,:,1)
      if(present(stat))stat=s
   end subroutine rpart_xpred_survival

   subroutine rpart_xpred_survival_startstop(model,x,start,stop,status,groups,pred,cp,weights,offset,stat)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::x(:,:),start(:),stop(:),status(:)
      integer,intent(in)::groups(:)
      real(dp),intent(out)::pred(:,:)
      real(dp),intent(in),optional::cp(:),weights(:),offset(:)
      integer,intent(out),optional::stat
      real(dp),allocatable::newtime(:),wt(:),cpuse(:),full(:,:,:)
      integer::s
      if(model%method/=RPART_EXP.or.size(start)/=size(x,1).or.size(stop)/=size(x,1).or. &
         size(status)/=size(x,1))then
         if(present(stat))stat=1;return
      end if
      allocate(newtime(size(stop)))
      call rpart_exp_transform_startstop(start,stop,status,newtime,offset,s)
      if(s/=0)then;if(present(stat))stat=10+s;return;end if
      call prepare_xpred(model,x,groups,pred,cp,weights,wt,cpuse,s)
      if(s/=0)then;if(present(stat))stat=s;return;end if
      allocate(full(size(x,1),size(cpuse),model%nresp))
      call xpred_core(model,x,newtime,status,wt,groups,cpuse,full,s)
      if(s==0)pred=full(:,:,1)
      if(present(stat))stat=s
   end subroutine rpart_xpred_survival_startstop

   subroutine rpart_xpred_full(model,x,y1,y2,groups,cp,pred,weights,stat)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::x(:,:),y1(:),y2(:),cp(:)
      integer,intent(in)::groups(:)
      real(dp),intent(out)::pred(:,:,:)
      real(dp),intent(in),optional::weights(:)
      integer,intent(out),optional::stat
      real(dp),allocatable::wt(:)
      integer::s
      if(size(y1)/=size(x,1).or.size(y2)/=size(x,1).or.size(pred,1)/=size(x,1).or. &
         size(pred,2)/=size(cp).or.size(pred,3)/=model%nresp)then
         if(present(stat))stat=1;return
      end if
      call prepare_weights(size(x,1),weights,wt,s)
      if(s==0)call xpred_core(model,x,y1,y2,wt,groups,cp,pred,s)
      if(present(stat))stat=s
   end subroutine rpart_xpred_full

   subroutine prepare_xpred(model,x,groups,pred,cp,weights,wt,cpuse,stat)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::groups(:)
      real(dp),intent(out)::pred(:,:)
      real(dp),intent(in),optional::cp(:),weights(:)
      real(dp),allocatable,intent(out)::wt(:),cpuse(:)
      integer,intent(out)::stat
      stat=0
      if(size(x,2)/=model%nvar.or.size(groups)/=size(x,1))then;stat=2;return;end if
      if(present(cp))then
         allocate(cpuse(size(cp)));cpuse=cp
      else
         call rpart_default_xpred_cp(model,cpuse)
      end if
      if(size(pred,1)/=size(x,1).or.size(pred,2)/=size(cpuse))then;stat=3;return;end if
      if(size(cpuse)==0.or.any(cpuse<0.0_dp))then;stat=4;return;end if
      call prepare_weights(size(x,1),weights,wt,stat)
   end subroutine prepare_xpred

   subroutine prepare_weights(n,weights,wt,stat)
      integer,intent(in)::n
      real(dp),intent(in),optional::weights(:)
      real(dp),allocatable,intent(out)::wt(:)
      integer,intent(out)::stat
      stat=0;allocate(wt(n));wt=1.0_dp
      if(present(weights))then
         if(size(weights)/=n.or.any(weights<0.0_dp))then;stat=5;return;end if
         wt=weights
      end if
      if(sum(wt)<=0.0_dp)stat=6
   end subroutine prepare_weights

   subroutine xpred_core(model,x,y1,y2,wt,groups,cp,pred,stat)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::x(:,:),y1(:),y2(:),wt(:),cp(:)
      integer,intent(in)::groups(:)
      real(dp),intent(out)::pred(:,:,:)
      integer,intent(out)::stat
      type(rpart_model)::sub
      integer,allocatable::train(:)
      real(dp),allocatable::resp(:)
      real(dp)::trainwt,cpscale,totalwt
      integer::f,i,j,k,nfold,ntr
      stat=0
      if(size(groups)/=size(x,1).or.any(groups<1))then;stat=7;return;end if
      nfold=maxval(groups)
      if(nfold<2)then;stat=8;return;end if
      do f=1,nfold
         if(count(groups==f)==0)then;stat=9;return;end if
      end do
      pred=0.0_dp;totalwt=sum(wt)
      do f=1,nfold
         ntr=count(groups/=f)
         if(ntr==0)then;stat=10;return;end if
         allocate(train(ntr));k=0
         do i=1,size(groups)
            if(groups(i)/=f)then;k=k+1;train(k)=i;end if
         end do
         call rpart_fit_subset_from_template(model,x,y1,y2,wt,train,sub)
         trainwt=sum(wt(train))
         if(sub%root_risk>0.0_dp)then
            cpscale=model%root_risk*(trainwt/totalwt)/sub%root_risk
         else
            cpscale=1.0_dp
         end if
         do i=1,size(groups)
            if(groups(i)/=f)cycle
            do j=1,size(cp)
               call rpart_predict_one(sub,x(i,:),cp(j)*cpscale,resp)
               pred(i,j,1:size(resp))=resp
            end do
         end do
         deallocate(train)
      end do
   end subroutine xpred_core

end module rpart_xpred

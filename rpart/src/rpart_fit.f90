module rpart_fit
   use rpart_kinds, only : dp, i8
   use rpart_types
   use rpart_utils, only : validate_control, shuffle_int
   use rpart_methods, only : init_classification, init_poisson, evaluate_node, point_error
   use rpart_tree, only : grow_tree, fix_complexity, compute_final_counts
   use rpart_cp, only : build_cp_table, compute_variable_importance
   use rpart_predict, only : rpart_predict_values, rpart_predict_where, rpart_predict_one
   use rpart_survival, only : rpart_exp_transform_right, rpart_exp_transform_startstop
   implicit none
   private
   public :: rpart_fit_regression, rpart_fit_classification, rpart_fit_poisson
   public :: rpart_fit_survival, rpart_fit_survival_startstop
   public :: rpart_fit_subset_from_template

contains

   subroutine rpart_fit_regression(x,y,model,control,weights,ncat,cost,xgroups,seed,stat)
      real(dp),intent(in)::x(:,:),y(:)
      type(rpart_model),intent(out)::model
      type(rpart_control),intent(in),optional::control
      real(dp),intent(in),optional::weights(:),cost(:)
      integer,intent(in),optional::ncat(:),xgroups(:)
      integer(i8),intent(in),optional::seed
      integer,intent(out),optional::stat
      real(dp),allocatable::wt(:),y2(:)
      integer::s
      call setup_common(x,model,control,weights,ncat,cost,wt,s)
      if(s/=0.or.size(y)/=size(x,1))then;if(s==0)s=10;if(present(stat))stat=s;return;end if
      model%method=RPART_ANOVA;model%nresp=1
      allocate(y2(size(y)));y2=0.0_dp
      call fit_core(model,x,y,y2,wt)
      call maybe_xval(model,x,y,y2,wt,xgroups,seed)
      if(present(stat))stat=0
   end subroutine rpart_fit_regression

   subroutine rpart_fit_classification(x,y,model,control,weights,ncat,cost,prior,loss,split_rule,xgroups,seed,stat)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::y(:)
      type(rpart_model),intent(out)::model
      type(rpart_control),intent(in),optional::control
      real(dp),intent(in),optional::weights(:),cost(:),prior(:),loss(:,:)
      integer,intent(in),optional::ncat(:),split_rule,xgroups(:)
      integer(i8),intent(in),optional::seed
      integer,intent(out),optional::stat
      real(dp),allocatable::wt(:),y1(:),y2(:)
      integer::s
      call setup_common(x,model,control,weights,ncat,cost,wt,s)
      if(s/=0.or.size(y)/=size(x,1))then;if(s==0)s=10;if(present(stat))stat=s;return;end if
      model%method=RPART_CLASS
      call init_classification(model,y,wt,prior,loss,split_rule,s)
      if(s/=0)then;if(present(stat))stat=20+s;return;end if
      allocate(y1(size(y)),y2(size(y)));y1=real(y,dp);y2=0.0_dp
      call fit_core(model,x,y1,y2,wt)
      call maybe_xval(model,x,y1,y2,wt,xgroups,seed)
      if(present(stat))stat=0
   end subroutine rpart_fit_classification

   subroutine rpart_fit_poisson(x,time,event,model,control,weights,ncat,cost,shrink,error_method,xgroups,seed,stat)
      real(dp),intent(in)::x(:,:),time(:),event(:)
      type(rpart_model),intent(out)::model
      type(rpart_control),intent(in),optional::control
      real(dp),intent(in),optional::weights(:),cost(:),shrink
      integer,intent(in),optional::ncat(:),error_method,xgroups(:)
      integer(i8),intent(in),optional::seed
      integer,intent(out),optional::stat
      real(dp),allocatable::wt(:)
      integer::s
      call setup_common(x,model,control,weights,ncat,cost,wt,s)
      if(s/=0.or.size(time)/=size(x,1).or.size(event)/=size(x,1))then
         if(s==0)s=10;if(present(stat))stat=s;return
      end if
      model%method=RPART_POISSON
      call init_poisson(model,time,event,wt,shrink,error_method,s)
      if(s/=0)then;if(present(stat))stat=30+s;return;end if
      call fit_core(model,x,time,event,wt)
      call maybe_xval(model,x,time,event,wt,xgroups,seed)
      if(present(stat))stat=0
   end subroutine rpart_fit_poisson

   subroutine rpart_fit_survival(x,time,status,model,control,weights,ncat,cost,shrink,error_method,offset,xgroups,seed,stat)
      real(dp),intent(in)::x(:,:),time(:),status(:)
      type(rpart_model),intent(out)::model
      type(rpart_control),intent(in),optional::control
      real(dp),intent(in),optional::weights(:),cost(:),shrink,offset(:)
      integer,intent(in),optional::ncat(:),error_method,xgroups(:)
      integer(i8),intent(in),optional::seed
      integer,intent(out),optional::stat
      real(dp),allocatable::wt(:),newtime(:)
      integer::s
      call setup_common(x,model,control,weights,ncat,cost,wt,s)
      if(s/=0.or.size(time)/=size(x,1).or.size(status)/=size(x,1))then
         if(s==0)s=10;if(present(stat))stat=s;return
      end if
      allocate(newtime(size(time)))
      call rpart_exp_transform_right(time,status,newtime,offset,s)
      if(s/=0)then;if(present(stat))stat=40+s;return;end if
      model%method=RPART_EXP
      call init_poisson(model,newtime,status,wt,shrink,error_method,s)
      if(s/=0)then;if(present(stat))stat=50+s;return;end if
      model%method=RPART_EXP
      call fit_core(model,x,newtime,status,wt)
      call maybe_xval(model,x,newtime,status,wt,xgroups,seed)
      if(present(stat))stat=0
   end subroutine rpart_fit_survival

   subroutine rpart_fit_survival_startstop(x,start,stop,status,model,control,weights,ncat,cost,shrink,error_method, &
                                            offset,xgroups,seed,stat)
      real(dp),intent(in)::x(:,:),start(:),stop(:),status(:)
      type(rpart_model),intent(out)::model
      type(rpart_control),intent(in),optional::control
      real(dp),intent(in),optional::weights(:),cost(:),shrink,offset(:)
      integer,intent(in),optional::ncat(:),error_method,xgroups(:)
      integer(i8),intent(in),optional::seed
      integer,intent(out),optional::stat
      real(dp),allocatable::wt(:),newtime(:)
      integer::s
      call setup_common(x,model,control,weights,ncat,cost,wt,s)
      if(s/=0.or.size(start)/=size(x,1).or.size(stop)/=size(x,1).or.size(status)/=size(x,1))then
         if(s==0)s=10;if(present(stat))stat=s;return
      end if
      allocate(newtime(size(stop)))
      call rpart_exp_transform_startstop(start,stop,status,newtime,offset,s)
      if(s/=0)then;if(present(stat))stat=40+s;return;end if
      model%method=RPART_EXP
      call init_poisson(model,newtime,status,wt,shrink,error_method,s)
      if(s/=0)then;if(present(stat))stat=50+s;return;end if
      model%method=RPART_EXP
      call fit_core(model,x,newtime,status,wt)
      call maybe_xval(model,x,newtime,status,wt,xgroups,seed)
      if(present(stat))stat=0
   end subroutine rpart_fit_survival_startstop

   subroutine setup_common(x,model,control,weights,ncat,cost,wt,stat)
      real(dp),intent(in)::x(:,:)
      type(rpart_model),intent(out)::model
      type(rpart_control),intent(in),optional::control
      real(dp),intent(in),optional::weights(:),cost(:)
      integer,intent(in),optional::ncat(:)
      real(dp),allocatable,intent(out)::wt(:)
      integer,intent(out)::stat
      integer::j,n
      n=size(x,1);stat=0
      if(n<1.or.size(x,2)<1)then;stat=1;return;end if
      model%nvar=size(x,2)
      if(present(control))model%control=control
      call validate_control(model%control,stat);if(stat/=0)return
      allocate(wt(n))
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            stat = 2
            return
         end if
         wt = weights
      else
         wt = 1.0_dp
      end if
      if(sum(wt)<=0.0_dp)then;stat=3;return;end if
      allocate(model%ncat(model%nvar));model%ncat=0
      if(present(ncat))then;if(size(ncat)/=model%nvar.or.any(ncat<0))then;stat=4;return;end if;model%ncat=ncat;end if
      allocate(model%vcost(model%nvar));model%vcost=1.0_dp
      if(present(cost))then;if(size(cost)/=model%nvar.or.any(cost<=0.0_dp))then;stat=5;return;end if;model%vcost=cost;end if
      do j=1,model%nvar
         if(model%ncat(j)==1)model%ncat(j)=0
      end do
      model%total_weight=sum(wt)
   end subroutine setup_common

   subroutine fit_core(model,x,y1,y2,wt)
      type(rpart_model),intent(inout)::model
      real(dp),intent(in)::x(:,:),y1(:),y2(:),wt(:)
      integer,allocatable::obs(:)
      real(dp),allocatable::resp(:)
      real(dp)::risk,sumrisk
      integer::i,nsplit
      allocate(obs(size(x,1)));obs=[(i,i=1,size(x,1))]
      call evaluate_node(model,y1,y2,obs,wt,resp,risk)
      model%root_risk=risk
      allocate(model%root);model%root%id=1;model%root%complexity=risk
      call grow_tree(model,x,y1,y2,wt,obs,model%root,0,sumrisk,nsplit)
      call fix_complexity(model%root,model%root%complexity)
      call build_cp_table(model)
      call compute_variable_importance(model)
      call compute_final_counts(model,x)
      allocate(model%where(size(x,1)),model%fitted(size(x,1)))
      call rpart_predict_where(model,x,model%where,0.0_dp)
      call rpart_predict_values(model,x,model%fitted,0.0_dp)
   end subroutine fit_core

   subroutine maybe_xval(model,x,y1,y2,wt,xgroups,seed)
      type(rpart_model),intent(inout)::model
      real(dp),intent(in)::x(:,:),y1(:),y2(:),wt(:)
      integer,intent(in),optional::xgroups(:)
      integer(i8),intent(in),optional::seed
      integer,allocatable::grp(:),idx(:)
      integer::i,nfold
      integer(i8)::state
      if(model%control%xval<=1.or.size(model%cptable)==0)return
      allocate(grp(size(x,1)))
      if(present(xgroups))then
         if(size(xgroups)/=size(x,1))return
         grp=xgroups;nfold=maxval(grp)
      else
         nfold=model%control%xval;allocate(idx(size(x,1)));idx=[(i,i=1,size(x,1))]
         state=123456789_i8;if(present(seed))state=seed
         call shuffle_int(idx,state)
         do i=1,size(idx);grp(idx(i))=1+mod(i-1,nfold);end do
      end if
      if(nfold>1)call cross_validate(model,x,y1,y2,wt,grp,nfold)
   end subroutine maybe_xval

   subroutine cross_validate(model,x,y1,y2,wt,grp,nfold)
      type(rpart_model),intent(inout)::model
      real(dp),intent(in)::x(:,:),y1(:),y2(:),wt(:)
      integer,intent(in)::grp(:),nfold
      type(rpart_model)::sub
      integer,allocatable::train(:)
      real(dp),allocatable::cpcmp(:),xr(:),xs(:),resp(:)
      real(dp)::err,totalwt,cpscale,trainwt
      integer::f,i,j,k,ntr
      allocate(cpcmp(size(model%cptable)),xr(size(model%cptable)),xs(size(model%cptable)))
      cpcmp(1)=10.0_dp*model%cptable(1)%cp
      do j=2,size(cpcmp);cpcmp(j)=sqrt(max(0.0_dp,model%cptable(j-1)%cp*model%cptable(j)%cp));end do
      xr=0.0_dp;xs=0.0_dp;totalwt=sum(wt)
      do f=1,nfold
         ntr=count(grp/=f);if(ntr==0)cycle;allocate(train(ntr));k=0
         do i=1,size(grp);if(grp(i)/=f)then;k=k+1;train(k)=i;end if;end do
         call rpart_fit_subset_from_template(model,x,y1,y2,wt,train,sub)
         trainwt = sum(wt(train))
         if (sub%root_risk > 0.0_dp) then
            cpscale = model%root_risk*(trainwt/totalwt)/sub%root_risk
         else
            cpscale = 1.0_dp
         end if
         do i=1,size(grp)
            if(grp(i)/=f)cycle
            do j=1,size(cpcmp)
               call rpart_predict_one(sub,x(i,:),cpcmp(j)*cpscale,resp)
               err=point_error(model,y1(i),y2(i),resp(1))
               xr(j)=xr(j)+wt(i)*err;xs(j)=xs(j)+wt(i)*err*err
            end do
         end do
         deallocate(train)
      end do
      do j=1,size(model%cptable)
         if(model%root_risk>0.0_dp)then
            model%cptable(j)%xerror=xr(j)/model%root_risk
            model%cptable(j)%xstd=sqrt(max(0.0_dp,xs(j)-xr(j)*xr(j)/totalwt))/model%root_risk
         end if
      end do
   end subroutine cross_validate

   subroutine rpart_fit_subset_from_template(template,x,y1,y2,wt,obs,sub)
      type(rpart_model),intent(in)::template
      real(dp),intent(in)::x(:,:),y1(:),y2(:),wt(:)
      integer,intent(in)::obs(:)
      type(rpart_model),intent(out)::sub
      real(dp),allocatable::resp(:)
      real(dp)::risk,sumrisk,d,t
      integer::nsplit
      sub%method=template%method;sub%nvar=template%nvar;sub%nclass=template%nclass;sub%nresp=template%nresp
      sub%poisson_method=template%poisson_method;sub%split_rule=template%split_rule;sub%poisson_shrink=template%poisson_shrink
      sub%control=template%control;sub%control%xval=0
      allocate(sub%ncat(template%nvar),sub%vcost(template%nvar));sub%ncat=template%ncat;sub%vcost=template%vcost
      if(template%method==RPART_CLASS)then
         allocate(sub%prior(template%nclass),sub%altered_prior(template%nclass),sub%class_freq(template%nclass), &
                  sub%loss(template%nclass,template%nclass))
         sub%prior=template%prior;sub%altered_prior=template%altered_prior;sub%class_freq=template%class_freq;sub%loss=template%loss
      else if(template%method==RPART_POISSON.or.template%method==RPART_EXP)then
         if(sub%poisson_shrink>0.0_dp)then
            d=sum(y2(obs)*wt(obs));t=sum(y1(obs)*wt(obs));sub%poisson_alpha=1.0_dp/(sub%poisson_shrink**2)
            if(d>0.0_dp.and.t>0.0_dp)sub%poisson_beta=sub%poisson_alpha/(d/t)
         end if
      end if
      sub%total_weight = sum(wt(obs))
      call evaluate_node(sub,y1,y2,obs,wt,resp,risk)
      sub%root_risk = risk
      if (risk > 0.0_dp .and. template%total_weight > 0.0_dp) then
         sub%control%cp = template%control%cp * template%root_risk * &
                          (sub%total_weight/template%total_weight) / risk
      end if
      allocate(sub%root)
      sub%root%id = 1
      sub%root%complexity = risk
      call grow_tree(sub,x,y1,y2,wt,obs,sub%root,0,sumrisk,nsplit);call fix_complexity(sub%root,sub%root%complexity)
   end subroutine rpart_fit_subset_from_template

end module rpart_fit

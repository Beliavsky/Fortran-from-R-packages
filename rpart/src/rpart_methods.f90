module rpart_methods
   use rpart_kinds, only : dp
   use rpart_types
   use rpart_utils, only : finite_dp, argsort_real
   implicit none
   private
   public :: init_classification, init_poisson, evaluate_node, point_error
   public :: choose_split

contains

   subroutine init_classification(model, yclass, wt, prior_in, loss_in, split_rule, stat)
      type(rpart_model), intent(inout) :: model
      integer, intent(in) :: yclass(:)
      real(dp), intent(in) :: wt(:)
      real(dp), intent(in), optional :: prior_in(:), loss_in(:,:)
      integer, intent(in), optional :: split_rule
      integer, intent(out), optional :: stat
      integer :: i, s
      real(dp) :: sw, denom
      s = 0
      if (size(yclass) /= size(wt) .or. size(yclass) == 0) then
         s = 1; if (present(stat)) stat=s; return
      end if
      model%nclass = maxval(yclass)
      if (model%nclass < 2 .or. any(yclass < 1)) then
         s = 2; if (present(stat)) stat=s; return
      end if
      model%nresp = model%nclass + 2
      model%split_rule = RPART_GINI
      if (present(split_rule)) model%split_rule = split_rule
      if (model%split_rule /= RPART_GINI .and. model%split_rule /= RPART_INFORMATION) then
         s = 3; if (present(stat)) stat=s; return
      end if
      allocate(model%class_freq(model%nclass), model%prior(model%nclass), &
               model%altered_prior(model%nclass), model%loss(model%nclass,model%nclass))
      model%class_freq = 0.0_dp
      do i = 1, size(yclass)
         model%class_freq(yclass(i)) = model%class_freq(yclass(i)) + wt(i)
      end do
      sw = sum(wt)
      if (sw <= 0.0_dp) then
         s = 4; if (present(stat)) stat=s; return
      end if
      model%class_freq = model%class_freq / sw
      if (present(prior_in)) then
         if (size(prior_in) /= model%nclass .or. any(prior_in < 0.0_dp) .or. &
             abs(sum(prior_in)-1.0_dp) > 1.0e-10_dp) then
            s=5; if (present(stat)) stat=s; return
         end if
         model%prior = prior_in
      else
         model%prior = model%class_freq
      end if
      if (present(loss_in)) then
         if (size(loss_in,1) /= model%nclass .or. size(loss_in,2) /= model%nclass .or. &
             any(loss_in < 0.0_dp)) then
            s=6; if (present(stat)) stat=s; return
         end if
         model%loss = loss_in
      else
         model%loss = 1.0_dp
         do i=1,model%nclass
            model%loss(i,i)=0.0_dp
         end do
      end if
      do i=1,model%nclass
         if (abs(model%loss(i,i)) > 1.0e-14_dp .or. sum(model%loss(i,:)) <= 0.0_dp) then
            s=7; if (present(stat)) stat=s; return
         end if
      end do
      denom = 0.0_dp
      do i=1,model%nclass
         denom = denom + model%prior(i)*sum(model%loss(i,:))
      end do
      do i=1,model%nclass
         if (model%class_freq(i) > 0.0_dp) then
            model%altered_prior(i) = model%prior(i)*sum(model%loss(i,:)) / &
                                      (max(denom,tiny(1.0_dp))*model%class_freq(i))
         else
            model%altered_prior(i) = 0.0_dp
         end if
      end do
      ! model%prior stores the user prior; risk rescaling is prior/frequency on demand.
      if (present(stat)) stat=s
   end subroutine init_classification

   subroutine init_poisson(model, time, event, wt, shrink, error_method, stat)
      type(rpart_model), intent(inout) :: model
      real(dp), intent(in) :: time(:), event(:), wt(:)
      real(dp), intent(in), optional :: shrink
      integer, intent(in), optional :: error_method
      integer, intent(out), optional :: stat
      real(dp) :: cv, total_event, total_time
      integer :: s
      s=0
      if (size(time)/=size(event) .or. size(time)/=size(wt) .or. any(time<=0.0_dp) .or. any(event<0.0_dp)) then
         s=1; if(present(stat)) stat=s; return
      end if
      model%poisson_shrink=1.0_dp
      if(present(shrink)) model%poisson_shrink=shrink
      if(model%poisson_shrink<0.0_dp) then
         s=2; if(present(stat)) stat=s; return
      end if
      model%poisson_method=RPART_DEVIANCE
      if(present(error_method)) model%poisson_method=error_method
      if(model%poisson_method/=RPART_DEVIANCE .and. model%poisson_method/=RPART_SQRT) then
         s=3; if(present(stat)) stat=s; return
      end if
      total_event=dot_product(event,wt)
      total_time=dot_product(time,wt)
      cv=model%poisson_shrink
      if(cv<=0.0_dp .or. total_event<=0.0_dp) then
         model%poisson_alpha=0.0_dp
         model%poisson_beta=0.0_dp
      else
         model%poisson_alpha=1.0_dp/(cv*cv)
         model%poisson_beta=model%poisson_alpha/(total_event/total_time)
      end if
      model%nresp=2
      if(present(stat)) stat=s
   end subroutine init_poisson

   subroutine evaluate_node(model, y1, y2, obs, wt, response, risk)
      type(rpart_model), intent(in) :: model
      real(dp), intent(in) :: y1(:), y2(:), wt(:)
      integer, intent(in) :: obs(:)
      real(dp), allocatable, intent(out) :: response(:)
      real(dp), intent(out) :: risk
      integer :: i, j, pred, k
      real(dp) :: sw, mu, d, t, lambda, dev, tmp, nodeprob, rr
      real(dp), allocatable :: cnt(:)
      select case(model%method)
      case(RPART_ANOVA)
         allocate(response(1)); sw=0.0_dp; mu=0.0_dp
         do k=1,size(obs)
            i=obs(k); sw=sw+wt(i); mu=mu+wt(i)*y1(i)
         end do
         if(sw>0.0_dp) mu=mu/sw
         risk=0.0_dp
         do k=1,size(obs)
            i=obs(k); risk=risk+wt(i)*(y1(i)-mu)**2
         end do
         response(1)=mu
      case(RPART_CLASS)
         allocate(response(model%nclass+2),cnt(model%nclass)); cnt=0.0_dp
         do k=1,size(obs)
            i=obs(k); j=nint(y1(i)); cnt(j)=cnt(j)+wt(i)
         end do
         pred=1; risk=huge(1.0_dp)
         do i=1,model%nclass
            rr=0.0_dp
            do j=1,model%nclass
               if(model%class_freq(j)>0.0_dp) rr=rr+cnt(j)*model%loss(j,i)*model%prior(j)/model%class_freq(j)
            end do
            if(rr<risk) then
               risk=rr; pred=i
            end if
         end do
         nodeprob=0.0_dp
         do j=1,model%nclass
            if(model%class_freq(j)>0.0_dp) nodeprob=nodeprob+cnt(j)*model%prior(j)/model%class_freq(j)
         end do
         response(1)=real(pred,dp)
         response(2:1+model%nclass)=cnt
         response(model%nclass+2)=nodeprob
      case(RPART_POISSON,RPART_EXP)
         allocate(response(2)); d=0.0_dp; t=0.0_dp
         do k=1,size(obs)
            i=obs(k); d=d+y2(i)*wt(i); t=t+y1(i)*wt(i)
         end do
         lambda=(d+model%poisson_alpha)/(t+model%poisson_beta)
         dev=0.0_dp
         do k=1,size(obs)
            i=obs(k); tmp=y2(i)
            dev=dev-(lambda*y1(i)-tmp)*wt(i)
            if(tmp>0.0_dp) dev=dev+tmp*log(lambda*y1(i)/tmp)*wt(i)
         end do
         response(1)=lambda; response(2)=d; risk=-2.0_dp*dev
      end select
   end subroutine evaluate_node

   real(dp) function point_error(model, y1, y2, pred) result(err)
      type(rpart_model), intent(in) :: model
      real(dp), intent(in) :: y1, y2, pred
      integer :: actual, pclass
      real(dp) :: tmp, dev
      select case(model%method)
      case(RPART_ANOVA)
         err=(y1-pred)**2
      case(RPART_CLASS)
         actual=nint(y1); pclass=nint(pred)
         if(model%class_freq(actual)>0.0_dp) then
            err=model%prior(actual)/model%class_freq(actual)*model%loss(actual,pclass)
         else
            err=0.0_dp
         end if
      case(RPART_POISSON,RPART_EXP)
         if(model%poisson_method==RPART_DEVIANCE) then
            tmp=y2; dev=tmp-pred*y1
            if(tmp>0.0_dp) dev=dev+tmp*log(pred*y1/tmp)
            err=-2.0_dp*dev
         else
            err=(sqrt(y2)-sqrt(max(0.0_dp,pred*y1)))**2
         end if
      case default
         err=0.0_dp
      end select
   end function point_error

   subroutine choose_split(model, xcol, y1, y2, obs, wt, ncat, myrisk, split)
      type(rpart_model), intent(in) :: model
      real(dp), intent(in) :: xcol(:), y1(:), y2(:), wt(:)
      integer, intent(in) :: obs(:), ncat
      real(dp), intent(in) :: myrisk
      type(rpart_split), intent(out) :: split
      integer, allocatable :: keep(:), ord(:)
      real(dp), allocatable :: xx(:), yy1(:), yy2(:), ww(:)
      integer :: i,k,n
      split%improve=0.0_dp; split%ncat=ncat; split%count=0
      n=0
      allocate(keep(size(obs)))
      do k=1,size(obs)
         i=obs(k)
         if(finite_dp(xcol(i)) .and. wt(i)>0.0_dp) then
            n=n+1; keep(n)=i
         end if
      end do
      split%count=n
      if(n<2*model%control%minbucket) return
      allocate(xx(n),yy1(n),yy2(n),ww(n))
      do k=1,n
         i=keep(k); xx(k)=xcol(i); yy1(k)=y1(i); yy2(k)=y2(i); ww(k)=wt(i)
      end do
      if(ncat==0) then
         call argsort_real(xx,ord)
         xx=xx(ord); yy1=yy1(ord); yy2=yy2(ord); ww=ww(ord)
         if(xx(1)>=xx(n)) return
      end if
      select case(model%method)
      case(RPART_ANOVA)
         call split_anova(xx,yy1,ww,ncat,model%control%minbucket,myrisk,split)
      case(RPART_CLASS)
         call split_class(model,xx,yy1,ww,ncat,model%control%minbucket,split)
      case(RPART_POISSON,RPART_EXP)
         call split_poisson(xx,yy1,yy2,ww,ncat,model%control%minbucket,split)
      end select
   end subroutine choose_split

   subroutine split_anova(x,y,wt,ncat,edge,myrisk,split)
      real(dp), intent(in) :: x(:),y(:),wt(:),myrisk
      integer,intent(in)::ncat,edge
      type(rpart_split),intent(inout)::split
      integer::i,j,n,left_n,right_n,where,nc
      real(dp)::right_wt,left_wt,right_sum,left_sum,grand,best,tmp
      real(dp),allocatable::sums(:),wts(:),means(:)
      integer,allocatable::counts(:),cats(:),cs(:)
      n=size(x); right_wt=sum(wt); grand=dot_product(y,wt)/right_wt
      if(ncat==0) then
         left_wt=0.0_dp; left_sum=0.0_dp; right_sum=0.0_dp
         left_n=0; right_n=n; best=0.0_dp; where=0
         do i=1,n-1
            left_wt=left_wt+wt(i); right_wt=right_wt-wt(i); left_n=left_n+1; right_n=right_n-1
            tmp=(y(i)-grand)*wt(i); left_sum=left_sum+tmp; right_sum=right_sum-tmp
            if(right_n>=edge .and. left_n>=edge .and. x(i+1)>x(i)) then
               tmp=left_sum*left_sum/left_wt+right_sum*right_sum/right_wt
               if(tmp>best) then
                  best=tmp; where=i
                  if(left_sum<right_sum) then; split%direction=RPART_LEFT; else; split%direction=RPART_RIGHT; end if
               end if
            end if
         end do
         if(myrisk>0.0_dp) split%improve=best/myrisk
         if(where>0) split%spoint=(x(where)+x(where+1))/2.0_dp
      else
         allocate(sums(ncat),wts(ncat),means(ncat),counts(ncat),cats(ncat),cs(ncat))
         sums=0.0_dp;wts=0.0_dp;counts=0;cs=RPART_RIGHT
         do i=1,n
            j=nint(x(i)); if(j<1.or.j>ncat) cycle
            counts(j)=counts(j)+1;wts(j)=wts(j)+wt(i);sums(j)=sums(j)+(y(i)-grand)*wt(i)
         end do
         do j=1,ncat
            if(counts(j)>0) then;means(j)=sums(j)/wts(j);else;means(j)=0.0_dp;cs(j)=0;end if
         end do
         call order_present_by_value(counts,means,cats,nc)
         left_wt=0.0_dp; left_sum=0.0_dp; right_wt=sum(wts); right_sum=0.0_dp
         left_n=0;right_n=sum(counts);best=0.0_dp
         allocate(split%csplit(ncat));split%csplit=0
         do i=1,max(0,nc-1)
            j=cats(i);cs(j)=RPART_LEFT;left_n=left_n+counts(j);right_n=right_n-counts(j)
            left_wt=left_wt+wts(j);right_wt=right_wt-wts(j);left_sum=left_sum+sums(j);right_sum=right_sum-sums(j)
            if(left_n>=edge.and.right_n>=edge) then
               tmp=left_sum*left_sum/left_wt+right_sum*right_sum/right_wt
               if(tmp>best) then
                  best=tmp
                  if(left_sum/left_wt>right_sum/right_wt) then;split%csplit=-cs;else;split%csplit=cs;end if
               end if
            end if
         end do
         if(myrisk>0.0_dp) split%improve=best/myrisk
      end if
   end subroutine split_anova

   subroutine split_class(model,x,y,wt,ncat,edge,split)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::x(:),y(:),wt(:)
      integer,intent(in)::ncat,edge
      type(rpart_split),intent(inout)::split
      real(dp),allocatable::left(:),right(:),ccnt(:,:),awt(:),rate(:)
      integer,allocatable::countn(:),tsplit(:),order(:),state(:)
      real(dp)::lwt,rwt,total_ss,best,tmp,lmean,rmean
      integer::i,j,k,n,ltot,rtot,where,nc,changed
      n=size(x);allocate(left(model%nclass),right(model%nclass));left=0.0_dp;right=0.0_dp
      lwt=0.0_dp;rwt=0.0_dp;ltot=0;rtot=0
      do i=1,n
         j=nint(y(i));rwt=rwt+model%altered_prior(j)*wt(i);right(j)=right(j)+wt(i);rtot=rtot+1
      end do
      total_ss=node_impurity(model,right,rwt);best=total_ss
      if(ncat==0) then
         where=0
         do i=1,n-1
            j=nint(y(i));rwt=rwt-model%altered_prior(j)*wt(i);lwt=lwt+model%altered_prior(j)*wt(i)
            rtot=rtot-1;ltot=ltot+1;right(j)=right(j)-wt(i);left(j)=left(j)+wt(i)
            if(ltot>=edge.and.rtot>=edge.and.x(i+1)>x(i)) then
               tmp=node_impurity(model,left,lwt)+node_impurity(model,right,rwt)
               if(tmp<best) then
                  best=tmp;where=i
                  call class_means(model,left,lwt,right,rwt,lmean,rmean)
                  if(lmean<rmean) then;split%direction=RPART_LEFT;else;split%direction=RPART_RIGHT;end if
               end if
            end if
         end do
         split%improve=total_ss-best
         if(where>0) split%spoint=(x(where)+x(where+1))/2.0_dp
      else
         allocate(ccnt(model%nclass,ncat),awt(ncat),rate(ncat),countn(ncat),tsplit(ncat),state(ncat))
         ccnt=0.0_dp;awt=0.0_dp;rate=0.0_dp;countn=0;tsplit=RPART_RIGHT
         do i=1,n
            j=nint(y(i));k=nint(x(i));if(k<1.or.k>ncat)cycle
            awt(k)=awt(k)+model%altered_prior(j)*wt(i);countn(k)=countn(k)+1;ccnt(j,k)=ccnt(j,k)+wt(i)
         end do
         do i=1,ncat
            if(awt(i)<=0.0_dp) tsplit(i)=0
            if(model%nclass==2.and.awt(i)>0.0_dp) rate(i)=ccnt(1,i)/awt(i)
         end do
         allocate(split%csplit(ncat));split%csplit=0
         if(model%nclass==2) then
            call order_present_by_value(countn,rate,order,nc)
            left=0.0_dp;lwt=0.0_dp;ltot=0;right=0.0_dp;rwt=0.0_dp;rtot=0
            do i=1,ncat
               right=right+ccnt(:,i);rwt=rwt+awt(i);rtot=rtot+countn(i)
            end do
            do i=1,max(0,nc-1)
               k=order(i);tsplit(k)=RPART_LEFT;left=left+ccnt(:,k);right=right-ccnt(:,k)
               lwt=lwt+awt(k);rwt=rwt-awt(k);ltot=ltot+countn(k);rtot=rtot-countn(k)
               if(ltot>=edge.and.rtot>=edge) then
                  tmp=node_impurity(model,left,lwt)+node_impurity(model,right,rwt)
                  if(tmp<best) then
                     best=tmp;call class_means(model,left,lwt,right,rwt,lmean,rmean)
                     if(lmean<rmean) then;split%csplit=tsplit;else;split%csplit=-tsplit;end if
                  end if
               end if
            end do
         else
            ! Exact upstream gray-code enumeration for multiclass categorical splits.
            state=0
            do i=1,ncat
               if(countn(i)>0) state(i)=1
            end do
            left=0.0_dp;lwt=0.0_dp;ltot=0;right=0.0_dp;rwt=0.0_dp;rtot=0
            do i=1,ncat
               right=right+ccnt(:,i);rwt=rwt+awt(i);rtot=rtot+countn(i)
            end do
            do
               changed=gray_next(state)
               if(changed==0) exit
               k=changed
               if(tsplit(k)==RPART_LEFT) then
                  tsplit(k)=RPART_RIGHT;left=left-ccnt(:,k);right=right+ccnt(:,k);lwt=lwt-awt(k);rwt=rwt+awt(k)
                  ltot=ltot-countn(k);rtot=rtot+countn(k)
               else
                  tsplit(k)=RPART_LEFT;left=left+ccnt(:,k);right=right-ccnt(:,k);lwt=lwt+awt(k);rwt=rwt-awt(k)
                  ltot=ltot+countn(k);rtot=rtot-countn(k)
               end if
               if(ltot>=edge.and.rtot>=edge) then
                  tmp=node_impurity(model,left,lwt)+node_impurity(model,right,rwt)
                  if(tmp<best) then
                     best=tmp;call class_means(model,left,lwt,right,rwt,lmean,rmean)
                     if(lmean<rmean) then;split%csplit=tsplit;else;split%csplit=-tsplit;end if
                  end if
               end if
            end do
         end if
         split%improve=total_ss-best
      end if
   end subroutine split_class

   real(dp) function node_impurity(model,cnt,totwt) result(v)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::cnt(:),totwt
      integer::j
      real(dp)::p
      v=0.0_dp;if(totwt<=0.0_dp)return
      do j=1,model%nclass
         p=model%altered_prior(j)*cnt(j)/totwt
         if(model%split_rule==RPART_INFORMATION) then
            if(p>0.0_dp)v=v-totwt*p*log(p)
         else
            v=v+totwt*p*(1.0_dp-p)
         end if
      end do
   end function node_impurity

   subroutine class_means(model,left,lwt,right,rwt,lmean,rmean)
      type(rpart_model),intent(in)::model
      real(dp),intent(in)::left(:),right(:),lwt,rwt
      real(dp),intent(out)::lmean,rmean
      integer::j;real(dp)::p
      lmean=0.0_dp;rmean=0.0_dp
      if(lwt>0.0_dp)then
         do j=1,model%nclass;p=model%altered_prior(j)*left(j)/lwt;lmean=lmean+p*real(j-1,dp);end do
      end if
      if(rwt>0.0_dp)then
         do j=1,model%nclass;p=model%altered_prior(j)*right(j)/rwt;rmean=rmean+p*real(j-1,dp);end do
      end if
   end subroutine class_means

   subroutine split_poisson(x,time,event,wt,ncat,edge,split)
      real(dp),intent(in)::x(:),time(:),event(:),wt(:)
      integer,intent(in)::ncat,edge
      type(rpart_split),intent(inout)::split
      integer::i,j,n,left_n,right_n,where,nc
      real(dp)::ld,rd,lt,rt,dev,best,tmp,lam1,lam2
      real(dp),allocatable::death(:),wtime(:),rate(:)
      integer,allocatable::countn(:),order(:)
      n=size(x);rd=dot_product(event,wt);rt=dot_product(time,wt);right_n=n
      if(rd<=0.0_dp.or.rt<=0.0_dp)return
      lam2=rd/rt;dev=rd*log(lam2);best=dev
      if(ncat==0)then
         ld=0.0_dp;lt=0.0_dp;left_n=0;where=0
         do i=1,n-1
            ld=ld+event(i)*wt(i);rd=rd-event(i)*wt(i);lt=lt+time(i)*wt(i);rt=rt-time(i)*wt(i)
            left_n=left_n+1;right_n=right_n-1
            if(left_n>=edge.and.right_n>=edge.and.x(i+1)>x(i))then
               lam1=ld/lt;lam2=rd/rt;tmp=0.0_dp
               if(lam1>0.0_dp)tmp=tmp+ld*log(lam1);if(lam2>0.0_dp)tmp=tmp+rd*log(lam2)
               if (tmp > best) then
                  best = tmp; where = i
                  if (lam1 < lam2) then
                     split%direction = RPART_LEFT
                  else
                     split%direction = RPART_RIGHT
                  end if
               end if
            end if
         end do
         split%improve=-2.0_dp*(dev-best);if(where>0)split%spoint=(x(where)+x(where+1))/2.0_dp
      else
         allocate(death(ncat),wtime(ncat),rate(ncat),countn(ncat));death=0.0_dp;wtime=0.0_dp;rate=0.0_dp;countn=0
         do i = 1, n
            j = nint(x(i))
            if (j < 1 .or. j > ncat) cycle
            countn(j) = countn(j) + 1
            death(j) = death(j) + event(i)*wt(i)
            wtime(j) = wtime(j) + time(i)*wt(i)
         end do
         do j=1,ncat;if(countn(j)>0)rate(j)=death(j)/wtime(j);end do
         call order_present_by_value(countn,rate,order,nc)
         allocate(split%csplit(ncat));split%csplit=0;ld=0.0_dp;lt=0.0_dp;left_n=0;where=0
         rd=sum(death);rt=sum(wtime);right_n=sum(countn)
         do i=1,max(0,nc-1)
            j=order(i);left_n=left_n+countn(j);right_n=right_n-countn(j);lt=lt+wtime(j);rt=rt-wtime(j);ld=ld+death(j);rd=rd-death(j)
            if(left_n>=edge.and.right_n>=edge)then
               lam1=ld/lt;lam2=rd/rt;tmp=0.0_dp;if(lam1>0.0_dp)tmp=tmp+ld*log(lam1);if(lam2>0.0_dp)tmp=tmp+rd*log(lam2)
               if (tmp > best) then
                  best = tmp; where = i
                  if (lam1 < lam2) then
                     split%direction = RPART_LEFT
                  else
                     split%direction = RPART_RIGHT
                  end if
               end if
            end if
         end do
         split%improve=-2.0_dp*(dev-best)
         if(where>0)then
            do i=1,where;split%csplit(order(i))=split%direction;end do
            do i=where+1,nc;split%csplit(order(i))=-split%direction;end do
         end if
      end if
   end subroutine split_poisson

   subroutine order_present_by_value(count,val,order,nc)
      integer,intent(in)::count(:)
      real(dp),intent(in)::val(:)
      integer,allocatable,intent(out)::order(:)
      integer,intent(out)::nc
      integer::i,j,t
      nc = 0
      do i = 1, size(count)
         if (count(i) > 0) nc = nc + 1
      end do
      allocate(order(nc));j=0
      do i=1,size(count);if(count(i)>0)then;j=j+1;order(j)=i;end if;end do
      do i=2,nc
         t=order(i);j=i-1
         do while(j>=1)
            if(val(order(j))<=val(t))exit
            order(j+1)=order(j);j=j-1
         end do
         order(j+1)=t
      end do
   end subroutine order_present_by_value

   integer function gray_next(state) result(changed)
      integer,intent(inout)::state(:)
      integer::i
      changed=0
      do i=1,size(state)-1
         if(state(i)==1)then;state(i)=2;changed=i;return
         else if(state(i)==2)then;state(i)=1
         end if
      end do
   end function gray_next

end module rpart_methods

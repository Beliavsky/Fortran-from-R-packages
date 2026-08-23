module rfast_repeated
   use rfast_special, only : dp, pi, f_cdf, student_t_cdf
   use rfast_arrays, only : colmeans
   implicit none
   private

   type, public :: variance_components_result
      real(dp), allocatable :: info(:,:)
      real(dp), allocatable :: ranef(:,:)
      integer :: status = 0
   end type variance_components_result

   type, public :: random_intercept_result
      real(dp), allocatable :: beta(:), se_beta(:), ranef(:)
      real(dp) :: sigma_tau = 0.0_dp
      real(dp) :: sigma_error = 0.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: deviance = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      real(dp) :: mean = 0.0_dp
      integer :: iterations = 0
      integer :: status = 0
   end type random_intercept_result

   public :: rm_anova, rm_anovas, rm_lines, varcomps_mom
   public :: varcomps_mle, colvarcomps_mle, varcomps_mle_balanced
   public :: rint_mle, rint_reg, rint_regbx, colrint_regbx, rint_regs

contains

   pure real(dp) function f_sf(x,df1,df2) result(p)
      real(dp),intent(in)::x,df1,df2
      p=max(0.0_dp,1.0_dp-f_cdf(x,df1,df2))
   end function f_sf

   real(dp) function f_quantile(prob,df1,df2) result(q)
      real(dp),intent(in)::prob,df1,df2
      real(dp)::lo,hi,mid,p
      integer::it
      if(prob<=0.0_dp)then;q=0.0_dp;return;else if(prob>=1.0_dp)then;q=huge(1.0_dp);return;end if
      lo=0.0_dp;hi=1.0_dp
      do while(f_cdf(hi,df1,df2)<prob.and.hi<1e12_dp);hi=2.0_dp*hi;end do
      do it=1,200
         mid=0.5_dp*(lo+hi);p=f_cdf(mid,df1,df2)
         if(p<prob)then;lo=mid;else;hi=mid;end if
         if(abs(hi-lo)<=1e-12_dp*max(1.0_dp,mid))exit
      end do
      q=0.5_dp*(lo+hi)
   end function f_quantile

   function rm_anova(y) result(out)
      real(dp),intent(in)::y(:,:)
      real(dp)::out(2)
      integer::n,d,i,j
      real(dp)::grand,sst,ssr,mst,msr,rowm(size(y,1)),colm(size(y,2))
      n=size(y,1);d=size(y,2);rowm=sum(y,dim=2)/real(d,dp);colm=sum(y,dim=1)/real(n,dp);grand=sum(rowm)/real(n,dp)
      sst=real(n,dp)*sum((colm-grand)**2);ssr=0.0_dp
      do i=1,n;do j=1,d;ssr=ssr+(y(i,j)-rowm(i)-colm(j)+grand)**2;end do;end do
      mst=sst/real(d-1,dp);msr=ssr/real((d-1)*(n-1),dp);out(1)=mst/msr;out(2)=f_sf(out(1),real(d-1,dp),real((d-1)*(n-1),dp))
   end function rm_anova

   function rm_anovas(y) result(out)
      real(dp),intent(in)::y(:,:,:)
      real(dp)::out(size(y,3),2)
      integer::k
      do k=1,size(y,3);out(k,:)=rm_anova(y(:,:,k));end do
   end function rm_anovas

   function rm_lines(y,x) result(out)
      real(dp),intent(in)::y(:,:,:),x(:)
      real(dp)::out(size(y,3),2)
      integer::n,d,p,i,k
      real(dp)::z(size(x)),den,be(size(y,1)),mb,sb,t
      n=size(y,1);d=size(y,2);p=size(y,3)
      if(size(x)/=d)then;out=huge(1.0_dp);return;end if
      z=x-sum(x)/real(d,dp);den=sum(z*z);z=z/den
      do k=1,p
         do i=1,n;be(i)=dot_product(z,y(i,:,k));end do
         mb=sum(be)/real(n,dp);sb=sqrt(sum((be-mb)**2)/real(n-1,dp));t=sqrt(real(n,dp))*mb/max(tiny(1.0_dp),sb)
         out(k,1)=t;out(k,2)=min(1.0_dp,2.0_dp*max(0.0_dp,1.0_dp-student_t_cdf(abs(t),real(n-1,dp))))
      end do
   end function rm_lines

   subroutine group_layout(id,groups,index,counts)
      integer,intent(in)::id(:)
      integer,allocatable,intent(out)::groups(:),index(:),counts(:)
      integer,allocatable::tmp(:)
      integer::i,j,k,key,ng
      allocate(tmp(size(id)));tmp=id
      do i=2,size(tmp);key=tmp(i);j=i-1;do while(j>=1);if(tmp(j)<=key)exit;tmp(j+1)=tmp(j);j=j-1;end do;tmp(j+1)=key;end do
      ng=1;do i=2,size(tmp);if(tmp(i)/=tmp(i-1))ng=ng+1;end do
      allocate(groups(ng),counts(ng),index(size(id)));groups(1)=tmp(1);k=1
      do i=2,size(tmp);if(tmp(i)/=tmp(i-1))then;k=k+1;groups(k)=tmp(i);end if;end do
      counts=0
      do i=1,size(id)
         do j=1,ng;if(id(i)==groups(j))then;index(i)=j;counts(j)=counts(j)+1;exit;end if;end do
      end do
   end subroutine group_layout

   function varcomps_mom(x,id) result(res)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::id(:)
      real(dp)::res(size(x,2),5)
      integer,allocatable::groups(:),idx(:),cnt(:)
      real(dp),allocatable::gs(:,:),a(:),b(:),sx2(:),mst(:),mse(:),fa(:),rv(:),ratio(:)
      real(dp)::ql,qu,l,u
      integer::n,p,k,i,j
      if(size(id)/=size(x,1))then;res=huge(1.0_dp);return;end if
      call group_layout(id,groups,idx,cnt);n=size(x,1);p=size(x,2);k=size(groups)
      allocate(gs(k,p),a(p),b(p),sx2(p),mst(p),mse(p),fa(p),rv(p),ratio(p));gs=0.0_dp
      do i=1,n;gs(idx(i),:)=gs(idx(i),:)+x(i,:);end do
      sx2=sum(x*x,dim=1);a=0.0_dp
      do i=1,k;a=a+gs(i,:)**2/real(cnt(i),dp);end do
      b=sum(gs,dim=1)**2/real(n,dp);mst=(a-b)/real(k-1,dp);mse=(sx2-a)/real(n-k,dp);fa=mst/mse
      rv=max(0.0_dp,(mst-mse)/real(k,dp));ratio=rv/(rv+mse)
      ql=f_quantile(0.025_dp,real(n-k,dp),real(k-1,dp));qu=f_quantile(0.975_dp,real(n-k,dp),real(k-1,dp))
      do j=1,p
         l=(fa(j)*ql-1.0_dp)/real(k,dp);u=(fa(j)*qu-1.0_dp)/real(k,dp)
         res(j,:)=[rv(j),mse(j),ratio(j),l/(1.0_dp+l),u/(1.0_dp+u)]
      end do
   end function varcomps_mom

   pure real(dp) function variance_objective(d,n,ni,ni2hi2,sse) result(val)
      real(dp), intent(in) :: d, n, ni(:), ni2hi2(:), sse
      real(dp) :: rem
      rem = sse - d*sum(ni2hi2/(1.0_dp + ni*d))
      if (rem <= 0.0_dp) then
         val = huge(1.0_dp)
      else
         val = sum(log(1.0_dp + ni*d)) + n*log(rem)
      end if
   end function variance_objective

   subroutine gold_rat3(n,ni,ni2,sse,hi2,tol,dval,fval)
      real(dp), intent(in) :: n, ni(:), ni2(:), sse, hi2(:), tol
      real(dp), intent(out) :: dval, fval
      real(dp), parameter :: ratio = 0.618033988749895_dp
      real(dp) :: a,b,x1,x2,f1,f2,width
      real(dp), allocatable :: ni2hi2(:)
      allocate(ni2hi2(size(ni)))
      ni2hi2 = ni2*hi2
      a=0.0_dp; b=50.0_dp
      x1=b-ratio*b; x2=ratio*b
      f1=variance_objective(x1,n,ni,ni2hi2,sse)
      f2=variance_objective(x2,n,ni,ni2hi2,sse)
      width=b-a
      do while(abs(width)>tol)
         if(f2>f1)then
            b=x2; width=b-a; x2=x1; f2=f1
            x1=b-ratio*width
            f1=variance_objective(x1,n,ni,ni2hi2,sse)
         else
            a=x1; width=b-a; x1=x2; f1=f2
            x2=a+ratio*width
            f2=variance_objective(x2,n,ni,ni2hi2,sse)
         end if
      end do
      dval=0.5_dp*(x1+x2)
      fval=0.5_dp*(f1+f2)
   end subroutine gold_rat3

   function colvarcomps_mle(x,id,want_ranef,tol,maxiters) result(res)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: id(:)
      logical, intent(in), optional :: want_ranef
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiters
      type(variance_components_result) :: res
      integer, allocatable :: groups(:),idx(:),cnt(:)
      real(dp), allocatable :: ni(:),ni2(:),gs(:,:),gm(:),oneplus(:),hi2(:)
      integer :: n,p,k,i,j,it,mit
      real(dp) :: eps,sxy,b1,b2,sse,dval,fval,down,sigma
      logical :: wr
      wr=.false.; if(present(want_ranef))wr=want_ranef
      eps=1.0e-8_dp; if(present(tol))eps=tol
      mit=100; if(present(maxiters))mit=maxiters
      if(size(id)/=size(x,1) .or. size(id)==0)then;res%status=1;return;end if
      call group_layout(id,groups,idx,cnt)
      if(all(cnt==cnt(1)))then
         res=varcomps_mle_balanced(x,id,wr)
         return
      end if
      n=size(x,1); p=size(x,2); k=size(groups)
      allocate(ni(k),ni2(k),gs(k,p),gm(k),oneplus(k),hi2(k))
      ni=real(cnt,dp); ni2=ni*ni; gs=0.0_dp
      do i=1,n
         gs(idx(i),:)=gs(idx(i),:)+x(i,:)
      end do
      allocate(res%info(p,3)); res%info=0.0_dp
      if(wr)then;allocate(res%ranef(k,p));res%ranef=0.0_dp;end if
      do j=1,p
         sxy=sum(x(:,j)); gm=gs(:,j)/ni; b1=sxy/real(n,dp)
         sse=sum((x(:,j)-b1)**2); hi2=(gm-b1)**2
         call gold_rat3(real(n,dp),ni,ni2,sse,hi2,eps,dval,fval)
         oneplus=1.0_dp+ni*dval
         down=real(n,dp)-dval*sum(ni2/oneplus)
         if(abs(down)<=tiny(1.0_dp))then;res%status=max(res%status,2);cycle;end if
         b2=(sxy-dval*sum((ni*gs(:,j))/oneplus))/down
         it=2
         do while(it<mit .and. abs(b2-b1)>eps)
            it=it+1; b1=b2
            sse=sum((x(:,j)-b1)**2); hi2=(gm-b1)**2
            call gold_rat3(real(n,dp),ni,ni2,sse,hi2,eps,dval,fval)
            oneplus=1.0_dp+ni*dval
            down=real(n,dp)-dval*sum(ni2/oneplus)
            if(abs(down)<=tiny(1.0_dp))exit
            b2=(sxy-dval*sum((ni*gs(:,j))/oneplus))/down
         end do
         if(abs(down)<=tiny(1.0_dp))then;res%status=max(res%status,2);cycle;end if
         sigma=sse/real(n,dp)
         res%info(j,2)=sigma/(1.0_dp+dval)
         res%info(j,1)=sigma-res%info(j,2)
         res%info(j,3)=-0.5_dp*(fval+real(n,dp)*(1.83787706640935_dp-log(real(n,dp))+1.0_dp))
         if(wr)res%ranef(:,j)=(gm-b2)*(dval*ni/oneplus)
      end do
   end function colvarcomps_mle

   function varcomps_mle(x,id,want_ranef,tol,maxiters) result(res)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: id(:)
      logical, intent(in), optional :: want_ranef
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiters
      type(variance_components_result) :: res
      real(dp), allocatable :: xx(:,:)
      type(variance_components_result) :: tmp
      logical :: wr
      real(dp) :: eps
      integer :: mit
      wr=.false.;if(present(want_ranef))wr=want_ranef
      eps=1.0e-8_dp;if(present(tol))eps=tol
      mit=100;if(present(maxiters))mit=maxiters
      allocate(xx(size(x),1));xx(:,1)=x
      tmp=colvarcomps_mle(xx,id,wr,eps,mit)
      res=tmp
   end function varcomps_mle

   function varcomps_mle_balanced(x,id,want_ranef) result(res)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::id(:)
      logical,intent(in),optional::want_ranef
      type(variance_components_result)::res
      integer,allocatable::groups(:),idx(:),cnt(:)
      real(dp),allocatable::gm(:,:),mx(:),ex(:,:),seid(:,:),co(:,:),sml(:),dml(:),tau(:),sx(:),sx2(:),com(:,:)
      integer::n,p,k,d,i,j
      real(dp)::f
      logical::wr
      wr=.false.;if(present(want_ranef))wr=want_ranef
      if(size(id)/=size(x,1))then;res%status=1;return;end if
      call group_layout(id,groups,idx,cnt);n=size(x,1);p=size(x,2);k=size(groups)
      if(any(cnt/=cnt(1)))then;res%status=2;return;end if
      d=cnt(1);f=1.0_dp-1.0_dp/real(d,dp)
      allocate(gm(k,p),mx(p),ex(n,p),seid(k,p),co(k,p),sml(p),dml(p),tau(p),sx(p),sx2(p),com(k,p))
      gm=0.0_dp;do i=1,n;gm(idx(i),:)=gm(idx(i),:)+x(i,:);end do;gm=gm/real(d,dp);mx=sum(gm,dim=1)/real(k,dp)
      ex=x-spread(mx,1,n);seid=0.0_dp
      do i=1,n;seid(idx(i),:)=seid(idx(i),:)+ex(i,:)**2;end do
      do i=1,k;co(i,:)=(gm(i,:)-mx)**2;end do
      sx2=sum(co,dim=1);sml=sum(seid-real(d,dp)*co,dim=1)/real(n,dp)/f;dml=sx2/real(k,dp)/sml-1.0_dp/real(d,dp);tau=dml*sml
      do j=1,p;if(dml(j)<0.0_dp)then;sml(j)=sml(j)+tau(j);tau(j)=0.0_dp;end if;end do
      sx=sum(seid,dim=1);sx2=sx2*real(d*d,dp);allocate(res%info(p,3))
      do j=1,p
         res%info(j,1)=tau(j);res%info(j,2)=sml(j)
         res%info(j,3)=-0.5_dp*(real(n,dp)*log(sml(j))+real(k,dp)*log(1.0_dp+real(d,dp)*tau(j)/sml(j)) &
             +sx(j)/sml(j)-tau(j)/(sml(j)**2+real(d,dp)*tau(j)*sml(j))*sx2(j))-0.5_dp*real(n,dp)*log(2.0_dp*pi)
      end do
      if(wr)then
         com=0.0_dp;do i=1,n;com(idx(i),:)=com(idx(i),:)+ex(i,:);end do;com=com/real(d,dp);allocate(res%ranef(k,p))
         do j=1,p;res%ranef(:,j)=tau(j)/(tau(j)+sml(j)/real(d,dp))*com(:,j);end do
      end if
   end function varcomps_mle_balanced


   function rint_mle(x,id,want_ranef,tol,maxiters) result(out)
      real(dp),intent(in)::x(:)
      integer,intent(in)::id(:)
      logical,intent(in),optional::want_ranef
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiters
      type(random_intercept_result)::out
      integer,allocatable::groups(:),idx(:),cnt(:)
      real(dp),allocatable::ni(:),ni2(:),gs(:),gm(:),hi2(:),oneplus(:)
      real(dp)::eps,sx,b1,b2,sse,dval,fval,down,sigma
      integer::n,k,i,it,mi
      logical::wr
      if(size(x)==0.or.size(id)/=size(x))then;out%status=1;return;end if
      wr=.false.;if(present(want_ranef))wr=want_ranef
      eps=1.0e-9_dp;if(present(tol))eps=tol;mi=100;if(present(maxiters))mi=maxiters
      call group_layout(id,groups,idx,cnt);n=size(x);k=size(groups)
      allocate(ni(k),ni2(k),gs(k),gm(k),hi2(k),oneplus(k));ni=real(cnt,dp);ni2=ni*ni;gs=0.0_dp
      do i=1,n;gs(idx(i))=gs(idx(i))+x(i);end do
      gm=gs/ni;sx=sum(x);b1=sx/real(n,dp);sse=sum((x-b1)**2);hi2=(gm-b1)**2
      call gold_rat3(real(n,dp),ni,ni2,sse,hi2,eps,dval,fval)
      oneplus=1.0_dp+ni*dval;down=real(n,dp)-dval*sum(ni2/oneplus)
      if(abs(down)<=tiny(1.0_dp))then;out%status=2;return;end if
      b2=(sx-dval*sum(ni*gs/oneplus))/down;it=2
      do while(it<mi.and.abs(b2-b1)>eps)
         it=it+1;b1=b2;sse=sum((x-b1)**2);hi2=(gm-b1)**2
         call gold_rat3(real(n,dp),ni,ni2,sse,hi2,eps,dval,fval)
         oneplus=1.0_dp+ni*dval;down=real(n,dp)-dval*sum(ni2/oneplus)
         if(abs(down)<=tiny(1.0_dp))then;out%status=2;return;end if
         b2=(sx-dval*sum(ni*gs/oneplus))/down
      end do
      sigma=sse/real(n,dp);out%sigma_error=sigma/(1.0_dp+dval);out%sigma_tau=sigma-out%sigma_error
      out%mean=b2;out%iterations=it
      out%loglik=-0.5_dp*(fval+real(n,dp)*(1.83787706640935_dp-log(real(n,dp))+1.0_dp))
      out%deviance=-2.0_dp*out%loglik;out%bic=out%deviance+3.0_dp*log(real(n,dp))
      allocate(out%beta(1),out%se_beta(1));out%beta=[b2]
      out%se_beta(1)=sqrt(max(0.0_dp,out%sigma_error/down))
      if(wr)then;allocate(out%ranef(k));out%ranef=(gm-b2)*(dval*ni/oneplus);end if
   end function rint_mle

   function rint_reg(y,x,id,tol,want_ranef,maxiters) result(out)
      real(dp),intent(in)::y(:),x(:,:)
      integer,intent(in)::id(:)
      real(dp),intent(in),optional::tol
      logical,intent(in),optional::want_ranef
      integer,intent(in),optional::maxiters
      type(random_intercept_result)::out
      integer,allocatable::groups(:),idx(:),cnt(:)
      real(dp),allocatable::xx(:,:),ni(:),ni2(:),sxg(:,:),syg(:),mx(:,:),my(:),b1(:),b2(:)
      real(dp),allocatable::resid(:),hi2(:),oneplus(:),a(:,:),rhs(:),step(:),ainv(:,:)
      real(dp)::eps,sse,dval,fval,sigma
      integer::n,p,k,i,j,it,mi,info
      logical::wr
      if(size(y)==0.or.size(x,1)/=size(y).or.size(id)/=size(y))then;out%status=1;return;end if
      wr=.false.;if(present(want_ranef))wr=want_ranef
      eps=1.0e-8_dp;if(present(tol))eps=tol;mi=100;if(present(maxiters))mi=maxiters
      n=size(y);p=size(x,2)+1;allocate(xx(n,p));xx(:,1)=1.0_dp;xx(:,2:)=x
      call group_layout(id,groups,idx,cnt);k=size(groups);allocate(ni(k),ni2(k),sxg(k,p),syg(k),mx(k,p),my(k))
      ni=real(cnt,dp);ni2=ni*ni;sxg=0.0_dp;syg=0.0_dp
      do i=1,n;sxg(idx(i),:)=sxg(idx(i),:)+xx(i,:);syg(idx(i))=syg(idx(i))+y(i);end do
      do i=1,k;mx(i,:)=sxg(i,:)/ni(i);my(i)=syg(i)/ni(i);end do
      allocate(b1(p),b2(p),resid(n),hi2(k),oneplus(k),a(p,p),rhs(p),step(p),ainv(p,p))
      a=matmul(transpose(xx),xx);rhs=matmul(transpose(xx),y);call solve_linear_local(a,rhs,b1,info)
      if(info/=0)then;out%status=info;return;end if
      do it=1,mi
         resid=y-matmul(xx,b1);sse=sum(resid*resid);hi2=(my-matmul(mx,b1))**2
         call gold_rat3(real(n,dp),ni,ni2,sse,hi2,eps,dval,fval);oneplus=1.0_dp+ni*dval
         a=matmul(transpose(xx),xx)-dval*matmul(transpose(sxg/spread(oneplus,2,p)),sxg)
         rhs=matmul(transpose(xx),y)-dval*matmul(transpose(sxg),syg/oneplus)
         call solve_linear_local(a,rhs,b2,info);if(info/=0)then;out%status=info;return;end if
         if(sum(abs(b2-b1))<=eps)exit;b1=b2
      end do
      b1=b2;resid=y-matmul(xx,b1);sse=sum(resid*resid);hi2=(my-matmul(mx,b1))**2
      call gold_rat3(real(n,dp),ni,ni2,sse,hi2,eps,dval,fval);oneplus=1.0_dp+ni*dval
      sigma=(sse-dval*sum(ni2*hi2/oneplus))/real(n,dp);out%sigma_error=sigma;out%sigma_tau=dval*sigma
      out%loglik=-0.5_dp*fval-0.5_dp*real(n,dp)*(1.83787706640935_dp-log(real(n,dp))+1.0_dp)
      out%deviance=-2.0_dp*out%loglik;out%bic=out%deviance+real(p+2,dp)*log(real(n,dp));out%iterations=it
      allocate(out%beta(p),out%se_beta(p));out%beta=b1
      a=matmul(transpose(xx),xx)-dval*matmul(transpose(sxg/spread(oneplus,2,p)),sxg)
      call inverse_matrix_local(a,ainv,info)
      if(info==0)then
         do j=1,p;out%se_beta(j)=sqrt(max(0.0_dp,ainv(j,j)*sigma));end do
      else;out%se_beta=huge(1.0_dp);end if
      if(wr)then
         allocate(out%ranef(k));out%ranef=dval*ni/oneplus
         do i=1,k;out%ranef(i)=out%ranef(i)*sum(pack(resid,idx==i))/ni(i);end do
      end if
   end function rint_reg

   function rint_regbx(y,x,id) result(out)
      real(dp),intent(in)::y(:),x(:,:)
      integer,intent(in)::id(:)
      type(random_intercept_result)::out
      integer,allocatable::groups(:),idx(:),cnt(:)
      real(dp),allocatable::xx(:,:),a(:,:),rhs(:),b(:),resid(:),seid(:),myid(:),com(:)
      real(dp)::sml,dml,tau,sz,sz2,f
      integer::n,p,k,d,i,info
      if(size(y)==0.or.size(x,1)/=size(y).or.size(id)/=size(y))then;out%status=1;return;end if
      call group_layout(id,groups,idx,cnt);if(any(cnt/=cnt(1)))then;out%status=2;return;end if
      n=size(y);p=size(x,2)+1;k=size(groups);d=cnt(1);f=1.0_dp-1.0_dp/real(d,dp)
      allocate(xx(n,p));xx(:,1)=1.0_dp;xx(:,2:)=x;allocate(a(p,p),rhs(p),b(p),resid(n),seid(k),myid(k),com(k))
      a=matmul(transpose(xx),xx);rhs=matmul(transpose(xx),y);call solve_linear_local(a,rhs,b,info)
      if(info/=0)then;out%status=info;return;end if
      resid=y-matmul(xx,b);seid=0.0_dp;myid=0.0_dp
      do i=1,n;seid(idx(i))=seid(idx(i))+resid(i)**2;myid(idx(i))=myid(idx(i))+y(i);end do
      myid=myid/real(d,dp);com=(myid-sum(myid)/real(k,dp))**2
      sml=sum(seid-real(d,dp)*com)/real(n,dp)/f;dml=sum(com)/real(k,dp)/sml-1.0_dp/real(d,dp);tau=dml*sml
      if(tau<0.0_dp)then;sml=sml+tau;tau=0.0_dp;end if
      com=0.0_dp;do i=1,n;com(idx(i))=com(idx(i))+resid(i);end do
      sz=sum(seid);sz2=sum(com*com);out%sigma_tau=tau;out%sigma_error=sml
      out%loglik=-0.5_dp*(real(n,dp)*log(sml)+real(k,dp)*log(1.0_dp+real(d,dp)*tau/sml)+sz/sml &
          -tau/(sml*sml+real(d,dp)*tau*sml)*sz2)-0.5_dp*real(n,dp)*log(2.0_dp*pi)
      out%deviance=-2.0_dp*out%loglik;out%bic=out%deviance+real(p+2,dp)*log(real(n,dp));allocate(out%beta(p));out%beta=b
      allocate(out%ranef(k));out%ranef=tau/(tau+sml/real(d,dp))*com/real(d,dp);out%iterations=1
   end function rint_regbx

   function colrint_regbx(y,x,id) result(out)
      real(dp),intent(in)::y(:,:),x(:,:)
      integer,intent(in)::id(:)
      type(random_intercept_result),allocatable::out(:)
      integer::j
      allocate(out(size(y,2)))
      do j=1,size(y,2);out(j)=rint_regbx(y(:,j),x,id);end do
   end function colrint_regbx

   function rint_regs(y,x,id,tol,maxiters,logged) result(out)
      real(dp),intent(in)::y(:),x(:,:)
      integer,intent(in)::id(:)
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiters
      logical,intent(in),optional::logged
      real(dp)::out(size(x,2),2)
      type(random_intercept_result)::fit
      real(dp),allocatable::one(:,:)
      real(dp)::stat,pv
      integer::j,n
      logical::lg
      n=size(y);lg=.false.;if(present(logged))lg=logged;allocate(one(n,1))
      do j=1,size(x,2)
         one(:,1)=x(:,j);fit=rint_reg(y,one,id,tol,.false.,maxiters)
         if(fit%status/=0.or..not.allocated(fit%beta).or..not.allocated(fit%se_beta))then
            out(j,:)=huge(1.0_dp);cycle
         end if
         stat=(fit%beta(2)/max(tiny(1.0_dp),fit%se_beta(2)))**2
         pv=max(tiny(1.0_dp),1.0_dp-f_cdf(stat,1.0_dp,real(max(1,n-4),dp)))
         out(j,:)=[stat,merge(log(pv),pv,lg)]
      end do
   end function rint_regs

   subroutine solve_linear_local(a,b,x,info)
      use rfast_linalg, only : solve_linear
      real(dp),intent(in)::a(:,:),b(:);real(dp),intent(out)::x(:);integer,intent(out)::info
      call solve_linear(a,b,x,info)
   end subroutine solve_linear_local

   subroutine inverse_matrix_local(a,ainv,info)
      use rfast_linalg, only : inverse_matrix
      real(dp),intent(in)::a(:,:);real(dp),intent(out)::ainv(:,:);integer,intent(out)::info
      call inverse_matrix(a,ainv,info)
   end subroutine inverse_matrix_local

end module rfast_repeated

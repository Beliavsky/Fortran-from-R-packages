module rebayes_repeated
   use rebayes_kinds, only : dp, pi
   use rebayes_math, only : gamma_pdf, normal_pdf, student_pdf, safe_log
   use rebayes_kw, only : kw_result, kw_control, kw_fit
   use rebayes_mixtures, only : mixture_fit, bivariate_mixture_fit, gvmix
   implicit none
   private
   public :: group_stats, wgvmix, wglvmix, wlvmix, wtlvmix
contains
   subroutine group_stats(y,id,w,t,s,m,wsum,sumlogw)
      real(dp),intent(in)::y(:),w(:)
      integer,intent(in)::id(:)
      real(dp),allocatable,intent(out)::t(:),s(:),m(:),wsum(:),sumlogw(:)
      integer,allocatable::uid(:),cnt(:)
      integer::i,j,ng,n
      logical :: found
      real(dp),allocatable::sy(:),sy2(:)
      n=size(y);if(size(id)/=n.or.size(w)/=n)error stop "group_stats: sizes"
      allocate(uid(n));ng=0
      do i=1,n
         found=.false.
         do j=1,ng
            if(uid(j)==id(i))then
               found=.true.; exit
            end if
         end do
         if(.not.found)then
            ng=ng+1;uid(ng)=id(i)
         end if
      end do
      allocate(t(ng),s(ng),m(ng),wsum(ng),sumlogw(ng),cnt(ng),sy(ng),sy2(ng))
      t=0;s=0;m=0;wsum=0;sumlogw=0;cnt=0;sy=0;sy2=0
      do i=1,n
         j=1;do while(uid(j)/=id(i));j=j+1;end do
         cnt(j)=cnt(j)+1;wsum(j)=wsum(j)+w(i);sy(j)=sy(j)+w(i)*y(i);sy2(j)=sy2(j)+w(i)*y(i)*y(i)
         sumlogw(j)=sumlogw(j)+log(max(w(i),tiny(1.0_dp)))
      end do
      do j=1,ng
         m(j)=real(cnt(j),dp);t(j)=sy(j)/wsum(j)
         if(cnt(j)>1)then;s(j)=(sy2(j)-t(j)*t(j)*wsum(j))/real(cnt(j)-1,dp);else;s(j)=0;end if
      end do
   end subroutine group_stats

   subroutine wgvmix(y,id,w,v,fit,control)
      real(dp),intent(in)::y(:),w(:),v(:);integer,intent(in)::id(:)
      type(mixture_fit),intent(out)::fit;type(kw_control),intent(in),optional::control
      real(dp),allocatable::t(:),s(:),m(:),ws(:),slw(:),a(:,:),ww(:),d(:)
      type(kw_result)::kw;integer::i,j,n
      call group_stats(y,id,w,t,s,m,ws,slw);n=size(s)
      allocate(a(n,size(v)),ww(n),d(size(v)));ww=1.0_dp/real(n,dp);d=1.0_dp
      do j=1,size(v);do i=1,n
         a(i,j)=gamma_pdf(s(i),0.5_dp*(m(i)-1.0_dp),v(j)/(0.5_dp*(m(i)-1.0_dp)))
      end do;end do
      call kw_fit(a,d,ww,kw,control)
      allocate(fit%grid(size(v)),fit%mass(size(v)),fit%g(n),fit%posterior_mean(n))
      fit%grid=v;fit%mass=kw%f;fit%g=kw%g;fit%posterior_mean=matmul(a,kw%f*v)/max(kw%g,tiny(1.0_dp))
      fit%loglik=sum(safe_log(kw%g));fit%status=kw%status;fit%iterations=kw%iterations;fit%kkt_gap=kw%kkt_gap
   end subroutine wgvmix

   subroutine wglvmix(y,id,w,u,v,fit,control)
      real(dp),intent(in)::y(:),w(:),u(:),v(:);integer,intent(in)::id(:)
      type(bivariate_mixture_fit),intent(out)::fit;type(kw_control),intent(in),optional::control
      real(dp),allocatable::t(:),s(:),m(:),ws(:),slw(:),a(:,:),ww(:),d(:),ug(:),vg(:)
      type(kw_result)::kw;real(dp)::r,z,logk
      integer::i,j,l,c,n
      call group_stats(y,id,w,t,s,m,ws,slw);n=size(s)
      allocate(a(n,size(u)*size(v)),ww(n),d(size(u)*size(v)),ug(size(u)*size(v)),vg(size(u)*size(v)))
      ww=1.0_dp/real(n,dp);d=1.0_dp;c=0
      do l=1,size(v);do j=1,size(u);c=c+1;ug(c)=u(j);vg(c)=v(l)
         do i=1,n
            r=0.5_dp*(m(i)-1.0_dp);z=(t(i)-u(j))*sqrt(ws(i)/v(l))
            a(i,c)=gamma_pdf(s(i),r,v(l)/r)*normal_pdf(z,0.0_dp,1.0_dp)
         end do
      end do;end do
      call kw_fit(a,d,ww,kw,control)
      allocate(fit%u(size(u)),fit%v(size(v)),fit%mass(size(d)))
      allocate(fit%g(n),fit%post_u(n),fit%post_v(n),fit%post_product(n),fit%a(n,size(d)))
      fit%u=u;fit%v=v;fit%mass=kw%f;fit%g=kw%g;fit%a=a
      fit%post_u=matmul(a,kw%f*ug)/max(kw%g,tiny(1.0_dp));fit%post_v=matmul(a,kw%f*vg)/max(kw%g,tiny(1.0_dp))
      fit%post_product=matmul(a,kw%f*ug*vg)/max(kw%g,tiny(1.0_dp))
      logk=0.0_dp
      do i=1,n
         r=0.5_dp*(m(i)-1.0_dp)
         logk=logk+log_gamma(r)-r*log(r)-0.5_dp*log(ws(i))-r*log(2.0_dp*pi) &
            -(r-1.0_dp)*log(max(s(i),tiny(1.0_dp)))+0.5_dp*slw(i)
      end do
      fit%loglik=sum(safe_log(kw%g))+logk;fit%status=kw%status
   end subroutine wglvmix

   subroutine wlvmix(y,id,w,u,v,fit,maxit,tol,control)
      real(dp),intent(in)::y(:),w(:),u(:),v(:);integer,intent(in)::id(:)
      type(bivariate_mixture_fit),intent(out)::fit
      integer,intent(in),optional::maxit;real(dp),intent(in),optional::tol
      type(kw_control),intent(in),optional::control
      real(dp),allocatable::t(:),s(:),m(:),ws(:),slw(:),a3(:,:,:),b(:,:),c(:,:),ww(:),du(:),dv(:)
      real(dp),allocatable::fu(:),fv(:),g(:),au(:,:),av(:,:),ug(:),vg(:),flat(:,:)
      type(kw_result)::ku,kv;real(dp)::r,sd,ll,llold,eps
      integer::i,j,l,it,n,ni,cidx
      call group_stats(y,id,w,t,s,m,ws,slw);n=size(s);ni=2;if(present(maxit))ni=maxit;eps=1e-4_dp;if(present(tol))eps=tol
      allocate(a3(n,size(v),size(u)),ww(n),du(size(u)),dv(size(v)),fu(size(u)),fv(size(v)))
      ww=1.0_dp/real(n,dp);du=1;dv=1
      do j=1,size(u);do l=1,size(v);do i=1,n
         r=0.5_dp*(m(i)-1.0_dp);sd=sqrt(v(l)/ws(i))
         a3(i,l,j)=gamma_pdf(s(i),r,v(l)/r)*normal_pdf(t(i),u(j),sd)
      end do;end do;end do
      allocate(c(n,size(v)))
      do l=1,size(v);do i=1,n;c(i,l)=gamma_pdf(s(i),0.5_dp*(m(i)-1),v(l)/(0.5_dp*(m(i)-1)));end do;end do
      call kw_fit(c,dv,ww,kv,control);fv=kv%f
      allocate(b(n,size(u)))
      do i=1,n;do j=1,size(u);b(i,j)=sum(a3(i,:,j)*fv);end do;end do
      call kw_fit(b,du,ww,ku,control);fu=ku%f;ll=sum(safe_log(matmul(b,fu)));llold=-huge(1.0_dp)
      do it=2,ni
         if(ll-llold<=eps)exit;llold=ll
         do i=1,n;do l=1,size(v);c(i,l)=sum(a3(i,l,:)*fu);end do;end do
         call kw_fit(c,dv,ww,kv,control);fv=kv%f
         do i=1,n;do j=1,size(u);b(i,j)=sum(a3(i,:,j)*fv);end do;end do
         call kw_fit(b,du,ww,ku,control);fu=ku%f;ll=sum(safe_log(matmul(b,fu)))
      end do
      allocate(au(n,size(u)),av(n,size(v)),g(n),ug(size(u)*size(v)),vg(size(u)*size(v)),flat(n,size(u)*size(v)))
      do i=1,n;do j=1,size(u);au(i,j)=sum(a3(i,:,j)*fv);end do;do l=1,size(v);av(i,l)=sum(a3(i,l,:)*fu);end do;end do
      g=matmul(au,fu);cidx=0
      do l=1,size(v);do j=1,size(u);cidx=cidx+1;ug(cidx)=u(j);vg(cidx)=v(l);flat(:,cidx)=a3(:,l,j);end do;end do
      allocate(fit%u(size(u)),fit%v(size(v)),fit%mass(size(u)*size(v)))
      allocate(fit%g(n),fit%post_u(n),fit%post_v(n),fit%post_product(n),fit%a(n,size(u)*size(v)))
      fit%u=u;fit%v=v;fit%g=g;fit%a=flat;cidx=0
      do l=1,size(v);do j=1,size(u);cidx=cidx+1;fit%mass(cidx)=fu(j)*fv(l);end do;end do
      fit%post_u=matmul(au,fu*u)/max(g,tiny(1.0_dp))
      fit%post_v=matmul(av,fv*v)/max(g,tiny(1.0_dp))
      fit%post_product=matmul(flat,fit%mass*ug*vg)/max(g,tiny(1.0_dp));fit%loglik=ll
      fit%status=max(ku%status,kv%status)
   end subroutine wlvmix

   subroutine wtlvmix(y,id,w,u,v,fit,control)
      real(dp),intent(in)::y(:),w(:),u(:),v(:);integer,intent(in)::id(:)
      type(bivariate_mixture_fit),intent(out)::fit;type(kw_control),intent(in),optional::control
      real(dp),allocatable::t(:),s(:),m(:),ws(:),slw(:),au0(:,:),ww(:),du(:),fu(:),fv(:),a3(:,:,:),flat(:,:)
      real(dp),allocatable::aum(:,:),avm(:,:),g(:),ug(:),vg(:)
      type(kw_result)::ku;type(mixture_fit)::vf;real(dp)::r,sd,logk
      integer::i,j,l,c,n
      call group_stats(y,id,w,t,s,m,ws,slw);n=size(s);allocate(au0(n,size(u)),ww(n),du(size(u)));ww=1.0_dp/real(n,dp);du=1
      do j=1,size(u);do i=1,n
         sd=sqrt(s(i)/ws(i));au0(i,j)=student_pdf((t(i)-u(j))/sd,m(i)-1.0_dp)/sd
      end do;end do
      call gvmix(s,m,v,vf,ww,control);fv=vf%mass;call kw_fit(au0,du,ww,ku,control);fu=ku%f
      allocate(a3(n,size(v),size(u)))
      do j=1,size(u);do l=1,size(v);do i=1,n
         r=0.5_dp*(m(i)-1.0_dp);sd=sqrt(v(l)/ws(i));a3(i,l,j)=gamma_pdf(s(i),r,v(l)/r)*normal_pdf(t(i),u(j),sd)
      end do;end do;end do
      allocate(aum(n,size(u)),avm(n,size(v)),g(n),flat(n,size(u)*size(v)),ug(size(u)*size(v)),vg(size(u)*size(v)))
      do i=1,n;do j=1,size(u);aum(i,j)=sum(a3(i,:,j)*fv);end do;do l=1,size(v);avm(i,l)=sum(a3(i,l,:)*fu);end do;end do
      g=matmul(aum,fu);c=0
      do l=1,size(v);do j=1,size(u);c=c+1;flat(:,c)=a3(:,l,j);ug(c)=u(j);vg(c)=v(l);end do;end do
      allocate(fit%u(size(u)),fit%v(size(v)),fit%mass(size(u)*size(v)))
      allocate(fit%g(n),fit%post_u(n),fit%post_v(n),fit%post_product(n),fit%a(n,size(u)*size(v)))
      fit%u=u;fit%v=v;fit%g=g;fit%a=flat;c=0
      do l=1,size(v);do j=1,size(u);c=c+1;fit%mass(c)=fu(j)*fv(l);end do;end do
      fit%post_u=matmul(aum,fu*u)/max(g,tiny(1.0_dp));fit%post_v=matmul(avm,fv*v)/max(g,tiny(1.0_dp))
      fit%post_product=matmul(flat,fit%mass*ug*vg)/max(g,tiny(1.0_dp))
      logk=0
      do i=1,n;r=0.5_dp*(m(i)-1);logk=logk+log_gamma(r)-r*log(r)-0.5_dp*log(ws(i))-r*log(2*pi) &
         -(r-1)*log(max(s(i),tiny(1.0_dp)))+0.5_dp*slw(i);end do
      fit%loglik=sum(safe_log(g))+logk;fit%status=max(ku%status,vf%status)
   end subroutine wtlvmix
end module rebayes_repeated

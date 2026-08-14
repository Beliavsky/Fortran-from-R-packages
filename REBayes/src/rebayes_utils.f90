module rebayes_utils
   use rebayes_kinds, only : dp
   use rebayes_math, only : normal_pdf, sample_discrete, normalize_prob
   implicit none
   private
   public :: kw_smooth, kw_quantiles, bw_kw, kw_random_sample
   public :: posterior_summary, lfdr_1d, thresh_fdr, l1_step_distance
   public :: kw2_marginal_quantiles, bw_kw2, finv
   abstract interface
      function scalar_fun(x) result(y)
         import dp
         real(dp),intent(in)::x
         real(dp)::y
      end function scalar_fun
   end interface
contains
   subroutine kw_smooth(x,mass,bw,k,density)
      real(dp),intent(in)::x(:),mass(:),bw
      integer,intent(in)::k
      real(dp),intent(out)::density(:)
      integer::i,j,n
      real(dp)::t,c,dx
      n=size(x);if(size(mass)/=n.or.size(density)/=n)error stop "kw_smooth: size"
      density=0.0_dp
      do i=1,n
         do j=1,n
            t=(x(j)-x(i))/bw
            select case(k)
            case(1)
               c=normal_pdf(t,0.0_dp,1.0_dp)
            case(2)
               if(abs(t)<1.0_dp)then;c=15.0_dp/16.0_dp*(1.0_dp-t*t)**2;else;c=0.0_dp;end if
            case default
               if(abs(t)<1.0_dp)then;c=35.0_dp/32.0_dp*(1.0_dp-t*t)**3;else;c=0.0_dp;end if
            end select
            density(i)=density(i)+mass(j)*c/bw
         end do
      end do
      if(n>1)then
         dx=x(2)-x(1)
         if(abs(dx)>tiny(1.0_dp))density=density/max(sum(density)*abs(dx),tiny(1.0_dp))
      end if
   end subroutine kw_smooth

   subroutine kw_quantiles(x,mass,q,out)
      real(dp),intent(in)::x(:),mass(:),q(:)
      real(dp),intent(out)::out(:)
      real(dp),allocatable::p(:)
      integer::i,j,n
      n=size(x);if(size(mass)/=n.or.size(out)/=size(q))error stop "kw_quantiles: size"
      allocate(p(n));p=max(mass,0.0_dp);call normalize_prob(p)
      do i=2,n;p(i)=p(i)+p(i-1);end do
      do j=1,size(q)
         out(j)=x(n)
         do i=1,n
            if(p(i)>=q(j))then;out(j)=x(i);exit;end if
         end do
      end do
   end subroutine kw_quantiles

   real(dp) function bw_kw(x,mass,k,minbw) result(bw)
      real(dp),intent(in)::x(:),mass(:)
      real(dp),intent(in),optional::k,minbw
      real(dp)::kk,mb,med(1),q(1)
      kk=1.0_dp;if(present(k))kk=k;mb=0.1_dp;if(present(minbw))mb=minbw
      q=0.5_dp;call kw_quantiles(x,mass,q,med)
      bw=max(kk*sum(abs(x-med(1))*mass),mb)
   end function bw_kw

   subroutine kw_random_sample(x,mass,out)
      real(dp),intent(in)::x(:),mass(:)
      real(dp),intent(out)::out(:)
      real(dp),allocatable::p(:)
      integer::i
      allocate(p(size(mass)));p=max(mass,0.0_dp);call normalize_prob(p)
      do i=1,size(out);out(i)=x(sample_discrete(p));end do
   end subroutine kw_random_sample

   subroutine posterior_summary(a,grid,mass,loss,out)
      real(dp),intent(in)::a(:,:),grid(:),mass(:),loss
      real(dp),intent(out)::out(:)
      real(dp),allocatable::post(:),den(:)
      integer::i,j,n,m,idx
      n=size(a,1);m=size(a,2)
      if(size(grid)/=m.or.size(mass)/=m.or.size(out)/=n)error stop "posterior_summary: size"
      allocate(post(m),den(n));den=matmul(a,mass)
      do i=1,n
         post=a(i,:)*mass/max(den(i),tiny(1.0_dp));call normalize_prob(post)
         if(abs(loss-2.0_dp)<100*epsilon(1.0_dp))then
            out(i)=sum(post*grid)
         else if(abs(loss)<100*epsilon(1.0_dp))then
            idx=maxloc(post,dim=1);out(i)=grid(idx)
         else if(loss>0.0_dp.and.loss<=1.0_dp)then
            idx=1
            do j=2,m;post(j)=post(j)+post(j-1);end do
            do j=1,m;if(post(j)>=merge(0.5_dp,loss,abs(loss-1.0_dp)<100*epsilon(1.0_dp)))then;idx=j;exit;end if;end do
            out(i)=grid(idx)
         else
            error stop "posterior_summary: unsupported loss"
         end if
      end do
   end subroutine posterior_summary

   subroutine lfdr_1d(a,grid,mass,cnull,right_tail,out)
      real(dp),intent(in)::a(:,:),grid(:),mass(:),cnull
      logical,intent(in)::right_tail
      real(dp),intent(out)::out(:)
      real(dp),allocatable::nullmass(:),den(:)
      integer::j
      allocate(nullmass(size(mass)),den(size(a,1)));nullmass=mass
      do j=1,size(grid)
         if(right_tail)then
            if(grid(j)>=cnull)nullmass(j)=0.0_dp
         else
            if(grid(j)<cnull)nullmass(j)=0.0_dp
         end if
      end do
      den=matmul(a,mass)
      out=1.0_dp-matmul(a,nullmass)/max(den,tiny(1.0_dp))
   end subroutine lfdr_1d

   real(dp) function thresh_fdr(lambda,stat,v) result(ans)
      real(dp),intent(in)::lambda,stat(:),v(:)
      real(dp)::num,den
      integer::i
      if(size(stat)/=size(v))error stop "thresh_fdr: size"
      num=0.0_dp;den=0.0_dp
      do i=1,size(stat)
         if(stat(i)>lambda)then;num=num+1.0_dp-v(i);den=den+1.0_dp;end if
      end do
      if(den>0.0_dp)then;ans=num/den;else;ans=0.0_dp;end if
   end function thresh_fdr

   real(dp) function l1_step_distance(xf,ff,xg,fg) result(ans)
      real(dp),intent(in)::xf(:),ff(:),xg(:),fg(:)
      real(dp),allocatable::x(:)
      real(dp)::fval,gval
      integer::i,j,n,k
      allocate(x(size(xf)+size(xg)));x=[xf,xg];call sort_unique(x,n)
      ans=0.0_dp
      do i=1,n-1
         fval=0.0_dp;gval=0.0_dp
         do j=1,size(xf);if(xf(j)<=x(i))fval=ff(j);end do
         do k=1,size(xg);if(xg(k)<=x(i))gval=fg(k);end do
         ans=ans+abs(fval-gval)*(x(i+1)-x(i))
      end do
   end function l1_step_distance

   subroutine sort_unique(x,nout)
      real(dp),intent(inout)::x(:);integer,intent(out)::nout
      real(dp)::key;integer::i,j,k
      do i=2,size(x);key=x(i);j=i-1;do while(j>=1);if(x(j)<=key)exit;x(j+1)=x(j);j=j-1;end do;x(j+1)=key;end do
      k=1
      do i=2,size(x)
         if(abs(x(i)-x(k))>100*epsilon(1.0_dp)*max(1.0_dp,abs(x(i))))then;k=k+1;x(k)=x(i);end if
      end do
      nout=k
   end subroutine sort_unique

   subroutine finv(target,f,interval,root,status,tol,max_iter)
      real(dp),intent(in)::target,interval(2)
      procedure(scalar_fun)::f
      real(dp),intent(out)::root
      integer,intent(out)::status
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::max_iter
      real(dp)::a,b,fa,fb,mid,fm,tt,width
      integer::it,nit,k
      tt=1.0e-10_dp;if(present(tol))tt=tol;nit=100;if(present(max_iter))nit=max_iter
      a=interval(1);b=interval(2);fa=f(a)-target;fb=f(b)-target
      do k=1,20
         if(fa*fb<=0.0_dp)exit
         width=b-a
         a=a-width;b=b+width;fa=f(a)-target;fb=f(b)-target
      end do
      if(fa*fb>0.0_dp)then;root=0.5_dp*(a+b);status=1;return;end if
      do it=1,nit
         mid=0.5_dp*(a+b);fm=f(mid)-target
         if(abs(fm)<tt.or.abs(b-a)<tt*(1.0_dp+abs(mid)))exit
         if(fa*fm<=0.0_dp)then;b=mid;fb=fm;else;a=mid;fa=fm;end if
      end do
      root=0.5_dp*(a+b);status=merge(0,1,it<=nit)
   end subroutine finv

   subroutine kw2_marginal_quantiles(u,v,mass,q,out)
      real(dp),intent(in)::u(:),v(:),mass(:),q(:)
      real(dp),intent(out)::out(:)
      real(dp),allocatable::marg(:)
      integer::i,j,c
      if(size(mass)/=size(u)*size(v))error stop "kw2_marginal_quantiles: mass"
      allocate(marg(size(u)));marg=0.0_dp;c=0
      do j=1,size(v);do i=1,size(u);c=c+1;marg(i)=marg(i)+mass(c);end do;end do
      call kw_quantiles(u,marg,q,out)
   end subroutine kw2_marginal_quantiles

   subroutine bw_kw2(u,v,mass,bw)
      real(dp),intent(in)::u(:),v(:),mass(:);real(dp),intent(out)::bw(2)
      real(dp),allocatable::mu(:),mv(:);integer::i,j,c
      if(size(mass)/=size(u)*size(v))error stop "bw_kw2: mass"
      allocate(mu(size(u)),mv(size(v)));mu=0.0_dp;mv=0.0_dp;c=0
      do j=1,size(v);do i=1,size(u);c=c+1;mu(i)=mu(i)+mass(c);mv(j)=mv(j)+mass(c);end do;end do
      bw(1)=bw_kw(u,mu);bw(2)=bw_kw(v,mv)
   end subroutine bw_kw2
end module rebayes_utils

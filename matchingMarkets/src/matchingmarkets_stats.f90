module matchingmarkets_stats
   use matchingmarkets_kinds, only : dp
   use matchingmarkets_types, only : khb_result_t
   use matchingmarkets_utils, only : normal_cdf, normal_pdf
   implicit none
   private
   public :: ols_fit, probit_fit, khb, inverse_mills_ratio, normal_cdf, normal_pdf

   interface
      subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
         integer,intent(in)::n,nrhs,lda,ldb
         integer,intent(out)::ipiv(*)
         integer,intent(out)::info
         double precision,intent(inout)::a(lda,*),b(ldb,*)
      end subroutine dgesv
   end interface
contains
   subroutine ols_fit(x,y,beta,vcov,resid,ok)
      real(dp),intent(in)::x(:,:),y(:)
      real(dp),allocatable,intent(out)::beta(:),vcov(:,:),resid(:)
      logical,intent(out)::ok
      real(dp),allocatable::xtx(:,:),rhs(:,:),inv(:,:),eye(:,:)
      integer,allocatable::ipiv(:)
      integer::n,p,info,j
      n=size(x,1);p=size(x,2);ok=.false.
      if(size(y)/=n) return
      allocate(xtx(p,p),rhs(p,1),ipiv(p));xtx=matmul(transpose(x),x);rhs(:,1)=matmul(transpose(x),y)
      call dgesv(p,1,xtx,p,ipiv,rhs,p,info);if(info/=0)return
      allocate(beta(p));beta=rhs(:,1);allocate(resid(n));resid=y-matmul(x,beta)
      allocate(inv(p,p),eye(p,p));eye=0.0_dp;do j=1,p;eye(j,j)=1.0_dp;end do
      xtx=matmul(transpose(x),x);call dgesv(p,p,xtx,p,ipiv,eye,p,info);if(info/=0)return
      inv=eye;allocate(vcov(p,p));vcov=inv*sum(resid*resid)/real(max(1,n-p),dp);ok=.true.
   end subroutine ols_fit

   subroutine probit_fit(x,y,beta,vcov,converged,max_iter,tol)
      real(dp),intent(in)::x(:,:),y(:)
      real(dp),allocatable,intent(out)::beta(:),vcov(:,:)
      logical,intent(out)::converged
      integer,intent(in),optional::max_iter
      real(dp),intent(in),optional::tol
      integer::n,p,it,mit,i,j,info
      real(dp)::eps,eta,phi,pv,w,z,stepnorm
      real(dp),allocatable::h(:,:),g(:,:),bnew(:),eye(:,:)
      integer,allocatable::ipiv(:)
      n=size(x,1);p=size(x,2);mit=100;if(present(max_iter))mit=max_iter
      eps=1.0e-9_dp;if(present(tol))eps=tol
      allocate(beta(p));beta=0.0_dp;allocate(ipiv(p));converged=.false.
      do it=1,mit
         allocate(h(p,p),g(p,1));h=0.0_dp;g=0.0_dp
         do i=1,n
            eta=dot_product(x(i,:),beta);phi=max(normal_pdf(eta),tiny(1.0_dp));pv=normal_cdf(eta)
            pv=max(1.0e-12_dp,min(1.0_dp-1.0e-12_dp,pv))
            w=phi*phi/(pv*(1.0_dp-pv))
            z=(y(i)-pv)*phi/(pv*(1.0_dp-pv))
            do j=1,p
               g(j,1)=g(j,1)+x(i,j)*z
            end do
            h=h+w*outer(x(i,:),x(i,:))
         end do
         call dgesv(p,1,h,p,ipiv,g,p,info);if(info/=0)exit
         allocate(bnew(p));bnew=beta+g(:,1);stepnorm=maxval(abs(bnew-beta));beta=bnew
         deallocate(h,g,bnew)
         if(stepnorm<eps)then;converged=.true.;exit;end if
      end do
      if(allocated(h)) deallocate(h)
      if(allocated(g)) deallocate(g)
      allocate(h(p,p));h=0.0_dp
      do i=1,n
         eta=dot_product(x(i,:),beta);phi=max(normal_pdf(eta),tiny(1.0_dp));pv=normal_cdf(eta)
         pv=max(1.0e-12_dp,min(1.0_dp-1.0e-12_dp,pv));w=phi*phi/(pv*(1.0_dp-pv))
         h=h+w*outer(x(i,:),x(i,:))
      end do
      allocate(eye(p,p));eye=0.0_dp;do j=1,p;eye(j,j)=1.0_dp;end do
      call dgesv(p,p,h,p,ipiv,eye,p,info);allocate(vcov(p,p));if(info==0)then;vcov=eye;else;vcov=0.0_dp;end if
   end subroutine probit_fit

   function khb(x,y,z_index) result(res)
      real(dp),intent(in)::x(:,:),y(:)
      integer,intent(in)::z_index
      type(khb_result_t)::res
      integer::n,p,q,j,k
      real(dp),allocatable::xr(:,:),xa(:,:),bf(:),vf(:,:),br(:),vr(:,:),ba(:),va(:,:),ra(:),z(:),xfs(:,:),bfs(:),vfs(:,:)
      logical::ok,cf,cr,cfs
      real(dp)::num,den,zstat
      n=size(x,1);p=size(x,2);if(z_index<1.or.z_index>p)error stop 'khb: bad z_index'
      q=p-1;allocate(xr(n,q));k=0
      do j=1,p;if(j==z_index)cycle;k=k+1;xr(:,k)=x(:,j);end do
      call probit_fit(x,y,bf,vf,cf);call probit_fit(xr,y,br,vr,cr)
      allocate(xa(n,q));xa=xr;z=x(:,z_index);call ols_fit(xa,z,ba,va,ra,ok)
      allocate(xfs(n,q+1));xfs(:,1:q)=xr;xfs(:,q+1)=ra
      call probit_fit(xfs,y,bfs,vfs,cfs)
      allocate(res%p_value(q),res%reduced_coef(q),res%full_coef(q),res%rescaled_coef(q))
      res%reduced_coef=br;res%rescaled_coef=bfs(:q);k=0
      do j=1,p;if(j==z_index)cycle;k=k+1;res%full_coef(k)=bf(j);end do
      k=0
      do j=1,p
         if(j==z_index)cycle;k=k+1
         num=bfs(k)-bf(j)
         den=sqrt(max(1.0e-30_dp,bf(z_index)**2*va(k,k)**2+ba(k)**2*vf(z_index,z_index)**2))
         zstat=num/den;res%p_value(k)=1.0_dp-normal_cdf(zstat)
      end do
      res%converged=cf.and.cr.and.cfs.and.ok
   end function khb

   function inverse_mills_ratio(x,binary_one) result(v)
      real(dp),intent(in)::x(:)
      logical,intent(in),optional::binary_one
      real(dp),allocatable::v(:)
      logical::one;integer::i
      one=.true.;if(present(binary_one))one=binary_one;allocate(v(size(x)))
      do i=1,size(x)
         if(one)then;v(i)=normal_pdf(x(i))/max(normal_cdf(x(i)),1.0e-14_dp)
         else;v(i)=-normal_pdf(x(i))/max(1.0_dp-normal_cdf(x(i)),1.0e-14_dp);end if
      end do
   end function inverse_mills_ratio

   pure function outer(a,b) result(c)
      real(dp),intent(in)::a(:),b(:);real(dp)::c(size(a),size(b));integer::i,j
      do j=1,size(b);do i=1,size(a);c(i,j)=a(i)*b(j);end do;end do
   end function outer
end module matchingmarkets_stats

! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_covariance
   use robustbase_kinds, only: dp, huge_penalty
   use robustbase_sort, only: median, sort_real_with_index
   use robustbase_scale, only: qn_scale, mad_scale
   use robustbase_linalg, only: covariance_matrix, symmetric_eigen, invert_symmetric
   use robustbase_probability, only: chi_square_quantile
   implicit none
   private
   public :: robust_cov_result, comedian, cov_comedian, cov_comed, cov_gk, cov_ogk, cov_mcd, &
             robust_mahalanobis, adjusted_outlyingness
   type :: robust_cov_result
      real(dp), allocatable :: center(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: raw_center(:)
      real(dp), allocatable :: raw_covariance(:,:)
      real(dp), allocatable :: distances(:)
      logical, allocatable :: weights(:)
      integer :: h = 0
      integer :: iterations = 0
      logical :: converged = .false.
   end type
contains
   function comedian(x,y) result(c)
      real(dp),intent(in)::x(:),y(:)
      real(dp)::c,mx,my
      if(size(x)/=size(y)) error stop "comedian: size mismatch"
      mx=median(x);my=median(y);c=median((x-mx)*(y-my))
   end function

   subroutine cov_comedian(x,center,cov)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(out)::center(:),cov(:,:)
      integer::p,j,k
      p=size(x,2)
      if(size(center)/=p .or. any(shape(cov)/=[p,p])) error stop "cov_comedian: size mismatch"
      do j=1,p;center(j)=median(x(:,j));end do
      do j=1,p
         do k=j,p
            cov(j,k)=comedian(x(:,j),x(:,k));cov(k,j)=cov(j,k)
         end do
      end do
   end subroutine

   subroutine cov_comed(x,result,n_iter,reweight)
      real(dp),intent(in)::x(:,:)
      type(robust_cov_result),intent(out)::result
      integer,intent(in),optional::n_iter
      logical,intent(in),optional::reweight
      integer::n,p,iter,ni,j,k,info,nw
      real(dp),allocatable::rho(:,:),vals(:),u(:,:),q(:,:),z(:,:),mads(:),medz(:),invc(:,:),d2(:)
      real(dp)::cut
      logical::rw
      n=size(x,1);p=size(x,2);ni=2;if(present(n_iter))ni=max(0,n_iter);rw=.true.;if(present(reweight))rw=reweight
      allocate(rho(p,p),vals(p),u(p,p),q(p,p),z(n,p),mads(p),medz(p),invc(p,p),d2(n))
      do j=1,p
         mads(j)=mad_scale(x(:,j));if(mads(j)<=1.0e-14_dp)mads(j)=1.0_dp
      end do
      do j=1,p
         do k=j,p
            rho(j,k)=1.4826_dp**2*comedian(x(:,j),x(:,k))/(mads(j)*mads(k));rho(k,j)=rho(j,k)
         end do
      end do
      do iter=0,ni
         call symmetric_eigen(rho,vals,u,info);if(info/=0)error stop 'cov_comed: eigen failure'
         q=spread(mads,2,p)*u
         z=matmul(x,spread(1.0_dp/mads,2,p)*u)
         do j=1,p
            mads(j)=mad_scale(z(:,j));if(mads(j)<=1.0e-14_dp)mads(j)=1.0_dp
         end do
         rho=matmul(q*spread(mads*mads,1,p),transpose(q))
      end do
      allocate(result%center(p),result%covariance(p,p),result%raw_center(p),result%raw_covariance(p,p),result%distances(n),result%weights(n))
      do j=1,p;medz(j)=median(z(:,j));end do
      result%raw_center=matmul(q,medz);result%raw_covariance=rho
      call invert_symmetric(rho,invc,info,ridge=1.0e-12_dp)
      do j=1,n;d2(j)=dot_product(x(j,:)-result%raw_center,matmul(invc,x(j,:)-result%raw_center));end do
      cut=chi_square_quantile(0.975_dp,real(p,dp));result%weights=d2<=cut;result%distances=sqrt(max(d2,0.0_dp));result%center=result%raw_center;result%covariance=result%raw_covariance
      if(rw)then
         nw=count(result%weights)
         if(nw>p)then
            do j=1,p;result%center(j)=sum(pack(x(:,j),result%weights))/real(nw,dp);end do
            call covariance_matrix(pack_rows_local(x,result%weights),result%center,result%covariance)
         end if
      end if
      result%h=count(result%weights);result%iterations=ni;result%converged=.true.
   contains
      function pack_rows_local(a,mask) result(b)
         real(dp),intent(in)::a(:,:);logical,intent(in)::mask(:)
         real(dp),allocatable::b(:,:);integer::jj
         allocate(b(count(mask),size(a,2)))
         do jj=1,size(a,2);b(:,jj)=pack(a(:,jj),mask);end do
      end function
   end subroutine cov_comed

   function cov_gk(x,y) result(c)
      real(dp),intent(in)::x(:),y(:)
      real(dp)::c,sx,sy,sp,sm
      real(dp),allocatable::u(:),v(:)
      if(size(x)/=size(y)) error stop "cov_gk: size mismatch"
      sx=qn_scale(x);sy=qn_scale(y)
      if(sx<=0.0_dp .or. sy<=0.0_dp) then;c=0.0_dp;return;end if
      u=(x-median(x))/sx;v=(y-median(y))/sy
      sp=qn_scale(u+v,finite_correction=.false.)
      sm=qn_scale(u-v,finite_correction=.false.)
      c=0.25_dp*(sp*sp-sm*sm)*sx*sy
   end function

   subroutine cov_ogk(x,result,reweight)
      real(dp),intent(in)::x(:,:)
      type(robust_cov_result),intent(out)::result
      logical,intent(in),optional::reweight
      integer::n,p,j,k,info,nw
      real(dp),allocatable::z(:,:),g(:,:),vals(:),vecs(:,:),proj(:,:),sc(:),cproj(:),invc(:,:),dist(:)
      real(dp)::cut
      logical::rw
      n=size(x,1);p=size(x,2);rw=.true.;if(present(reweight))rw=reweight
      allocate(result%center(p),result%covariance(p,p),result%raw_center(p),result%raw_covariance(p,p),result%distances(n),result%weights(n))
      allocate(z(n,p),g(p,p),vals(p),vecs(p,p),proj(n,p),sc(p),cproj(p),invc(p,p),dist(n))
      do j=1,p
         result%raw_center(j)=median(x(:,j));sc(j)=qn_scale(x(:,j))
         if(sc(j)<=sqrt(epsilon(1.0_dp))) sc(j)=1.0_dp
         z(:,j)=(x(:,j)-result%raw_center(j))/sc(j)
      end do
      g=0.0_dp
      do j=1,p
         g(j,j)=1.0_dp
         do k=j+1,p
            g(j,k)=cov_gk(z(:,j),z(:,k));g(k,j)=g(j,k)
         end do
      end do
      call symmetric_eigen(g,vals,vecs,info)
      if(info/=0) error stop "cov_ogk: eigen failure"
      proj=matmul(z,vecs)
      do j=1,p;cproj(j)=median(proj(:,j));vals(j)=max(qn_scale(proj(:,j))**2,1.0e-12_dp);end do
      result%raw_covariance=0.0_dp
      do j=1,p
         do k=1,p
            result%raw_covariance(j,k)=sc(j)*sc(k)*sum(vecs(j,:)*vals*vecs(k,:))
         end do
      end do
      result%raw_center=result%raw_center+sc*matmul(vecs,cproj)
      result%center=result%raw_center;result%covariance=result%raw_covariance
      call invert_symmetric(result%raw_covariance,invc,info,ridge=1.0e-12_dp)
      do j=1,n
         dist(j)=dot_product(x(j,:)-result%raw_center,matmul(invc,x(j,:)-result%raw_center))
      end do
      result%distances=sqrt(max(dist,0.0_dp));cut=chi_square_quantile(0.975_dp,real(p,dp));result%weights=dist<=cut
      if(rw) then
         nw=count(result%weights)
         if(nw>p) then
            do j=1,p;result%center(j)=sum(pack(x(:,j),result%weights))/real(nw,dp);end do
            call covariance_matrix(pack_rows(x,result%weights),result%center,result%covariance)
         end if
      end if
      result%h=count(result%weights);result%converged=.true.;result%iterations=1
   contains
      function pack_rows(a,mask) result(b)
         real(dp),intent(in)::a(:,:);logical,intent(in)::mask(:)
         real(dp),allocatable::b(:,:);integer::jj
         allocate(b(count(mask),size(a,2)))
         do jj=1,size(a,2);b(:,jj)=pack(a(:,jj),mask);end do
      end function
   end subroutine

   subroutine cov_mcd(x,result,alpha,n_starts,max_csteps)
      real(dp),intent(in)::x(:,:)
      type(robust_cov_result),intent(out)::result
      real(dp),intent(in),optional::alpha
      integer,intent(in),optional::n_starts,max_csteps
      integer::n,p,h,ns,mc,s,i,j,info,it,nw
      integer,allocatable::subset(:),idx(:),best_subset(:)
      real(dp),allocatable::center(:),cov(:,:),invc(:,:),d(:),ds(:),vals(:),vecs(:,:),bestc(:),bestcov(:,:)
      real(dp)::a,bestdet,det,cut
      logical,allocatable::mask(:)
      n=size(x,1);p=size(x,2);a=0.75_dp;if(present(alpha))a=alpha
      h=max(p+1,min(n,int(floor(a*real(n,dp)))));ns=100;if(present(n_starts))ns=n_starts;mc=30;if(present(max_csteps))mc=max_csteps
      allocate(subset(h),idx(n),best_subset(h),center(p),cov(p,p),invc(p,p),d(n),ds(n),vals(p),vecs(p,p),bestc(p),bestcov(p,p),mask(n))
      bestdet=huge_penalty;best_subset=0
      do s=1,ns
         call random_subset(n,p+1,subset(1:p+1))
         call subset_stats(x,subset(1:p+1),center,cov)
         do it=1,mc
            call invert_symmetric(cov,invc,info,ridge=1.0e-10_dp);if(info/=0)exit
            do i=1,n;d(i)=dot_product(x(i,:)-center,matmul(invc,x(i,:)-center));end do
            ds=d;call sort_real_with_index(ds,idx);subset=idx(1:h)
            call subset_stats(x,subset,center,cov)
         end do
         call symmetric_eigen(cov,vals,vecs,info)
         if(info==0 .and. all(vals>1.0e-14_dp)) then
            det=product(vals)
            if(det<bestdet) then;bestdet=det;best_subset=subset;bestc=center;bestcov=cov;end if
         end if
      end do
      if(any(best_subset==0)) then
         best_subset=[(i,i=1,h)];call subset_stats(x,best_subset,bestc,bestcov)
      end if
      allocate(result%center(p),result%covariance(p,p),result%raw_center(p),result%raw_covariance(p,p),result%distances(n),result%weights(n))
      result%raw_center=bestc;result%raw_covariance=bestcov
      call invert_symmetric(bestcov,invc,info,ridge=1.0e-10_dp)
      do i=1,n;d(i)=dot_product(x(i,:)-bestc,matmul(invc,x(i,:)-bestc));end do
      cut=chi_square_quantile(0.975_dp,real(p,dp));mask=d<=cut;nw=count(mask)
      result%weights=mask;result%distances=sqrt(max(d,0.0_dp));result%center=bestc;result%covariance=bestcov
      if(nw>p) then
         do j=1,p;result%center(j)=sum(pack(x(:,j),mask))/real(nw,dp);end do
         call covariance_matrix(pack_rows(x,mask),result%center,result%covariance)
      end if
      result%h=h;result%iterations=mc;result%converged=.true.
   contains
      subroutine random_subset(nn,kk,out)
         integer,intent(in)::nn,kk;integer,intent(out)::out(:)
         integer::ii,cand
         real(dp)::rr
         if(size(out)/=kk) error stop "random_subset size"
         ii=0
         do while(ii<kk)
            call random_number(rr);cand=1+int(rr*real(nn,dp));cand=min(nn,cand)
            if(ii==0 .or. .not.any(out(1:ii)==cand)) then;ii=ii+1;out(ii)=cand;end if
         end do
      end subroutine
      subroutine subset_stats(a,ind,c,cv)
         real(dp),intent(in)::a(:,:);integer,intent(in)::ind(:);real(dp),intent(out)::c(:),cv(:,:)
         real(dp),allocatable::b(:,:);integer::jj
         allocate(b(size(ind),size(a,2)))
         do jj=1,size(a,2);b(:,jj)=a(ind,jj);c(jj)=sum(b(:,jj))/real(size(ind),dp);end do
         call covariance_matrix(b,c,cv)
      end subroutine
      function pack_rows(a,m) result(b)
         real(dp),intent(in)::a(:,:);logical,intent(in)::m(:);real(dp),allocatable::b(:,:);integer::jj
         allocate(b(count(m),size(a,2)));do jj=1,size(a,2);b(:,jj)=pack(a(:,jj),m);end do
      end function
   end subroutine

   subroutine robust_mahalanobis(x,center,cov,distances)
      real(dp),intent(in)::x(:,:),center(:),cov(:,:)
      real(dp),intent(out)::distances(:)
      real(dp),allocatable::inv(:,:)
      integer::i,info
      allocate(inv(size(cov,1),size(cov,2)));call invert_symmetric(cov,inv,info,ridge=1.0e-12_dp)
      if(info/=0) error stop "robust_mahalanobis: inverse failure"
      do i=1,size(x,1);distances(i)=sqrt(max(0.0_dp,dot_product(x(i,:)-center,matmul(inv,x(i,:)-center))));end do
   end subroutine

   subroutine adjusted_outlyingness(x,outlyingness,center,n_directions)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(out)::outlyingness(:),center(:)
      integer,intent(in),optional::n_directions
      integer::n,p,nd,d,i,j
      real(dp),allocatable::dir(:),proj(:)
      real(dp)::normd,med,sl,sr,z,u
      n=size(x,1);p=size(x,2);nd=max(n,100);if(present(n_directions))nd=n_directions
      if(size(outlyingness)/=n .or. size(center)/=p) error stop "adjusted_outlyingness: size mismatch"
      do j=1,p;center(j)=median(x(:,j));end do
      allocate(dir(p),proj(n));outlyingness=0.0_dp
      do d=1,nd
         if(d<=n) then;dir=x(d,:)-center;else;call random_number(dir);dir=2.0_dp*dir-1.0_dp;end if
         normd=sqrt(sum(dir*dir));if(normd<=1.0e-14_dp)cycle;dir=dir/normd;proj=matmul(x-spread(center,1,n),dir);med=median(proj)
         sl=qn_scale(pack(proj,proj<=med));sr=qn_scale(pack(proj,proj>=med));if(sl<=0.0_dp)sl=mad_scale(proj);if(sr<=0.0_dp)sr=mad_scale(proj)
         do i=1,n
            if(proj(i)>=med) then;u=max(sr,1.0e-12_dp);else;u=max(sl,1.0e-12_dp);end if
            z=abs(proj(i)-med)/u;outlyingness(i)=max(outlyingness(i),z)
         end do
      end do
   end subroutine
end module robustbase_covariance

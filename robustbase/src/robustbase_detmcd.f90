! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_detmcd
   use robustbase_kinds, only: dp, huge_penalty
   use robustbase_sort, only: median, sort_real_with_index
   use robustbase_scale, only: qn_scale
   use robustbase_probability, only: normal_quantile, chi_square_quantile
   use robustbase_linalg, only: covariance_matrix, symmetric_eigen, invert_symmetric
   use robustbase_covariance, only: robust_cov_result, cov_ogk
   implicit none
   private
   public :: detmcd_result, cov_detmcd

   type :: detmcd_result
      type(robust_cov_result) :: estimate
      integer, allocatable :: best_subset(:)
      integer, allocatable :: initial_orderings(:,:)
      integer, allocatable :: csteps(:)
      integer :: best_start = 0
      real(dp) :: log_determinant = huge_penalty
      logical :: exact_fit = .false.
      real(dp), allocatable :: hyperplane_coefficients(:)
      real(dp), allocatable :: hyperplane_center(:)
      integer :: points_on_hyperplane = 0
   end type detmcd_result
contains
   subroutine cov_detmcd(x, result, alpha, max_csteps, save_orderings)
      real(dp), intent(in) :: x(:,:)
      type(detmcd_result), intent(out) :: result
      real(dp), intent(in), optional :: alpha
      integer, intent(in), optional :: max_csteps
      logical, intent(in), optional :: save_orderings
      real(dp), allocatable :: z(:,:), loc(:), sc(:), orient(:,:,:), center(:), cov(:,:), best_center(:), best_cov(:,:), &
                              invcov(:,:), dist(:), sorted_dist(:), eigvals(:), eigvecs(:,:), xw(:,:), &
                              exact_center(:), exact_cov(:,:), exact_normal(:), projection(:)
      integer, allocatable :: orderings(:,:), subset(:), previous(:), idx(:), best_subset(:)
      logical, allocatable :: mask(:)
      real(dp) :: a, objective, best_objective, cutoff
      logical :: save_ord, has_exact
      integer :: n, p, h, mc, j, k, s, it, info, nw, exact_count, kmin, candidate_count

      n=size(x,1)
      p=size(x,2)
      if(n<max(2,p+1) .or. p<1) error stop "cov_detmcd: invalid dimensions"
      a=0.75_dp
      if(present(alpha)) a=alpha
      if(a<0.5_dp .or. a>1.0_dp) error stop "cov_detmcd: alpha must be in [0.5,1]"
      h=max(p+1,min(n,h_alpha_n(a,n,p)))
      mc=200
      if(present(max_csteps)) mc=max(1,max_csteps)
      save_ord=.true.
      if(present(save_orderings)) save_ord=save_orderings

      allocate(z(n,p),loc(p),sc(p),orient(p,p,6),orderings(n,6),subset(h),previous(h),idx(n),best_subset(h), &
               center(p),cov(p,p),best_center(p),best_cov(p,p),invcov(p,p),dist(n),sorted_dist(n),eigvals(p),eigvecs(p,p),mask(n),result%csteps(6), &
               exact_center(p),exact_cov(p,p),exact_normal(p),projection(n))
      do j=1,p
         loc(j)=median(x(:,j))
         sc(j)=qn_scale(x(:,j))
         if(sc(j)<=sqrt(epsilon(1.0_dp))) sc(j)=1.0_dp
         z(:,j)=(x(:,j)-loc(j))/sc(j)
      end do
      call six_orientations(z,orient)
      do s=1,6
         call initial_order(z,orient(:,:,s),orderings(:,s))
      end do

      best_objective=huge_penalty
      best_subset=0
      best_center=0.0_dp
      best_cov=0.0_dp
      result%csteps=0
      result%exact_fit=.false.
      has_exact=.false.;exact_count=0;exact_center=0.0_dp;exact_cov=0.0_dp;exact_normal=0.0_dp
      do s=1,6
         subset=orderings(1:h,s)
         previous=0
         objective=huge_penalty
         do it=1,mc
            call subset_stats(z,subset,center,cov)
            call symmetric_eigen(cov,eigvals,eigvecs,info)
            if(info/=0 .or. any(eigvals<=1.0e-14_dp)) then
               objective=huge_penalty
               if(info==0)then
                  kmin=minloc(eigvals,dim=1)
                  projection=matmul(z-spread(center,1,n),eigvecs(:,kmin))
                  candidate_count=count(abs(projection)<=1.0e-8_dp*(1.0_dp+maxval(abs(z))))
                  if(candidate_count>=h .and. candidate_count>exact_count)then
                     has_exact=.true.;exact_count=candidate_count;exact_center=center;exact_cov=cov;exact_normal=eigvecs(:,kmin)
                  end if
               end if
               exit
            end if
            objective=sum(log(eigvals))
            call invert_symmetric(cov,invcov,info,ridge=0.0_dp)
            if(info/=0) exit
            do j=1,n
               dist(j)=dot_product(z(j,:)-center,matmul(invcov,z(j,:)-center))
            end do
            sorted_dist=dist
            call sort_real_with_index(sorted_dist,idx)
            previous=subset
            subset=idx(1:h)
            if(all(subset==previous)) exit
         end do
         result%csteps(s)=min(it,mc)
         if(objective<best_objective) then
            best_objective=objective
            best_subset=subset
            best_center=center
            best_cov=cov
            result%best_start=s
         end if
      end do
      if(has_exact)then
         result%exact_fit=.true.
         best_center=exact_center;best_cov=exact_cov;result%best_start=0
         projection=abs(matmul(z-spread(exact_center,1,n),exact_normal))
         sorted_dist=projection;call sort_real_with_index(sorted_dist,idx);best_subset=idx(1:h)
         best_objective=-huge_penalty
      else if(any(best_subset==0))then
         error stop "cov_detmcd: no nonsingular initial subset"
      end if

      allocate(result%best_subset(h))
      result%best_subset=best_subset
      if(save_ord) then
         allocate(result%initial_orderings(n,6))
         result%initial_orderings=orderings
      else
         allocate(result%initial_orderings(0,0))
      end if
      if(result%exact_fit)then
         result%log_determinant=-huge_penalty
      else
         result%log_determinant=best_objective+2.0_dp*sum(log(sc))
      end if
      allocate(result%estimate%center(p),result%estimate%covariance(p,p),result%estimate%raw_center(p), &
               result%estimate%raw_covariance(p,p),result%estimate%distances(n),result%estimate%weights(n))
      result%estimate%raw_center=loc+sc*best_center
      do j=1,p
         do k=1,p
            result%estimate%raw_covariance(j,k)=best_cov(j,k)*sc(j)*sc(k)
         end do
      end do
      result%estimate%center=result%estimate%raw_center
      result%estimate%covariance=result%estimate%raw_covariance
      allocate(result%hyperplane_coefficients(p),result%hyperplane_center(p))
      result%hyperplane_coefficients=0.0_dp;result%hyperplane_center=result%estimate%raw_center;result%points_on_hyperplane=0
      if(result%exact_fit)then
         result%hyperplane_coefficients=exact_normal/sc
         result%hyperplane_coefficients=result%hyperplane_coefficients/max(sqrt(sum(result%hyperplane_coefficients**2)),1.0e-14_dp)
         result%hyperplane_center=result%estimate%raw_center
         projection=matmul(x-spread(result%hyperplane_center,1,n),result%hyperplane_coefficients)
         mask=abs(projection)<=1.0e-7_dp*(1.0_dp+maxval(abs(x)))
         result%points_on_hyperplane=count(mask)
         result%estimate%weights=mask
         result%estimate%distances=abs(projection)
      else
         call invert_symmetric(result%estimate%raw_covariance,invcov,info,ridge=1.0e-12_dp)
         do j=1,n
            dist(j)=dot_product(x(j,:)-result%estimate%raw_center,matmul(invcov,x(j,:)-result%estimate%raw_center))
         end do
         cutoff=chi_square_quantile(0.975_dp,real(p,dp))
         mask=dist<=cutoff
         result%estimate%weights=mask
         result%estimate%distances=sqrt(max(dist,0.0_dp))
         nw=count(mask)
         if(nw>p) then
            do j=1,p
               result%estimate%center(j)=sum(pack(x(:,j),mask))/real(nw,dp)
            end do
            allocate(xw(nw,p))
            do j=1,p
               xw(:,j)=pack(x(:,j),mask)
            end do
            call covariance_matrix(xw,result%estimate%center,result%estimate%covariance)
         end if
      end if
      result%estimate%h=h
      result%estimate%iterations=maxval(result%csteps)
      result%estimate%converged=.true.
   end subroutine cov_detmcd

   integer function h_alpha_n(alpha,n,p) result(h)
      real(dp),intent(in)::alpha
      integer,intent(in)::n,p
      integer::n2
      n2=(n+p+1)/2
      h=int(floor(real(2*n2-n,dp)+2.0_dp*real(n-n2,dp)*alpha))
   end function h_alpha_n

   subroutine six_orientations(z,orient)
      real(dp),intent(in)::z(:,:)
      real(dp),intent(out)::orient(:,:,:)
      real(dp),allocatable::work(:,:),cor(:,:),values(:),vectors(:,:),ranks(:,:),norms(:),sorted(:),sub(:,:),center(:),cov(:,:)
      type(robust_cov_result)::ogk
      integer,allocatable::idx(:)
      integer::n,p,j,info,half
      n=size(z,1);p=size(z,2)
      allocate(work(n,p),cor(p,p),values(p),vectors(p,p),ranks(n,p),norms(n),sorted(n),idx(n),center(p),cov(p,p))

      work=tanh(z)
      call correlation(work,cor)
      call symmetric_eigen(cor,values,vectors,info)
      if(info/=0) error stop "cov_detmcd: tanh start failed"
      orient(:,:,1)=vectors

      call column_ranks(z,ranks)
      call correlation(ranks,cor)
      call symmetric_eigen(cor,values,vectors,info)
      if(info/=0) error stop "cov_detmcd: Spearman start failed"
      orient(:,:,2)=vectors

      do j=1,p
         work(:,j)=normal_scores(ranks(:,j),n)
      end do
      call correlation(work,cor)
      call symmetric_eigen(cor,values,vectors,info)
      if(info/=0) error stop "cov_detmcd: normal-score start failed"
      orient(:,:,3)=vectors

      do j=1,n
         norms(j)=sqrt(sum(z(j,:)**2))
         if(norms(j)>epsilon(1.0_dp)) then
            work(j,:)=z(j,:)/norms(j)
         else
            work(j,:)=z(j,:)
         end if
      end do
      cor=matmul(transpose(work),work)
      call symmetric_eigen(cor,values,vectors,info)
      if(info/=0) error stop "cov_detmcd: spatial-sign start failed"
      orient(:,:,4)=vectors

      sorted=norms
      call sort_real_with_index(sorted,idx)
      half=(n+1)/2
      allocate(sub(half,p))
      sub=z(idx(1:half),:)
      do j=1,p
         center(j)=sum(sub(:,j))/real(half,dp)
      end do
      call covariance_matrix(sub,center,cov)
      call symmetric_eigen(cov,values,vectors,info)
      if(info/=0) error stop "cov_detmcd: BACON start failed"
      orient(:,:,5)=vectors

      call cov_ogk(z,ogk,reweight=.false.)
      call symmetric_eigen(ogk%raw_covariance,values,vectors,info)
      if(info/=0) error stop "cov_detmcd: OGK start failed"
      orient(:,:,6)=vectors
   end subroutine six_orientations

   subroutine initial_order(z,p_matrix,ordering)
      real(dp),intent(in)::z(:,:),p_matrix(:,:)
      integer,intent(out)::ordering(:)
      real(dp),allocatable::projected(:,:),lambda(:),sqrtcov(:,:),sqrtinv(:,:),transformed(:,:),estloc(:),score(:,:),distance(:),sorted(:)
      integer::n,p,i,j,k
      n=size(z,1);p=size(z,2)
      allocate(projected(n,p),lambda(p),sqrtcov(p,p),sqrtinv(p,p),transformed(n,p),estloc(p),score(n,p),distance(n),sorted(n))
      projected=matmul(z,p_matrix)
      do j=1,p
         lambda(j)=qn_scale(projected(:,j))
         if(lambda(j)<=sqrt(epsilon(1.0_dp))) lambda(j)=1.0_dp
      end do
      sqrtcov=0.0_dp;sqrtinv=0.0_dp
      do k=1,p
         do i=1,p
            do j=1,p
               sqrtcov(i,j)=sqrtcov(i,j)+p_matrix(i,k)*lambda(k)*p_matrix(j,k)
               sqrtinv(i,j)=sqrtinv(i,j)+p_matrix(i,k)*p_matrix(j,k)/lambda(k)
            end do
         end do
      end do
      transformed=matmul(z,sqrtinv)
      do j=1,p
         estloc(j)=median(transformed(:,j))
      end do
      estloc=matmul(estloc,sqrtcov)
      score=matmul(z-spread(estloc,1,n),p_matrix)
      do i=1,n
         distance(i)=sum((score(i,:)/lambda)**2)
      end do
      sorted=distance
      call sort_real_with_index(sorted,ordering)
   end subroutine initial_order

   subroutine subset_stats(x,subset,center,cov)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::subset(:)
      real(dp),intent(out)::center(:),cov(:,:)
      real(dp),allocatable::xs(:,:)
      integer::j,m
      m=size(subset)
      allocate(xs(m,size(x,2)))
      xs=x(subset,:)
      do j=1,size(x,2)
         center(j)=sum(xs(:,j))/real(m,dp)
      end do
      call covariance_matrix(xs,center,cov)
   end subroutine subset_stats

   subroutine correlation(x,cor)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(out)::cor(:,:)
      real(dp),allocatable::center(:),cov(:,:),sd(:)
      integer::p,j,k
      p=size(x,2)
      allocate(center(p),cov(p,p),sd(p))
      do j=1,p
         center(j)=sum(x(:,j))/real(size(x,1),dp)
      end do
      call covariance_matrix(x,center,cov)
      do j=1,p
         sd(j)=sqrt(max(cov(j,j),epsilon(1.0_dp)))
      end do
      do j=1,p
         do k=1,p
            cor(j,k)=cov(j,k)/(sd(j)*sd(k))
         end do
      end do
      do j=1,p
         cor(j,j)=1.0_dp
      end do
   end subroutine correlation

   subroutine column_ranks(x,ranks)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(out)::ranks(:,:)
      real(dp),allocatable::sorted(:)
      integer,allocatable::idx(:)
      integer::n,p,j,left,right,k
      n=size(x,1);p=size(x,2)
      allocate(sorted(n),idx(n))
      do j=1,p
         sorted=x(:,j)
         call sort_real_with_index(sorted,idx)
         left=1
         do while(left<=n)
            right=left
            do while(right<n)
               if(abs(sorted(right+1)-sorted(left))>epsilon(1.0_dp)*max(1.0_dp,abs(sorted(left)))) exit
               right=right+1
            end do
            do k=left,right
               ranks(idx(k),j)=0.5_dp*real(left+right,dp)
            end do
            left=right+1
         end do
      end do
   end subroutine column_ranks

   elemental function normal_scores(rank,n) result(score)
      real(dp),intent(in)::rank
      integer,intent(in)::n
      real(dp)::score,prob
      prob=(rank-1.0_dp/3.0_dp)/(real(n,dp)+1.0_dp/3.0_dp)
      score=normal_quantile(max(1.0e-12_dp,min(1.0_dp-1.0e-12_dp,prob)))
   end function normal_scores
end module robustbase_detmcd
